from stellargraph.mapper import CorruptedGenerator, HinSAGENodeGenerator
from sklearn.preprocessing import StandardScaler
from stellargraph import StellarGraph
from stellargraph.utils import plot_history
from stellargraph.layer import DeepGraphInfomax, HinSAGE

import pandas as pd
import os
import random
import numpy as np
import tensorflow as tf

from keras.api._v2.keras import Model
from keras.api._v2.keras.optimizers import Adam
from keras.api._v2.keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.api._v2.keras import backend as K

# Reproducibility setup
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
tf.random.set_seed(SEED)
os.environ['PYTHONHASHSEED'] = str(SEED)
os.environ['TF_DETERMINISTIC_OPS'] = '1'
tf.config.experimental.enable_op_determinism()


def run_deep_graph_infomax(base_model, generator, epochs, model_save_path):
    """
    Run Deep Graph Infomax using a given base model and generator.
    Save the trained embedding model for future explanation.
    """
    lncRNA_nodes = LPPI_graph.nodes(node_type='lncRNA')

    # Corrupted generator for DGI
    corrupted_generator = CorruptedGenerator(generator)
    lncRNA_gen = corrupted_generator.flow(lncRNA_nodes, shuffle=True, seed=SEED)
    infomax = DeepGraphInfomax(base_model, corrupted_generator)

    # Build training model
    x_in, x_out = infomax.in_out_tensors()
    model = Model(inputs=x_in, outputs=x_out)

    # Callbacks for early stopping and checkpointing
    es = EarlyStopping(monitor="loss", patience=8, mode='min', verbose=1)
    checkpoint_path = "best_weights_tmp.h5"
    checkpoint = ModelCheckpoint(filepath=checkpoint_path, monitor="loss", save_best_only=True,
                                 save_weights_only=True, mode="min", verbose=0)

    model.compile(loss=tf.nn.sigmoid_cross_entropy_with_logits, optimizer=Adam(learning_rate=1e-2))
    history = model.fit(
        lncRNA_gen,
        epochs=epochs,
        verbose=2,
        callbacks=[es, checkpoint],
        workers=1,
        use_multiprocessing=False
    )
    plot_history(history)

    # Load best weights
    model.load_weights(checkpoint_path)

    # Build the embedding model (encoder only)
    x_emb_in, x_emb_out = base_model.in_out_tensors()
    emb_model = Model(inputs=x_emb_in, outputs=x_emb_out)

    # Save model for future use
    emb_model.save(model_save_path)
    print(f"Saved embedding model to: {model_save_path}")

    # Generate embeddings
    node_gen = generator.flow(lncRNA_nodes, targets=None, shuffle=False)
    embedding = emb_model.predict(node_gen)

    if os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)

    return embedding


def run_deterministic_hinsage(LPPI_graph, layer_sizes, num_samples, epochs, model_save_path, seed=42):
    """
    Wrapper to initialize generator and model and run DGI.
    """
    K.clear_session()
    np.random.seed(seed)
    tf.random.set_seed(seed)
    random.seed(seed)

    generator = HinSAGENodeGenerator(
        LPPI_graph, batch_size=32, num_samples=num_samples,
        head_node_type="lncRNA", seed=seed
    )

    model = HinSAGE(
        layer_sizes=layer_sizes,
        activations=["LeakyReLU", "linear"],
        generator=generator
    )
    return run_deep_graph_infomax(model, generator, epochs=epochs, model_save_path=model_save_path)


# Main loop for each tissue
human_tissue = ['heart', 'lung', 'stomach']
mouse_tissue = ['heart', 'lung', 'brain']
species = 'human'

for tissue in human_tissue:
    print(f"\nProcessing tissue: {tissue}")

    # Load data
    proteins = pd.read_csv(f"../Annotate/{species}/transformed_protein_annotation.csv")
    lncRNAs = pd.read_csv(f"../Annotate/{species}/valid_{tissue}_annotation.csv")
    LPPI = pd.read_csv(f'../Annotate/{species}/weighted_valid_inter.csv')

    # Preprocess feature data
    proteins.set_index(proteins['protein_ID'], inplace=True)
    proteins = proteins.drop('protein_ID', axis=1)
    lncRNAs.set_index(lncRNAs['lncRNA_ID'], inplace=True)
    lncRNAs = lncRNAs.drop('lncRNA_ID', axis=1)

    # Scale lncRNA features
    scaler_lncRNAs = StandardScaler()
    lncRNAs_scaled = scaler_lncRNAs.fit_transform(lncRNAs)

    # Build heterogeneous graph
    LPPI_graph = StellarGraph(
        {"lncRNA": pd.DataFrame(lncRNAs_scaled, index=lncRNAs.index),
         "protein": pd.DataFrame(proteins, index=proteins.index)},
        LPPI
    )
    print(LPPI_graph.info())

    # Output paths
    tissue_tag = f"{species}_{tissue}"
    emb_save_path = f"./{species}/lncRNA_embeddings_{tissue}.csv"
    model_save_path = f"./{species}/trained_hinsage_model_{tissue}.h5"

    # Train and save embedding model
    node_embedding = run_deterministic_hinsage(
        LPPI_graph,
        layer_sizes=[64, 256],
        num_samples=[20, 15],
        epochs=1000,
        model_save_path=model_save_path,
        seed=42
    )

    # Save node embeddings
    embedding_df = pd.DataFrame(node_embedding, index=LPPI_graph.nodes(node_type='lncRNA'))
    embedding_df.to_csv(emb_save_path, header=None)
    print(f"Saved embeddings to: {emb_save_path}")
