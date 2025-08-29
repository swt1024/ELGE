from stellargraph.mapper import  CorruptedGenerator, HinSAGENodeGenerator
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

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
tf.random.set_seed(SEED)
os.environ['PYTHONHASHSEED'] = str(SEED)
os.environ['TF_DETERMINISTIC_OPS'] = '1' 
tf.config.experimental.enable_op_determinism()

def run_deep_graph_infomax(
    base_model, generator, epochs, reorder=lambda sequence, subjects: subjects
):
    lncRNA_nodes = LPPI_graph.nodes(node_type='lncRNA')

    corrupted_generator = CorruptedGenerator(generator)
    lncRNA_gen = corrupted_generator.flow(lncRNA_nodes, shuffle=True, seed=SEED)
    infomax = DeepGraphInfomax(base_model, corrupted_generator)

    x_in, x_out = infomax.in_out_tensors()

    model = Model(inputs=x_in, outputs=x_out)

    es = EarlyStopping(monitor="loss", 
                   min_delta=0, 
                   patience=8, 
                   verbose=1, 
                   mode='min')

    checkpoint_path = "best_weights_tmp.h5"
    checkpoint = ModelCheckpoint(filepath=checkpoint_path,
                                monitor="loss",
                                save_best_only=True,
                                save_weights_only=True,
                                mode="min",verbose=0,)

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
    model.load_weights(checkpoint_path)

    x_emb_in, x_emb_out = base_model.in_out_tensors()

    emb_model = Model(inputs=x_emb_in, outputs=x_emb_out)
    node_gen = generator.flow(lncRNA_nodes, targets=None, shuffle=False)
    embedding = emb_model.predict(node_gen)

    if os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)

    return embedding

def run_deterministic_hinsage(LPPI_graph, layer_sizes, num_samples, epochs, seed=42):
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
    return run_deep_graph_infomax(model, generator, epochs=epochs)


layer_size_1=[16,32,64]
layer_size_2=[32,64,128,256]

samples_num_1=[5,10,15,20]
samples_num_2=[10,15,20,25]


# Tuning hyperparameters
proteins = pd.read_csv(f"../Annotate/human/transformed_protein_annotation.csv")
lncRNAs = pd.read_csv(f"../Annotate/human/valid_heart_annotation.csv")
LPPI = pd.read_csv(f'../Annotate/human/weighted_valid_inter.csv')

proteins.set_index(proteins['protein_ID'], inplace=True)
proteins = proteins.drop('protein_ID', axis=1)
lncRNAs.set_index(lncRNAs['lncRNA_ID'], inplace=True)
lncRNAs = lncRNAs.drop('lncRNA_ID', axis=1)

scaler_lncRNAs = StandardScaler()
lncRNAs_scaled = scaler_lncRNAs.fit_transform(lncRNAs)

# create StellarGraph Object
LPPI_graph = StellarGraph({"lncRNA": pd.DataFrame(lncRNAs_scaled, index=lncRNAs.index),
                        "protein": pd.DataFrame(proteins, index=proteins.index)},
                        LPPI)
print(LPPI_graph.info())

for i in layer_size_1:
    for j in layer_size_2:
        print(f"Train begin: {i} for first layer,{j} for second layer.")
        node_embedding = run_deterministic_hinsage(LPPI_graph, [i,j], [10,15], epochs=1000, seed=42)

        embedding_df = pd.DataFrame(node_embedding, index=LPPI_graph.nodes(node_type='lncRNA'))
        embedding_df.to_csv(f'human/layer_size/lncRNA_embeddings_heart_{i}_{j}.csv', header=None)

# for i in samples_num_1:
#     for j in samples_num_2:
#         print(f"Train begin: {i} for first layer,{j} for second layer.")
#         node_embedding = run_deterministic_hinsage(LPPI_graph, [64,256], [i,j], epochs=1000, seed=42)

#         embedding_df = pd.DataFrame(node_embedding, index=LPPI_graph.nodes(node_type='lncRNA'))
#         embedding_df.to_csv(f'mouse/samples_num/lncRNA_embeddings_heart_{i}_{j}.csv', header=None)