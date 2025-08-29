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


# ========== Customizable Feature Selection Section ==========

# Define your feature columns for lncRNA and protein here (edit as needed)
lncRNA_feature_cols = []     # e.g., ['expression', 'conservation', ...]
protein_feature_cols = []    # e.g., ['degree', 'essentiality', ...]

# Define ablation configs: name → {'lncRNA': cols, 'protein': cols}
feature_config = {
    'lppi_only': {
		'lncRNA':None,
		'protein': None,
	},
    '+seq': {
        'lncRNA': ['CpG_count','CpG_islands','GC_content (%)','Length'],     
		'protein': None,     
    },
    '+seq': {
        'lncRNA': ['CpG_count','CpG_islands','GC_content (%)','Length',
                   'phyloP_mean_score','phyloP_max_score','phastCons_mean_score','phastCons_max_score'], 
		'protein': None,               
    },
    '+epi': {
        'lncRNA': ['CpG_count','CpG_islands','GC_content (%)','Length',
                   'phyloP_mean_score','phyloP_max_score','phastCons_mean_score','phastCons_max_score',
                   'DHS_peak_counts','DHS_max_signalValue','H3K4me3_peak_counts','H3K4me3_max_signalValue',
                   'H3K9me3_peak_counts','H3K9me3_max_signalValue','H3K27ac_peak_counts','H3K27ac_max_signalValue',
                   'H3K36me3_peak_counts','H3K36me3_max_signalValue','H3K27me3_peak_counts','H3K27me3_max_signalValue'], 
		'protein': None,         
    }
}

# ========== Main Loop for Ablation ==========

species = 'human'
tissue = 'heart'

if species == 'human':
	samples_num = [20,15]
	layer_sizes = [64,256]	
else:
	samples_num = [10,15]
	layer_sizes = [64,64]


# Load original annotation tables
proteins = pd.read_csv(f"../Annotate/{species}/transformed_protein_annotation.csv")
lncRNAs = pd.read_csv(f"../Annotate/{species}/valid_{tissue}_annotation.csv")
LPPI = pd.read_csv(f'../Annotate/{species}/weighted_valid_inter.csv')

proteins.set_index('protein_ID', inplace=True)
lncRNAs.set_index('lncRNA_ID', inplace=True)

for ablation_name, config in feature_config.items():
	# ---- Prepare lncRNA node features ----
	if not config['lncRNA']:
		# Use zeros as features (no features)
		lnc_feat = pd.DataFrame(1, index=lncRNAs.index, columns=["f0"])
	else:
		lnc_feat = lncRNAs[config['lncRNA']].copy()
		scaler_lncRNAs = StandardScaler()
		lnc_feat = pd.DataFrame(
			scaler_lncRNAs.fit_transform(lnc_feat),
			index=lnc_feat.index,
			columns=lnc_feat.columns
		)

	# ---- Prepare protein node features ----
	if not config['protein']:
		prot_feat = pd.DataFrame(1, index=proteins.index, columns=["f0"])
	else:
		prot_feat = proteins[config['protein']].copy()
		scaler_proteins = StandardScaler()
		prot_feat = pd.DataFrame(
			scaler_proteins.fit_transform(prot_feat),
			index=prot_feat.index,
			columns=prot_feat.columns
		)

	# ---- Build graph ----
	LPPI_graph = StellarGraph(
		{"lncRNA": lnc_feat, "protein": prot_feat}, LPPI
	)
	print(f"[{tissue}][{ablation_name}] graph info:")
	print(LPPI_graph.info())

	# ---- Train & Save embeddings ----
	model_save_path = f"./{species}/{tissue}_{ablation_name}_weights.h5"

	node_embedding = run_deterministic_hinsage(
		LPPI_graph,
		layer_sizes=layer_sizes,
		num_samples=samples_num,
		epochs=1000,
		model_save_path=model_save_path,
		seed=42
	)

	# Save node embeddings
	emb_save_path = f"./{species}/lncRNA_embeddings_{tissue}_{ablation_name}.csv"
	embedding_df = pd.DataFrame(node_embedding, index=LPPI_graph.nodes(node_type='lncRNA'))
	embedding_df.to_csv(emb_save_path, header=None)
	print(f"Saved embeddings to: {emb_save_path}")

