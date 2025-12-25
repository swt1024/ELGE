#!/bin/bash

# Usage: ./run_model.sh <species>

SPECIES=$1

# Define default values for tissues and model parameters based on species
if [ "$SPECIES" == "human" ]; then
    # For human, set tissues and model parameters
    TISSUES="heart lung stomach"
    LAYER_SIZES="64 256"
    SAMPLES_NUM="20 25"
elif [ "$SPECIES" == "mouse" ]; then
    # For mouse, set tissues and model parameters
    TISSUES="heart lung brain"
    LAYER_SIZES="32 64"
    SAMPLES_NUM="10 15"
else
    echo "Invalid species: $SPECIES. Choose 'human' or 'mouse'."
    exit 1
fi

# Loop over each tissue in the tissue list
for TISSUE in $TISSUES; do
    echo "Processing tissue: $TISSUE"

    # Run the Python script for each tissue, processing only one LPPI network at a time
    python train.py --layer_sizes $LAYER_SIZES \
      --samples_num $SAMPLES_NUM \
      --lncRNA_nodes_file "../annotate/$SPECIES/valid_${TISSUE}_annotation.csv" \
	  --protein_nodes_file "../annotate/$SPECIES/transformed_protein_annotation.csv" \
      --lppi_file "../annotate/$SPECIES/weighted_valid_inter_${TISSUE}.csv" \
      --embedding_save_path "./$SPECIES/unweighted/lncRNA_embeddings_${TISSUE}.csv"
      
    echo "Finished processing tissue: $TISSUE"
done
