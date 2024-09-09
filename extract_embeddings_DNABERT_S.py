# extract_embeddings.py
import sys
import numpy as np
from transformers import AutoTokenizer, AutoModel
import torch
import pdb
from tqdm import tqdm

def extract_embeddings(sequences, model_name='zhihan1996/DNABERT-S'):
    tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True)

    # Ensure model is on CPU
    model.cpu()

    embeddings = []
    for seq in tqdm(sequences):
        inputs = tokenizer(seq, return_tensors='pt')['input_ids']
        inputs = inputs.to('cpu')

        # Assuming you are not receiving hidden states, just the final output
        outputs = model(inputs)  # No output_hidden_states=True if it doesn't work

        # Assuming outputs[0] contains the token-wise embeddings for the last layer
        last_hidden_states = outputs[0]  # [1, sequence_length, 768]

        # Average the embeddings across the sequence length to get a single vector
        embedding = torch.mean(last_hidden_states, dim=1).squeeze().detach().numpy()  # Mean pooling
        embeddings.append(embedding)

    return np.array(embeddings)

def main():
    if len(sys.argv) < 3:
        print("Usage: python extract_embeddings.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Read sequences from the input file
    sequences = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                sequences.append(line)

    # Extract embeddings
    embeddings = extract_embeddings(sequences, 'zhihan1996/DNABERT-S')

    # Save embeddings to a file
    np.save(output_file, embeddings)

if __name__ == "__main__":
    main()
