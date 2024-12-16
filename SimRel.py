import pandas as pd
from goatools import obo_parser
from goatools.semantic import TermCounts, resnik_sim, get_info_content
from collections import defaultdict
import numpy as np
from itertools import islice

# Load GO ontology
godag = obo_parser.GODag("go-basic.obo")
print(f"Total GO terms loaded into godag: {len(godag)}")

# Load annotations (all peptides)
annotations_df = pd.read_csv("all_go_terms.tsv", sep='\t', names=["ProteinID", "Known_GO_Terms"], skiprows=1)
annotations = defaultdict(list)

# Parse annotations to create a dictionary {ProteinID: [GO_Terms]}
for _, row in annotations_df.iterrows():
    protein_id = row["ProteinID"]
    go_terms = row["Known_GO_Terms"].split(",")
    annotations[protein_id].extend(go_terms)



# Create TermCounts object for annotation-based frequency calculation (all data)
term_counts = TermCounts(godag, annotations)




for go_term in term_counts.gocnts.keys():
    freq = term_counts.get_term_freq(go_term)
    print(f"GO Term: {go_term}, Frequency: {freq}")


# Load predicted GO terms (696 peptides)
df_predicted = pd.read_csv("go_terms_DeepFRI.tsv", sep='\t', names=["ProteinID", "Predicted_GO_Terms"], skiprows=1)



# Merge known and predicted datasets (only the 696 peptides)
df = pd.merge(annotations_df, df_predicted, on="ProteinID", how="inner")



# Parse GO terms in both datasets
df["Known_GO_Terms"] = df["Known_GO_Terms"].apply(lambda x: [term.strip() for term in x.split(",")])
df["Predicted_GO_Terms"] = df["Predicted_GO_Terms"].apply(lambda x: [term.strip() for term in x.split(",")])


# Prepare to store unknown predicted GO terms
unknown_terms = []

def resolve_obsolete_term(term, terms):
    if term in godag and godag[term].is_obsolete:
        replacements = godag[term].replaced_by or godag[term].consider
        if replacements:
            replacement = replacements[0]
            if replacement not in terms:
                return replacement
        return None
    return term

def filter_terms_by_namespace(terms, namespace):
 
    filtered_terms = []
    for term in terms:
        resolved_term = resolve_obsolete_term(term, terms)
        if resolved_term and resolved_term in godag and godag[resolved_term].namespace == namespace:
            filtered_terms.append(resolved_term)
    return filtered_terms


def calculate_bma_similarity(known_terms, predicted_terms, term_counts, protein_id):
    """
    Calculate Best Match Average (BMA) similarity using GOATOOLS.
    """
   

    if not known_terms:
        # No known terms; log predictions as unknown
        if predicted_terms:
            for term in predicted_terms:
                if term in godag:
                    unknown_terms.append([
                        term, godag[term].name, godag[term].namespace, protein_id
                    ])
        return "NA"

    if not predicted_terms:
        # No predictions; similarity score is 0
        return 0

    sim_matrix = []
    for known in known_terms:
        row = []
        for predicted in predicted_terms:
            if known in godag and predicted in godag:
                lca_ic = resnik_sim(known, predicted, godag, term_counts)
                max_ic = max(get_info_content(known, term_counts), get_info_content(predicted, term_counts))
                sim_score = lca_ic / max_ic if max_ic > 0 else 0
                row.append(sim_score)
            else:
                row.append(0)
        sim_matrix.append(row)
            
    if not sim_matrix:
        return 0

    sim_matrix = np.array(sim_matrix)
    row_max = np.max(sim_matrix, axis=1).sum() if len(known_terms) > 0 else 0
    col_max = np.max(sim_matrix, axis=0).sum() if len(predicted_terms) > 0 else 0
    bma_score = (row_max + col_max) / (len(known_terms) + len(predicted_terms))
   
    return bma_score


# Calculate similarities for each protein
results = []
for _, row in df.iterrows():
    protein_id = row["ProteinID"]
    known_terms = row["Known_GO_Terms"]
    predicted_terms = row["Predicted_GO_Terms"]

    # Separate GO terms by namespace
    known_bp = filter_terms_by_namespace(known_terms, "biological_process")
    known_cc = filter_terms_by_namespace(known_terms, "cellular_component")
    known_mf = filter_terms_by_namespace(known_terms, "molecular_function")
    
    predicted_bp = filter_terms_by_namespace(predicted_terms, "biological_process")
    predicted_cc = filter_terms_by_namespace(predicted_terms, "cellular_component")
    predicted_mf = filter_terms_by_namespace(predicted_terms, "molecular_function")
    
    
    # Calculate BMA similarities for each ontology
    bp_score = calculate_bma_similarity(known_bp, predicted_bp, term_counts, protein_id)
    cc_score = calculate_bma_similarity(known_cc, predicted_cc, term_counts, protein_id)
    mf_score = calculate_bma_similarity(known_mf, predicted_mf, term_counts, protein_id)

    # Handle missing ontology scores
    bp_score = bp_score if bp_score is not None else "NA"
    cc_score = cc_score if cc_score is not None else "NA"
    mf_score = mf_score if mf_score is not None else "NA"

    # Calculate total score as the average of the three ontologies (if available)
    ontology_scores = [s for s in [bp_score, cc_score, mf_score] if s != "NA"]
    total_score = np.mean(ontology_scores) if ontology_scores else "NA"

    results.append([protein_id, mf_score, cc_score, bp_score, total_score])
    
  


# Save results to output file
result_df = pd.DataFrame(results, columns=["ProteinID", "MF_Score", "CC_Score", "BP_Score", "Total_Score"])
result_df.to_csv("simrel_scores_output.tsv", sep="\t", index=False)
print("Similarity scores calculated and saved to simrel_scores_output.tsv")

# Save unknown predicted terms
unknown_terms_df = pd.DataFrame(unknown_terms, columns=[
    "Predicted_GO_Term_ID", "Predicted_GO_Term_Name", "Subontology", "ProteinID"
])
unknown_terms_df.to_csv("unknown_predicted_terms.tsv", sep="\t", index=False)
print("Unknown predicted terms saved to unknown_predicted_terms.tsv")
