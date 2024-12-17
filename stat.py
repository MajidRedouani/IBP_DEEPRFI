import pandas as pd
from scipy.stats import mannwhitneyu

def compare_unmatched_tools(file_tool1, file_tool2, output_file):
    """
    Compares two tools' performance using Mann-Whitney U test for unpaired data
    and determines which tool performed better.

    Args:
        file_tool1 (str): Filename for Tool 1's performance scores.
        file_tool2 (str): Filename for Tool 2's performance scores.
        output_file (str): Name of the file to save the results.
    """
    # Load the performance data
    data_tool1 = pd.read_csv(file_tool1, sep="\t")
    data_tool2 = pd.read_csv(file_tool2, sep="\t")

    # Strip leading/trailing whitespace from column names
    data_tool1.columns = data_tool1.columns.str.strip()
    data_tool2.columns = data_tool2.columns.str.strip()

    # Initialize results
    results = []

    # List of sub-ontologies to compare
    subontologies = ["MF_Score", "CC_Score", "BP_Score", "Total_Score"]

    for ontology in subontologies:
        if ontology in data_tool1.columns and ontology in data_tool2.columns:
            # Extract scores for the sub-ontology
            scores_tool1 = pd.to_numeric(data_tool1[ontology], errors="coerce").dropna()
            scores_tool2 = pd.to_numeric(data_tool2[ontology], errors="coerce").dropna()

            # Perform Mann-Whitney U test
            if not scores_tool1.empty and not scores_tool2.empty:
                stat, p_value = mannwhitneyu(scores_tool1, scores_tool2, alternative="two-sided")

                # Calculate medians
                median_tool1 = scores_tool1.median()
                median_tool2 = scores_tool2.median()

                
                # Determine which tool performed better
                if p_value > 0.05:  # If the result is not statistically significant
                    better_tool = "Tie"
                elif median_tool1 > median_tool2:
                    better_tool = file_tool1
                else:
                    better_tool = file_tool2

                results.append({
                    "Ontology": ontology,
                    "MannWhitney_Stat": stat,
                    "P_Value": p_value,
                    "Tool_With_Better_Performance": better_tool
                })
            else:
                print(f"No valid data for {ontology}.")
                results.append({
                    "Ontology": ontology,
                    "MannWhitney_Stat": None,
                    "P_Value": None,
                    "Tool_With_Better_Performance": "No Data"
                })
        else:
            print(f"{ontology} not found in both files.")
            results.append({
                "Ontology": ontology,
                "MannWhitney_Stat": None,
                "P_Value": None,
                "Tool_With_Better_Performance": "Not Found"
            })

    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)

    # Save results to a file
    results_df.to_csv(output_file, sep="\t", index=False)

    print(f"Comparison results saved to {output_file}")


# Example usage
file_tool1 = "simrel_scores_output_blast2go.tsv"
file_tool2 = "simrel_scores_output_mDeepFRI.tsv"
file_tool3 = "simrel_scores_output_OGDeepFRI.tsv"
output_file = "comparison_blast2go_vs_mDeepFRI_unmatched.tsv"
output_file2 = "comparison_blast2go_vs_OGDeepFRI_unmatched.tsv"
output_file3 = "comparison_mDeepFRI_vs_OGDeepFRI_unmatched.tsv"


compare_unmatched_tools(file_tool1, file_tool2, output_file)
compare_unmatched_tools(file_tool1, file_tool3, output_file2)
compare_unmatched_tools(file_tool2, file_tool3, output_file3)
