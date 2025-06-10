import gradio as gr
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from pprint import pprint


def find_protospacer_by_grna_id(grna_id: str):
    """
    Search Benchling database by gRNA ID and retrieve the protospacer sequence by the grna_id.

    Args:
        grna_id (str): The gRNA ID to search for

    Returns:
        str: The protospacer sequence
    """
    # fake data to mimic data retrieval from Benchling database
    grna_data = {
        "TAR206": "GGTGCTAGCCTTGCGTTCCG",
        "TAR207": "CTGGCCGAAGCACCCGAGCA",
        "TAR208": "TGCGGAAACCTTCTAGGGTG",
        "TAR209": "GCGGAAACCTTCTAGGGTGT",
        "TAR210": "TCAAGCACCCACACCCTAGA",
    }

    return grna_data.get(grna_id, "No protospacer found")


def initiate_aligner(aligner_type: str):
    """
    Initiate a Biopython pairwise aligner with custom penalty score matrix.
    """
    if aligner_type == "protospacer":
        alphabet = "ACGTN"
        aln_array = np.zeros([5, 5])
        np.fill_diagonal(aln_array, 10)
        aln_array[np.where(aln_array == 0)] = -10
        aln_array[:, -1] = 0
        aln_array[-1, :] = 0

        custom_aln_matrix = substitution_matrices.Array(
            alphabet=alphabet,
            data=aln_array
        )

        aligner = PairwiseAligner()

        aligner.substitution_matrix = custom_aln_matrix
        aligner.mode = "global"
        aligner.open_gap_score = -15
        aligner.extend_gap_score = -20
        # aligner.target_end_gap_score = -10
        # aligner.query_end_gap_score = -10
        # aligner.target_left_open_gap_score = -10
        # aligner.query_left_open_gap_score = -10

    return aligner


def compute_mismatch_between_protospacer_and_off_target(
    off_target_sequence: str,
    protospacer_sequence: str
) -> dict:
    """
    Align query off-target sequence to protospacer sequence. Return the mismatch annotation.

    Args:
        off_target_sequence (str): The off-target sequence to align
        protospacer_sequence (str): The protospacer sequence to align to

    Returns:
        dict: Annotation of the mismatch between the off-target sequence and protospacer sequence, and the aligned sequences.
    """
    aligner = initiate_aligner("protospacer")
    alignments = aligner.align(protospacer_sequence, off_target_sequence)

    if len(alignments) > 0:
        best_alignment = alignments[0]

        aligned_protospacer_seq = best_alignment[0]
        aligned_off_target_seq = best_alignment[1]

        mismatch_annotation = {
            "aligned_protospacer_seq": aligned_protospacer_seq,
            "aligned_off_target_seq": aligned_off_target_seq,
            "number_of_mismatch": 0,
            "number_of_dna_bulge": 0,
            "number_of_rna_bulge": 0,
        }

        for dna_base, rna_base in zip(aligned_protospacer_seq, aligned_off_target_seq):
            if dna_base != rna_base:
                if dna_base == "-":
                    mismatch_annotation["number_of_rna_bulge"] += 1
                elif rna_base == "-":
                    mismatch_annotation["number_of_dna_bulge"] += 1
                else:
                    mismatch_annotation["number_of_mismatch"] += 1

        return mismatch_annotation


mismatch_demo = gr.Interface(
    fn=compute_mismatch_between_protospacer_and_off_target,
    inputs=[gr.Textbox("off-target sequence"),
            gr.Textbox("protospacer sequence")],
    outputs="text",
    title="Mismatch between protospacer and off-target",
    description="Enter protospacer and off-target sequence to compute mismatch between them."
)


grna_demo = gr.Interface(
    fn=find_protospacer_by_grna_id,
    inputs=gr.Textbox(label="Enter gRNA ID", placeholder="e.g., TAR206"),
    outputs="text",
    description="Get the protospacer sequence by gRNA ID",
)

# Create tabbed interface
demo = gr.TabbedInterface(
    [grna_demo, mismatch_demo], ["gRNA", "Mismatch"], title="gRNA Information Hub"
)

if __name__ == "__main__":
    demo.launch(mcp_server=True)

    # ontarget = "GGTGCTAGCCTTGCGTTCCG"
    # offtarget = "GGTGCTAGCGTTGGTTCCG"
    # mismatch_annotation = compute_mismatch_between_protospacer_and_off_target(
    #     protospacer_sequence=ontarget,
    #     off_target_sequence=offtarget
    # )

    # pprint(mismatch_annotation)
