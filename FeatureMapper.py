import argparse
import os
import math
from Bio import SeqIO, Align
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
def parse_arguments():
    parser = argparse.ArgumentParser(description="Annotate query sequence based on reference gene exons and introns")
    parser.add_argument("--reference_folder", required=False,default="None", help="folder with reference genes in GenBank format")
    parser.add_argument("--input_gbk", required=False, default= "None", help="Reference gene in GenBank format")
    parser.add_argument("--query", required=True, help="Query sequence in FASTA format")
    parser.add_argument("--organism", required=True, help="Organism name")
    parser.add_argument("--isolate", required=True, help="Isolate information")
    parser.add_argument("--output", required=True, help="Output file path")
    parser.add_argument("--gene", required=True, help="Gene name")
    parser.add_argument("--format", default="tbl", choices=["tbl", "gbk"], help="Output format (default: tbl)")
    return parser.parse_args()

def align_sequences(reference, query):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -5
    aligner.extend_gap_score = 0
    aligner.match_score = 1
    aligner.mismatch_score = -1
    alignments = aligner.align(reference, query)
    if not alignments:
        raise ValueError("No alignments were found. Check the sequences or adjust alignment settings.")
    return alignments[0]

def map_features_to_query(reference_record, query_record, alignment, gene_name, mismatch_tolerance=0.1):
    annotated_features = []
    ref_to_query_map = alignment.aligned

    numero = 1
    for feature in reference_record.features:
        
        if feature.type == "exon":
            feature_length = feature.location.end - feature.location.start
            max_mismatches = int(feature_length * mismatch_tolerance)
            mapped = False
            
            for (ref_start, ref_end), (query_start, query_end) in zip(ref_to_query_map[0], ref_to_query_map[1]):
                overlap_start = max(feature.location.start, ref_start)
                overlap_end = min(feature.location.end, ref_end)
                overlap_length = overlap_end - overlap_start

                if overlap_length >= feature_length - max_mismatches:
                    # Calculate adjusted start and end considering possible mismatches
                    query_feature_start = int(query_start + (overlap_start - ref_start))
                    query_feature_end = int(query_start + (overlap_end - ref_start))
                    new_feature = SeqFeature(FeatureLocation(query_feature_start, query_feature_end, strand=feature.location.strand),
                                             type=feature.type,
                                             qualifiers={"gene": gene_name, "number": numero})
                    annotated_features.append(new_feature)
                    mapped = True
                    print(f"Feature {feature.type} at {feature.location} mapped to {new_feature.location}")
                    break

            if not mapped:
                print(f"Warning: Feature {feature.type} at {feature.location} could not be mapped. Mismatches may exceed tolerance.")
            else:
                numero += 1
            
    return annotated_features

def write_tbl(annotations, output_path):
    with open(output_path, "w") as tbl_file:
        tbl_file.write(">Feature " + annotations['seqid'] + "\n")

        # First, write the gene feature
        for feature in annotations['features']:
            if feature.type == "gene":
                start = '<' + str(feature.location.start + 1) if 'partial_start' in feature.qualifiers else str(feature.location.start + 1)
                end = '>' + str(feature.location.end) if 'partial_end' in feature.qualifiers else str(feature.location.end)
                tbl_file.write(f"{start}\t{end}\t{feature.type}\n")
                for key, value in feature.qualifiers.items():
                    if key not in ["partial_start", "partial_end"]:
                        tbl_file.write(f"\t\t\t{key}\t{value}\n")

        # Handle compound CDS features separately
        for feature in annotations['features']:
            if feature.type == "CDS" and isinstance(feature.location, CompoundLocation):
                # Write CDS header line first
                first_part = feature.location.parts[0]
                start = '<' + str(first_part.start + 1) if 'partial_start' in feature.qualifiers else str(first_part.start + 1)
                end = '>' + str(first_part.end) if 'partial_end' in feature.qualifiers else str(first_part.end)
                tbl_file.write(f"{start}\t{end}\t{feature.type}\n")

                # Write remaining parts without the type
                for part in feature.location.parts[1:]:
                    start = '<' + str(part.start + 1) if 'partial_start' in feature.qualifiers else str(part.start + 1)
                    end = '>' + str(part.end) if 'partial_end' in feature.qualifiers else str(part.end)
                    tbl_file.write(f"{start}\t{end}\n")

                # Write qualifiers for the CDS feature once
                tbl_file.write(f"\t\t\tgene\t{feature.qualifiers['gene']}\n")
                tbl_file.write(f"\t\t\tproduct\t{feature.qualifiers['product']}\n")
                if 'translation' in feature.qualifiers:
                    tbl_file.write(f"\t\t\ttranslation\t{feature.qualifiers['translation'][0]}\n")

        # Write all non-CDS exon features next
        for feature in annotations['features']:
            if feature.type == "exon" and not isinstance(feature.location, CompoundLocation):
                start = '<' + str(feature.location.start + 1) if 'partial_start' in feature.qualifiers else str(feature.location.start + 1)
                end = '>' + str(feature.location.end) if 'partial_end' in feature.qualifiers else str(feature.location.end)
                tbl_file.write(f"{start}\t{end}\t{feature.type}\n")
                for key, value in feature.qualifiers.items():
                    if key not in ["partial_start", "partial_end"]:
                        tbl_file.write(f"\t\t\t{key}\t{value}\n")

        # Write any remaining non-source, non-gene, non-CDS, non-exon features
        for feature in annotations['features']:
            if feature.type not in ["source", "gene", "CDS", "exon"]:
                start = '<' + str(feature.location.start + 1) if 'partial_start' in feature.qualifiers else str(feature.location.start + 1)
                end = '>' + str(feature.location.end) if 'partial_end' in feature.qualifiers else str(feature.location.end)
                tbl_file.write(f"{start}\t{end}\t{feature.type}\n")
                for key, value in feature.qualifiers.items():
                    if key not in ["partial_start", "partial_end"]:
                        tbl_file.write(f"\t\t\t{key}\t{value}\n")


def create_features(reference_record, query_record, organism, isolate, gene_name, alignment):
    annotated_features = [SeqFeature(FeatureLocation(0, len(query_record.seq)), type="source", qualifiers={
        "organism": organism, "mol_type": "genomic DNA", "isolate": isolate})]

    gene_feature = SeqFeature(FeatureLocation(0, len(query_record.seq)), type="gene", qualifiers={"gene": gene_name})
    annotated_features.append(gene_feature)

    mapped_features = map_features_to_query(reference_record, query_record, alignment, gene_name)

    exon_locations = [feature.location for feature in mapped_features if feature.type == "exon"]
    if exon_locations:
        compound_location = CompoundLocation(exon_locations)
        cds_feature = SeqFeature(location=compound_location, type="CDS", qualifiers={"gene": gene_name, "product": gene_name})
        annotated_features.append(cds_feature)

    annotated_features.extend(mapped_features)

    return annotated_features

def write_genbank(query_record, output_path):
    with open(output_path, "w") as output_handle:
        SeqIO.write(query_record, output_handle, "genbank")

def read_genbank(input_file):
    # Load the GenBank file
    records = list(SeqIO.parse(input_file, "genbank"))
    
    # Process each record
    for record in records:
        new_features = []
        for feature in record.features:
            if feature.type == "CDS" and isinstance(feature.location, CompoundLocation):
                parts = feature.location.parts
                for index, part in enumerate(parts):
                    # Create an exon feature
                    exon = SeqFeature(
                        location=FeatureLocation(start=part.start, end=part.end, strand=part.strand),
                        type="exon",
                        qualifiers={
                            "number": str(index + 1),
                            "parent_cds": feature.qualifiers.get("locus_tag", ["unknown"])[0]
                        }
                    )
                    new_features.append(exon)
            new_features.append(feature)
        record.features = new_features
    
    return records

def main():
    args = parse_arguments()
    if args.reference_folder != "None" and args.input_gbk != "None":
        raise Exception("Please provide either a reference folder or a reference gene, not both.")
    elif args.reference_folder == "None" and args.input_gbk == "None":
        raise Exception("Please provide a reference folder or a reference gene.")
    
    reference_records = []
    if args.reference_folder != "None":
        for reference_file in os.listdir(args.reference_folder):
            reference_path = os.path.join(args.reference_folder, reference_file)
            reference_records.extend(read_genbank(reference_path))  # Extend the list with records from each file
    elif args.input_gbk != "None":
        reference_records.extend(read_genbank(args.input_gbk))  # Extend the list with records from the single file
    
    query_record = SeqIO.read(args.query, "fasta")
    if len(query_record) !=1:
        raise ValueError("Query should contain only one sequence.")
    else:
        query_record_id = query_record[0].id

    try:
        max_score = -math.inf
        best_reference = None
        best_alignment = None
        for reference_record in reference_records:
            temp_alignment = align_sequences(reference_record.seq, query_record.seq)
            if temp_alignment.score > max_score:
                max_score = temp_alignment.score
                best_alignment = temp_alignment
                best_reference = reference_record
        if best_reference is None or best_alignment is None:
            raise ValueError("No suitable alignment found.")

        print(f"Best alignment found with {best_reference.id} {best_reference.description} (score: {max_score})")
        print(best_alignment)

        annotated_features = create_features(best_reference, query_record, args.organism, args.isolate, args.gene, best_alignment)
        query_record.features = annotated_features
        query_record.id = query_record_id
        query_record.name = "Query"
        query_record.description = f"{args.organism} partial gene, exons and introns annotated"
        query_record.annotations["molecule_type"] = "genomic DNA"

        if args.format == "tbl":
            write_tbl({"seqid": query_record.id, "features": query_record.features}, args.output)
        else:
            write_genbank(query_record, args.output)

    except ValueError as e:
        print(e)
        return

if __name__ == "__main__":
    main()
