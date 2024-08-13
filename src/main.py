import argparse
import os
import math
from Bio import SeqIO, Align
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, BeforePosition, AfterPosition
from tqdm import tqdm

def parse_arguments():
    parser = argparse.ArgumentParser(description="Annotate query sequence based on reference gene exons and introns")
    parser.add_argument("--reference_folder", required=False, default="None", help="Folder with reference genes in GenBank format")
    parser.add_argument("--input_gbk", required=False, default="None", help="Reference gene in GenBank format")
    parser.add_argument("--query", required=True, help="Query sequence in FASTA format")
    parser.add_argument("--gene", required=False, default=None, help="Gene name")
    parser.add_argument("--product", required=False, default=None, help="CDS product")
    parser.add_argument("--output", required=True, help="Output file path")
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

def map_features_to_query(reference_record, query_record, alignment, gene_name, mismatch_tolerance=0.):
    annotated_features = []
    ref_to_query_map = alignment.aligned

    numero = 1
    for feature in reference_record.features:
        
        if feature.type == "exon" or (feature.type == "CDS" and not isinstance(feature.location, CompoundLocation)):
            feature_length = feature.location.end - feature.location.start
            max_mismatches = int(feature_length * mismatch_tolerance)
            mapped = False
            
            for (ref_start, ref_end), (query_start, query_end) in zip(ref_to_query_map[0], ref_to_query_map[1]):
                overlap_start = max(feature.location.start, ref_start)
                overlap_end = min(feature.location.end, ref_end)
                overlap_length = overlap_end - overlap_start

                if overlap_length >= feature_length - max_mismatches:
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


def create_features(reference_record, query_record, gene_name, alignment, product):
    annotated_features = []

    if gene_name is None:
        for feature in reference_record.features:
            if feature.type == "gene":
                gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
                break

    if product is None:
        for feature in reference_record.features:
            if feature.type == "CDS":
                product = feature.qualifiers.get("product", ["unknown"])[0]
                break

    # Map features to query and determine fuzzy positions
    mapped_features = map_features_to_query(reference_record, query_record, alignment, gene_name)
    mapped_exon_locations = [feature.location for feature in mapped_features if feature.type == "exon" or feature.type == "CDS"]

    # Determine if first and last exons are missing
    first_exon_missing = (len(mapped_exon_locations) > 0 and mapped_exon_locations[0].start > 0)
    last_exon_missing = (len(mapped_exon_locations) > 0 and mapped_exon_locations[-1].end < len(query_record.seq))

    # Handle the case where no exons are mapped
    if not mapped_exon_locations:
        gene_start = BeforePosition(0) if isinstance(reference_record.features[0].location.start, BeforePosition) else 0
        gene_end = AfterPosition(len(query_record.seq)) if isinstance(reference_record.features[-1].location.end, AfterPosition) else len(query_record.seq)
        gene_feature = SeqFeature(FeatureLocation(gene_start, gene_end), type="gene", qualifiers={"gene": gene_name})
        annotated_features.append(gene_feature)
        
        cds_feature = SeqFeature(FeatureLocation(gene_start, gene_end), type="CDS", qualifiers={"gene": gene_name, "product": product})
        annotated_features.append(cds_feature)
        
    else:
        # Add gene feature
        reference_gene = next((f for f in reference_record.features if f.type == "gene"), None)
        gene_has_partial_start = isinstance(reference_gene.location.start, BeforePosition) if reference_gene else False
        gene_has_partial_end = isinstance(reference_gene.location.end, AfterPosition) if reference_gene else False

        # Determine fuzzy positions for gene based on mapping
        gene_start = BeforePosition(mapped_exon_locations[0].start) if gene_has_partial_start or mapped_exon_locations[0].start > 0 else mapped_exon_locations[0].start
        gene_end = AfterPosition(mapped_exon_locations[-1].end) if gene_has_partial_end or mapped_exon_locations[-1].end < len(query_record.seq) else mapped_exon_locations[-1].end
        gene_feature = SeqFeature(FeatureLocation(gene_start, gene_end), type="gene", qualifiers={"gene": gene_name})
        annotated_features.append(gene_feature)

        # Handle fuzzy positions for the CDS based on alignment and mapping
        if first_exon_missing:
            mapped_exon_locations[0] = FeatureLocation(BeforePosition(mapped_exon_locations[0].start), mapped_exon_locations[0].end, strand=mapped_exon_locations[0].strand)
        if last_exon_missing:
            mapped_exon_locations[-1] = FeatureLocation(mapped_exon_locations[-1].start, AfterPosition(mapped_exon_locations[-1].end), strand=mapped_exon_locations[-1].strand)

        compound_location = CompoundLocation(mapped_exon_locations) if len(mapped_exon_locations) > 1 else mapped_exon_locations[0]
        qualifiers = {"gene": gene_name, "product": product}
        cds_feature = SeqFeature(location=compound_location, type="CDS", qualifiers=qualifiers)
        annotated_features.append(cds_feature)

    annotated_features.extend(mapped_features)

    return annotated_features



def write_tbl(annotations, output_path):
    with open(output_path, "w") as tbl_file:
        tbl_file.write(">Feature " + annotations['seqid'] + "\n")

        # Write gene features first
        for feature in annotations['features']:
            if feature.type == "gene":
                start = str(feature.location.start + 1)
                end = str(feature.location.end)
                tbl_file.write(f"{start}\t{end}\t{feature.type}\n")
                for key, value in feature.qualifiers.items():
                    tbl_file.write(f"\t\t\t{key}\t{value}\n")

        # Handle CDS features
        written_cds = False
        for feature in annotations['features']:
            if feature.type == "CDS" and not written_cds:
                if isinstance(feature.location, CompoundLocation):
                    # Handle compound locations
                    for i, part in enumerate(feature.location.parts):
                        start = str(part.start + 1)
                        end = str(part.end)
                        if i == 0:
                            tbl_file.write(f"{start}\t{end}\tCDS\n")
                        else:
                            tbl_file.write(f"{start}\t{end}\n")
                else:
                    # Handle non-compound CDS
                    start = str(feature.location.start + 1)
                    end = str(feature.location.end)
                    tbl_file.write(f"{start}\t{end}\tCDS\n")

                # Write the CDS qualifiers
                tbl_file.write(f"\t\t\tgene\t{feature.qualifiers['gene']}\n")
                if 'product' in feature.qualifiers:
                    tbl_file.write(f"\t\t\tproduct\t{feature.qualifiers['product']}\n")
                if 'translation' in feature.qualifiers:
                    tbl_file.write(f"\t\t\ttranslation\t{feature.qualifiers['translation'][0]}\n")

                written_cds = True

        # Write exon features if present
        for feature in annotations['features']:
            if feature.type == "exon":
                start = str(feature.location.start + 1)
                end = str(feature.location.end)
                tbl_file.write(f"{start}\t{end}\t{feature.type}\n")
                for key, value in feature.qualifiers.items():
                    tbl_file.write(f"\t\t\t{key}\t{value}\n")



def read_genbank(input_file):
    records = list(SeqIO.parse(input_file, "genbank"))
    for record in records:
        new_features = []
        for feature in record.features:
            if feature.type == "CDS" and isinstance(feature.location, CompoundLocation):
                parts = feature.location.parts
                for index, part in enumerate(parts):
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

def read_fasta(fastafile):
    total_size = os.path.getsize(fastafile)
    with open(fastafile) as f, tqdm(total=total_size, desc="Reading FASTA file", unit="B", unit_scale=True, unit_divisor=1024) as pbar:
        total_records = 0
        for line in f:
            pbar.update(len(line))
            if line.startswith(">"):
                total_records += 1
    
    with tqdm(total=total_records, desc="Parsing FASTA file", unit=" Records") as pbar:
        records = []
        for record in SeqIO.parse(fastafile, "fasta"):
            records.append(record)
            pbar.update(1)
    
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
            reference_records.extend(read_genbank(reference_path))
    elif args.input_gbk != "None":
        reference_records.extend(read_genbank(args.input_gbk))
    
    query_record = SeqIO.read(args.query, "fasta")
    query_record_parsed = read_fasta(args.query)
    if len(query_record_parsed) != 1:
        raise ValueError("Query should contain only one sequence.")
    else:
        query_record_id = query_record_parsed[0].id

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

        annotated_features = create_features(best_reference, query_record, args.gene, best_alignment, args.product)
        query_record.features = annotated_features
        query_record.id = query_record_id
        query_record.name = "Query"
        query_record.description = f"Partial gene, exons and introns annotated"
        query_record.annotations["molecule_type"] = "genomic DNA"

        write_tbl({"seqid": query_record.id, "features": query_record.features}, args.output)

    except ValueError as e:
        print(e)
        return

if __name__ == "__main__":
    main()
