import sys


def variant_key(vcf_line):
    arr = vcf_line.split()
    if len(arr) < 5:
        return None
    chromosome = arr[0].replace("chr", "")
    position = arr[1]
    ref = arr[3]
    alt = arr[4]
    return "%s:%s:%s>%s" % (chromosome, position, ref, alt)


def build_whitelist_dict(vcf_whitelist_path):
    whitelist_dict = {}
    with open(vcf_whitelist_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            key = variant_key(line)
            whitelist_dict[key] = line.strip()
    return whitelist_dict


def intersect_whitelist(whitelist_dict, vcf_path):
    variants = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            key = variant_key(line)
            if key in whitelist_dict:
                variants.append(whitelist_dict.get(key))

    print("%d whitelist variants found in VCF." % len(variants))
    counter = 0
    for line in variants:
        print(line)
        counter += 1
        if counter >= 100:
            print("... and %s more..." % (len(variants) - counter))


def main(vcf_whitelist_path, vcf_path):
    whitelist_dict = build_whitelist_dict(vcf_whitelist_path)
    print("%s variants in the white list." % len(whitelist_dict.keys()))
    intersect_whitelist(whitelist_dict, vcf_path)


if __name__ == "__main__":
    white_list = sys.argv[1]
    vcf_file = sys.argv[2]
    main(white_list, vcf_file)


"""
gsutil cp gs://analysis_results/17927/annovar__dna/annovar__dna_Annovar.hg38_multianno.vcf .
gsutil cp gs://davelab_db/vcf/vcf_intersect.py .
gsutil cp gs://davelab_db/vcf/whitelist.vcf .
python vcf_intersect.py "whitelist.vcf" "annovar__dna_Annovar.hg38_multianno.vcf"
"""
