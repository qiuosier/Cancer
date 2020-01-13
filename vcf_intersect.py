import sys


class VCFFile:
    def __init__(self, file_path):
        self.file_path = file_path


if __name__ == "__main__":
    white_list = sys.argv[1]
    vcf_file = sys.argv[2]
    white_dict = {}
    with open(white_list) as f:
        for line in f:
            if not line.startswith("chr"):
                continue
            arr = line.split()
            chromosome = arr[0]
            position = arr[1]
            ref = arr[3]
            alt = arr[4]
            key = "%s%s%s%s" % (chromosome, position, ref, alt)
            white_dict[key] = line

    total = 0
    # print("In the white list:")
    # for key in white_dict.keys():
    #     count = len(white_dict[key])
    #     total += count
    #     print("%s has %s variants." % (key, count))
    print("%s variants in the white list." % total)

    variants = []
    with open(vcf_file) as f:
        for line in f:
            if not line.startswith("chr"):
                continue
            arr = line.split()
            chromosome = arr[0]
            position = arr[1]
            ref = arr[3]
            alt = arr[4]
            key = "%s%s%s%s" % (chromosome, position, ref, alt)
            if white_dict.get(key):
                variants.append(white_dict.get(key))

    print("%d variants found in VCF." % len(variants))
    if len(variants) < 100:
        for line in variants:
            print(line)

"""
gsutil cp gs://analysis_results/17927/annovar__dna/annovar__dna_Annovar.hg38_multianno.vcf .
gsutil cp gs://davelab_db/vcf/vcf_intersect.py .
gsutil cp gs://davelab_db/vcf/whitelist.vcf .
python vcf_intersect.py "whitelist.vcf" "annovar__dna_Annovar.hg38_multianno.vcf"
"""
