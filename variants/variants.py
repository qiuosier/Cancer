import csv
import copy
import sys
from commons.Aries.storage import StorageFile


class VariantCallList:
    def __init__(self, uri):
        self.uri = uri
        self.file_obj = None
        self.header_line = None
        self.read_headers()

    def variant_key(self, line):
        """Returns a unique identifier for a variant represented by
        a line or an object containing all information of a variant in the file.
        The format of the identifier will be:
        "Chromosome:Position:Ref>Alt"
        Chromosome should be a number, X or Y, without the "chr"
        This unique identifier should be consistent for all sub-classes.
        The same variant should have the same identifier regardless of the file format.

        Args:
            line: A text line or an object containing all information of a variant

        """
        raise NotImplementedError

    def build_index(self):
        return dict()

    def read_headers(self):
        raise NotImplementedError

    def write_headers(self, *args, **kwargs):
        raise NotImplementedError

    # Handles the open and close of the file
    def __enter__(self):
        self.file_obj = open(self.uri, 'r')
        return self.file_obj

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file_obj.close()


class VCFFile(VariantCallList):

    FILTER_COLUMN = 6

    def __init__(self, uri):
        self.meta = dict()
        self.headers = []
        self.__count = None
        super().__init__(uri)

    def variants(self):
        with self as f:
            for line in f:
                if line.startswith("#"):
                    continue
                key = self.variant_key(line)
                val = line
                yield key, val

    @staticmethod
    def read_header(header_line):
        arr = header_line.split("=", 1)
        key = arr[0]
        val = arr[1]
        if val.startswith("<") and val.startswith(">"):
            pairs = val[1:-1].split(",")
            val = dict()
            for pair in pairs:
                arr = pair.split("=", 1)
                k = arr[0]
                v = arr[1]
                val[k] = v
        return key, val

    def read_headers(self):
        with self as f:
            for line in f:
                if not line.startswith("#"):
                    break
                if line.startswith("##"):
                    key, val = self.read_header(line)
                    val_list = self.meta.get(key, [])
                    val_list.append(val)
                    self.meta[key] = val_list
                    continue
                if line.startswith("#"):
                    self.header_line = line

    def write_headers(self, to_file, meta=None):
        """

        Args:
            to_file:
            meta:

        Returns:

        """
        if meta is None:
            meta = self.meta
        for key, val_list in meta.items():
            for val in val_list:
                if isinstance(val, dict):
                    value = "<" + ",".join(["%s=%s" % (k, v) for k, v in val.items()]) + ">"
                else:
                    value = val
                to_file.write("%s=%s" % (key, value))

    def variant_key(self, line):
        arr = line.split()
        if len(arr) < 5:
            return None
        chromosome = arr[0].replace("chr", "")
        position = arr[1]
        ref = arr[3]
        alt = arr[4]
        return "%s:%s:%s>%s" % (chromosome, position, ref, alt)

    def build_index(self):
        whitelist_dict = {}
        with self as f:
            for line in f:
                if line.startswith("#"):
                    continue
                key = self.variant_key(line)
                whitelist_dict[key] = line.strip()
        return whitelist_dict

    def apply_filter(self, output_vcf_path, filter_id, filter_description, filter_func, passed_only=False):
        meta = copy.deepcopy(self.meta)
        filter_list = meta.get("FILTER", [])
        filter_list.append({
            "ID": filter_id,
            "Description": filter_description
        })
        with open(output_vcf_path, 'w') as output_vcf:
            self.write_headers(output_vcf, meta)
            counter = 0
            for key, variant in self.variants():
                counter += 1
                columns = variant.split("\t")
                passed = filter_func(key, variant)
                if not passed:
                    if passed_only:
                        continue
                    # Filter not passed
                    if columns[self.FILTER_COLUMN].strip() in ["PASS", ".", ""]:
                        columns[self.FILTER_COLUMN] = filter_id
                    else:
                        columns[self.FILTER_COLUMN] = columns[self.FILTER_COLUMN] + "," + filter_id
                output_vcf.write("\t".join(columns))
        self.__count = counter
        print("%d total variants in VCF." % counter)
        return VCFFile(output_vcf_path)

    def count(self):
        if self.__count is None:
            self.__count = sum(1 for _ in self.variants())
        return self.__count


class CSVFile(VariantCallList):
    header_keys = ["chr", "start", "end", "ref", "alt"]

    def __init__(self, uri, delimiter=","):
        self.delimiter = delimiter
        self.columns = dict()
        self.header_line = None
        super().__init__(uri)

    def variant_key(self, row):
        columns = self.columns
        return "%s:%s:%s>%s" % (
            str(row[columns.get("chr")]).replace("chr", ""),
            row[columns.get("start")],
            row[columns.get("ref")],
            row[columns.get("alt")],
        )

    def read_headers(self):
        with self as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=self.delimiter)
            for row in csv_reader:
                self.header_line = self.delimiter.join(row).strip()
                self.columns = self.column_index(row, self.header_keys)
                break
        for key in self.header_keys:
            if self.columns.get(key) is None:
                raise AttributeError("Column %s not found in %s." % (key, self.columns))

    def write_headers(self, to_file):
        to_file.write(self.header_line + "\n")

    @staticmethod
    def column_index(header_row, keys):
        columns = {}
        data = [str(v).lower() for v in header_row]
        for key in keys:
            columns[key] = data.index(key.lower())
        # print(columns)
        return columns

    def build_index(self):
        whitelist_dict = {
            "header": self.header_line
        }
        with self as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=self.delimiter)
            for i, row in enumerate(csv_reader, start=1):
                if i == 1:
                    continue
                if row[self.columns.get("start")] != row[self.columns.get("end")]:
                    print("Variant has different start and end")
                key = self.variant_key(row)
                whitelist_dict[key] = self.delimiter.join(row).strip()
        return whitelist_dict


class WhitelistFilter:
    def __init__(self, variant_calls):
        """

        Args:
            variant_calls (VariantCallList):
        """
        self.variant_calls = variant_calls
        self.index = variant_calls.build_index()
        self.in_whitelist = set()

    def filter_variant(self, key, variant):
        if key in self.index:
            self.in_whitelist.add(key)
            return True
        return False

    def print_passed(self):
        print(self.variant_calls.header_line)
        counter = 0
        for key in self.in_whitelist:
            print(self.index[key])
            counter += 1
            if counter >= 100:
                print("... and %s more..." % (len(self.in_whitelist) - counter))

    def save_passed(self, to_file_path):
        with open(to_file_path, "w") as f:
            self.variant_calls.write_headers(f)
            for key in self.in_whitelist:
                f.write(self.index[key] + "\n")


def main(whitelist_path, vcf_path, output_vcf_path=None, output_csv_path=None):
    print("Filtering VCF with %s" % whitelist_path)
    # whitelist_dict = build_whitelist_from_vcf(whitelist_path)
    if whitelist_path.endswith(".vcf"):
        whitelist = VCFFile(whitelist_path)
    elif whitelist_path.endswith(".csv"):
        whitelist = CSVFile(whitelist_path)
    elif whitelist_path.endswith(".tsv"):
        whitelist = CSVFile(whitelist_path, "\t")
    else:
        raise TypeError("Whitelist file type is not supported.")
    vcf = VCFFile(vcf_path)
    whitelist_filter = WhitelistFilter(whitelist)
    print("%s variants in the white list." % len(whitelist_filter.index.keys()))
    description = "In Whitelist: %s" % whitelist_path
    vcf.apply_filter(output_vcf_path, "Whitelist", description, whitelist_filter.filter_variant, passed_only=True)
    print("%d whitelist variants found in VCF." % len(whitelist_filter.in_whitelist))
    whitelist_filter.print_passed()
    if output_csv_path:
        whitelist_filter.save_passed(output_csv_path)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Example command:")
        print("python vcf_intersect.py whitelist.csv raw.vcf output.vcf output.tsv")
        print("output.csv is optional.")
    else:
        main(*sys.argv[1:])
