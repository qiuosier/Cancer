import csv
import copy
import sys
from commons.Aries.table import TableCSVFile


class VCFVariant:
    def __init__(self, chrom, pos, ref, alt, info=None, annotation=None, raw_data=None):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt

        self.raw_data = raw_data
        self.rs_id = None

        if annotation is None:
            self.annotation = dict()
        self.annotation = annotation

        if info is None:
            self.info = dict()
        self.info = info

    @property
    def key(self):
        return "%s:g.%s%s>%s" % (self.chrom, self.pos, self.ref, self.alt)

    def __str__(self):
        return self.key


class VariantsFile:
    def __init__(self, uri):
        self.uri = uri
        self.file_obj = None
        # header_line should NOT contain the ending line break
        self.header_line = None
        self.read_headers()

    def variant_key(self, line):
        """Returns a unique identifier for a variant represented by
        a line or an object containing all information of a variant in the file.
        The format of the identifier must be:
        "[Chromosome]:p.[Position][Ref]>[Alt]"
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

    @property
    def lines(self):
        with open(self.uri, 'r') as f:
            for line in f:
                yield line

    # Handles the open and close of the file
    def __enter__(self):
        self.file_obj = open(self.uri, 'r')
        return self.file_obj

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file_obj.close()


class VCFVariants(VariantsFile):

    FILTER_COLUMN = 6

    def __init__(self, uri):
        self.meta = dict()
        self.headers = []

        # Caches the number of total variants
        self.__count = None
        super().__init__(uri)

    def variants(self):
        for line in self.lines:
            if line.startswith("#"):
                continue
            key = self.variant_key(line)
            val = line
            yield key, val

    @staticmethod
    def parse_meta(meta_line):
        arr = meta_line.split("=", 1)
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
        for line in self.lines:
            if not line.startswith("#"):
                break
            if line.startswith("##"):
                key, val = self.parse_meta(line)
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
        # Write meta data
        if meta is None:
            meta = self.meta
        for key, val_list in meta.items():
            for val in val_list:
                if isinstance(val, dict):
                    value = "<" + ",".join(["%s=%s" % (k, v) for k, v in val.items()]) + ">"
                else:
                    value = val
                to_file.write("%s=%s" % (key, value))
        # Write header line
        to_file.write(self.header_line)

    def variant_key(self, line):
        arr = line.split()
        if len(arr) < 5:
            return None
        chromosome = arr[0].replace("chr", "")
        position = arr[1]
        ref = arr[3]
        alt = arr[4]
        return "%s:g.%s%s>%s" % (chromosome, position, ref, alt)

    def build_index(self):
        whitelist_dict = {}
        for line in self.lines:
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
        return VCFVariants(output_vcf_path)

    def count(self):
        if self.__count is None:
            self.__count = sum(1 for _ in self.variants())
        return self.__count


class CSVVariants(VariantsFile):
    header_keys = ["chr", "start", "end", "ref", "alt"]

    def __init__(self, uri, **kwargs):
        self.table = TableCSVFile(uri, **kwargs)
        self.delimiter = self.table.kwargs.get("delimiter")
        if not self.delimiter:
            self.delimiter = self.table.kwargs.get("dialect").delimiter
        self.columns = dict()
        super().__init__(uri)
        self.header_line = self.delimiter.join(["\"%s\"" % h for h in self.table.headers])

    def variant_key(self, row):
        columns = self.columns
        return "%s:g.%s%s>%s" % (
            str(row[columns.get("chr")]).replace("chr", ""),
            row[columns.get("start")],
            row[columns.get("ref")],
            row[columns.get("alt")],
        )

    def read_headers(self):
        headers = [str(h).lower() for h in self.table.headers]
        self.columns = {headers[i]: i for i in range(len(headers))}
        for key in self.header_keys:
            if key not in self.columns:
                raise AttributeError("Column %s not found in %s." % (key, self.table.headers))

    def write_headers(self, to_file):
        to_file.write(self.header_line + "\n")

    def build_index(self):
        whitelist_dict = {
            "header": self.header_line
        }
        for i, row in enumerate(self.table.rows):
            if i <= self.table.header_row:
                continue
            if row[self.columns.get("start")] != row[self.columns.get("end")]:
                print("Variant has different start and end")
            key = self.variant_key(row)
            whitelist_dict[key] = self.delimiter.join(["\"%s\"" % r for r in row]).strip()
        return whitelist_dict


class WhitelistFilter:
    def __init__(self, variant_calls):
        """

        Args:
            variant_calls (VariantsFile):
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
