import logging
from .Aries.storage import StorageFile
from .variants.variants import CSVFile
logger = logging.getLogger(__name__)


class VCF:
    def __init__(self, uri):
        self.uri = uri


class Variant:
    def __init__(self, vcf_line, annotation=None):
        self.line = vcf_line
        self.columns = self.line.split("\t")
        if annotation is None:
            self.annotation = dict()
        self.annotation = annotation

    @property
    def chromosome(self):
        return self.columns[0]

    @property
    def position(self):
        return self.columns[1]

    @property
    def rs_id(self):
        return self.columns[2]

    @property
    def ref(self):
        return self.columns[3]

    @property
    def alt(self):
        return self.columns[4]

    @property
    def info(self):
        data = dict()
        pairs = self.columns[7].split(";")
        for pair in pairs:
            arr = pair.split("=", 1)
            key = arr[0]
            val = arr[1]
            data[key] = val
        return data


class InMemoryVCF(VCF):
    def __init__(self, uri, annotation_uri):
        super().__init__(uri)
        self.content = StorageFile.init(uri).read()
        if isinstance(self.content, bytes):
            self.content = self.content.decode()
        self.content = self.content.split("\n")
        self.headers = []
        self.variants = []
        self.annotations = self.load_annotation(annotation_uri)
        for line in self.content:
            if not line:
                continue
            if line.startswith("#"):
                self.headers.append(line)
            else:
                key = self.variant_key(line)
                self.variants.append(Variant(line, self.annotations.get(key)))

    def group_by_gene(self):
        groups = dict()
        for v in self.variants:
            if not v.annotation:
                continue
            gene = v.annotation.get('Gene')
            v_list = groups.get(gene, [])
            v_list.append(v)
            groups[gene] = v_list
        logger.debug(groups)
        return groups

    def variant_key(self, line):
        arr = line.split()
        if len(arr) < 5:
            return None
        chromosome = arr[0].replace("chr", "")
        position = arr[1]
        ref = arr[3]
        alt = arr[4]
        return "%s:%s:%s>%s" % (chromosome, position, ref, alt)

    @property
    def count(self):
        return len(self.variants)

    @staticmethod
    def load_annotation(uri):
        new_index = dict()
        csv = CSVFile(uri, "\t")
        index = csv.build_index()
        headers = csv.header_line.split("\t")
        for k, v in index.items():
            arr = v.split("\t")
            data = dict()
            for i in range(len(arr)):
                if i < len(headers):
                    data[headers[i]] = arr[i]
            new_index[k] = data
        return new_index
