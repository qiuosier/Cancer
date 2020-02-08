import logging
from .Aries.storage import StorageFile
from .variants.files import CSVVariants
logger = logging.getLogger(__name__)


class VCF:
    def __init__(self, uri):
        self.uri = uri


class Variant:
    def __init__(self, vcf_line, annotation=None):
        self.line = vcf_line
        self.columns = self.line.split("\t")
        self.annotation = annotation
        if not self.annotation:
            self.annotation = dict()

    @property
    def chromosome(self):
        return str(self.columns[0]).upper().strip("CHR")

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
            if len(arr) < 2:
                continue
            key = arr[0]
            val = arr[1]
            data[key] = val
        return data

    def __str__(self):
        return "%s:g.%s%s>%s" % (self.chromosome, self.position, self.ref, self.alt)


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
        return "%s:g.%s%s>%s" % (chromosome, position, ref, alt)

    @property
    def count(self):
        return len(self.variants)

    @staticmethod
    def load_annotation(uri):
        logger.debug("Loading annotations from %s" % uri)
        new_index = dict()
        csv = CSVVariants(uri)
        index = csv.build_index()
        headers = csv.headers
        for k, v in index.items():
            arr = v.split(csv.delimiter)
            data = dict()
            for i in range(len(arr)):
                if i < len(headers):
                    data[headers[i]] = arr[i]
            new_index[k] = data
        return new_index
