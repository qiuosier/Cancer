from .Aries.storage import StorageFile


class VCF:
    def __init__(self, uri):
        self.uri = uri


class InMemoryVCF(VCF):
    def __init__(self, uri):
        super().__init__(uri)
        self.content = StorageFile.init(uri).read()
        if isinstance(self.content, bytes):
            self.content = self.content.decode()
        self.content = self.content.split("\n")
        self.headers = []
        self.variants = []
        for line in self.content:
            if not line:
                continue
            if line.startswith("#"):
                self.headers.append(line)
            else:
                self.variants.append(line)

    @property
    def count(self):
        return len(self.variants)
