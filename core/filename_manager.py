from os.path import join


class File:
    def __init__(self, tag, step, attributes):
        self.tag = tag
        self.step = step
        self.attributes = attributes

    def hasattrs(self, attributes):
        return attributes <= self.attributes

    @property
    def filename(self):
        self.filename = join(self.step.__name__, '{tag}.{hash}'.format(tag=self.tag, hash=hash(frozenset(self.attributes.items()))))
        return self.filename


class FilenameManager:
    def __init__(self, steps):
        self.files = {step: set() for step in steps}

    def __repr__(self):
        # for step in self['step']:
        raise NotImplemented

    def declare_output(self, tag, step, attributes=None):
        file = File(tag=tag, step=step, attributes=attributes)
        self.files[step].add(file)
        return file.filename


    def find_output(self, tag, step, attributes={}):
        return [f for f in self.files[step] if f.tag == tag and f.hasattrs(attributes)][0].filename
