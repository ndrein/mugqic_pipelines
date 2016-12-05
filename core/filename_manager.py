import traceback
from collections import namedtuple
from os.path import join


class File:
    # def __init__(self, tag, step, attributes):
    #     self.tag = tag
    #     self.step = step
    #     self.attributes = attributes
    def __init__(self, tag, readset_name=None):
        self.tag = tag
        self.readset_name = readset_name

    @property
    def filename(self):
        self.filename = join(self.tag.step, getattr(self, 'readset_name', ''), self.tag.name)
        return self.filename


class FilenameManager:
    def __init__(self, steps):
        self.files = {step: set() for step in steps}

    # def declare_output(self, tag, step, attributes=None):
    #     file = File(tag=tag, step=step, attributes=attributes)
    #     self.files[step].add(file)
    #     return file.filename
    def declare_output(self, tag, readset_name=None):
        file = File(tag=tag, readset_name=readset_name)
        self.files[tag.step].add(file)
        return file.filename

    def find_output(self, tag, readset_name=None):
        return [f for f in self.files[tag.step] if f.tag.name == tag.name and f.readset_name == readset_name][0].filename


FileTag = namedtuple('FileTag', ['step', 'name'])
