import traceback
from collections import namedtuple
from os.path import join


class File:
    def __init__(self, tag, step_name, readset_name=None):
        self.tag = tag
        self.step = step_name
        self.readset_name = readset_name

    @property
    def filename(self):
        return join(self.step, getattr(self, 'readset_name', ''), self.tag)


class FilenameManager:
    def __init__(self, steps):
        self.files = {step: set() for step in steps}

    def declare_output(self, tag, step, readset_name=None):
        file = File(tag, step.__name__, readset_name)
        self.files[step].add(file)
        return file.filename

    def find_output(self, tag, step, readset_name=None):
        try:
            return [f for f in self.files[step] if f.tag == tag and f.readset_name == readset_name][0].filename
        except IndexError:
            raise Exception('Unable to find output file:\n'
                            'Step: {step}\n'
                            'Tag: {tag}\n'.format(step=step.__name__, tag=tag) +
                            ('Readset: {readset}\n'.format(readset=readset_name) if readset_name else ''))

# FileTag = namedtuple('FileTag', ['step', 'name'])
