from os.path import join


class FilenameManager:
    def __init__(self, steps):
        self.files = {step: dict() for step in steps}

    def __repr__(self):
        # for step in self['step']:
        raise NotImplemented

    @staticmethod
    def make_filename(step, attributes):
        # TODO: use extension to make filename
        return join(step, attributes.get('readset_name', ''), hash(attributes))

    @staticmethod
    def is_contained_in(dict1, dict2):
        """
        Determine if dict1 is contained in dict2
        """
        return dict1.items() in dict2.items()

    def set(self, step, attributes=None, filename=''):
        self.files[step] = filename if filename is not None else self.make_filename(step, attributes)

    def get(self, step, attributes=None):
        # TODO: handle no match
        matching_files = [file for file in self.files[step] if self.is_contained_in(attributes, self.files[step])]
        assert(len(matching_files) == 1)
        return matching_files[0].filename

