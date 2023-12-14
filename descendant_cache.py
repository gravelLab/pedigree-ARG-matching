from abc import abstractmethod


class DescendantCache:

    def record_descendants(self, left_parent_id: int, right_parent_id: int, child_id: int):
        self.record_child_descendants(left_parent_id, child_id)
        self.record_child_descendants(right_parent_id, child_id)

    @abstractmethod
    def record_proband(self, proband: int):
        pass

    @abstractmethod
    def get_vertex_descendants(self, vertex_id: int):
        pass

    @abstractmethod
    def get_dictionary(self):
        pass

    @abstractmethod
    def record_child_descendants(self, parent_id, child_id):
        pass
