class Node:
    def __init__(self, val, payload=None):
        self.left = None
        self.right = None
        self.value = val
        self.payload = payload

class Tree:
    def __init__(self, need_payload=False):
        self.root = None
        self.need_payload = need_payload

    def delete_tree(self):
        self.root = None

    def is_empty(self):
        return self.root == None

    def print_tree(self):
        if not self.root == None:
            self._print_tree(self.root)

    def _print_tree(self, node):
        if not node == None:
            self._print_tree(node.left)
            print(str(node.value)+' ')
            self._print_tree(node.right)

    def find(self, val):
        if not self.root == None:
            return self._find(val, self.root)
        else:
            return None

    def _find(self, val, node):
        if val == node.value:
            return node
        elif(val < node.value and not node.left == None):
            self._find(val, node.left)
        elif(val > node.value and not node.right == None):
            self._find(val, node.right)

    def add(self, new_node):
        if self.root == None:
            self.root = new_node
        else:
            self._add(new_node, self.root)

    def _add(self, new_node, node):
        if new_node.value < node.value:
            if not node.left == None:
                self._add(new_node, node.left)
            else:
                node.left = new_node
        else:
            if not node.right == None:
                self._add(new_node, node.right)
            else:
                node.right = new_node

    def in_order_traversal(self):
        if self.root == None:
            return None
        result = []
        self._in_order_traversal(result, self.root)
        return result

    def _in_order_traversal(self, result, node):
        if not node.left == None:
            self._in_order_traversal(result, node.left)
        if self.need_payload == True:
            result.append((node.value, node.payload))
        else:
            result.append(node.value)
        if not node.right == None:
            self._in_order_traversal(result, node.right)
