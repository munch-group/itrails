"""Fill in missing branch lengths of a Newick tree, assuming the tree is ultrametric.

Usage:
    python fill_ultrametric.py "(A,(B,(E:3,D:3):5):5);"
    echo "(A,(B,(E:3,D:3):5):5);" | python fill_ultrametric.py

or from Python:
    from fill_ultrametric import fill_ultrametric_newick
    fill_ultrametric_newick("(A,(B,(E:3,D:3):5):5);")   # -> '(A:13,(B:8,(E:3,D:3):5):5);'

Every node has a height (distance down to the tips; leaves are at 0) and every
non-root node has a branch length, tied together by
    height(parent) = height(node) + length(node).
Known lengths are propagated through these equations until nothing changes.
If the given lengths contradict each other, or leave some length undetermined,
a ValueError is raised naming the offending nodes.
"""

import math
import re
import sys


class Node:
    __slots__ = ("name", "length", "children", "parent", "height")

    def __init__(self, name=None, length=None):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None
        self.height = None

    def __repr__(self):
        return self.name or "(" + ",".join(repr(c) for c in self.children) + ")"


_TOKEN = re.compile(r"\s*(\(|\)|,|;|:|'[^']*'|[^\s()':,;\[]+)")
_COMMENT = re.compile(r"\[[^\]]*\]")


def parse_newick(text):
    """Minimal Newick parser: nested parentheses, optional labels (quoted or not),
    optional ':length', [comments] are stripped."""
    text = _COMMENT.sub("", text.strip())
    tokens = _TOKEN.findall(text)
    pos = 0

    def peek():
        return tokens[pos] if pos < len(tokens) else None

    def take():
        nonlocal pos
        if pos >= len(tokens):
            raise ValueError("unexpected end of Newick string")
        tok = tokens[pos]
        pos += 1
        return tok

    def node():
        n = Node()
        if peek() == "(":
            take()
            while True:
                child = node()
                child.parent = n
                n.children.append(child)
                if peek() == ",":
                    take()
                else:
                    break
            if take() != ")":
                raise ValueError("expected ')' in Newick string")
        if peek() not in (None, "(", ")", ",", ":", ";"):
            n.name = take().strip("'")
        if peek() == ":":
            take()
            n.length = float(take())
        return n

    root = node()
    if peek() != ";":
        raise ValueError("expected ';' at end of Newick string")
    return root


def _fmt(x):
    s = "%.12g" % x
    return s if s != "-0" else "0"


def to_newick(node):
    def rec(n):
        s = "(" + ",".join(rec(c) for c in n.children) + ")" if n.children else ""
        if n.name is not None:
            s += n.name
        if n.length is not None:
            s += ":" + _fmt(n.length)
        return s
    return rec(node) + ";"


def walk(node):
    yield node
    for c in node.children:
        yield from walk(c)


def fill_ultrametric(root, rel_tol=1e-9):
    """Fill in missing branch lengths in place, assuming an ultrametric tree.
    Returns the root. Raises ValueError on contradictions or underdetermined lengths."""
    nodes = list(walk(root))
    for n in nodes:
        if not n.children:
            n.height = 0.0

    def same(a, b):
        return math.isclose(a, b, rel_tol=rel_tol, abs_tol=1e-12)

    changed = True
    while changed:
        changed = False
        for n in nodes:
            p = n.parent
            if p is None:
                continue
            # child height + child branch -> parent height
            if n.height is not None and n.length is not None:
                h = n.height + n.length
                if p.height is None:
                    p.height = h
                    changed = True
                elif not same(p.height, h):
                    raise ValueError(
                        f"tree is not ultrametric: node {p!r} has height {_fmt(p.height)} "
                        f"from one child but {_fmt(h)} via child {n!r}"
                    )
            # parent height - child branch -> child height
            if n.height is None and p.height is not None and n.length is not None:
                n.height = p.height - n.length
                changed = True
            # parent height - child height -> child branch
            if n.length is None and p.height is not None and n.height is not None:
                n.length = p.height - n.height
                if n.length < 0 and not same(n.length, 0.0):
                    raise ValueError(
                        f"given lengths force a negative branch ({_fmt(n.length)}) above node {n!r}"
                    )
                changed = True

    missing = [n for n in nodes if n.parent is not None and n.length is None]
    if missing:
        raise ValueError(
            "branch lengths are underdetermined for: "
            + ", ".join(repr(n) for n in missing)
            + " (no path of known lengths ties the subtree to the tips)"
        )
    return root


def fill_ultrametric_newick(newick):
    return to_newick(fill_ultrametric(parse_newick(newick)))


if __name__ == "__main__":
    src = sys.argv[1:] or [line for line in sys.stdin.read().split("\n") if line.strip()]
    for s in src:
        print(fill_ultrametric_newick(s))
