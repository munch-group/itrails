"""Build a rooted binary tree from rooted three-leaf trees ("trios").

Two builders are provided.

build_aho(trios)
    The BUILD algorithm of Aho, Sagiv, Szymanski and Ullman (1981,
    https://doi.org/10.1137/0210030).  It uses all trios at once, so the result
    does not depend on their order.  It raises ValueError if the trios
    contradict each other, or if they leave some relationship undetermined (the
    error names the groups whose relationship is missing).  When it succeeds,
    the tree it returns is the only one that displays every trio.  Use this
    when each trio is a true relationship.

build_tree(trios) / add_trio(tree, trio)
    Grafts trios onto a growing tree one at a time.  Every node of a trio is
    matched to a node of the tree:

    * a leaf matches the tree leaf with the same key, if there is one;
    * an inner node matches the lowest common ancestor (LCA) of the matches of
      its two children, if both children have a match.

    A trio node without a match is new.  It is inserted directly above the match
    of its matched child, with a copy of the unmatched child as its sibling, i.e.
    as close to the leaves already in the tree as the trio allows.  A trio pins a
    new leaf to a single edge only in the tight cases (a new outgroup when the
    shared pair's LCA is the root, or a new sister to a shared leaf whose partner
    is its direct neighbour).  Otherwise the placement is a choice, so the result
    can depend on the order in which trios are added, and a later trio may even
    conflict with a placement that an earlier one left open.  build_tree's
    most_overlap_first option reduces this but does not remove it in general.

    Grafting never changes the relationships among leaves that are already in
    the tree, so earlier trios stay satisfied.  A trio whose leaves are all
    present is only checked; a contradiction raises ValueError.

node_map(tree, trios)
    Maps every node of a built tree to the trio nodes it is built from.  Both
    sides are numbered by pre-order position, the order in which from_tuples
    creates nodes; preorder() lists a tree's nodes in that order.
"""


class TreeNode:

    def __init__(self, key):
        self.key = key
        self.left = None
        self.right = None


def from_tuples(t):

    if type(t) is not tuple:
        return TreeNode(t)
    tree = TreeNode(t)
    tree.left = from_tuples(t[0])
    tree.right = from_tuples(t[1])
    return tree


def to_tuples(node):
    """Inverse of from_tuples."""
    if node.left is None:
        return node.key
    return (to_tuples(node.left), to_tuples(node.right))


def copy_tree(node):
    new = TreeNode(node.key)
    if node.left is not None:
        new.left = copy_tree(node.left)
        new.right = copy_tree(node.right)
    return new


def leaves(node):
    """Leaf nodes below node, left to right."""
    if node.left is None:
        return [node]
    return leaves(node.left) + leaves(node.right)


def leaf_keys(node):
    return {leaf.key for leaf in leaves(node)}


def preorder(node):
    """Nodes of the subtree at node in pre-order: the order from_tuples creates them."""
    if node.left is None:
        return [node]
    return [node] + preorder(node.left) + preorder(node.right)


def refresh_keys(node):
    """Reset every inner key to the nested-tuple form of its subtree, as from_tuples does."""
    if node.left is None:
        return node.key
    node.key = (refresh_keys(node.left), refresh_keys(node.right))
    return node.key


def triplet(trio):
    """The trio's leaf keys as (a, b, c), with a and b the closest pair."""
    keys = [leaf.key for leaf in leaves(trio)]
    if len(keys) != 3 or len(set(keys)) != 3:
        raise ValueError(f'not a trio with three distinct leaves: {to_tuples(trio)}')
    if trio.left.left is not None:      # closest pair on the left
        a, b, c = keys
    else:
        c, a, b = keys
    return a, b, c


# ---------------------------------------------------------------------------
# BUILD (Aho et al., 1981): all trios at once
# ---------------------------------------------------------------------------

def build_aho(trios):
    """The tree that displays every trio, built from all trios at once.

    Every trio's closest pair must sit in the same clade under the root, so the
    connected components of the graph joining those pairs are the root's clades;
    the same applies recursively within each clade.  If the trios contradict
    each other some clade cannot be split at all, and if they leave a
    relationship undetermined some clade splits into three or more parts; both
    raise ValueError.  Otherwise every clade splits in exactly two, which any
    tree displaying the trios must reproduce, so the result is the unique such
    tree.  Inner keys follow from_tuples's convention.
    """
    trios = list(trios)
    triplets = [triplet(t) for t in trios]
    order = list(dict.fromkeys(leaf.key for t in trios for leaf in leaves(t)))
    return from_tuples(_aho(order, triplets))


def _aho(order, triplets):
    """Nested-tuple tree on the leaves in `order` that displays every triplet."""
    if len(order) == 1:
        return order[0]
    # Aho graph: join the closest pair of every triplet.  Each connected
    # component must lie within one clade under the root, so the components
    # are those clades.
    rep = {leaf: leaf for leaf in order}

    def find(leaf):
        while rep[leaf] != leaf:
            leaf = rep[leaf]
        return leaf

    for a, b, _ in triplets:
        rep[find(a)] = find(b)
    clades = {}
    for leaf in order:
        clades.setdefault(find(leaf), []).append(leaf)
    clades = list(clades.values())
    if len(clades) == 1:
        raise ValueError('these trios contradict each other: '
                         + ', '.join(str(((a, b), c)) for a, b, c in triplets))
    subtrees = [_aho(clade, [t for t in triplets if set(t) <= set(clade)])
                for clade in clades]
    if len(clades) > 2:
        raise ValueError('the trios do not resolve the relationships among '
                         + ', '.join(map(str, subtrees)))
    return tuple(subtrees)


# ---------------------------------------------------------------------------
# Iterative grafting: one trio at a time
# ---------------------------------------------------------------------------

def add_trio(tree, trio):
    """Graft trio onto tree in place and return the root.

    The root is a new node when the trio's outgroup is added above the old root,
    so always keep working with the returned node.
    """
    trio_keys = triplet(trio)

    parent = {}    # tree node -> its parent (None for the root)
    leaf_of = {}   # leaf key -> tree leaf

    def register(node, par):
        parent[node] = par
        if node.left is None:
            leaf_of[node.key] = node
        else:
            register(node.left, node)
            register(node.right, node)

    register(tree, None)

    if not any(key in leaf_of for key in trio_keys):
        raise ValueError(f'trio {to_tuples(trio)} shares no leaf with the tree')

    root = tree

    def ancestors(node):
        chain = []
        while node is not None:
            chain.append(node)
            node = parent[node]
        return chain

    def lca(u, v):
        above_u = set(ancestors(u))
        return next(w for w in ancestors(v) if w in above_u)

    def graft(node):
        """Return the tree node matching this trio node, inserting what is missing.

        None means nothing below this trio node is in the tree yet.
        """
        nonlocal root
        if node.left is None:
            return leaf_of.get(node.key)
        lm, rm = graft(node.left), graft(node.right)
        if lm is not None and rm is not None:
            w = lca(lm, rm)
            if w is lm or w is rm:
                raise ValueError(f'trio {to_tuples(trio)} conflicts with tree {to_tuples(root)}')
            return w
        if lm is None and rm is None:
            return None
        # Exactly one child is in the tree, so this trio node is new: insert it
        # directly above the matched child, with a copy of the other child as sibling.
        new = TreeNode(None)
        if rm is None:                       # left child is in the tree
            old, fresh = lm, copy_tree(node.right)
            new.left, new.right = old, fresh
        else:                                # right child is in the tree
            old, fresh = rm, copy_tree(node.left)
            new.left, new.right = fresh, old
        par = parent[old]
        if par is None:
            root = new
        elif par.left is old:
            par.left = new
        else:
            par.right = new
        parent[new] = par
        parent[old] = new
        register(fresh, new)
        return new

    graft(trio)
    refresh_keys(root)
    return root


def build_tree(trios, most_overlap_first=False):
    """Build a tree from trios by adding them one at a time.

    By default trios are added in the given order.  With most_overlap_first the
    trio sharing the most leaves with the current tree is added next, which
    keeps arbitrary placements (trios sharing a single leaf) to a minimum.
    """
    todo = list(trios)
    tree = copy_tree(todo.pop(0))
    while todo:
        i = 0
        if most_overlap_first:
            present = leaf_keys(tree)
            i = max(range(len(todo)), key=lambda i: len(leaf_keys(todo[i]) & present))
        tree = add_trio(tree, todo.pop(i))
    return tree


def displays(tree, trio):
    """True if tree induces the trio's topology on the trio's three leaves."""
    a, b, c = triplet(trio)
    pair = {a, b}

    def has_clade_separating(node):
        below = leaf_keys(node)
        if pair <= below and c not in below:
            return True
        return node.left is not None and (
            has_clade_separating(node.left) or has_clade_separating(node.right))

    return {a, b, c} <= leaf_keys(tree) and has_clade_separating(tree)


# ---------------------------------------------------------------------------
# Provenance: which trio nodes a tree node is built from
# ---------------------------------------------------------------------------

def node_map(tree, trios):
    """Map each node of tree to the trio nodes it is built from.

    Returns {i: [(t, j), ...]}: tree node number i is built from node number j
    of trio number t.  Nodes are numbered by position in a pre-order traversal,
    the order from_tuples creates them; preorder() lists them in that order.

    A trio leaf corresponds to the tree leaf with the same key, and a trio
    inner node to the lowest common ancestor in tree of the leaves below it.
    Every trio must be displayed by tree (ValueError otherwise); that is what
    makes the two inner nodes of a trio land on distinct tree nodes.
    """
    nodes = preorder(tree)
    index = {node: i for i, node in enumerate(nodes)}
    clade = {node: leaf_keys(node) for node in nodes}

    def lca(keys):
        """Deepest tree node with all of keys below it."""
        node = tree
        while node.left is not None:
            child = node.left if keys <= clade[node.left] else node.right
            if not keys <= clade[child]:
                return node
            node = child
        return node

    mapping = {i: [] for i in range(len(nodes))}
    for t, trio in enumerate(trios):
        if not displays(tree, trio):
            raise ValueError(f'trio {to_tuples(trio)} is not displayed by the tree')
        for j, node in enumerate(preorder(trio)):
            mapping[index[lca(leaf_keys(node))]].append((t, j))
    return mapping


if __name__ == '__main__':
    from itertools import permutations

    trio1 = from_tuples((('H', 'C'), 'G'))
    trio2 = from_tuples((('H', 'G'), 'O'))
    trio3 = from_tuples((('H', 'N'), 'C'))
    trio4 = from_tuples(('H', ('O', 'X')))
    trios = [trio1, trio2, trio3, trio4]
    expected = (((('H', 'N'), 'C'), 'G'), ('O', 'X'))

    def canonical(t):
        """Nested tuples with a fixed child order, to compare trees up to mirroring."""
        if type(t) is not tuple:
            return t
        return tuple(sorted((canonical(t[0]), canonical(t[1])), key=str))

    # ----- BUILD -----
    tree = build_aho(trios)
    print('build_aho:', to_tuples(tree))
    assert to_tuples(tree) == expected
    assert all(displays(tree, t) for t in trios)
    assert tree.key == expected                       # inner keys follow from_tuples's convention
    assert all(canonical(to_tuples(build_aho(order))) == canonical(expected)
               for order in permutations(trios))

    # Brute force: the expected tree is the only rooted binary tree on these
    # six species that displays all four trios.
    def rooted_binary_trees(keys):
        if len(keys) == 1:
            yield keys[0]
            return
        for t in rooted_binary_trees(keys[1:]):
            yield from insert_everywhere(keys[0], t)

    def insert_everywhere(leaf, t):
        yield (leaf, t)
        if type(t) is tuple:
            yield from ((l, t[1]) for l in insert_everywhere(leaf, t[0]))
            yield from ((t[0], r) for r in insert_everywhere(leaf, t[1]))

    candidates = list(rooted_binary_trees(['H', 'C', 'G', 'O', 'N', 'X']))
    consistent = [t for t in candidates if all(displays(from_tuples(t), trio) for trio in trios)]
    assert len(candidates) == 945 and len(consistent) == 1
    assert canonical(consistent[0]) == canonical(expected)
    print(f'brute force: {len(consistent)} of {len(candidates)} rooted binary trees displays all four trios')

    for bad, why in [(trios + [from_tuples((('H', 'G'), 'C'))], 'contradiction'),
                     ([trio1, from_tuples((('O', 'X'), 'H'))], 'unresolved')]:
        try:
            build_aho(bad)
        except ValueError as e:
            print(f'build_aho rejected ({why}):', e)
        else:
            raise AssertionError(why)

    # ----- iterative grafting -----
    tree = copy_tree(trio1)
    print('\nstart:', to_tuples(tree))
    for trio in trios[1:]:
        tree = add_trio(tree, trio)
        print('+', to_tuples(trio), '->', to_tuples(tree))

    assert to_tuples(tree) == expected
    assert all(displays(tree, t) for t in trios)
    assert tree.key == expected
    assert to_tuples(build_tree(trios)) == expected

    # A trio whose leaves are all present is only checked.
    assert to_tuples(add_trio(tree, from_tuples((('C', 'N'), 'G')))) == expected
    for bad in [(('H', 'G'), 'C'), ('H', ('C', 'N')), (('A', 'B'), 'D')]:
        try:
            add_trio(tree, from_tuples(bad))
        except ValueError as e:
            print('add_trio rejected:', e)
        else:
            raise AssertionError(bad)

    # Order matters when a trio is added while sharing only one leaf.
    ok = failed = other = 0
    for order in permutations(trios):
        try:
            result = to_tuples(build_tree(order))
        except ValueError:
            failed += 1
        else:
            ok += result == expected
            other += result != expected
        assert to_tuples(build_tree(order, most_overlap_first=True)) == expected
    print(f'in given order: {ok} of 24 orders give the expected tree, '
          f'{other} give another tree, {failed} hit a conflict')
    print('most_overlap_first=True: all 24 orders give the expected tree')

    # ----- node map -----
    tree = build_aho(trios)
    mapping = node_map(tree, trios)
    print('\nnode map: tree node [i] <- (trio t, node j) ...')
    for i, node in enumerate(preorder(tree)):
        origins = ', '.join(f'({t}, {j}) {to_tuples(preorder(trios[t])[j])}' for t, j in mapping[i])
        print(f'  [{i}] {to_tuples(node)} <- {origins}')
    assert mapping == {
        0: [(1, 0), (3, 0)], 1: [(0, 0), (1, 1)], 2: [(0, 1), (2, 0)], 3: [(2, 1)],
        4: [(0, 2), (1, 2), (2, 2), (3, 1)], 5: [(2, 3)], 6: [(0, 3), (2, 4)],
        7: [(0, 4), (1, 3)], 8: [(3, 2)], 9: [(1, 4), (3, 3)], 10: [(3, 4)]}
    # every trio node is mapped exactly once, every tree node is built from something,
    # and the map is the same for the iteratively built tree
    assert sorted(pair for pairs in mapping.values() for pair in pairs) == [
        (t, j) for t in range(4) for j in range(5)]
    assert all(mapping.values())
    assert node_map(build_tree(trios), trios) == mapping
    try:
        node_map(tree, trios + [from_tuples((('H', 'G'), 'C'))])
    except ValueError as e:
        print('node_map rejected:', e)
    else:
        raise AssertionError

    # from_tuples creates nodes in pre-order
    created = []
    original_init = TreeNode.__init__

    def recording_init(self, key):
        original_init(self, key)
        created.append(self)

    TreeNode.__init__ = recording_init
    try:
        recorded = from_tuples(expected)
    finally:
        TreeNode.__init__ = original_init
    assert created == preorder(recorded)
    print('from_tuples creates nodes in pre-order: confirmed')
