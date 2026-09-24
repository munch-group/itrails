"""Turn an msprime Demography (or its ``asdict()`` output) into a Newick string.

Population-split events define the tree: every PopulationSplit event makes the
``ancestral`` population the parent of each ``derived`` population, and the
event ``time`` is the height of the ancestral node above the present. Leaf
populations sit at their ``default_sampling_time`` (0 if None). Branch lengths
are height differences, in generations.

Only PopulationSplit events are used. Admixture events cannot be represented as
a tree and raise a ValueError; other events (size changes, migration changes)
are ignored.

    >>> demography_to_newick(spec)
    '((hg38:272982.0827,panTro5:272982.0827)hg38_panTro5:21456.7795,gorGor5:294438.8622)hg38_panTro5_gorGor5;'
"""


def demography_to_newick(demography, internal_labels=True, fmt="%.10g"):
    spec = demography.asdict() if hasattr(demography, "asdict") else demography
    populations = [p["name"] for p in spec["populations"]]
    leaf_time = {
        p["name"]: p.get("default_sampling_time") or 0.0 for p in spec["populations"]
    }

    children, height = {}, {}
    for ev in spec.get("events", []):
        cls = ev.get("__class__", "")
        if cls.endswith("Admixture"):
            raise ValueError("demography contains an Admixture event; not a tree")
        if not cls.endswith("PopulationSplit"):
            continue
        anc = ev["ancestral"]
        if anc in children:
            raise ValueError(f"population {anc!r} is ancestral in more than one split")
        children[anc] = list(ev["derived"])
        height[anc] = ev["time"]

    parent = {d: a for a, ds in children.items() for d in ds}
    for name in populations:
        if name not in children:
            height[name] = leaf_time[name]
    roots = [name for name in populations if name not in parent]
    if len(roots) != 1:
        raise ValueError(f"expected exactly one root population, found {roots}")

    def rec(name):
        s = ""
        if name in children:
            s = "(" + ",".join(rec(c) for c in children[name]) + ")"
            if internal_labels:
                s += name
        else:
            s = name
        if name in parent:
            length = height[parent[name]] - height[name]
            if length < 0:
                raise ValueError(f"negative branch above {name!r}: split times are not ordered")
            s += ":" + fmt % length
        return s

    return rec(roots[0]) + ";"


if __name__ == "__main__":
    import ast, sys
    print(demography_to_newick(ast.literal_eval(sys.stdin.read())))
