# %%
from itertools import combinations_with_replacement
# import hashlib

###
### uv run fasta2smi -i stupid_plan/peps/ligands.fasta -o stupid_plan/peps/smiles.txt
###


def main(filename="stupid_plan/peps/ligands.fasta"):
    lex = [
        "K",
        "E",
        "D",
        "R",
    ]

    l1 = list(combinations_with_replacement(lex, 1))
    l2 = list(combinations_with_replacement(lex, 2))
    l3 = list(combinations_with_replacement(lex, 3))
    l_all = l1 + l2 + l3

    # %%
    with open(filename, "w") as f:
        for (
            i,
            j,
        ) in enumerate(l_all):
            p = "".join(j)
            # hex_hash = hashlib.md5(p.encode("utf-8")).hexdigest()
            header = f"P_{p}_{str(len(p))}_{str(len(set(p)))}"
            f.write(f">{header}\n")
            f.write(f"{p}\n")

        # with open("stupid_plan/mols/all.csv", "r") as f2:
        #     lines = f2.readlines()
        #     for line in lines[0:]:
        #         fields = line.strip().split(",")


if __name__ == "__main__":
    main()
    print("Done writing ligands to file.")
