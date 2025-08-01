# %%
def main():
    with open("stupid_plan/mols/all_mols.txt", "w") as f0:
        with open("stupid_plan/peps/smiles.txt", "r") as f1:
            i = 0
            for line in f1:
                fields = line.strip().split(": ")
                header = (
                    fields[0]
                    .replace(",", "")
                    .replace(":", "")
                    .replace("-", "_")
                    .strip()
                )
                # hex_hash = hashlib.md5(fields[0].encode("utf-8")).hexdigest() # unused
                out = f"'S{i}_{header}'" + ",'" + fields[1] + "'\n"
                i += 1
                print(out)
                f0.write(out)

            with open("stupid_plan/mols/all.txt", "r") as f2:
                for line in f2:
                    fields = line.strip().split(",")
                    out = f"'S{i}_{fields[1]}'" + "," + fields[2] + "'\n"
                    i += 1
                    print(out)
                    f0.write(out)


# %%
if __name__ == "__main__":
    main()
    print("Done writing all_mols.txt.")
