import pandas as pd

def count_samples(roh, count):
    samples = []
    for i in range(1, 278):
        counts = str(roh.count("S" + str(i) + "\n"))
        if counts == count:
            samples.append("S" + str(i))
    return samples

def write_output(df, pool, filename):
    output = df.query('POOL == @pool')
    output[["CHR", "BP1", "BP2"]].to_csv(filename, sep=' ', index=False, header=False)

ovi_roh = ""
botat_roh = ""
evansia_roh = ""
evansib_roh = ""
evansic_roh = ""
gambiense_roh = ""
gambienseII_roh = ""
brucei_roh = ""
rhodesiense_roh = ""

with open("data/plink/plink.hom.overlap.csv") as file:
    for line in file:
        fields = line.rstrip("\n").split(",")
        if "T.b.equiperdumtypeOVI" in fields:
            ovi_roh += str(fields[0]) + "\n"
        if "T.b.equiperdumtypeBoTat" in fields:
            botat_roh += str(fields[0]) + "\n"
        if "T.b.evansitypeA" in fields:
            evansia_roh += str(fields[0]) + "\n"
        if "T.b.evansitypeB" in fields:
            evansib_roh += str(fields[0]) + "\n"
        if "T.b.evansitypeC" in fields:
            evansic_roh += str(fields[0]) + "\n"
        if "T.b.gambiense" in fields:
            gambiense_roh += str(fields[0]) + "\n"
        if "T.b.gambienseII" in fields:
            gambienseII_roh += str(fields[0]) + "\n"
        if "T.b.brucei" in fields:
            brucei_roh += str(fields[0]) + "\n"
        if "T.b.rhodesiense" in fields:
            rhodesiense_roh += str(fields[0]) + "\n"

ovi = count_samples(ovi_roh, "4")
botat = count_samples(botat_roh, "2")
evansia = count_samples(evansia_roh, "28")
evansib = count_samples(evansib_roh, "2")
evansic = count_samples(evansic_roh, "1")
gambiense = count_samples(gambiense_roh, "20")
gambienseII = count_samples(gambienseII_roh, "5")
brucei = count_samples(brucei_roh, "12")
rhodesiense = count_samples(rhodesiense_roh, "9")

df = pd.read_csv('data/plink/plink.hom.overlap.csv')
union = df.loc[df['FID'] == "UNION"]

union['BP1'] = union['BP1'].astype(str).str.replace('.0', '')
union['BP2'] = union['BP2'].astype(str).str.replace('.0', '')

write_output(union, ovi, "data/circos/equi_ovi_roh.txt")
write_output(union, botat, "data/circos/equi_botat_roh.txt")
write_output(union, evansia, "data/circos/evansi_a_roh.txt")
write_output(union, evansib, "data/circos/evansi_b_roh.txt")
write_output(union, evansic, "data/circos/evansi_c_roh.txt")
write_output(union, gambiense, "data/circos/gambiense_out_roh.txt")