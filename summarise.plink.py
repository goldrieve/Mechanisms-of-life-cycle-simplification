import pandas as pd

ovi_roh = ""
botat_roh = ""
evansia_roh = ""
evansib_roh = ""
evansic_roh = ""
gambiense_roh = ""
gambienseII_roh = ""
brucei_roh = ""
rhodesiense_roh = ""

for line in open("data/plink/plink.hom.overlap.csv"):
    fields = line.rstrip("\n").split(",")
    if "T.b.equiperdumtypeOVI" in fields:
        ovi_roh = ovi_roh + str(fields[0]) + "\n"
    if "T.b.equiperdumtypeBoTat" in fields:
        botat_roh = botat_roh + str(fields[0]) + "\n"
    if "T.b.evansitypeA" in fields:
        evansia_roh = evansia_roh + str(fields[0]) + "\n"
    if "T.b.evansitypeB" in fields:
        evansib_roh = evansib_roh + str(fields[0]) + "\n"
    if "T.b.evansitypeC" in fields:
        evansic_roh = evansic_roh + str(fields[0]) + "\n"
    if "T.b.gambiense" in fields:
        gambiense_roh = gambiense_roh + str(fields[0]) + "\n"
    if "T.b.gambienseII" in fields:
        gambienseII_roh = gambienseII_roh + str(fields[0]) + "\n"
    if "T.b.brucei" in fields:
        brucei_roh = brucei_roh + str(fields[0]) + "\n"
    if "T.b.rhodesiense" in fields:
        rhodesiense_roh = rhodesiense_roh + str(fields[0]) + "\n"

ovi = []
for i in range (1,278):
    counts = str(ovi_roh.count("S" + str(i) + "\n"))
    if counts == "4":
        ovi = ovi + ["S" + str(i)]

botat = []
for i in range (1,278):
    counts = str(botat_roh.count("S" + str(i) + "\n"))
    if counts == "2":
        botat = botat + ["S" + str(i)]

evansia = []
for i in range (1,278):
    counts = str(evansia_roh.count("S" + str(i) + "\n"))
    if counts == "28":
        evansia = evansia + ["S" + str(i)]

evansib = []
for i in range (1,278):
    counts = str(evansib_roh.count("S" + str(i) + "\n"))
    if counts == "2":
        evansib = evansib + ["S" + str(i)]

evansic = []
for i in range (1,278):
    counts = str(evansic_roh.count("S" + str(i) + "\n"))
    if counts == "1":
        evansic = evansic + ["S" + str(i)]

gambiense = []
for i in range (1,278):
    counts = str(gambiense_roh.count("S" + str(i) + "\n"))
    if counts == "20":
        gambiense = gambiense + ["S" + str(i)]

gambienseII = []
for i in range (1,278):
    counts = str(gambienseII_roh.count("S" + str(i) + "\n"))
    if counts == "5":
        gambienseII = gambienseII + ["S" + str(i)]

brucei = []
for i in range (1,278):
    counts = str(brucei_roh.count("S" + str(i) + "\n"))
    if counts == "12":
        brucei = brucei + ["S" + str(i)]

rhodesiense = []
for i in range (1,278):
    counts = str(rhodesiense_roh.count("S" + str(i) + "\n"))
    if counts == "9":
        rhodesiense = rhodesiense + ["S" + str(i)]

df = pd.read_csv('data/plink/plink.hom.overlap.csv')
union = df.loc[df['FID'] == "UNION"]

union['BP1'] = union['BP1'].astype(str).str.replace('.0', '')
union['BP2'] = union['BP2'].astype(str).str.replace('.0', '')

ovi_out = union.query('POOL == @ovi')
botat_out = union.query('POOL == @botat')
evansia_out = union.query('POOL == @evansia')
evansib_out = union.query('POOL == @evansib')
evansic_out = union.query('POOL == @evansic')
gambiense_out = union.query('POOL == @gambiense')

ovi_out[["CHR", "BP1", "BP2"]].to_csv("data/circos/equi_ovi_roh.txt", sep=' ', index=False, header=False)
botat_out[["CHR", "BP1", "BP2"]].to_csv("data/circos/equi_botat_roh.txt", sep=' ', index=False, header=False)
evansia_out[["CHR", "BP1", "BP2"]].to_csv("data/circos/evansi_a_roh.txt", sep=' ', index=False, header=False)
evansib_out[["CHR", "BP1", "BP2"]].to_csv("data/circos/evansi_b_roh.txt", sep=' ', index=False, header=False)
evansic_out[["CHR", "BP1", "BP2"]].to_csv("data/circos/evansi_c_roh.txt", sep=' ', index=False, header=False)
gambiense_out[["CHR", "BP1", "BP2"]].to_csv("data/circos/gambiense_out_roh.txt", sep=' ', index=False, header=False)
