install_github(repo="DreamBox7", username="weiyi-bitw", ref="master")
require(DreamBox7)
require(synapseClient)
require(fields)

synapseLogin()


#==================================================
#
# Figure 1: Top three genes of the CIN attractor metagene in pancan 12
#
#==================================================


gene1 = "CENPA"
gene2 = "DLGAP5"
gene3 = "MELK"

fig_cin <- pancan12Scatter(gene1, gene2, gene3)

#==================================================
#
# Figure 2: Top three genes of the MES attractor metagene in pancan 12
#
#==================================================


gene1 = "COL5A2"
gene2 = "SPARC"
gene3 = "VCAN"

fig_mes <- pancan12Scatter(gene1, gene2, gene3)

#==================================================
#
# Figure 3: Top three genes of the LYM attractor metagene in pancan 12
#
#==================================================

gene1 = "PTPRC"
gene2 = "CD53"
gene3 = "LCP2"

fig_lym <- pancan12Scatter(gene1, gene2, gene3)

#==================================================
#
# Figure 4: Top three genes of the chr8q24.3 attractor metagene in pancan 12
#
#==================================================

gene1 = "EXOSC4"
gene2 = "PUF60"
gene3 = "BOP1"

fig_8q24_3 <- pancan12Scatter(gene1, gene2, gene3)

