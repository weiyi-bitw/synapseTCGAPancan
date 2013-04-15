install.packages(c("devtools", "survival", "rms", "fields"))

source("http://depot.sagebase.org/CRAN.R")
pkgInstall("synapseClient", "staging")

require(devtools)
install_github("devtools", ref="devtools-1.1")
install_github("rGithubClient", "brian-bot", ref="rGithubClient-0.7")
install_github("cafr", "weiyi-bitw", ref="master")
install_github("DreamBox7", "weiyi-bitw", ref="master")

