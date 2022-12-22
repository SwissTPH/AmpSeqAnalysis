#######################
# All the commands needed to install HaplotypR
#
# created 21.10.2022
# monica.golumbeanu@unibas.ch
#######################

# Installing ShortRead. 
# This installation will take some time especially if yout decide to update all packages!
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead")

# Installing devtools
install.packages("devtools")

# Installing git2r
install.packages("git2r")

# Installing Rswarm, Rvsearch and HaplotypR
library(devtools)
library(git2r)
path = file.path(tempfile(pattern="Rswarm-"), "Rswarm")
dir.create(path, recursive=TRUE)
repo = clone("https://github.com/lerch-a/Rswarm.git", path)
clone("https://github.com/torognes/swarm.git", file.path(path, "src", "swarm"))
install(path)

path = file.path(tempfile(pattern="Rvsearch-"), "Rvsearch")
dir.create(path, recursive=TRUE)
repo = clone("https://github.com/lerch-a/Rvsearch.git", path)
clone("https://github.com/torognes/vsearch.git", file.path(path, "src", "vsearch"))
install(path)

# This seems to throw an error so we don't run it
# detach("package:HaplotypR", unload=TRUE)

devtools::install_github("lerch-a/HaplotypR", force = TRUE)
