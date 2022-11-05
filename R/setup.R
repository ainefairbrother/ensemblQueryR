library(devtools)
library(usethis)

usethis::use_git()

gitcreds::gitcreds_set("https://github.com/ainefairbrother/ensemblQueryR.git")

# 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

usethis::create_from_github(repo_spec = "https://github.com/ainefairbrother/ensemblQueryR.git", fork=F)
