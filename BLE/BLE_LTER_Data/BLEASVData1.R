library(readr)
BLE_ASVs_csv <- read_csv("~/Desktop/SIO/Classes/2023-24/Winter2024/SIOB278_MarineMicrobialSeminar/GITHUB/BLE/BLE_ASVs.csv.gz")
View(BLE_ASVs_csv)

BLE_ASVs_flipped <- t(BLE_ASVs_csv)
