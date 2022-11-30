# importing library

library(dplyr)
library(dslabs)
library(readr)
library(ggplot2)
library(esquisse)

# importing data

df <- read.csv2("base de dados.csv")
str(df)

# adjusting data

# transforming data from logical to numeric. fixing the errors

df$ExactMolWt=as.numeric(df$ExactMolWt)
df$MolLogp=as.numeric(df$MolLogp)
df$TPSA=as.numeric(df$TPSA)

str(df)

# calculating 

# mean valors

names(df)

mean(df$NumAtoms)
mean(df$ExactMolWt)
mean(df$NumRotableBonds)
mean(df$MolLogp)
mean(df$RingCount)
mean(df$NumHAcceptors)
mean(df$TPSA)
mean(df$NumHDonors)

# max valors

max(df$NumAtoms)
max(df$ExactMolWt)
max(df$NumRotableBonds)
max(df$MolLogp)
max(df$RingCount)
max(df$NumHAcceptors)
max(df$TPSA)
max(df$NumHDonors)

# min valors
min(df$NumAtoms)
min(df$ExactMolWt)
min(df$NumRotableBonds)
min(df$MolLogp)
min(df$RingCount)
min(df$NumHAcceptors)
min(df$TPSA)
min(df$NumHDonors)

# making data frame

MMM <- data.frame(Descriptors = c("NumAtoms", "ExactMolWt", "NumRotableBonds",
                        "MolLogP", "RingCount", "NumHAcceptors",
                        "TPSA", "NumHdonors"),
                 Max = c(max(df$NumAtoms),
                         max(df$ExactMolWt),
                         max(df$NumRotableBonds),
                         max(df$MolLogp),
                         max(df$RingCount),
                         max(df$NumHAcceptors),
                         max(df$TPSA),
                         max(df$NumHDonors)),
                 Mean = c(mean(df$NumAtoms),
                          mean(df$ExactMolWt),
                          mean(df$NumRotableBonds),
                          mean(df$MolLogp),
                          mean(df$RingCount),
                          mean(df$NumHAcceptors),
                          mean(df$TPSA),
                          mean(df$NumHDonors)),
                 Min = c(min(df$NumAtoms),
                         min(df$ExactMolWt),
                         min(df$NumRotableBonds),
                         min(df$MolLogp),
                         min(df$RingCount),
                         min(df$NumHAcceptors),
                         min(df$TPSA),
                         min(df$NumHDonors)),
                 stringsAsFactors = FALSE)

# max molecular formula valors

df$MolecularFormula[which.max(df$NumAtoms)]
df$MolecularFormula[which.max(df$ExactMolWt)]
df$MolecularFormula[which.max(df$NumRotableBonds)]
df$MolecularFormula[which.max(df$MolLogp)]
df$MolecularFormula[which.max(df$RingCount)]
df$MolecularFormula[which.max(df$NumHAcceptors)]
df$MolecularFormula[which.max(df$TPSA)]
df$MolecularFormula[which.max(df$NumHDonors)]

# making max molecular formula data frame valors

MAXMF <- data.frame(Descriptors = c("NumAtoms", "ExactMolWt", "NumRotableBonds",
                                       "MolLogP", "RingCount", "NumHAcceptors",
                                       "TPSA", "NumHdonors"),
                  MolecularFormula = c(df$MolecularFormula[which.max(df$NumAtoms)],
                          df$MolecularFormula[which.max(df$ExactMolWt)],
                          df$MolecularFormula[which.max(df$NumRotableBonds)],
                          df$MolecularFormula[which.max(df$MolLogp)],
                          df$MolecularFormula[which.max(df$RingCount)],
                          df$MolecularFormula[which.max(df$NumHAcceptors)],
                          df$MolecularFormula[which.max(df$TPSA)],
                          df$MolecularFormula[which.max(df$NumHDonors)]),
                  Max = c(max(df$NumAtoms),
                          max(df$ExactMolWt),
                          max(df$NumRotableBonds),
                          max(df$MolLogp),
                          max(df$RingCount),
                          max(df$NumHAcceptors),
                          max(df$TPSA),
                          max(df$NumHDonors)),
                  stringsAsFactors = FALSE)

# min molecular formula valors

df$MolecularFormula[which.min(df$NumAtoms)]
df$MolecularFormula[which.min(df$ExactMolWt)]
df$MolecularFormula[which.min(df$NumRotableBonds)]
df$MolecularFormula[which.min(df$MolLogp)]
df$MolecularFormula[which.min(df$RingCount)]
df$MolecularFormula[which.min(df$NumHAcceptors)]
df$MolecularFormula[which.min(df$TPSA)]
df$MolecularFormula[which.min(df$NumHDonors)]

# making max molecular formula data frame valors

MINMF <- data.frame(Descriptors = c("NumAtoms", "ExactMolWt", "NumRotableBonds",
                                  "MolLogP", "RingCount", "NumHAcceptors",
                                  "TPSA", "NumHdonors"),
                  MolecularFormula = c(df$MolecularFormula[which.min(df$NumAtoms)],
                                       df$MolecularFormula[which.min(df$ExactMolWt)],
                                       df$MolecularFormula[which.min(df$NumRotableBonds)],
                                       df$MolecularFormula[which.min(df$MolLogp)],
                                       df$MolecularFormula[which.min(df$RingCount)],
                                       df$MolecularFormula[which.min(df$NumHAcceptors)],
                                       df$MolecularFormula[which.min(df$TPSA)],
                                       df$MolecularFormula[which.min(df$NumHDonors)]),
                  Min = c(min(df$NumAtoms),
                          min(df$ExactMolWt),
                          min(df$NumRotableBonds),
                          min(df$MolLogp),
                          min(df$RingCount),
                          min(df$NumHAcceptors),
                          min(df$TPSA),
                          min(df$NumHDonors)),
                  stringsAsFactors = FALSE)


# ploting datas

# ploting ExactMolWt x NumAtoms
ggplot(df) +
  aes(x = ExactMolWt, y = NumAtoms, colour = ExactMolWt) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  labs(title = "Plot", subtitle = "ExactMolWt x NumAtoms") +
  theme_bw()

# ploting ExactMolWt x RingCount
ggplot(df) +
  aes(x = ExactMolWt, y = RingCount, colour = ExactMolWt) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  labs(title = "Plot", subtitle = "ExactMolWt x RingCount") +
  theme_bw()

# ploting ExacMolWt x TPSA
ggplot(df) +
  aes(x = ExactMolWt, y = TPSA, colour = ExactMolWt) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  labs(title = "Plot ", subtitle = "ExactMolWt x TPSA") +
  theme_bw()

# ploting ExactMolWt x NumHAcceptors
ggplot(df) +
  aes(x = ExactMolWt, y = NumHAcceptors, colour = ExactMolWt) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  labs(
    title = "Plot",
    subtitle = "ExactMolWt x NumHAcceptors"
  ) +
  theme_bw()

# ploting ExactMolWt x NumHDonors
ggplot(df) +
  aes(x = ExactMolWt, y = NumHDonors, colour = ExactMolWt) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  labs(title = "Plot", subtitle = "ExactMolWt x NumHDonors") +
  theme_bw()

# descriptors histogram

# NumAtoms
ggplot(df) +
  aes(x = NumAtoms) +
  geom_histogram(bins = 30L, fill = "#0C4C8A") +
  labs(title = "Histogram", subtitle = "NumAtoms") +
  theme_bw()

# ExactMolWt
ggplot(df) +
  aes(x = ExactMolWt) +
  geom_histogram(bins = 30L, fill = "#0C4C8A") +
  labs(title = "Histogram", subtitle = "ExactMolWt") +
  theme_bw()

# NumRotableBonds
ggplot(df) +
  aes(x = NumRotableBonds) +
  geom_histogram(bins = 30L, fill = "#0C4C8A") +
  labs(title = "Histogram", subtitle = "NumRotableBonds") +
  theme_bw()

# MolLogP
ggplot(df) +
  aes(x = MolLogp) +
  geom_histogram(bins = 30L, fill = "#0C4C8A") +
  labs(title = "Histogram", subtitle = "MolLogP") +
  theme_bw()

# RingCount
ggplot(df) +
  aes(x = RingCount) +
  geom_histogram(bins = 30L, fill = "#0C4C8A") +
  labs(title = "Histogram", subtitle = "RingCount") +
  theme_bw()

# NumHAcceptors
ggplot(df) +
  aes(x = NumHAcceptors) +
  geom_histogram(bins = 30L, fill = "#0C4C8A") +
  labs(title = "Histogram", subtitle = "NumHAcceptors") +
  theme_bw()

# TPSA

ggplot(df) +
  aes(x = TPSA) +
  geom_histogram(bins = 30L, fill = "#0C4C8A") +
  labs(title = "Histogram", subtitle = "TPSA") +
  theme_bw()

# NumHDonors
ggplot(df) +
  aes(x = NumHDonors) +
  geom_histogram(bins = 30L, fill = "#0C4C8A") +
  labs(title = "Histogram", subtitle = "NumHDonors") +
  theme_bw()

