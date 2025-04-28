

# Adaptation.Rdata is for the number of reactions added to the models, including its subsystem.
# re2.Rdata is for all of the reactions in supermodel, showing that this reaction exists in the E. coli model or not, including its subsystem information.
# NS_E.coli.Rdata is the total number of simulations for respective BOF and medium type.

library(ggplot2)
##for E.coli 
load("~/DATA/E.coli/Randomization simulation/randomization_ecoli.Rdata")
##for E.coli model simulation
load("~/DATA/E.coli/Randomization simulation/Normal OBJ/randomization.Rdata")


# Extract the model names
model_names <- colnames(re2)[-(1:2)]
# Function to get number of simulations for models where the reaction is absent
get_simulation_count <- function(row) {
  absent_models <- model_names[row[model_names] == 0]
  total_simulations <- sum(NS_E.coli$nSimulation[NS_E.coli$model %in% absent_models])
  return(total_simulations)
}
# Apply the function to each row
re2$nSimulations.without.reaction <- apply(re2, 1, get_simulation_count)
# This will give you the desired dataframe containing, for each reaction, 
# the subsystem and the total number of simulations run with models 
# where the reaction was absent.
simulations <- re2[, c("Reaction", "Subsystem", "nSimulations.without.reaction")]
hist(re2$nSimulations.without.reaction) # most reactions are missing in almost all simulations (n=27806)

# remove inconsistency in the model names:
model_names <- intersect(model_names, colnames(Adaptation)[-(1:2)])
# remove reactions without any simulations 
simulations <- subset(simulations, nSimulations.without.reaction > 0)
# Function to get sum of gains for each reaction
get_gains_count <- function(reaction_name) {
  # If the reaction is not in the Adaptation dataframe, return 0
  if (!reaction_name %in% Adaptation$reaction) return(0)
  
  # Otherwise, sum the gains for that reaction across all models
  gains_for_reaction <- Adaptation[Adaptation$reaction == reaction_name, model_names]
  total_gains <- sum(as.numeric(gains_for_reaction))
  return(total_gains)
}

# Apply the function to each reaction in reactions_df
# - Now reactions_df will have a new column "total.gains" with the sum of gains for each reaction:
simulations$total.gains <- sapply(simulations$Reaction, get_gains_count)
simulations <- simulations[!is.na(simulations$Subsystem), ]

# Do NOT remove reactions that were never gained, as these are irrelevant - that would ruin the statistics for enrichment!
# simulations <- subset(simulations, total.gains > 0)

# function to calculate the mean fraction of gains per simulation across the genes in each subsystem
# - the variable "subsystems" is ONLY used at the end to make sure that the results are sorted consistently
gains.perSimulation <- function(simulations, subsystems) {
  # gains per simulation over all reactions and subsystems:
  simulations.total <- sum(simulations$nSimulations.without.reaction)
  gains.total <- sum(simulations$total.gains)  # overall mean = gains.total/simulations.total
  gains.perSimulation.total <- gains.total/simulations.total
  # gains per simulation, averaged over all reactions in each subsystem:
  gains.subsystem <- aggregate(total.gains ~ Subsystem, data=simulations, sum) # this is now sum(total.gains) instead of mean(gains.perSimulation)
  simulations.subsystem <- aggregate(nSimulations.without.reaction ~ Subsystem, data=simulations, sum)
  results <- merge(gains.subsystem, simulations.subsystem, by="Subsystem")
  rownames(results) <- results$Subsystem
  # count number of reactions in subsystem:
  simulations$count <- 1
  tmp <- aggregate(count ~ Subsystem, data=simulations, sum)
  results$n.reactions <- tmp$count
  
  results$gains.perSimulation <- results$total.gains / results$nSimulations.without.reaction
  
  # Calculate the gains per simulation for the complement of each subsystem S
  results$gains.perSimulation.notInSubsys <- (gains.total - results$total.gains) / (simulations.total - results$nSimulations.without.reaction)
  
  # calculate odds ratio: 
  results$OR <- results$gains.perSimulation / results$gains.perSimulation.notInSubsys
  
  # Order results based on input subsystems and return
  ordered_results <- results[subsystems, ]
  return(ordered_results)
}

# for debugging: generate some mock data with known results:
if (FALSE) {
  subsystems <- c('A','B','C','D','E')
  simulations <- data.frame(Reaction=1:100, Subsystem=rep(subsystems, times=20), 
                            nSimulations.without.reaction=rep(c(1000,1000,1000,1000,100000), times=20), 
                            total.gains=rep(c(10,0,0,1,10), times=20))
}

# calculate the number of reactions in each subsystem:
reaction.counts <- table(simulations$Subsystem)
simulations$gains.perSimulation <- simulations$total.gains / simulations$nSimulations.without.reaction

# These are the observed results:
subsystems <- unique(simulations$Subsystem)  # subsystems.nonZero <- unique(simulations[simulations$total.gains > 0, "Subsystem"])
results.observed <- gains.perSimulation(simulations, subsystems)
# gains.perSimulation.observed <- results.observed$gains.perSimulation
# OR.observed <- results.observed$OR
subsystems.nonZero <- unique(simulations[simulations$total.gains > 0, "Subsystem"])
subsystems.nonZero <- as.character(subsystems.nonZero$Subsystem)
###### Randomizations #########
set.seed(100000)
R <- 100000 # use 10000 or more for the final calculation (if it doesn't take too long)
# Initialize two dataframes to store R simulation results for each subsystem
# gains.perSimulation.perReact.rand <- data.frame(matrix(nrow = length(subsystems), ncol = R))
# rownames(gains.perSimulation.perReact.rand) <- subsystems
# OR.rand <- data.frame(matrix(nrow = length(subsystems), ncol = R))
gains.perSimulation.rand <- data.frame(row.names=subsystems)
OR.rand <- data.frame(row.names=subsystems)
for (i in 1:R) {
  simulations.rand <- simulations
 # simulations.rand$Subsystem <- sample(subsystems.nonZero,length(simulations.rand$Subsystem),replace = T)
 simulations.rand$Subsystem <- sample(simulations.rand$Subsystem)
  results.i <- gains.perSimulation(simulations.rand, subsystems)
  gains.perSimulation.rand[ , i] <- results.i$gains.perSimulation
  OR.rand[ , i] <- results.i$OR
}

# calculate p-values:
results <- results.observed # data.frame(gains.perSimulation=as.numeric(gains.perSimulation.observed), OR=as.numeric(OR.observed))
# (the as.numeric(...) should not be necessary, but otherwise we get a strange bug...)
rownames(results) <- subsystems
results <- results[ , -1]
# results$gains.perSimulation.rand <- apply(gains.perSimulation.rand, 1, mean, na.rm=TRUE)
# results$OR.rand.median <- apply(OR.rand, 1, median, na.rm=TRUE)
# OR.rand1 <- OR.rand
# OR.rand1[OR.rand1==0] <- 1e-6   # remove zeros before taking the log
# results$log10.OR.rand <- apply(log10(OR.rand1), 1, median, na.rm=TRUE)
n.valid.rand <- apply(!is.na(gains.perSimulation.rand), 1, sum)
results$p <- apply(gains.perSimulation.rand >= results$gains.perSimulation, 1, sum, na.rm=TRUE) / n.valid.rand
results[n.valid.rand < R/2, "p"] <- NA    # remove p-value if there are not enough valid randomizations
results$p.adj.BH <- p.adjust(results$p, method="BH")

# remove those with 0 observed gains - these cannot be enriched:
results1 <- subset(results, gains.perSimulation > 0)
results.sorted <- results1[order(results1$p), ]
ggplot(results.sorted, aes(x=p)) + geom_histogram()
# ggplot(results.sorted, aes(x = p)) + stat_ecdf(geom = "step", color="blue") 
# ggplot(results.sorted, aes(x = p.adj.BH)) + stat_ecdf(geom = "step", color="blue") 

# the true OR is nicely distributed around 1 (Gaussian on log-scale)
ggplot(results.sorted, aes(x=OR)) + geom_histogram() + scale_x_log10() + geom_vline(xintercept = 1, color = "red")
# the random OR should show the same (maybe more narrow) distribution around 1, but it doesn't 
# - the reason is that for the true data, we ignore those with 0 gains (OR=0), but for the random data we can't
# -- so OR=0 for most models
# ggplot(results.sorted, aes(x=log10.OR.rand)) + geom_histogram()
# ggplot(results, aes(x=OR.rand.median)) + geom_histogram() + geom_vline(xintercept = 1, color = "red")#+ scale_x_log10()
# expected gains per simulation (averaged over all reactions):
expect <- mean(simulations$gains.perSimulation)
# ggplot(results.sorted, aes(x=gains.perSimulation.rand)) + geom_histogram() + geom_vline(xintercept = expect, color = "red") # scale_x_log10() 

# significant results:
results.significant <- results.sorted[results.sorted$p.adj.BH < 0.05, ]
results.significant
head(results.sorted, n=20)
plot(results.sorted$p.adj.BH)
sum(results.sorted$p.adj.BH < 0.05)
write.table(results.sorted, file="~/DATA/E.coli/Randomization simulation/randomization_results.tab", sep="\t")
#write.csv2(results.sorted, file="randomization.csv")

# for display:
results.short <- results.sorted
rownames(results.short) <- 1:nrow(results.short)
head(results.short, 20)

