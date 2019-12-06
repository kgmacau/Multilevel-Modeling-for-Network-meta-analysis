### Frequentist method
# Multiple treatments comparisons using frequentist ranking probability method
library(netmeta)

# Senn data
data(Senn2013)
net1 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data = Senn2013, sm = "MD",
                comb.random = FALSE)

nr1 <- netrank(net1)
nr1




# AI algorithm data
AI<-read.csv("C:/Users/212697818/Desktop/AI.csv")
net2 <- netmeta(TE, seTE, treat1, treat2, studlab,
                data = AI, sm = "MD",
                comb.random = FALSE)

nr2 <- netrank(net2)
nr2


### Bayesian method
# Multiple treatments comparisons using bayesian ranking probability method
library(gemtc)
install.packages("rjags")
library(rjags)

reAI<-read.csv("C:/Users/212697818/Desktop/reAI.csv")
network <- mtc.network(reAI)
summary(network)
plot(network)
model <- mtc.model(network, 
                   linearModel = "random",
                   n.chain = 4)

mcmc1 <- mtc.run(model, n.adapt = 50, n.iter = 1000, thin = 10)
mcmc2 <- mtc.run(model, n.adapt = 5000, n.iter = 100000, thin = 10)

gelman.plot(mcmc1)
gelman.plot(mcmc2)


rank <- rank.probability(mcmc2, preferredDirection = -1)
plot(rank, beside=TRUE, cex.names=0.5)

# Reference
# https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/bayesian-network-meta-analysis.html