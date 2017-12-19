library(bnlearn)
library(igraph)
library(Matrix)
library(ggplot2)
library(gRain)
library(lattice)
library(RColorBrewer)
library(gridExtra)
library(Rgraphviz)
library(parallel)
library(SID)


# For some reason these variables need to be instantiated outside 
# the function in the global scope.
a.cur.lv = "a"
b.cur.lv = "b"
c.cur.lv = "c"


# This calculates the KL-Divergence from the
# probability distribution associated with the original BN to the 
# one reconstructed from the data. It also calculates the 
# SID between the original struture and the reconstructed one.
# These are calculated "times.to.repeat" times, and then averaged,
# where a fresh BN is reconstructed each time on the same number of
# samples.
# It also finds (exactly) the probabilities P(D|A) and P(D|~A), 
# however only once as it can be quite expensive to do.
KLMean = function(bn,
                  num.samples.to.gen,
                  times.to.repeat = 1,
                  test.funct = "x2",
                  score.funct = "bde",
                  param.learn.method = "mle"){
  # Holds running sum of graph edit distances
  sumSID = 0
  klsum = 0
  
  # Create an adjacency matrix for the original BN
  bn.am = as(bnlearn::as.graphNEL(bn), "matrix")
  
  
  for (i in 1:times.to.repeat){
    # Randomly generate 'num.samples.to.gen' samples from the original BN.
    gen.samples = rbn(bn, num.samples.to.gen)

    # Learn BN structure from those samples using the restrict-maximise heuristic.
    # Semi-parametric X2 test used for independence (restrict phase).
    # Inter-Associative Markov Blanket (IAMB) used to find the markov blanket of each node, as it is 
    # fast and reasonably accurate and easy to understand.
    # Bayesian Information Criterion score used in the maximise phase - good score of how well the
    # network fits the data, penalises for having too many edges.
    # Hill-climbing algorithm used to explore restricted search-space.
    # Denis & Scutari has good explanation of rsmax2.
    cust = rsmax2(gen.samples, restrict="iamb", maximize = "hc")

    # Learn BN parameters with simple Maximum Likelihood Estimate method.
    bn.lnd.mle = bn.fit(cust, data = gen.samples, method = "mle")
    
    # Generate adjacency matrix for learned network
    lrnd.am = as(bnlearn::as.graphNEL(bn.lnd.mle), "matrix")

    # Try using Structural Intervention Distance rather than GED
    sumSID = sumSID + structIntervDist(bn.am, lrnd.am)$sid
  }
  
  # Take the average SID
  t = sumSID/times.to.repeat

  return(t)
}

# Specify the network structure
dag.test = model2network("[A][B|A][C|A][D|B:C]")
graphviz.plot(dag.test)


# Label levels
A.lv = c("a", "na")
B.lv = c("b", "nb")
C.lv = c("c", "nc")
D.lv = c("d", "nd")

# Define probabilities
A.true    = 0.40
B.A.true  = 0.70
B.nA.true = 0.50
C.A.true  = 0.50
C.nA.true = 0.60

# Assign to arrays
A.prob = array(c(A.true, 1 - A.true), dim = 2, dimnames = list(A = A.lv))
B.prob = array(c(B.A.true, 1 - B.A.true, B.nA.true, 1-B.nA.true), dim = c(2, 2), dimnames = list(B = B.lv, A = A.lv))
C.prob = array(c(C.A.true, 1 - C.A.true, C.nA.true, 1-C.nA.true), dim = c(2, 2), dimnames = list(C = C.lv, A = A.lv))
D.prob = array(c(0.1, 0.9, 0.8, 0.2, 0.4, 0.6, 0.5, 0.5), dim=c(2,2,2), dimnames = list(D = D.lv, B = B.lv, C = C.lv))

# Create distribution & fit to BN
loc.dist = list(A = A.prob, B = B.prob, C = C.prob, D = D.prob)
bn.actual = custom.fit(dag.test, loc.dist)

# Get rid of the heinous default colours when 
# displaying the probability distributions for each node
myColours <- brewer.pal(6,"Greys")
my.settings <- list(
    superpose.polygon=list(col=myColours[1], border="transparent"),
    strip.background=list(col=myColours[1]),
    plot.polygon=list(col=myColours[2]),
    strip.border=list(col="black")
)
trellis.par.set(my.settings)

# Display the probability distributions of each node with barcharts.
prA = bn.fit.barchart(bn.actual$A, main = "P(A)")
prB = bn.fit.barchart(bn.actual$B, main = "P(B|A)")
prC = bn.fit.barchart(bn.actual$C, main = "P(C|A)")
prD = bn.fit.barchart(bn.actual$D, main = "P(D|B,C)")

# This where the data generation takes place. It takes a long time, 
# but saves the results in a text file.

# It goes through from i = min to max.number.of.observations, calling KLMean each time and 
# instructing it to generate "i" samples, and average the values out over 'repeat.times' times.
# This way we get to see how a BN generated from 1,2,3,4,...,max.number.of.observations samples
# compares to the original BN.

repeat.times = 5
min.number.of.observations = 50
max.num.observations = 1000

res.data = data.frame(seq(from = min.number.of.observations, to = max.num.observations), vector(mode = "numeric", length = max.num.observations - min.number.of.observations + 1))
colnames(res.data) = c("Samples", "Mean.SID")

c=1
for (i in seq(from = min.number.of.observations, to = max.num.observations)){
  print(i)
  res.data[c, 1] = i
  res.data[c, 2] = KLMean(bn.actual, i, times.to.repeat = repeat.times)
  c = c + 1
}


# Run below line if generating new data to stop log-problems. 
# Will lose some low-sample data points.
res.data = res.data[is.finite(rowSums(res.data)),]

plot(res.data[,2] ~ res.data[,1])








# 
# 
# 
# 
# 
# 
# 
# low.col = "black"
# high.col = "yellow"
# 
# 
# # Uncomment the below line if reading the data from a file rather than regenerating.
# res.data = read.table("newmetric1000samples25rep.txt", sep="\t")
# res.data = res.data[-c(4,5)]
# 
# plot(res.data)
# 
# # Plot 
# ggplot(res.data, aes(x=Samples, y=Mean.Distance, color=log(KL.Divergence))) + geom_point() + scale_color_gradient(low=low.col, high=high.col)
# 
# 
# # SID ~ Samples LM
# res.data.high = res.data[100:length(res.data[,1]),]
# plot(Mean.Distance ~ Samples, data = res.data.high, pch ="+", xlab = "Samples", ylab = "SID")
# SID.lm = lm(Mean.Distance ~ Samples, data = res.data.high)
# 
# # Diagnostic plots for SID~Samples Linear model
# par(mfrow=c(2,2))
# plot(SID.lm)
# par(mfrow=c(1,1))
# 
# # Polynomial Model?
# SID.poly.lm = lm(Mean.Distance ~ poly(Samples, 2), data = res.data.high)
# anova(SID.lm, SID.poly.lm)
# 
# # Diagnostic plots for SID~Samples Polynomial model
# par(mfrow=c(2,2))
# plot(SID.poly.lm)
# par(mfrow=c(1,1))
# 
# summary(SID.poly.lm)
# 
# # Draw quadratic through data
# # NOT WORKING???
# prd = data.frame(Samples = seq(from = range(res.data.high$Samples)[1], to = range(res.data.high$Samples)[2], length.out = length(res.data.high$Samples)))
# plot(Mean.Distance ~ Samples, data = res.data.high, pch ="+", xlab = "Samples", ylab = "SID")
# points(predict(SID.poly.lm, newdata=prd), type='l', col='red')
# 
# 
# 
# 
# # KL Divergence ~ Samples LM
# plot(KL.Divergence ~ Samples, data = res.data.high, pch ="+", xlab = "Samples", ylab = "Kullback–Leibler Divergence")
# kl.dist.lm = lm(`KL.Divergence` ~ Samples, data = res.data.high)
# abline(kl.dist.lm,lwd=2,col=2)
# 
# # Diagnostic plots for SID~Samples Linear model
# par(mfrow=c(2,2))
# plot(kl.dist.lm)
# par(mfrow=c(1,1))
# 
# summary(kl.dist.lm)
# 
# 
# 
# ggplot(res.data, aes(x=Samples, y=KL.Divergence, color=Mean.Distance)) + geom_point() + scale_color_gradient(low=low.col, high=high.col)
# 
# # Can you predict the sample size by how well it works?
# test.lm = lm(Samples ~ poly(Mean.Distance,2) * KL.Divergence, data = res.data.high)
# summary(test.lm)
# par(mfrow=c(2,2))
# plot(test.lm)
# par(mfrow=c(1,1))
# 
# 
# 
# 
# 
# 
# # # Print DAG for report
# # graphviz.plot(model2network("[S][R|S][G|S][W|R:G]"))
# # 
# # # Print interventional DAG for report
# graphviz.plot(model2network("[S][R|S][G][W|R:G]"), highlight=list(nodes=c("G"), fill="grey", col="black"))
# graphviz.plot(model2network("[S][R|S][G|S][W|R:G]"), highlight = list(arcs = c(from = "S", to = "G"), col="white"))
# 
# graphviz.plot(model2network("[S][R|S][G|S][W|R:G]"), highlight = list(nodes=c("G"), arcs = c(from = "S", to = "G"), col=c("white")))
# 
# 
# # Plot Markov blanket for chapter
# # This is the one from Scutari
# graphviz.plot(model2network("[A][B|A:E][C|A:D][D|A][E|D:F][F|G][G][H|F:I][I][J|H]"), highlight=list(nodes=c("A", "B", "D", "E", "F"), fill=c("grey90", "grey90", "grey90", "grey50", "grey90"), col=c("black", "black","black","black","black")), main="MB(E)")
# # This one is for the chapter
# graphviz.plot(model2network("[A][B|A][C|A][D|B:C][E|C][I|E][H|D:F][F][G|F][J|G:I][K|D:M][L|K:H][M|N:O][N|A][O|N]"), layout="fdp", highlight=list(nodes=c("B", "C", "D", "K", "H", "F", "M"), fill=c("grey90"), col=c("black")), main = "MB(D)")
# 
# 
# 
# # # Plot 2 node network for equivalence
# # graphviz.plot(model2network("[A][B|A]"), layout = "fdp")
# # graphviz.plot(model2network("[A|B][B]"))
# 
# 
# plot(res.data$Mean.Distance ~ res.data$Samples)
