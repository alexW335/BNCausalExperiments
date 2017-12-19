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


# This calculates the SID-Divergence from the
# probability distribution associated with the original BN to the 
# one reconstructed from the data. It also calculates the 
# SID between the original struture and the reconstructed one.
# These are calculated "times.to.repeat" times, and then averaged,
# where a fresh BN is reconstructed each time on the same number of
# samples.
# It also finds (exactly) the probabilities P(D|A) and P(D|~A), 
# however only once as it can be quite expensive to do.
SIDMean = function(bn,
                  num.samples.to.gen,
                  times.to.repeat = 1,
                  test.funct = "x2",
                  score.funct = "bde",
                  param.learn.method = "mle"){
  # Holds running sum of graph edit distances
  sumSID = 0
  sup = 0
  sdown = 0

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
    sidout = structIntervDist(bn.am, lrnd.am)
    
    sumSID = sumSID + sidout$sid
    sup    = sup    + sidout$sidUpperBound
    sdown  = sdown  + sidout$sidLowerBound
    
  }
  
  return(c(sumSID, sdown, sup)/times.to.repeat)
}

# Specify the network structure
dag.test = model2network("[A][B][C|A:B][D|B][E|B][G|E][H|C:D]")
graphviz.plot(dag.test)


# Label levels
A.lv = c("a", "-a")
B.lv = c("b", "-b")
C.lv = c("c", "-c")
D.lv = c("d", "-d")
E.lv = c("e", "-e")

G.lv = c("g", "-g")
H.lv = c("h", "-h")

# Define probabilities
A      = 0.3
B      = 0.7
EgB    = 0.2
EgnB   = 0.8
GgE    = 0.4
GgnE   = 0.6
DgB    = 0.3
DgnB   = 0.1
CgAB   = 0.1
CgAnB  = 0.7
CgnAB  = 0.8
CgnAnB = 0.4
HgCD   = 0.2
HgCnD  = 0.9
HgnCD  = 0.7
HgnCnD = 0.5


# Assign to arrays
A.prob = array(c(A, 1 - A), dim = 2, dimnames = list(A = A.lv))
B.prob = array(c(B, 1 - B), dim = 2, dimnames = list(B = B.lv))
E.prob = array(c(EgB, 1-EgB, EgnB, 1- EgnB), dim=c(2,2), dimnames = list(E = E.lv, B = B.lv))
G.prob = array(c(GgE, 1-GgE, GgnE, 1- GgnE), dim=c(2,2), dimnames = list(G = G.lv, E = E.lv))
D.prob = array(c(DgB, 1-DgB, DgnB, 1- DgnB), dim=c(2,2), dimnames = list(D = D.lv, B = B.lv))
C.prob = array(c(CgAB, 1-CgAB, CgnAB, 1-CgnAB, CgAnB, 1-CgAnB, CgnAnB, 1-CgnAnB), 
               dim=c(2,2,2), dimnames = list(C = C.lv, A = A.lv, B = B.lv))
H.prob = array(c(HgCD, 1-HgCD, HgnCD, 1-HgnCD, HgCnD, 1-HgCnD, HgnCnD, 1-HgnCnD), 
               dim=c(2,2,2), dimnames = list(H = H.lv, C = C.lv, D = D.lv))

# Create distribution & fit to BN
loc.dist = list(A = A.prob, B = B.prob, C = C.prob, D = D.prob, E = E.prob, G = G.prob, H = H.prob)
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
prB = bn.fit.barchart(bn.actual$B, main = "P(B)")
prC = bn.fit.barchart(bn.actual$C, main = "P(C|A,B)")
prD = bn.fit.barchart(bn.actual$D, main = "P(D|B)")
prE = bn.fit.barchart(bn.actual$E, main = "P(E|B)")
prG = bn.fit.barchart(bn.actual$G, main = "P(G|E)")
prH = bn.fit.barchart(bn.actual$H, main = "P(H|C,D)")

# This where the data generation takes place. It takes a long time, 
# but saves the results in a text file.

# It goes through from i = min to max.number.of.observations, calling SIDMean each time and 
# instructing it to generate "i" samples, and average the values out over 'repeat.times' times.
# This way we get to see how a BN generated from 1,2,3,4,...,max.number.of.observations samples
# compares to the original BN.

repeat.times = 10
min.number.of.observations = 1
max.num.observations = 1000

res.data = data.frame(seq(from = min.number.of.observations, to = max.num.observations), vector(mode = "numeric", length = max.num.observations - min.number.of.observations + 1),
                      vector(mode = "numeric", length = max.num.observations - min.number.of.observations + 1), vector(mode = "numeric", length = max.num.observations - min.number.of.observations + 1))
colnames(res.data) = c("Samples", "Mean SID", "Lower", "Upper")

c=1
for (i in seq(from = min.number.of.observations, to = max.num.observations)){
  print(i)
  res.data[c, 1] = i
  s = SIDMean(bn.actual, i, times.to.repeat = repeat.times)
  res.data[c, 2] = s[1]
  res.data[c, 3] = s[2]
  res.data[c, 4] = s[3]
  
  c = c + 1
}


# Run below line if generating new data to stop log-problems. 
# Will lose some low-sample data points.
res.data = res.data[is.finite(rowSums(res.data)),]

plot(res.data[,2] ~ res.data[,1], pch ="+", xlab = "Samples", ylab = "SID")

# ggplot(res.data, aes('Samples', 'Mean SID')) + geom_ribbon(aes(ymin = res.data[,3], ymax = res.data[,4]), fill = "grey70") + geom_line()




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
# ggplot(res.data, aes(x=Samples, y=Mean.Distance, color=log(SID.Divergence))) + geom_point() + scale_color_gradient(low=low.col, high=high.col)
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
# # SID Divergence ~ Samples LM
# plot(SID.Divergence ~ Samples, data = res.data.high, pch ="+", xlab = "Samples", ylab = "Kullback–Leibler Divergence")
# SID.dist.lm = lm(`SID.Divergence` ~ Samples, data = res.data.high)
# abline(SID.dist.lm,lwd=2,col=2)
# 
# # Diagnostic plots for SID~Samples Linear model
# par(mfrow=c(2,2))
# plot(SID.dist.lm)
# par(mfrow=c(1,1))
# 
# summary(SID.dist.lm)
# 
# 
# 
# ggplot(res.data, aes(x=Samples, y=SID.Divergence, color=Mean.Distance)) + geom_point() + scale_color_gradient(low=low.col, high=high.col)
# 
# # Can you predict the sample size by how well it works?
# test.lm = lm(Samples ~ poly(Mean.Distance,2) * SID.Divergence, data = res.data.high)
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
