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
library(parallel)
library(rlist)
library(gam)
library(stats)


# install.packages("SID", lib = "H:/My Documents/libs")
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rgraphviz", lib = "H:/My Documents/libs")
# biocLite("RBGL", lib = "H:/My Documents/libs")

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
SIDMean = function(bn,
                   num.samples.to.gen,
                   times.to.repeat = 1,
                   test.funct = "x2",
                   score.funct = "bde",
                   param.learn.method = "mle"){
    # Holds running sum of graph edit distances
    sdev = c()
    
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
        
        sdev[i] = sidout$sid
        
    }
    return(c(mean(sdev), sd(sdev)))
}

meanGraph = function(bn,
                   num.samples.to.gen,
                   times.to.repeat = 1,
                   test.funct = "x2",
                   score.funct = "bde",
                   param.learn.method = "mle"){

    # Randomly generate 'num.samples.to.gen' samples from the original BN.
    gen.samples = rbn(bn, num.samples.to.gen)
    
    bn.am = as(bnlearn::as.graphNEL(bn), "matrix")
    
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
 
    return(lrnd.am)
}


# Specify the network structure
dag.test = model2network("[A][B][C|A:B][D|B][E|B][G|E][H|C:D]")
# graphviz.plot(dag.test)

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

# # Get rid of the heinous default colours when 
# # displaying the probability distributions for each node
# myColours <- brewer.pal(6,"Greys")
# my.settings <- list(
#     superpose.polygon=list(col=myColours[1], border="transparent"),
#     strip.background=list(col=myColours[1]),
#     plot.polygon=list(col=myColours[2]),
#     strip.border=list(col="black")
# )
# trellis.par.set(my.settings)
# 
# # Display the probability distributions of each node with barcharts.
# prA = bn.fit.barchart(bn.actual$A, main = "P(A)")
# prB = bn.fit.barchart(bn.actual$B, main = "P(B)")
# prC = bn.fit.barchart(bn.actual$C, main = "P(C|A,B)")
# prD = bn.fit.barchart(bn.actual$D, main = "P(D|B)")
# prE = bn.fit.barchart(bn.actual$E, main = "P(E|B)")
# prG = bn.fit.barchart(bn.actual$G, main = "P(G|E)")
# prH = bn.fit.barchart(bn.actual$H, main = "P(H|C,D)")


# It goes through from i = min to max.number.of.observations, calling SIDMean each time and 
# instructing it to generate "i" samples, and average the values out over 'repeat.times' times.
# This way we get to see how a BN generated from 1,2,3,4,...,max.number.of.observations samples
# compares to the original BN.
repeat.times = 50
min.number.of.observations = 1
max.num.observations = 2500
step.size = 25


res.data = data.frame(seq(from = min.number.of.observations, to = max.num.observations), vector(mode = "numeric", length = max.num.observations - min.number.of.observations + 1),
                      vector(mode = "numeric", length = max.num.observations - min.number.of.observations + 1), vector(mode = "numeric", length = max.num.observations - min.number.of.observations + 1))
colnames(res.data) = c("Samples", "Mean SID", "Lower", "Upper")

sfunct = function(i){
    v = SIDMean(bn.actual, i, repeat.times)
    return(c(v[1], v[2]))
}

storeda = seq(from = min.number.of.observations, to = max.num.observations, by=step.size)

# Calculate the number of cores
detectCores()
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, list("SIDMean", "bn.actual", "repeat.times", "rbn", "rsmax2",
                       "bn.fit", "structIntervDist"))

res.data = parLapply(cl, storeda, fun=sfunct)

sdevs = lapply(seq(1,length(storeda)), function(i) res.data[[i]][2])
sdevs = unlist(sdevs)
res.data = lapply(seq(length(storeda)), function(i) res.data[[i]][1])
res.data = unlist(res.data)
stopCluster(cl)

# write.table(pts, 'datagen.txt')
# pts = read.table("datagen.txt")


# res.data = res.data[is.finite(rowSums(res.data))]

pts = data.frame(Samples=storeda, SID=res.data, SDev=sdevs)

plot(SID ~ Samples, pch =".", xlab = "Samples", ylab = "SID", main = "SID ~ Samples", data = pts)
plot(log(SID) ~ Samples, pch ="+", xlab = "Samples", ylab = "log(SID)", main = "log(SID) ~ Samples 500:2000", data = pts)

tail(pts)


h = ggplot(pts, aes(x=Samples, y=SID)) + geom_line() + geom_ribbon(aes(ymin = SID - 1.96*SDev, ymax = SID + 1.96*SDev, x = Samples, fill="band"), alpha = 0.3) 
print(h)    



#################### AVERAGE GRAPH ####################
num.to.average = 1000
num.samples = 2505

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl, list("meanGraph", "rbn", "rsmax2", "bn.fit", "bn.actual", "num.samples"))
graphs = parLapply(cl, seq(1, num.to.average), fun=function(i) meanGraph(bn.actual, num.samples))
stopCluster(cl)

g = Reduce('+', graphs)/num.to.average
g
net = graph_from_adjacency_matrix(g, mode="directed",weighted=TRUE)

plot(net,vertex.label=V(net)$name, 
     edge.color=rgb(E(net)$weight, 0, 1-E(net)$weight, (E(net)$weight + 1)/2),
     arrow.color=rgb(E(net)$weight, 0, 1-E(net)$weight, (E(net)$weight + 1)/2), 
     edge.label=NA, edge.width=3, vertex.color="white", vertex.size=25, autocurve.edges=T,
     label.color="black", layout = matrix(c(0,2,1,2,3,2,1,2,2,1,1,1,0,0), nrow=7, ncol=2))



################## MODEL FITTING ################## 
tlm = lm(log(SID) ~ Samples, data = pts[500:2000,])
plm = lm(log(SID) ~ poly(Samples, 5), data = pts[500:2000,])
anova(tlm, plm)

summary(plm)

par(mfrow=c(2,2))
plot(plm)
par(mfrow=c(1,1))

# Try a GAM
gm = gam(SID ~ s(Samples), data = pts)
gm$coefficients
summary(gm)

plot(SID ~ Samples, pch =".", xlab = "Samples", ylab = "SID", main = "SID ~ Samples", data = pts)
t = predict.gam(gm, newdata=data.frame(Samples = pts$Samples))
lines(t ~ pts$Samples, col='red', type='l')

par(mfrow=c(2,2))
gam.check(gm)
par(mfrow=c(1,1))

