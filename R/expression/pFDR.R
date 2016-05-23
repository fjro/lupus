attach(genes)
require(ggplot2)

res <- apply(genes[,-c(1,96)], 2, function(c) t.test(c)$statistic)
plot(res)

#ranks the individuals by gene expression level
#group is a vector of gene expression values
#level is the percentile to take e.g. 0.25 or 0.5
#upper is a boolean value: True for upper percentiles, False for lower percentiles
# Steps 1, 3 and 4
rank_individuals <- function(group, level, upper)
{
  count <- length(group)
  if (upper)
  {    
    start_index <- count - (count * level) +1
    sort(group)[start_index:count]
  }
  else
  {
    end_index <- count * level
    sort(group)[1:end_index]
  }
}

#calculate the maximum absolute t-statistic for each group (upper and lower 25%, 50%,75%)
#correspondes to Steps 2 and 5
max_t <- function(patients, controls)
{
  t_statistics <- rep(0,6)
  levels <- c(0.25,0.5,0.75)
  expression <- c(T,F) 
  index <- 1
  for (i in 1:length(expression))
  {
    for (j in 1:length(levels))
    {
      p <- rank_individuals(patients, levels[j], expression[i])
      c <- rank_individuals(controls, levels[j], expression[i])
      t_statistics[index] <- t.test(p,c)$statistic
      index <- index + 1
    }
  }
  max(abs(t_statistics))
}

#Resamples the patients and controls and calculates the maximum t-statistic
#data is the original data.frame
#gene is the name of the gene to test
# Step 6
resample_and_test <- function(data, gene)
{
  subjects <-sample(data[,gene],nrow(data))
  p=subjects[0:53]
  c=subjects[54:420]
  max_t(p, c)
}



#calculates a single p-value
#data is the original data frame
#patients is the original patient or treatment group
#controls is the original control group
#gene is the name of the gene to test
#n is the number of iterations
#Steps 7 and 8
p_value=function(data, patients, controls, gene, n)
{
  critical_p <- max_t(patients[,gene], controls[, gene])
  count <- 0
  for (i in 1:n)
  {
    if (resample_and_test(data, gene) > critical_p)
    {
      count <- count + 1
    }      
  }
  count/n
}

#calculates the p-values for each gene
#data is the original data.frame
#n is the number of iterations e.g. 10, 100, 1000 etc.
#Step 9
all_p_values <- function(data, n)
{
  genes <- colnames(data)[2:ncol(data)]
  patients <- data[diseased=="RED",]
  controls <- data[diseased=="GREEN",]
  p_values <- rep(0, length(genes))
  for(i in 1:length(genes))
  {
    time <- system.time(p_values[i] <- p_value(data, patients, controls, genes[i], n))
    cat(sprintf("Completed:%s in %s\"\n",  genes[i], round(time["elapsed"], 2)))
  }
  p_values
}


#p_values the p-values for each gene
#q the value to control the FDR procedure
#returns a list containing the Critical value and it's index
#Step 10
fdr <- function(p_values, q)
{
  p_values <- sort(p_values)
  p_critical <- 0
  index <- 0
  for(i in length(p_values):1)
  {
    p_critical <- (i/length(p_values))*q
    if(p_values[i] <= p_critical)
    {
      index <- i
      break
    }
    else
    {
      p_critical <-0
    }
  }
  list(Critical = p_critical, Index=index)
}

#pcer: the per comparison error rate
pcer <- function(x, genes)
{
  patients <- data[diseased=="RED",]
  controls <- data[diseased=="GREEN",]
  p <- rep(0, length(genes))
  
  for (i in 1:length(genes))
  {
    p[i] <- t.test(patients[,genes[i]],controls[,genes[i]])$p.value
  }
  p
}

genes[which(pcer(data,genes)<=0.05)] 

nrow(genes)

#calculate the p-values, find the critical value using FDF and then show the genes where the H0 is rejected
#10 iterations
n10 <- all_p_values(genes[,-96], 10)

p10 <- fdr(n10, 0.05)
#find the hypotheses to reject
geneNames <- colnames(genes)[-c(1,96)]
geneNames[which(n10 < p10$Critical)]
barplot(n10)
barplot(n10,names=geneNames,horiz=T)
require(ggplot2)
results = data.frame(cbind(data.frame(n10),geneNames))
typeof(results)
ggplot(results, aes(y = n10, x = geneNames, order=n10)) + geom_bar(stat="identity") +coord_flip()
n10
n100 <- all_p_values(genes[,-96], 100)
p100 <- fdr(n100, 0.05)
results = cbind(data.frame(n10),geneNames)
ggplot(results, aes(y = n100, x = geneNames)) + geom_bar(stat="identity") +coord_flip()

n1000 <- all_p_values(genes[,-96], 1000)
p1000 <- fdr(n1000, 0.01)
genes[which(n1000 < p1000$Critical)]
#X1038_s_at"

#careful, takes about 20 minutes
n10000 <- all_p_values(data[,-96], 10000)
p10000 <- fdr(n10000, 0.001)
genes[which(n10000 < p10000$Critical)]
#"X1038_s_at" "X1039_s_at"

n100000 <- all_p_values(genes[,-96], 100000)
p100000 <- fdr(n100000[,-1], 0.10)
p100000$Index
length(fdr(n100000[,1:2],0.1)$n100000)

geneNames[which(n100000$x < p100000$Critical)]
#"X1038_s_at" "X1039_s_at"
write.csv(n100000,"n100000.csv")
n100000 = read.csv("n100000.csv")
results = cbind(data.frame(n100000),geneNames)
colnames(results) = c('X','pvalue','gene')
ggplot(results, aes(y = pvalue, x = gene)) + geom_bar(stat="identity") 
ggplot(results, aes(y=pvalue,x=gene))+ scale_x_discrete(limits = results$gene)+ geom_bar(stat = "identity",position="dodge")+ theme(axis.text.x=element_text(angle=270, vjust = 0.5))


# use Storey's pFDR library
#http://www.inside-r.org/packages/cran/WGCNA/docs/qvalue
# to install uncomment the lines below, the library is not on CRAN
#source("http://bioconductor.org/biocLite.R") 
#biocLite("WGCNA") 
#install.packages("impute")
#install.packages("WGCNA")
require(WGCNA)

#read in the p-values from the permutation tests
pvalues <- read.csv("pvalues.csv")

#pcer
genes[which(pvalues$n1000000 <=0.05)] 

#fwer
genes[which(pvalues$n1000000 <=0.05/length(genes))] 

#test using pFDR at q=0.05
genes[which(qvalue(pvalues$n1000000)$qvalues<=0.05)]
#"X1038_s_at"

#test using pFDR at q=0.10
genes[which(qvalue(pvalues$n1000000)$qvalues<=0.10)]
#"X1038_s_at"