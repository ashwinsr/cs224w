library(fastSOM)
library(dplyr)
library(reshape)
library(frequencyConnectedness)
library(BigVAR)
library(foreign)
library(readxl)
library(lubridate)
library(zoo)
library(scales)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(viridis)

## Testing GVD -- swap ordering of the variables and see if it's the same

## Calculate inter-minute lagged difference in "close" price
dat <- read.csv("~/Downloads/^VIX-equity-1-minute-20160210.csv")
dat <- dat[,colnames(dat) %in% c("quote_datetime", "close")]

## Number of minutes measured over in a day, for random generation of data
## later on
nmin <- nrow(dat)
bw_minute_changes <- dat %>% mutate(Delta = close - lag(close))

## Remove first two minutes' observations, because the close is 
## instantiated to 0 (why?) which yields an NA difference and a 
## difference equal to the entire current stock price
bw_minute_changes <- bw_minute_changes[3:nrow(bw_minute_changes),]
hist(bw_minute_changes$Delta)

## Looks sufficiently normal such that we'll treat the data generating process
## as Gaussian for the purpose of constructing a test set
mu <- mean(bw_minute_changes$Delta)
sigma <- sqrt(var(bw_minute_changes$Delta))

ndays <- 1000
impl_volatility_1 <- sapply(1:ndays, function(x) {
  intra_day <- abs(rnorm(nmin, mu, sigma))
  daily_vol <- sum(log(intra_day)^2)
  return(daily_vol)
}) %>% unlist()

## Load generic data from BigVAR package and test dependence of variable
## ordering
data("volatilities")
## Add a small bias to allow us to take logs of the volatilities 
volatilities <- (volatilities) %>% na.omit()
volatilities <- volatilities + 1

volatilities1 <- volatilities[,1:5] 
volatilities2 <- volatilities[,c(5, 1:4)]

bv_est <- function(data) {
  model <- constructModel(as.matrix(data), p = 12, struct = "Basic", gran = c(50, 50), verbose = FALSE)
  results <- cv.BigVAR(model)
  return(results)
}

## test1 <- bv_est(log(volatilities1))
## test2 <- bv_est(log(volatilities2))

test1 <- VAR(volatilities1, p = 12, type = "const")
test2 <- VAR(volatilities2, p = 12, type = "const")

genFEVD(test1, n.ahead = 100, no.corr = FALSE)
genFEVD(test2, n.ahead = 100, no.corr = FALSE)

####################################################################################################
####################################################################################################

## Arguments are a data.frame of intra-day tick data for an arbitrary number of assets where the
## following conditions hold: 
## 
## (1) The first column of the data is titled "Date" and is of type date() (converted using 
##     as.Date()). 
## (2) There is an equal number of observations for each asset, with each time series referring
##     to the same set of time periods. 
## (3) Intra-day ticks are observed using the same periods (1 min., 5 min., etc.) across 
##     days. 
##
## And a parameter p of "look-back" for the vector auto-regression. 

constructSpilloverMatrix <- function(dat, p) {
  firstDifferences <- sapply(2:ncol(dat), function(x) {
    name <- colnames(dat)[x]
    deltaName <- paste("Delta_", name, sep = "")
    diff <- dat %>% group_by(Date) %>% mutate(!!deltaName := !!as.name(name) - lag(!!as.name(name)))
    return (diff[,ncol(diff)])
  }) %>% data.frame()
  firstDifferences$Date <- dat$Date
  firstDifferences <- firstDifferences[,c(ncol(firstDifferences), 1:(ncol(firstDifferences)-1))]
  firstDifferences <- firstDifferences %>% na.omit()
  
  ## Add 1 to each value to eliminate problems with zero change. Is this cool, given that first 
  ## differences are additive and not ratios? 
  constant <- abs(min(firstDifferences[,2:ncol(firstDifferences)])) + 0.1
  firstDifferences[,2:ncol(firstDifferences)] <- log(abs(firstDifferences[,2:ncol(firstDifferences)] + constant))^2
  dailyRealizedVol <- aggregate(. ~ Date, firstDifferences, sum)
  rownames(dailyRealizedVol) <- dailyRealizedVol$Date
  dailyRealizedVol$Date <- NULL
  dat <- zoo(dailyRealizedVol, order.by = rownames(dailyRealizedVol))
  
  ## Arbitrarily set p <- 7 (b/c VAR converges) -- probably has to do with # of
  ## free parameters one needs to estimate given sample of ~97 days... 
  est <- VAR(dailyRealizedVol, p, type = "const")
  fevd <- genFEVD(est, n.ahead = 100, no.corr = FALSE) %>% data.frame()
  spilloverMatrix <- fevd*100
  
  colnames(spilloverMatrix) <- str_replace_all(colnames(dailyRealizedVol), "Delta_", "")
  rownames(spilloverMatrix) <- NULL
  
  write.csv(spilloverMatrix, "~/Downloads/Bank_Network_GVD.csv")
  return (spilloverMatrix)
}

####################################################################################################
####################################################################################################

dat <- read_xlsx("~/Downloads/CS224W_140Days_JPM-WFC-BAC-C-USB-BK-PNC-AXP-GS-MS-FNM-FRE-AIG-2.xlsx", 
                 sheet = "Assets")
dat$Date <- as.Date(dat$Date)
spilloverMatrix <- constructSpilloverMatrix(dat, p = 1)

####################################################################################################
####################################################################################################

dat <- read_xlsx("~/Downloads/CS224_CDS_2006-2018.xlsx", sheet = "Processed")
dat$Date <- as.Date(dat$Date)

## Fill values to interpolate missing data from nearest neighbors, first going 
## "down" (preserving the distantness of known data) and then up
dat <- dat %>% fill(Venezuela:Vietnam, .direction = "down")
dat <- dat %>% fill(Venezuela:Vietnam, .direction = "up")
dat <- data.frame(dat)
rownames(dat) <- dat$Date

## Drop non-numeric columns (again -- no data...)
dat <- dat[,!(colnames(dat) %in% c("Iceland", "Morocco"))]

dat$Date <- NULL
dat <- zoo(dat, order.by = rownames(dat))

## Estimate for different horizon parameters up to p = 14 
for (p in 1:14) {
  est <- VAR(dat, p, type = "const")
  fevd <- genFEVD(est, n.ahead = 100, no.corr = FALSE) %>% data.frame()
  spilloverMatrix <- fevd*100
  
  colnames(spilloverMatrix) <- colnames(dat)
  rownames(spilloverMatrix) <- NULL
  
  filename <- paste("~/Downloads/CDS_Network_GVD_p", p, ".csv", sep = "")
  write.csv(spilloverMatrix, filename, row.names = FALSE)
}

write.csv(dat, "~/Downloads/CS224_CDS_2006-2018_filled.csv")

####################################################################################################
####################################################################################################

dat <- read_xlsx("~/Downloads/CS 224W CDS Data/CS224_CDS_2006-2018.xlsx", 
                 sheet = "Processed")
dat$Date <- as.Date(dat$Date)

## Fill values to interpolate missing data from nearest neighbors, first going 
## "down" (preserving the distantness of known data) and then up
dat <- dat %>% fill(Venezuela:Vietnam, .direction = "down")
dat <- dat %>% fill(Venezuela:Vietnam, .direction = "up")
dat <- data.frame(dat)
dat_sd <- dat %>% group_by(week = cut(Date, "week")) %>% 
  summarise_all(funs(sd))
dat_sd$Date <- NULL
dat_sd <- data.frame(dat_sd)
rownames(dat_sd) <- dat_sd$week
dat_sd$week <- NULL
dat_sd <- dat_sd[,!(colnames(dat_sd) %in% c("Iceland", "Morocco"))]
dat_sd <- dat_sd %>% na.omit()

dat_sd <- zoo(dat_sd, order.by = rownames(dat_sd))

## Estimate for different horizon parameters up to p = 14 
for (p in 1:14) {
  est <- VAR(dat_sd, p, type = "const")
  fevd <- genFEVD(est, n.ahead = 100, no.corr = FALSE) %>% data.frame()
  spilloverMatrix <- fevd*100
  
  colnames(spilloverMatrix) <- colnames(dat_sd)
  rownames(spilloverMatrix) <- NULL
  
  filename <- paste("~/Downloads/WeeklySTD_CDS_Network_GVD_p", p, ".csv", sep = "")
  write.csv(spilloverMatrix, filename, row.names = FALSE)
}

write.csv(dat_sd, "~/Downloads/CS224_CDS_2006-2018_filled_weeklySTD.csv")

####################################################################################################
####################################################################################################

scores <- read.csv("~/Downloads/CS224W_BumpChart_p.csv")
rownames(scores) <- scores$Country
scores$Country <- NULL

## Arguments are a data.frame of PageRank scores for each entity and a filename
## representing the save destination of the resultant image (with no extension, 
## e.g. simply "testImage").
## The scores data object should have rownames corresponding to each entity and
## column names corresponding to each value of "p" (just a number, or a name of
## the form X1, X2, etc.). 
pageRankBumpChart <- function(scores, filename) {
  ranks <- scores %>% mutate_all(funs(dense_rank(desc(.))))
  ranks$Country <- rownames(scores)
  ranks <- melt(ranks, id = c("Country"))
  ranks$variable <- str_replace_all(ranks$variable, "X", "") %>% as.character() %>% as.numeric()
  colnames(ranks) <- c("Country", "PValue", "Rank")
  ranks$Label <- ifelse(
    ((ranks$PValue == min(ranks$PValue)) | (ranks$PValue == max(ranks$PValue))), 
    paste(ranks$Country, " (", ranks$Rank, ")", sep = ""), NA)
  
  marginWidth <- 1
  breaks <- sort(unique(ranks$PValue))
  labels <- as.character(sort(unique(ranks$PValue)))
  
  ranks$Rank <- ranks$Rank * 2
  
  plt <- ggplot(ranks, aes(x = PValue, y = Rank)) + 
    geom_line(aes(color = Country)) + 
    geom_point(shape = 21, stroke = 1, fill = "white", aes(color = Country)) + 
    geom_label_repel(data = ranks[ranks$PValue == min(ranks$PValue),], 
                     aes(label = Label), size = 2, nudge_x = -1, direction = "x",
                     segment.size = 0, label.size = NA, box.padding = 0.1, force = 0) + 
    geom_label_repel(data = ranks[ranks$PValue == max(ranks$PValue),], 
                     aes(label = Label), size = 2, nudge_x = +1, direction = "x",
                     segment.size = 0, label.size = NA, box.padding = 0.1, force = 0) + 
    theme_minimal() + 
    theme(axis.title.x = element_text(face = "bold"), 
          axis.title.y = element_text(face = "bold"), 
          panel.grid.major.y = element_blank(), 
          panel.grid.major.x = element_line(size = 0.5),
          panel.grid.minor = element_blank(), 
          legend.position = "none", 
          axis.text.y = element_blank(), 
          title = element_text(face = "bold")) + 
    scale_x_continuous(lim = c(min(ranks$PValue)-marginWidth, max(ranks$PValue)+marginWidth), 
                       breaks = breaks,
                       labels = labels) + 
    xlab("Backwards-looking horizon; p in VAR(p)") + 
    ylab("PageRank rank output; importance rank") + 
    ggtitle("Sensitivity in PageRank rank output to VAR(p) horizon")
  save <- paste("~/Downloads/CS224W/final project/", filename, ".png", sep = "")
  ggsave(save, plt, device = "png")
  return(plt)
}

####################################################################################################
####################################################################################################

dat <- read_xlsx("~/Downloads/CS 224W CDS Data/CS224_CDS_2006-2018.xlsx", 
                 sheet = "Processed")
dat$Date <- as.Date(dat$Date)
not_na <- apply(dat, 2, function(x) {length(which(!is.na(x)))})
names(not_na) <- colnames(dat)
hist(not_na, breaks = 20)

for (i in c(1500, 2000, 2250, 2500, 2750, 3000)) {
  ncountries <- length(not_na[not_na > i]) - 1
  msg <- paste("# countries with >= ", i, " samples: ", ncountries, sep = "")
  print(msg)
}

# [1] "# countries with >= 1500 samples: 48"
# [1] "# countries with >= 2000 samples: 38"
# [1] "# countries with >= 2250 samples: 32"
# [1] "# countries with >= 2500 samples: 30"
# [1] "# countries with >= 2750 samples: 23"
# [1] "# countries with >= 3000 samples: 14"

## Let threshold be >= 2500 samples initially -- arbitrary but could --> more meaningful results

threshold <- 2500
keep <- names(not_na)[which(not_na >= threshold)]
not_na <- not_na[keep]
not_na <- sort(not_na)
not_na <- not_na[2:length(not_na)]

## Keep Greece, because it's consequential
keep <- names(not_na)
keep <- c(keep, "Greece")

dat <- dat[,colnames(dat) %in% keep]
ordering <- lapply(keep, function(x){
  return(match(x, colnames(dat)))
}) %>% unlist()
dat <- dat[,ordering]

dat <- dat %>% fill(Romania:Greece, .direction = "down")
dat <- dat %>% fill(Romania:Greece, .direction = "up")
dat <- data.frame(dat)
dat_sd <- dat %>% group_by(week = cut(Date, "week")) %>% 
  summarise_all(funs(sd))
dat_sd$Date <- NULL
dat_sd <- data.frame(dat_sd)
rownames(dat_sd) <- dat_sd$week
dat_sd$week <- NULL
dat_sd <- dat_sd[,!(colnames(dat_sd) %in% c("Iceland", "Morocco"))]
dat_sd <- dat_sd %>% na.omit()

dat_sd <- zoo(dat_sd, order.by = rownames(dat_sd))

## Estimate for different horizon parameters up to p = 14 
for (p in 1:14) {
  est <- VAR(dat_sd, p, type = "const")
  fevd <- genFEVD(est, n.ahead = 100, no.corr = FALSE) %>% data.frame()
  spilloverMatrix <- fevd*100
  
  colnames(spilloverMatrix) <- colnames(dat_sd)
  rownames(spilloverMatrix) <- NULL
  
  filename <- paste("~/Downloads/WeeklySTD_CDS_Network_GVD_p", p, "_keep2500_Greece.csv", sep = "")
  write.csv(spilloverMatrix, filename, row.names = FALSE)
}

write.csv(dat_sd, "~/Downloads/CS224_CDS_2006-2018_filled_weeklySTD_keep2500.csv")
