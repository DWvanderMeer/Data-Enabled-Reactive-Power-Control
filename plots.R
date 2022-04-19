# Author: Dennis van der Meer
# E-mail: dennis.van_der_meer[at]minesparis.psl.eu

# Create plots for the paper
# IEEE width limit 88.9 millimeters single column or 182 millimeters double column

library(ggplot2)
library(gtools)
library(ggridges)
library(dplyr)
library(ggfan)
plot.size = 8; line.size = 0.15; point.size = 0.6
ts <- seq(as.POSIXct("2009-01-01 00:00"), as.POSIXct("2009-01-01 23:59"),by = "5 min")

WORKING_DIR <- "~/Google Drive/My Drive/PhD-Thesis/My-papers/CopulaEDAs/Data-Enabled-Reactive-Power-Control/"
source(file.path(WORKING_DIR,"functions.R"))

################################################################################
# Figure 2
################################################################################
source(file.path(WORKING_DIR,"StartMatlabAndLoadData.R"))
PV_CS <- getVariable(matlab,'PV_CS')
PV_CS <- data.frame(Power=PV_CS$PV.CS,timeStep=ts)
GHIM_int <- getVariable(matlab, 'GHIM_int')
GHIM_int <- data.frame(GHIM_int$GHIM.int,timeStep=ts)
Load_P <- getVariable(matlab,'Load_P')
Load_P <- data.frame(Power=Load_P$Load.P,timeStep=ts)
acf_PVcs <- data.frame(ACF=acf(PV_CS$Power,plot = F)$acf[,1,1],Lag=seq(1,25,1))
acf_loadP <- data.frame(ACF=acf(Load_P$Power, plot = F)$acf[,1,1],Lag=seq(1,25,1))
acf_loadP$Profile <- "Load"
acf_PVcs$Profile <- "PV clear sky"
acf_PVint <- acf_PVint <- lapply(GHIM_int[,1:10],acf,plot=F)
lst = list()
for(i in 1:length(acf_PVint)){lst[[i]] <- acf_PVint[[i]]$acf[,1,1]}
acfs <- as.data.frame(do.call(cbind,lst))
colnames(acfs)[1:10] <- c(18, 48, 56, 66, 79, 83, 95, 250, 300, 450)
acfs$lag <- seq(1,nrow(acfs),1)
close(matlab)
Load_P$Profile <- "Load"
PV_CS$Profile <- "PV clear sky"
colnames(GHIM_int)[1:10] <- c(18, 48, 56, 66, 79, 83, 95, 250, 300, 450)
mydf <- GHIM_int %>%
  tidyr::pivot_longer(cols = -timeStep)

p1 <- ggplot() +
  geom_line(data = mydf, aes(x=timeStep,y=value,group=name,linetype = "PV cloudy"), size=line.size,alpha=0.8) + #
  geom_line(data = Load_P, aes(x=timeStep,y=Power, linetype = "Load"), size = line.size) +
  ylab("Power (-)") +
  xlab("Time (5 min.)") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.4,0,0), "lines"),
        text = element_text(family = "Times"),
        axis.text=element_text(family = "Times", size=plot.size),
        axis.title=element_text(family = "Times", size=plot.size),
        plot.title=element_text(family = "Times", size=plot.size),
        legend.title = element_text(family = "Times", size=plot.size),
        legend.position = "bottom",
        legend.text = element_text(family = "Times", size=plot.size),
        legend.box.margin = ggplot2::margin(0,-4,0,-10),
        legend.margin = ggplot2::margin(0,0,0,0),
        strip.text = element_text(family = "Times", size=plot.size),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_linetype_manual(name="",
                        values = c("Load"="dotted", "PV cloudy"="solid")) + # , "PV clear sky"="dashed"
  scale_x_datetime(date_breaks = "6 hour",
                   date_labels = "%H:%M") +
  scale_y_continuous(limits = c(0,1))

mydf_acfs <- acfs %>%
  tidyr::pivot_longer(cols = -lag)
p2 <- ggplot()+
  geom_line(data = mydf_acfs, aes(x=5*lag,y=value, group=name, linetype="PV cloudy"), size=line.size) +
  geom_line(data = acf_loadP, aes(x=5*Lag,y=ACF, linetype = "Load"), size = line.size) +
  scale_linetype_manual(name="",
                        values = c("Load"="dotted", "PV cloudy"="solid")) + # , "PV clear sky"="dashed"
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.4,0,0), "lines"),
        text = element_text(family = "Times"),
        axis.text=element_text(family = "Times", size=plot.size),
        axis.title=element_text(family = "Times", size=plot.size),
        plot.title=element_text(family = "Times", size=plot.size),
        legend.title = element_text(family = "Times", size=plot.size),
        legend.position = "bottom",
        legend.text = element_text(family = "Times", size=plot.size),
        legend.box.margin = ggplot2::margin(0,-4,0,-10),
        legend.margin = ggplot2::margin(0,0,0,0),
        strip.text = element_text(family = "Times", size=plot.size),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines")) +
  xlab("Lag (min.)") +
  ylab("ACF (-)")

legend <- cowplot::get_legend(p2)
p1 <- p1 + theme(legend.position="none")
p2 <- p2 + theme(legend.position="none")

p <- cowplot::plot_grid(p1,p2, align = "v", nrow = 2, rel_heights = c(1,1)) # 9.5,4,4,4
my.plot <- gridExtra::grid.arrange(p, legend, nrow = 2, ncol = 1, 
                        layout_matrix = rbind(1,2),heights = c(2.5, 0.3)) # widths = 2.7, 
ggsave(filename = "~/Desktop/ACF.pdf", plot = my.plot, 
       device = cairo_pdf, units = "cm", height = 5.5, width = 8.5)
################################################################################
# Figure 3
################################################################################
# Results are saved in a directory
# From: https://stackoverflow.com/questions/52490552/r-convert-nan-to-na
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

dirs <- list.dirs("~/Desktop/case2")
dirs <- dirs[nchar(dirs)<59 & nchar(dirs)>45] # This only works when the files are in dirs
dirs
ts <- seq(as.POSIXct("2009-01-01 00:00"), as.POSIXct("2009-01-01 23:59"),by = "5 min")

iter_res <- list()
j <- 1
for(dir in dirs){
  setwd(dir)
  if(file.exists("Results.RData")){load("Results.RData")}
  lst <- list()
  tmp_case <- strsplit(dir,"/")[[1]]
  case <- tmp_case[length(tmp_case)]
  # print(case)
  tmp_model <- strsplit(case, split="[[:upper:]]")[[1]]
  model <- tmp_model[1]
  # print(model)
  tmp_run <- strsplit(case, split="[[:lower:]]")[[1]]
  run <- tmp_run[length(tmp_run)]
  # print(run)
  if(startsWith(case,"ga") | startsWith(case,"Unity") | startsWith(case,"Volt")){ # 
    next
  } else {
    busVoltage <- read.table("busVoltage.txt")
    busVoltage$Run <- run
    busVoltage$Model <- model
  }
  # if(startsWith(case,"Unity") | startsWith(case,"Volt")){
  #   # busVoltage$Run <- case
  #   # busVoltage$Model <- case
  #   next
  # }
  busVoltage$Time <- ts
  iter_res[[j]] <- busVoltage
  j <- j+1
}
iter_res <- do.call(rbind,iter_res)
mydf <- tidyr::pivot_longer(iter_res, cols = -c("Run","Model","Time"))
busVoltage_VV <- read.table("~/Desktop/case2/voltVar/busVoltage.txt")
busVoltage_VV <- as.data.frame(t(apply(busVoltage_VV, 1, quantile, probs=c(0.01,0.99))))
colnames(busVoltage_VV) <- c(0.01,0.99)
busVoltage_VV$Time <- ts
busVoltage_VV <- tidyr::pivot_longer(busVoltage_VV, cols=-"Time")
busVoltage_PF <- read.table("~/Desktop/case2/UnityPF/busVoltage.txt")
busVoltage_PF <- as.data.frame(t(apply(busVoltage_PF, 1, quantile, probs=c(0.01,0.99))))
colnames(busVoltage_PF) <- c(0.01,0.99)
busVoltage_PF$Time <- ts
busVoltage_PF <- melt(busVoltage_PF,id.vars = "Time")
busVoltage_PF <- tidyr::pivot_longer(busVoltage_PF, cols = -"Time")
p3 <- ggplot() + 
  geom_interval(data=mydf %>% filter(Model!="volt") %>% filter(Model!="unity") %>% mutate(across(Model, .fns = toupper)), aes(x=Time,y=value,colour=Model, group=Model),
                intervals = 0.99, size=line.size) + # , show.legend = TRUE
  # geom_interval(data = busVoltage_PF, aes(x=Time,y=value,colour="pf1"),intervals = c(0.99),size=line.size) +
  geom_fan(data = busVoltage_PF, aes(x=Time,y=value,quantile=as.numeric(as.character(name)),fill="PF1"), alpha=0.25) + #
  # geom_interval(data = busVoltage_VV, aes(x=Time,y=value,colour="voltvar"),intervals = c(0.99),size=line.size) +
  geom_fan(data = busVoltage_VV, aes(x=Time,y=value,quantile=as.numeric(as.character(name)),fill="VoltVar"), alpha=0.25) + #
  facet_wrap(~Run) +
  xlab("Time (5 min.)") +
  ylab("Voltage (p.u.)") +
  # scale_colour_viridis_d(option = "plasma") +
  # scale_color_brewer(palette = "Dark2") +
  scale_colour_manual(values=c('#E41A1C','#377EB8','#4DAF4A')) + # 
  scale_fill_manual(values = c('#CCCCCC',"#999999")) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.4,0,0), "lines"),
        text = element_text(family = "Times"),
        axis.text=element_text(family = "Times", size=plot.size),
        axis.title=element_text(family = "Times", size=plot.size),
        plot.title=element_text(family = "Times", size=plot.size),
        legend.title = element_blank(),
        legend.text = element_text(family = "Times", size=plot.size),
        legend.box.margin = ggplot2::margin(-5,-4,0,-10),
        legend.margin = ggplot2::margin(-10,0,0,0),
        strip.text = element_text(family = "Times", size=plot.size),
        strip.text.x = element_text(margin = margin(0,0,0,0, "lines")),
        panel.spacing = unit(0, "lines"),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(linetype="none") + # Hide legend of interval
  scale_x_datetime(date_breaks = "6 hour",
                   date_labels = "%H:%M") 
ggsave(filename = "~/Desktop/case2Voltages.pdf", plot = p3,
       device = cairo_pdf, units = "cm", height = 5, width = 17)
################################################################################
# Figure 4
################################################################################
lstOfPVsizes = 8405*c(70, 140, 200, 50, 90, 160, 200, 140, 60, 80)
lstOfPVnodes <- c(18, 48, 56, 66, 79, 83, 95, 250, 300, 450)

# Find the occassions around noon where PSO could not satisfy the evaluation tolerance:
dir <- "~/Desktop/case2/psoPOP100TOL001"
setwd(dir)
optimalQ <- read.table("optimalQ.txt")
busVoltage <- read.table("busVoltage.txt")
bestEval <- read.table("bestEval.txt")
tmp <- bestEval[which(complete.cases(bestEval)==TRUE),]# 
tmp1 <- tmp[as.numeric(rownames(tmp)) %in% seq(121,160,1), ]
idx <- as.numeric(rownames(tmp1))

dirs <- list.dirs("~/Desktop/case2")
dirs <- dirs[nchar(dirs)<59 & nchar(dirs)>45] # This only works when the files are in dirs
dirs

iter_res <- list()
j <- 1
for(dir in dirs){
  setwd(dir)
  lst <- list()
  tmp_case <- strsplit(dir,"/")[[1]]
  case <- tmp_case[length(tmp_case)]
  tmp_model <- strsplit(case, split="[[:upper:]]")[[1]]
  model <- tmp_model[1]
  tmp_run <- strsplit(case, split="[[:lower:]]")[[1]]
  run <- tmp_run[length(tmp_run)]
  if(startsWith(case,"Unity") | startsWith(case,"ga")){
    next
  } else if(endsWith(case,"001") | endsWith(case,"Var")) {
    if(file.exists("optimalQ.txt")){optimalQ <- read.table("optimalQ.txt")} else {next}
    optimalQ <- optimalQ[idx,] # Take only times when PSO failed
    if(model=="volt"){
      optimalQ$Model <- "voltvar"
    } else {
      optimalQ$Model <- model
    }
  }
  iter_res[[j]] <- optimalQ
  j <- j+1
}
iter_res <- do.call(rbind,iter_res)
colnames(iter_res)[1:10] <- lstOfPVnodes

for(i in 1:nrow(iter_res)){
  iter_res[i,1:10] <- iter_res[i,1:10]/lstOfPVsizes
}

p4 <- iter_res %>%
  dplyr::filter(Model!="voltvar") %>%
  tidyr::pivot_longer(cols = -Model) %>%
  dplyr::mutate(name=factor(name, levels = c(18,56,48,250,66,79,300,450,83,95))) %>%
  ggplot(data = ., mapping = aes(x=name,y=value,colour=Model)) +
  geom_boxplot(outlier.size = point.size, size=line.size+0.05) +
  scale_colour_manual(values = c("aquamarine3","blue1","darkorange1")) +
  xlab("DER bus number") +
  ylab("Relative reactive power (-)") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.4,0,0), "lines"),
        text = element_text(family = "Times"),
        axis.text=element_text(family = "Times", size=plot.size),
        axis.title=element_text(family = "Times", size=plot.size),
        plot.title=element_text(family = "Times", size=plot.size),
        legend.title = element_blank(),
        legend.text = element_text(family = "Times", size=plot.size),
        legend.box.margin = ggplot2::margin(-5,-4,0,-10),
        legend.margin = ggplot2::margin(-10,0,0,0),
        strip.text = element_text(family = "Times", size=plot.size),
        strip.text.x = element_text(margin = margin(0.05,0,0.1,0, "lines")),
        panel.spacing = unit(0, "lines"))
ggsave(filename = "~/Desktop/relativeContributionwithPSO.pdf", plot = p4,
       device = cairo_pdf, units = "cm", height = 4, width = 17)
################################################################################
# Table I
################################################################################
# From: https://stackoverflow.com/questions/52490552/r-convert-nan-to-na
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

dirs <- list.dirs("~/Desktop/case2")
dirs <- dirs[nchar(dirs)<59 & nchar(dirs)>45] # This only works when the files are in dirs
pvNodes = c(19, 49, 57, 67, 80, 84, 96, 116, 118, 119) # Col index in result matrix, different from lstOfPVnodes

res_list <- list()
j <- 1
for(dir in dirs){
  setwd(dir)
  lst <- list()
  tmp_case <- strsplit(dir,"/")[[1]]
  case <- tmp_case[length(tmp_case)]
  tmp_model <- strsplit(case, split="[[:upper:]]")[[1]]
  model <- tmp_model[1]
  tmp_run <- strsplit(case, split="[[:lower:]]")[[1]]
  run <- tmp_run[length(tmp_run)]
  if(startsWith(case,"ga")){
    next
  } else if(startsWith(case,"pso")) {
    if(file.exists("bestEval.txt")){
      tmp_eval <- read.table("bestEval.txt")
      eva <- NULL
      tmp_eval[is.nan(tmp_eval)] <- NA
      for(i in 1:nrow(tmp_eval)){
        eva[[i]] <- sum(!is.na(tmp_eval[i,]))
      } 
    } else {next}
    if(file.exists("busVoltage.txt")){tmp_vlt <- read.table("busVoltage.txt")} else {next}
    tmp_res <- data.frame(Case=case_study,Run=run,Model=model,
                          muObj=round(mean(rowSums(abs(1-tmp_vlt[,pvNodes]))),digits=4),
                          sdObj=round(sd(rowSums(abs(1-tmp_vlt[,pvNodes]))),digits=3),
                          muItr=round(mean(do.call(c,eva)),digits = 2),
                          sdItr=round(sd(do.call(c,eva)),digits = 1),
                          sd=round(mean(apply(tmp_vlt[,pvNodes],1,sd)),digits=4))
    
  } else {
    if(file.exists("busVoltage.txt")){tmp_vlt <- read.table("busVoltage.txt")} else {next}
    if(file.exists("Results.RData")){
      load("Results.RData")
      eva <- NULL
      for(i in 1:length(results)){
        eva[[i]] <- results[[i]]$numGens
      }
    } else {next}
    tmp_res <- data.frame(Case=case_study,Run=run,Model=model,
                          muObj=round(mean(rowSums(abs(1-tmp_vlt[,pvNodes]))),digits=4),
                          sdObj=round(sd(rowSums(abs(1-tmp_vlt[,pvNodes]))),digits=3),
                          muItr=round(mean(do.call(c,eva)),digits = 2),
                          sdItr=round(sd(do.call(c,eva)),digits = 1),
                          sd=round(mean(apply(tmp_vlt[,pvNodes],1,sd)),digits=4)
    )
  }
  res_list[[j]] <- tmp_res
  j <- j+1
}
res <- do.call(rbind,res_list) %>%
  dplyr::select(-Case) %>%
  mutate_if(is.numeric, format, nsmall=1) %>%
  mutate(`Objective function (p.u.)` = paste(xtable::sanitize.numbers(format(res$muObj, scientific = TRUE, digits = 2),
                                                              type = "latex", math.style.exponents = TRUE), 
                                             "$\\pm$", 
                                             xtable::sanitize.numbers(format(res$sdObj, scientific = TRUE, digits = 2),
                                                              type = "latex", math.style.exponents = TRUE))) %>%
  mutate(`Number of iterations (-)` = paste(xtable::sanitize.numbers(format(res$muItr, scientific = FALSE, digits = 2),
                                                             type = "latex", math.style.exponents = FALSE), 
                                            "$\\pm$", 
                                            xtable::sanitize.numbers(format(res$sdItr, scientific = FALSE, digits = 2),
                                                             type = "latex", math.style.exponents = FALSE))) %>%
  mutate(`Standard deviation (p.u.)` = paste(xtable::sanitize.numbers(format(res$sd, scientific = TRUE, digits = 2),
                                                              type = "latex", math.style.exponents = TRUE))) %>%
  dplyr::select(-one_of(c("muObj", "sdObj", "muItr", "sdItr", "sd")))

print(xtable::xtable(res), include.rownames=FALSE, sanitize.text.function = function(x){x})
