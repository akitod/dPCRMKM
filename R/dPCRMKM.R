#'dPCRQC: Digital PCR quality control program
#' @importFrom XML xmlParse
#' @importFrom XML getNodeSet
#' @importFrom XML xmlGetAttr
#' @importFrom XML xmlValue
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ylim
#' @importFrom grDevices dev.off
#' @importFrom grDevices tiff
#' @importFrom graphics hist
#' @importFrom graphics plot
#' @importFrom stats median
#' @importFrom utils read.table
#' @importFrom utils write.table
#' 
#' @description Digital PCR quality control program. You can use only xml datas derived from the QuantStudio 3D Digital PCR instrument. Please insert the xml datas in the same folder and execute this program on the folder.
#' @author Akito Dobashi dobashi-jik@umin.ac.jp
#' @param quality_filter quality threshold (default : 0.6)
#' @export

dPCRQC <- function(quality_filter){
  if (is.na(quality_filter)){
    quality_filter <- 0.6
  }
  xmlfile <- dir()
  cnt <- grep(".xml", xmlfile)
  V1 <- NULL
  V2 <- NULL
  V3 <- NULL
  data_c <- NULL
  for ( i in 1:length(cnt)){
    data_a  <- xmlParse(xmlfile[cnt[i]])
    well    <- getNodeSet(data_a, "//chip/wells/well")
    quality <- sapply(well, function(x) xmlGetAttr(x, "quality"))
    fam    <- getNodeSet(data_a, "//chip/wells/well/fam")
    fam_a <- sapply(fam, function(x) value <- xmlValue(x))
    vic    <- getNodeSet(data_a, "//chip/wells/well/vic")
    vic_a <- sapply(vic, function(x) value <- xmlValue(x))
    data_b <- as.data.frame(cbind(as.numeric(as.character(quality)),
                                  as.numeric(as.character(fam_a)),
                                  as.numeric(as.character(vic_a))))
    data_c <- rbind(data_c, data_b)
  }
  data_d <- subset(data_c, V1 >= quality_filter)
  filter_pass <- length(data_d$V1) / length(data_c$V1)
  g <- ggplot(data_c, aes(x = V1))
  g <- g + geom_histogram(binwidth = 0.01) + xlab("quality") + ylab("count") +
    geom_vline(xintercept = quality_filter, colour = "red")
  tiff(filename = "dPCR_QC.tif", width = 1280, height = 960, units = "px", compression = "lzw");
  plot(g)
  dev.off()
  data_e <- hist(data_d$V2, breaks = seq(-5000, 20000, by = 20))
  FAM_COUNT1 <- which.max(data_e$counts)
  FAM_COUNT2 <- FAM_COUNT1 + 40
  FAM_MAX1 <- max(data_e$counts[c(FAM_COUNT2:length(data_e$counts))])
  FAM_COUNT3 <- max(match(FAM_MAX1, data_e$counts))
  FAM_MIN1 <- min(data_e$counts[c(FAM_COUNT1:FAM_COUNT3)])
  FAM_COUNT4 <- match(FAM_MIN1, data_e$counts)
  FAM_COUNT5 <- FAM_COUNT4[c(FAM_COUNT4 >= FAM_COUNT1) & c(FAM_COUNT4 <= FAM_COUNT3)]
  FAM_COUNT6 <- median(FAM_COUNT5)
  FAM_DIS <- (data_e$breaks[FAM_COUNT3] - data_e$breaks[FAM_COUNT1]) / 2
  FAM_BORDER <- data_e$breaks[FAM_COUNT6]
  data_f <- hist(data_d$V3, breaks = seq(-5000, 20000, by = 20))
  VIC_COUNT1 <- which.max(data_f$counts)
  VIC_COUNT2 <- VIC_COUNT1 + 40
  VIC_MAX1 <- max(data_f$counts[c(VIC_COUNT2:length(data_f$counts))])
  VIC_COUNT3 <- max(match(VIC_MAX1, data_f$counts))
  VIC_MIN1 <- min(data_f$counts[c(VIC_COUNT1:VIC_COUNT3)])
  VIC_COUNT4 <- match(VIC_MIN1, data_f$counts)
  VIC_COUNT5 <- VIC_COUNT4[c(VIC_COUNT4 >= VIC_COUNT1) & c(VIC_COUNT4 <= VIC_COUNT3)]
  VIC_COUNT6 <- median(VIC_COUNT5)
  VIC_DIS <- (data_f$breaks[VIC_COUNT3] - data_f$breaks[VIC_COUNT1]) / 2
  VIC_BORDER <- data_f$breaks[VIC_COUNT6]
  distance <- c(FAM_DIS, VIC_DIS)[which.max(c(FAM_DIS, VIC_DIS))]
  g <- ggplot(data_d, aes(x = V2))
  g <- g + geom_histogram(binwidth = 20) + xlab("FAM") + ylab("count") + xlim(-1000, 4000) + geom_vline(xintercept = FAM_BORDER, colour = "red")
  tiff(filename = "dPCR_FAM.tif", width = 1280, height = 960, units = "px", compression = "lzw");
  plot(g)
  dev.off()
  g <- ggplot(data_d, aes(x = V3))
  g <- g + geom_histogram(binwidth = 20) + xlab("VIC") + ylab("count") + xlim(-1000, 4000) + geom_vline(xintercept = VIC_BORDER, colour = "red")
  tiff(filename = "dPCR_VIC.tif", width = 1280, height = 960, units = "px", compression = "lzw");
  plot(g)
  dev.off()
  data_g <- cbind(c("quality_filter", "distance", "FAM_BORDER", "VIC_BORDER"), c(quality_filter, distance, FAM_BORDER, VIC_BORDER))
  write.table(data_g, "dPCR_KM.ini", sep = "\t", append = F, quote = F, row.names = F, col.names = F)
}

#'dPCRKM: Digital PCR mutation calculation program
#' @importFrom XML xmlParse
#' @importFrom XML getNodeSet
#' @importFrom XML xmlGetAttr
#' @importFrom XML xmlValue
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 ylim
#' @importFrom grDevices dev.off
#' @importFrom grDevices tiff
#' @importFrom graphics hist
#' @importFrom graphics plot
#' @importFrom stats median
#' @importFrom utils read.table
#' @importFrom utils write.table
#' 
#' @description Digital PCR mutation calculation program by modified k-means clustering algorithm. You can use only xml datas derived from the QuantStudio 3D Digital PCR instrument. Please insert the xml datas in the same folder and execute this program on the folder. In order to automatically change the setting value, I recommend execute after the dPCRQC program.
#' @author Akito Dobashi dobashi-jik@umin.ac.jp
#' @param quality_filter quality threshold (default : 0.6)
#' @export

dPCRKM <- function(quality_filter){
  if (is.na(quality_filter)){
    quality_filter <- 0.6
  }
  distance <- 1000
  FAM_BORDER <- 1000
  VIC_BORDER <- 1000
  CATEGORY <- NULL
  VIC <- NULL
  FAM <- NULL
  QUALITY <- NULL
  V1 <- NULL
  V2 <- NULL
  V3 <- NULL
  V4 <- NULL
  if (!file.exists("dPCR_KM.ini")){
    data_a <- cbind(c("quality_filter", "FAM_BORDER", "VIC_BORDER"), c(quality_filter, FAM_BORDER, VIC_BORDER))
    write.table(data_a, "dPCR_KM.ini", sep = "\t", append = F, quote = F, row.names = F, col.names = F)
  }
  data_b <- read.table("dPCR_KM.ini", header = F, sep = "\t")
  quality_filter <- subset(data_b, V1 == "quality_filter")[, 2]
  distance <- subset(data_b, V1 == "distance")[, 2]
  FAM_BORDER <- subset(data_b, V1 == "FAM_BORDER")[, 2]
  VIC_BORDER <- subset(data_b, V1 == "VIC_BORDER")[, 2]
  xmlfile <- dir()
  cnt <- grep(".xml", xmlfile)
  data_i <- c("SAMPLE_NAME", "NC", "WT", "MT", "MX", "Poisson", "WT_P", "MT_P", "MT_R", "-99%CI", "+99%CI", "WT_PT", "MT_PT", "MT_RT", "-99%CIT", "+99%CIT", "NN_VIC", "NN_FAM", "NP_VIC", "NP_FAM", "PN_VIC", "PN_FAM", "PP_VIC", "PP_FAM")
  for ( i in 1:length(cnt)){
    data_c  <- xmlParse(xmlfile[cnt[i]])
    well    <- getNodeSet(data_c, "//chip/wells/well")
    quality <- sapply(well, function(x) xmlGetAttr(x, "quality"))
    fam    <- getNodeSet(data_c, "//chip/wells/well/fam")
    fam_a <- sapply(fam, function(x) value <- xmlValue(x))
    vic    <- getNodeSet(data_c, "//chip/wells/well/vic")
    vic_a <- sapply(vic, function(x) value <- xmlValue(x))
    data_d <- as.data.frame(cbind(as.numeric(as.character(quality)),
                                  as.numeric(as.character(vic_a)),
                                  as.numeric(as.character(fam_a))))
    data_e <- subset(data_d, V1 >= quality_filter)
    data_epp <- rbind(c("QUALITY", "VIC", "FAM"), subset(data_e, V2 >= VIC_BORDER & V3 >= FAM_BORDER))
    data_epn <- rbind(c("QUALITY", "VIC", "FAM"), subset(data_e, V2 >= VIC_BORDER & V3 < FAM_BORDER))
    data_enp <- rbind(c("QUALITY", "VIC", "FAM"), subset(data_e, V2 < VIC_BORDER & V3 >= FAM_BORDER))
    data_enn <- rbind(c("QUALITY", "VIC", "FAM"), subset(data_e, V2 < VIC_BORDER & V3 < FAM_BORDER))
    colnames(data_epp ) <- c("QUALITY", "VIC", "FAM")
    colnames(data_epn ) <- c("QUALITY", "VIC", "FAM")
    colnames(data_enp ) <- c("QUALITY", "VIC", "FAM")
    colnames(data_enn ) <- c("QUALITY", "VIC", "FAM")
    data_epp <- subset(cbind(data_epp, 4), QUALITY != "QUALITY")
    data_epn <- subset(cbind(data_epn, 3), QUALITY != "QUALITY")
    data_enp <- subset(cbind(data_enp, 2), QUALITY != "QUALITY")
    data_enn <- subset(cbind(data_enn, 1), QUALITY != "QUALITY")
    data_g    <- rep(5, length(data_e$V1))
    data_e$V4 <- rep(6, length(data_e$V1))
    data_k <- data_g
    COUNT <- 1
    while (!identical(data_e$V4, data_g) & !identical(data_e$V4, data_k)){
      COUNT <- COUNT + 1
      data_k <- data_g
      data_g <- data_e$V4
      colnames(data_epp ) <- c("QUALITY", "VIC", "FAM", "CATEGORY")
      colnames(data_epn ) <- c("QUALITY", "VIC", "FAM", "CATEGORY")
      colnames(data_enp ) <- c("QUALITY", "VIC", "FAM", "CATEGORY")
      colnames(data_enn ) <- c("QUALITY", "VIC", "FAM", "CATEGORY")
      data_eppm <- c(mean(as.numeric(as.character(data_epp$VIC))), mean(as.numeric(as.character(data_epp$FAM))))
      data_epnm <- c(mean(as.numeric(as.character(data_epn$VIC))), mean(as.numeric(as.character(data_epn$FAM))))
      data_enpm <- c(mean(as.numeric(as.character(data_enp$VIC))), mean(as.numeric(as.character(data_enp$FAM))))
      data_ennm <- c(mean(as.numeric(as.character(data_enn$VIC))), mean(as.numeric(as.character(data_enn$FAM))))
      if (is.na(data_ennm[1])){
        data_ennm <- c(data_epnm[1] - 2 * distance, data_epnm[2])
      }
      if (is.na(data_eppm[1])){
        data_eppm <- c(data_ennm[1] + 2 * distance, data_ennm[2] + 2 * distance)
      }
      if (is.na(data_epnm[1])){
        data_epnm <- c(data_ennm[1] + 2 * distance, data_ennm[2])
      }
      if (is.na(data_enpm[1])){
        data_enpm <- c(data_ennm[1], data_ennm[2] + 2 * distance)
      }
      if (data_eppm[2] < FAM_BORDER - distance / 2){
        data_eppm[2] <- FAM_BORDER - distance / 2
      }
      if (data_enpm[2] < FAM_BORDER - distance / 2){
        data_enpm[2] <- FAM_BORDER - distance / 2
      }
      if (data_eppm[1] < VIC_BORDER - distance / 2){
        data_eppm[1] <- VIC_BORDER - distance / 2
      }
      if (data_epnm[1] < VIC_BORDER - distance / 2){
        data_epnm[1] <- VIC_BORDER - distance / 2
      }
      if (data_eppm[1] < data_enpm[1] + distance){
        data_eppm[1] <- data_enpm[1] + distance
      }
      if (data_eppm[1] < data_ennm[1] + distance){
        data_eppm[1] <- data_ennm[1] + distance
      }
      if (data_eppm[2] < data_epnm[2] + distance){
        data_eppm[2] <- data_epnm[2] + distance
      }
      if (data_eppm[2] < data_ennm[2] + distance){
        data_eppm[2] <- data_ennm[2] + distance
      }
      if (data_epnm[1] < data_enpm[1] + distance){
        data_epnm[1] <- data_enpm[1] + distance
      }
      if (data_epnm[1] < data_ennm[1] + distance){
        data_epnm[1] <- data_ennm[1] + distance
      }
      if (data_enpm[2] < data_epnm[2] + distance){
        data_enpm[2] <- data_epnm[2] + distance
      }
      if (data_enpm[2] < data_ennm[2] + distance){
        data_enpm[2] <- data_ennm[2] + distance
      }
      if (data_eppm[1] >= data_epnm[1] + distance / 2){
        data_eppm[1] <- data_epnm[1] + distance / 2
      }
      if (data_eppm[1] < data_epnm[1] - distance / 2){
        data_eppm[1] <- data_epnm[1] - distance / 2
      }
      if (data_enpm[1] >= data_ennm[1] + distance / 2){
        data_enpm[1] <- data_ennm[1] + distance / 2
      }
      if (data_enpm[1] < data_ennm[1] - distance / 2){
        data_enpm[1] <- data_ennm[1] - distance / 2
      }
      if (data_eppm[2] >= data_enpm[2] + distance / 2){
        data_eppm[2] <- data_enpm[2] + distance / 2
      }
      if (data_eppm[2] < data_enpm[2] - distance / 2){
        data_eppm[2] <- data_enpm[2] - distance / 2
      }
      if (data_epnm[2] >= data_ennm[2] + distance / 2){
        data_epnm[2] <- data_ennm[2] + distance / 2
      }
      if (data_epnm[2] < data_ennm[2] - distance / 2){
        data_epnm[2] <- data_ennm[2] - distance / 2
      }
      data_f <- NULL
      data_f$V1 <- (data_e$V2 - data_ennm[1]) ^ 2 + (data_e$V3 - data_ennm[2]) ^ 2
      data_f$V2 <- (data_e$V2 - data_enpm[1]) ^ 2 + (data_e$V3 - data_enpm[2]) ^ 2
      data_f$V3 <- (data_e$V2 - data_epnm[1]) ^ 2 + (data_e$V3 - data_epnm[2]) ^ 2
      data_f$V4 <- (data_e$V2 - data_eppm[1]) ^ 2 + (data_e$V3 - data_eppm[2]) ^ 2
      data_f <- as.data.frame(data_f)
      data_e$V4 <- apply(data_f, 1, which.min)
      data_epp <- subset(data_e, V4 == 4)
      data_epn <- subset(data_e, V4 == 3)
      data_enp <- subset(data_e, V4 == 2)
      data_enn <- subset(data_e, V4 == 1)
      if (COUNT == 1000){
        break
      }
    }
    colnames(data_e) <- c("QUALITY", "VIC", "FAM", "CATEGORY")
    data_e$CATEGORY <- gsub(1, "NN", data_e$CATEGORY)
    data_e$CATEGORY <- gsub(2, "NP", data_e$CATEGORY)
    data_e$CATEGORY <- gsub(3, "PN", data_e$CATEGORY)
    data_e$CATEGORY <- gsub(4, "PP", data_e$CATEGORY)
    FILE_NAME <- sub("xml", "tif", xmlfile[cnt[i]])
    g <- ggplot(data_e, aes(x = VIC, y = FAM))
    g <- g + geom_point(aes(colour = CATEGORY)) + xlab("VIC") + ylab("FAM") + xlim(-2000, 4000) + ylim(-2000, 4000)
    g <- g + scale_colour_manual(values = c(PP = "#35a16b", PN = "#66ccff", NP = "#ff99a0", NN = "#000000"))
    tiff(filename = FILE_NAME, width = 1280, height = 960, units = "px", compression = "lzw");
    plot(g)
    dev.off()
    OUT_FILE <- sub(".xml", "", xmlfile[cnt[i]])
    OUT_NN   <- length(data_enn$V1)
    OUT_PN   <- length(data_epn$V1)
    OUT_NP   <- length(data_enp$V1)
    OUT_PP   <- length(data_epp$V1)
    OUT_PS   <- -log(OUT_NN / (OUT_PP + OUT_PN + OUT_NP + OUT_NN))
    OUT_PSW  <- (OUT_PN + OUT_PP) * OUT_PS
    OUT_PSM  <- (OUT_NP + OUT_PP) * OUT_PS
    OUT_PSWT  <- OUT_PN * OUT_PS
    OUT_PSMT  <- OUT_NP * OUT_PS
    OUT_MR   <- OUT_PSM / (OUT_PSW + OUT_PSM)
    OUT_MRT  <- OUT_PSMT / (OUT_PSWT + OUT_PSMT)
    OUT_MMP  <- OUT_MR - 2.58 * ( (OUT_MR * (1 - OUT_MR) / (OUT_PSW + OUT_PSM)) ^ (1 / 2))
    OUT_MPP  <- OUT_MR + 2.58 * ( (OUT_MR * (1 - OUT_MR) / (OUT_PSW + OUT_PSM)) ^ (1 / 2))
    OUT_MMPT  <- OUT_MRT - 2.58 * ( (OUT_MRT * (1 - OUT_MRT) / (OUT_PSWT + OUT_PSMT)) ^ (1 / 2))
    OUT_MPPT  <- OUT_MRT + 2.58 * ( (OUT_MRT * (1 - OUT_MRT) / (OUT_PSWT + OUT_PSMT)) ^ (1 / 2))
    data_j <- c(OUT_FILE, OUT_NN, OUT_PN, OUT_NP, OUT_PP, OUT_PS, OUT_PSW, OUT_PSM, OUT_MR, OUT_MMP, OUT_MPP, OUT_PSWT, OUT_PSMT, OUT_MRT, OUT_MMPT, OUT_MPPT,
                data_ennm[1], data_ennm[2], data_enpm[1], data_enpm[2], data_epnm[1], data_epnm[2], data_eppm[1], data_eppm[2])
    data_i <- rbind(data_i, data_j)
  }
  write.table(data_i, "dPCR_OUT.txt", sep = "\t", append = F, quote = F, row.names = F, col.names = F)
}
