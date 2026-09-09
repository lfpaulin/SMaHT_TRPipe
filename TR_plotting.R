#!/usr/bin/env Rscript

options(width = 200)

tissue_colors <- list(
    BLOO = "#FF382EAA", 
    BUCC = "#75001FAA", 
    ESOP = "#FFB161AA", 
    COAS = "#C0B5E7AA", 
    COAD = "#9D4BBFAA", 
    LIVR = "#FF4AAAA9", 
    ADGL = "#EA7A00AA", 
    ADGR = "#EA7A00AA", 
    AORT = "#FAAB9BAA", 
    LUNG = "#1DDD30AA", 
    HART = "#DDD824AA", 
    TESL = "#AAF8B4AA", 
    TESR = "#AAF8B4AA", 
    OVAL = "#FFAFD4AA", 
    OVAR = "#FFAFD4AA", 
    FBRO = "#FFE98DAA", 
    SKSE = "#0EC1B8AA", 
    SKNE = "#91E5DBAA", 
    MUSC = "#AA9A41AA", 
    BRAI = "#BBE7FFFF", 
    BRFL = "#BBE7FFFF", 
    BRTL = "#76AEFFAA", 
    BRCE = "#1655C4AA", 
    BRHL = "#002A66AA", 
    BRHR = "#002A66AA", 
    BLOO_PB = "#FF382EAA", 
    BUCC_PB = "#75001FAA", 
    ESOP_PB = "#FFB161AA", 
    COAS_PB = "#C0B5E7AA", 
    COAD_PB = "#9D4BBFAA", 
    LIVR_PB = "#FF4AAAA9", 
    ADGL_PB = "#EA7A00AA", 
    ADGR_PB = "#EA7A00AA", 
    AORT_PB = "#FAAB9BAA", 
    LUNG_PB = "#1DDD30AA", 
    HART_PB = "#DDD824AA", 
    TESL_PB = "#AAF8B4AA", 
    TESR_PB = "#AAF8B4AA", 
    OVAL_PB = "#FFAFD4AA", 
    OVAR_PB = "#FFAFD4AA", 
    FBRO_PB = "#FFE98DAA", 
    SKSE_PB = "#0EC1B8AA", 
    SKNE_PB = "#91E5DBAA", 
    MUSC_PB = "#AA9A41AA", 
    BRAI_PB = "#BBE7FFFF", 
    BRFL_PB = "#BBE7FFFF", 
    BRTL_PB = "#76AEFFAA", 
    BRCE_PB = "#1655C4AA", 
    BRHL_PB = "#002A66AA", 
    BRHR_PB = "#002A66AA", 
    BLOO_PB2 = "#FF382EAA", 
    BUCC_PB2 = "#75001FAA", 
    ESOP_PB2 = "#FFB161AA", 
    COAS_PB2 = "#C0B5E7AA", 
    COAD_PB2 = "#9D4BBFAA", 
    LIVR_PB2 = "#FF4AAAA9", 
    ADGL_PB2 = "#EA7A00AA", 
    ADGR_PB2 = "#EA7A00AA", 
    AORT_PB2 = "#FAAB9BAA", 
    LUNG_PB2 = "#1DDD30AA", 
    HART_PB2 = "#DDD824AA", 
    TESL_PB2 = "#AAF8B4AA", 
    TESR_PB2 = "#AAF8B4AA", 
    OVAL_PB2 = "#FFAFD4AA", 
    OVAR_PB2 = "#FFAFD4AA", 
    FBRO_PB2 = "#FFE98DAA", 
    SKSE_PB2 = "#0EC1B8AA", 
    SKNE_PB2 = "#91E5DBAA", 
    MUSC_PB2 = "#AA9A41AA", 
    BRAI_PB2 = "#BBE7FFFF", 
    BRFL_PB2 = "#BBE7FFFF", 
    BRTL_PB2 = "#76AEFFAA", 
    BRCE_PB2 = "#1655C4AA", 
    BRHL_PB2 = "#002A66AA", 
    BRHR_PB2 = "#002A66AA", 
    BLOO_ONT = "#FF382EAA", 
    BUCC_ONT = "#75001FAA", 
    ESOP_ONT = "#FFB161AA", 
    COAS_ONT = "#C0B5E7AA", 
    COAD_ONT = "#9D4BBFAA", 
    LIVR_ONT = "#FF4AAAA9", 
    ADGL_ONT = "#EA7A00AA", 
    ADGR_ONT = "#EA7A00AA", 
    AORT_ONT = "#FAAB9BAA", 
    LUNG_ONT = "#1DDD30AA", 
    HART_ONT = "#DDD824AA", 
    TESL_ONT = "#AAF8B4AA", 
    TESR_ONT = "#AAF8B4AA", 
    OVAL_ONT = "#FFAFD4AA", 
    OVAR_ONT = "#FFAFD4AA", 
    FBRO_ONT = "#FFE98DAA", 
    SKSE_ONT = "#0EC1B8AA", 
    SKNE_ONT = "#91E5DBAA", 
    MUSC_ONT = "#AA9A41AA", 
    BRAI_ONT = "#BBE7FFFF", 
    BRFL_ONT = "#BBE7FFFF", 
    BRTL_ONT = "#76AEFFAA", 
    BRCE_ONT = "#1655C4AA", 
    BRHL_ONT = "#002A66AA", 
    BRHR_ONT = "#002A66AA", 
    BLOO_ONT2 = "#FF382EAA", 
    BUCC_ONT2 = "#75001FAA", 
    ESOP_ONT2 = "#FFB161AA", 
    COAS_ONT2 = "#C0B5E7AA", 
    COAD_ONT2 = "#9D4BBFAA", 
    LIVR_ONT2 = "#FF4AAAA9", 
    ADGL_ONT2 = "#EA7A00AA", 
    ADGR_ONT2 = "#EA7A00AA", 
    AORT_ONT2 = "#FAAB9BAA", 
    LUNG_ONT2 = "#1DDD30AA", 
    HART_ONT2 = "#DDD824AA", 
    TESL_ONT2 = "#AAF8B4AA", 
    TESR_ONT2 = "#AAF8B4AA", 
    OVAL_ONT2 = "#FFAFD4AA", 
    OVAR_ONT2 = "#FFAFD4AA", 
    FBRO_ONT2 = "#FFE98DAA", 
    SKSE_ONT2 = "#0EC1B8AA", 
    SKNE_ONT2 = "#91E5DBAA", 
    MUSC_ONT2 = "#AA9A41AA", 
    BRAI_ONT2 = "#BBE7FFFF", 
    BRFL_ONT2 = "#BBE7FFFF", 
    BRTL_ONT2 = "#76AEFFAA", 
    BRCE_ONT2 = "#1655C4AA", 
    BRHL_ONT2 = "#002A66AA", 
    BRHR_ONT2 = "#002A66AA"
)

require(argparser)
parser <- arg_parser("Boxplot from TSV")
parser <- add_argument(parser, "--tab", help = "path to TSV file")
parser <- add_argument(parser, "--out", help = "output file", default = "out.png")
parser <- add_argument(parser, "--donor", help = "donor name", default = "")
parser <- add_argument(parser, "--meta", help = "donor metadata, includes gene", default = "")
parser <- add_argument(parser, "--merget", help = "merged tissues", flag = TRUE)
parser <- add_argument(parser, "--dev", help = "debug output", flag = TRUE)
argv <- parse_args(parser)

if (argv$dev) {
    print(argv)
}

# counts
tab <- read.table(argv$tab, header=T)
tab[tab == -1] <- NA

hap <- read.table(paste(argv$tab, "2", sep=""), header=T)
hap[hap == -1] <- NA
hap[hap ==  0] <- 22  # full:15  empty:22  SQR
hap[hap ==  1] <- 21  # full:19  empty:21  CIR
hap[hap ==  2] <- 24  # full:17  empty:24  TRI

n_row <- nrow(tab)
n_col <- ncol(tab)
y_max <- max(tab[!is.na(tab)])
y_min <- min(tab[!is.na(tab)])
if (y_max <= 500){
    ylims <- c(0, as.integer(y_max*1.1))
} else if (y_max > 500 && y_min >= 0.5*y_max) {
    ylims <- c(as.integer(0.9*y_min), as.integer(y_max*1.1))
} else {
    ylims <- c(0, as.integer(y_max*1.1))
}

source(argv$meta)

# default: 5.1 4.1 4.1 2.1 (bottom, left, top, right)
custom_margins <- c(bottom=6, left=5, top=4, right=1)


plot_height <- 1000
plot_width <- 4000
dpi_res <- 150
x_spacing <- 0.1
if (argv$merget){
    x_spacing <- 0.07
    plot_width <- 2500
    dpi_res <- 200
}


png(argv$out, height=plot_height, width=plot_width, res=dpi_res)
par(mar=custom_margins)  
plot(NA, xlim=c(1, n_col), ylim=ylims, 
    main = paste(argv$donor, gene), cex.main = 2,
    xaxt = "n", xlab = NA, ylab = "Allele length", 
    cex.lab=2, cex.axis=1.2)

# axis(1, at = 1:n_col, labels = colnames(tab), las = 2)
y_offset <- y_max * 0.05
text(x = seq_along(colnames(tab)), y = par("usr")[3] - y_offset, labels = colnames(tab),
     srt = 45, adj = 1, xpd = TRUE, cex=1.2)
title(xlab = "Tissues", line = 4, cex.lab=2)


x_size <- nrow(tab)
for (index in seq(1, n_col)) {
    y <- tab[,index]
    x <- rnorm(x_size, index, x_spacing)
    pch_use <- hap[,index]
    # x <- rep(index, n_row)
    points(x, y, pch=pch_use, cex=1.5, col="#000", bg=tissue_colors[[colnames(tab)[index]]])
    abline(v=index, pwd=2, lty=3)

}
abline(h=trlen, pwd=2, lty=2)

dev.off()

# data <- read.table(argv$file, header = TRUE, sep = "\t")
# boxplot(data[[argv$col]], main = paste("Boxplot of", argv$col))

