#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rhdf5)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(scico)
  library(tictoc)
  library(Rfast)
  library(nlme)
  library(compositions)
})

# -----------------------------
# Simple CLI argument parsing
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) return(default)
  args[idx + 1]
}

input_dir  <- get_arg("--input_dir", ".")
output_dir <- get_arg("--output_dir", "results")

cat("Input dir: ", input_dir,  "\n")
cat("Output dir:", output_dir, "\n")

# Create output subdirectories
plots_root   <- file.path(output_dir, "plots")
qc_dir       <- file.path(plots_root, "segmentation_qc")
subsets_root <- file.path(plots_root, "subsets")
pca_dir      <- file.path(plots_root, "pca2")

dir.create(output_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plots_root,  recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir,      recursive = TRUE, showWarnings = FALSE)
dir.create(subsets_root,recursive = TRUE, showWarnings = FALSE)
dir.create(pca_dir,     recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Functions from the Rmd
# -----------------------------
he.rename_labels <- function(img, labels = c(
  "background","CMC","viable tumor","infiltrative tumor",
  "degraded viable tumor","degraded infiltrative tumor",
  "necrosis/tissue degradation","vasculature",
  "artifact","bubble","fold","fibre",
  "calcification","hemorrhage","debris"
)) {
  levels(img$segmentation) <- labels[as.numeric(levels(img$segmentation))]
  img
}

he.remove_nonbiological <- function(
  img,
  remove = c("background","CMC","artifact","bubble","fold","fibre","debris")
) {
  img$segmentation[img$segmentation %in% remove] <- NA
  img
}

he.montecarlo <- function(img, n = 1000, size = NULL, um = NULL,
                          ratio = 0.32249742002) {
  if (is.null(size))
    size <- um / ratio

  # Sample n centers within the tissue
  tictoc::tic("Sample centers")
  tissue  <- which(!is.na(img$segmentation))
  centers <- sample(tissue, n)
  tictoc::toc()

  # Order coordinates
  tictoc::tic("Sort coords")
  ox <- sort(img$pos$x, index.return = TRUE, method = "radix")
  oy <- sort(img$pos$y, index.return = TRUE, method = "radix")
  tictoc::toc()

  # Get pixels of the ROI
  tictoc::tic("Search ROIs")
  hx <- sapply(img$pos[centers, 1] + size/2,
               function(x) Rfast::binary_search(ox$x, x, TRUE))
  lx <- sapply(img$pos[centers, 1] - size/2,
               function(x) Rfast::binary_search(ox$x, x, TRUE))
  hy <- sapply(img$pos[centers, 2] + size/2,
               function(x) Rfast::binary_search(oy$x, x, TRUE))
  ly <- sapply(img$pos[centers, 2] - size/2,
               function(x) Rfast::binary_search(oy$x, x, TRUE))
  tictoc::toc()

  tictoc::tic("Build index sets")
  ix <- lapply(seq_along(centers), function(i) ox$ix[lx[i]:hx[i]])
  iy <- lapply(seq_along(centers), function(i) oy$ix[ly[i]:hy[i]])
  tictoc::toc()

  tictoc::tic("Intersect indices")
  rois <- lapply(seq_along(centers), function(i) intersect(ix[[i]], iy[[i]]))
  tictoc::toc()

  # ROI-level class tables
  tictoc::tic("Tabulate ROIs")
  mc <- as.data.frame(do.call(rbind,
                              lapply(rois,
                                     function(is) table(img$segmentation[is]))))
  tictoc::toc()

  mcnorm <- t(apply(mc, 1, function(x) x / sum(x)))
  list(mc = mc, mcnorm = mcnorm, centers = centers, rois = rois)
}

# -----------------------------
# Prepare metadata and load data
# -----------------------------
meta <- data.frame(
  file = list.files(input_dir,
                    pattern = "Probabilities.h5",
                    full.names = TRUE)
)

if (nrow(meta) == 0) {
  stop("No *_Probabilities.h5 files found in: ", input_dir)
}

meta$sample    <- gsub("_Probabilities.h5", "", basename(meta$file))
meta$biopsy    <- sapply(strsplit(meta$sample, "-"), function(x) x[[1]])
meta$replicate <- sapply(strsplit(meta$sample, "-"), function(x) x[[2]])
meta$csv       <- file.path(output_dir,
                            paste0(meta$sample, "_montecarlo.csv"))

write.table(
  meta,
  file.path(output_dir, "metadata.csv"),
  sep = ",",
  row.names = FALSE
)

cat("Found", nrow(meta), "samples\n")

# Load probabilities and build img objects
probs <- lapply(meta$file, function(x) h5read(x, "/exported_data"))
segmentation <- lapply(probs, function(x) apply(x, c(2, 3), which.max))
imgs <- lapply(segmentation, function(x) list(
  segmentation = factor(x),
  pos          = expand.grid(x = 1:nrow(x), y = 1:ncol(x))
))
imgs <- lapply(imgs, he.rename_labels)
imgs <- lapply(imgs, he.remove_nonbiological)

# -----------------------------
# Segmentation QC plots
# -----------------------------
labels <- c("viable tumor","infiltrative tumor","degraded viable tumor",
            "degraded infiltrative tumor","necrosis/tissue degradation",
            "vasculature","other")
seg_colors <- scico::scico(length(labels), palette = "managua")
seg_colors <- setNames(seg_colors, labels)

f <- 5  # downsampling factor

for (i in seq_along(imgs)) {
  df <- data.frame(imgs[[i]]$pos)
  df$i <- as.character(imgs[[i]]$segmentation)
  df$i[df$i %in% c("calcification","hemorrhage")] <- "other"
  df$i <- factor(df$i, levels = labels)
  df   <- df[df$x %% f == 0 & df$y %% f == 0, ]
  df$x <- df$x / f
  df$y <- df$y / f

  p <- ggplot(df, aes(x = x, y = -y, fill = i)) +
    geom_raster() +
    theme_void() +
    coord_fixed() +
    ggtitle(meta$sample[i]) +
    scale_fill_manual(values = seg_colors, na.value = "transparent")

  ggsave(
    file.path(qc_dir, paste0(meta$sample[i], ".tiff")),
    plot   = p,
    dpi    = 300,
    width  = 4,
    height = 4
  )
}

# -----------------------------
# Monte Carlo simulations
# -----------------------------
ums  <- c(200, 100, 50, 20)
sims <- lapply(ums, function(um) lapply(imgs, he.montecarlo, um = um))

# Save per-sample, per-um CSVs with normalized compositions
for (i_um in seq_along(sims)) {
  this_um <- ums[i_um]
  for (j in seq_along(imgs)) {
    out_csv <- gsub(".csv", paste0(this_um, "um.csv"), meta$csv[j])
    write.table(
      sims[[i_um]][[j]]$mcnorm,
      out_csv,
      sep       = ",",
      row.names = FALSE
    )
  }
}

# -----------------------------
# Example ROI ("Show some segments")
# -----------------------------
for (k in seq_along(sims)) {
  this_um <- ums[k]
  subset_dir <- file.path(subsets_root, paste0(this_um, "um"))
  dir.create(subset_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_along(imgs)) {
    for (j in 1:10) {
      is <- sims[[k]][[i]]$rois[[j]]
      df <- data.frame(imgs[[i]]$pos[is, ])
      df$i <- imgs[[i]]$segmentation[is]

      p <- ggplot(df, aes(x = x, y = -y, fill = i)) +
        geom_raster() +
        theme_void() +
        coord_fixed() +
        ggtitle(meta$sample[i]) +
        scale_fill_manual(values = seg_colors, na.value = "transparent") +
        theme(legend.position = "none")

      ggsave(
        file.path(subset_dir,
                  paste0(meta$sample[i], "_", j, ".tiff")),
        plot   = p,
        dpi    = 300,
        width  = 4,
        height = 4
      )
    }
  }
}

# -----------------------------
# PCA
# -----------------------------
data_all <- do.call(
  dplyr::bind_rows,
  lapply(sims, function(x)
    do.call(dplyr::bind_rows,
            lapply(x, function(y) as.data.frame(y$mcnorm))))
)
data_all[is.na(data_all)] <- 0

pca <- prcomp(data_all)

for (um in ums) {
  df <- data.frame(pca$x)
  # NOTE: these rep() patterns follow the original script
  df$size      <- rep(ums,         each = 7000)
  df$sample    <- rep(meta$sample, each = 1000)
  df$biopsy    <- rep(meta$biopsy, each = 1000)
  df$replicate <- rep(meta$replicate, each = 1000)

  ve   <- (pca$sdev)^2 / sum((pca$sdev)^2)
  xlab <- paste0("PC1 (", round(100 * ve[1], 2), "%)")
  ylab <- paste0("PC2 (", round(100 * ve[2], 2), "%)")

  colors <- c(
    colorRampPalette(c("#FFFFFF", "#9966CC","#000000"))(5)[2:4],
    colorRampPalette(c("#FFFFFF", "#FA8072","#000000"))(4)[2:3],
    colorRampPalette(c("#FFFFFF", "#008080","#000000"))(4)[2:3]
  )

  df_um <- subset(df, size == um)

  # scatter
  p1 <- ggplot(df_um[sample(1:nrow(df_um), nrow(df_um)), ],
               aes(x = PC1, y = PC2, color = biopsy)) +
    geom_point(alpha = 0.1) +
    theme_classic() +
    coord_fixed() +
    scale_color_manual(values = c("#9966CC","#FA8072","#008080")) +
    xlim(c(-1.5, 1.5)) + ylim(c(-1.5, 1.5)) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom") +
    xlab(xlab) + ylab(ylab)

  # density all biopsies
  p2 <- ggplot(df_um, aes(x = PC1, y = PC2, color = biopsy)) +
    geom_density2d(h = 0.5) +
    theme_classic() +
    coord_fixed() +
    scale_color_manual(values = c("#9966CC","#FA8072","#008080")) +
    xlim(c(-1.5, 1.5)) + ylim(c(-1.5, 1.5)) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom") +
    xlab(xlab) + ylab(ylab)

  ggsave(file.path(pca_dir, paste0("pca_all_", um, "um.tiff")),
         p2, dpi = 300, height = 4, width = 4)

  # per-biopsy density plots
  p_b1 <- ggplot(subset(df_um, biopsy == 1),
                 aes(x = PC1, y = PC2, color = sample)) +
    geom_density2d(h = 0.5) +
    theme_classic() + coord_fixed() +
    scale_color_manual(values = colors[1:3]) +
    xlim(c(-1.5, 1.5)) + ylim(c(-1.5, 1.5)) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom") +
    xlab(xlab) + ylab(ylab)

  ggsave(file.path(pca_dir, paste0("pca_1_", um, "um.tiff")),
         p_b1, dpi = 300, height = 4, width = 4)

  p_b2 <- ggplot(subset(df_um, biopsy == 2),
                 aes(x = PC1, y = PC2, color = sample)) +
    geom_density2d(h = 0.5) +
    theme_classic() + coord_fixed() +
    scale_color_manual(values = colors[4:5]) +
    xlim(c(-1.5, 1.5)) + ylim(c(-1.5, 1.5)) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom") +
    xlab(xlab) + ylab(ylab)

  ggsave(file.path(pca_dir, paste0("pca_2_", um, "um.tiff")),
         p_b2, dpi = 300, height = 4, width = 4)

  p_b3 <- ggplot(subset(df_um, biopsy == 3),
                 aes(x = PC1, y = PC2, color = sample)) +
    geom_density2d(h = 0.5) +
    theme_classic() + coord_fixed() +
    scale_color_manual(values = colors[6:7]) +
    xlim(c(-1.5, 1.5)) + ylim(c(-1.5, 1.5)) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom") +
    xlab(xlab) + ylab(ylab)

  ggsave(file.path(pca_dir, paste0("pca_3_", um, "um.tiff")),
         p_b3, dpi = 300, height = 4, width = 4)

  write.csv(
    df_um,
    file.path(pca_dir, paste0("pca_", um, "um.csv")),
    row.names = TRUE
  )

  # Loadings
  rot <- data.frame(pca$rotation)
  rot$tissue <- rownames(rot)
  rot$length <- apply(rot[, 1:2], 1,
                      function(x) sqrt(sum(x^2)))

  p_load <- ggplot(subset(rot, length > 0.05),
                   aes(x = 0, y = 0, xend = PC1, yend = PC2)) +
    geom_segment(arrow = arrow(length = unit(0.2, "cm")),
                 color = "#D4A017") +
    geom_text_repel(aes(x = PC1, y = PC2, label = tissue)) +
    theme_classic() +
    coord_fixed() +
    xlab(xlab) + ylab(ylab)

  ggsave(file.path(pca_dir, paste0("loadings_", um, "um.tiff")),
         p_load, dpi = 300, height = 4, width = 4)

  write.csv(
    rot,
    file.path(pca_dir, paste0("loadings_", um, "um.csv")),
    row.names = TRUE
  )
}

# -----------------------------
# Mixed-effects statistics
# -----------------------------
data <- data_all
data[is.na(data)] <- 0
colnames(data) <- c(
  "background","CMC","vt","it","dvt","dit","n","v",
  "artifact","bubble","fold","fibre","h","debris","c"
)
data$size      <- rep(ums,         each = 7000)
data$sample    <- rep(meta$sample, each = 1000)
data$biopsy    <- rep(meta$biopsy, each = 1000)
data$replicate <- rep(meta$replicate, each = 1000)

df_res <- expand.grid(
  pvalue    = NA_real_,
  comparison = c("1v2","1v3","2v3"),
  test       = c("viable_tumor","all","all_20um","all_corrected_20um")
)

# Viable tumor
out <- nlme::lme(
  fixed  = vt ~ biopsy*size + biopsy + size,
  data   = subset(data, biopsy %in% c(1, 2)),
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[1] <- coef(summary(out))[2, "p-value"]

out <- nlme::lme(
  fixed  = vt ~ biopsy*size + biopsy + size,
  data   = subset(data, biopsy %in% c(1, 3)),
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[2] <- coef(summary(out))[2, "p-value"]

out <- nlme::lme(
  fixed  = vt ~ biopsy*size + biopsy + size,
  data   = subset(data, biopsy %in% c(2, 3)),
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[3] <- coef(summary(out))[2, "p-value"]

# All classes
out <- nlme::lme(
  fixed  = n + vt + it + v + c ~ biopsy*size + biopsy + size,
  data   = subset(data, biopsy %in% c(1, 2)),
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[4] <- coef(summary(out))[2, "p-value"]

out <- nlme::lme(
  fixed  = n + vt + it + v + c ~ biopsy*size + biopsy + size,
  data   = subset(data, biopsy %in% c(1, 3)),
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[5] <- coef(summary(out))[2, "p-value"]

out <- nlme::lme(
  fixed  = n + vt + it + v + c ~ biopsy*size + biopsy + size,
  data   = subset(data, biopsy %in% c(2, 3)),
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[6] <- coef(summary(out))[2, "p-value"]

# All classes, 20um only
out <- nlme::lme(
  fixed  = n + vt + it + v + c ~ biopsy,
  data   = subset(data, biopsy %in% c(1, 2) & size == 20),
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[7] <- coef(summary(out))[2, "p-value"]

out <- nlme::lme(
  fixed  = n + vt + it + v + c ~ biopsy,
  data   = subset(data, biopsy %in% c(1, 3) & size == 20),
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[8] <- coef(summary(out))[2, "p-value"]

out <- nlme::lme(
  fixed  = n + vt + it + v + c ~ biopsy,
  data   = subset(data, biopsy %in% c(2, 3) & size == 20),
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[9] <- coef(summary(out))[2, "p-value"]

# All classes, CLR-corrected, 20um
tmp <- compositions::clr(data[, 1:15])
d20 <- subset(cbind(tmp, data[, 16:ncol(data)]),
              biopsy %in% c(1, 2) & size == 20)
out <- nlme::lme(
  fixed  = n + vt + it + v + c + dvt + dit + h ~ biopsy,
  data   = d20,
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[10] <- coef(summary(out))[2, "p-value"]

d20 <- subset(cbind(tmp, data[, 16:ncol(data)]),
              biopsy %in% c(1, 3) & size == 20)
out <- nlme::lme(
  fixed  = n + vt + it + v + c + dvt + dit + h ~ biopsy,
  data   = d20,
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[11] <- coef(summary(out))[2, "p-value"]

d20 <- subset(cbind(tmp, data[, 16:ncol(data)]),
              biopsy %in% c(2, 3) & size == 20)
out <- nlme::lme(
  fixed  = n + vt + it + v + c + dvt + dit + h ~ biopsy,
  data   = d20,
  random = ~1 | sample,
  method = "ML"
)
df_res$pvalue[12] <- coef(summary(out))[2, "p-value"]

write.csv(
  df_res,
  file.path(plots_root, "group_comparisons.csv"),
  row.names = FALSE
)

cat("Analysis complete.\n")
