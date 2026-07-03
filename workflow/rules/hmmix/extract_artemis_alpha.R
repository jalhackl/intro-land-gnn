# process artemis output, following https://github.com/MoiColl/HMMenhancements/blob/main/simulations/hmmix_simulated.ipynb
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

#---------------------------------------------------------------------------------------------------------------------
in_artemis <- read.table(infile, header = TRUE) # take Artemis file as input

# Select rows for alpha=0 or alpha=1
max_vals <- in_artemis[in_artemis$alpha == 0 | in_artemis$alpha == 1, ]

# Create min/max axis values as a data frame (not a matrix!)
minmaxaxis <- data.frame(
  max_x = max(max_vals[max_vals$alpha == 0, 'average_pointwise_posterior_probability']),
  min_x = min(max_vals[max_vals$alpha == 1, 'average_pointwise_posterior_probability']),
  max_y = max(max_vals[max_vals$alpha == 1, 'loglikelihood']),
  min_y = min(max_vals[max_vals$alpha == 0, 'loglikelihood'])
)

# extract best alpha without having to generate plot
in_artemis1 <- in_artemis %>%
  # add min, max vals
  mutate(
    max_x = minmaxaxis$max_x,
    min_x = minmaxaxis$min_x,
    max_y = minmaxaxis$max_y,
    min_y = minmaxaxis$min_y
  ) %>%
  # calc distances following https://github.com/MoiColl/HMMenhancements/blob/main/simulations/hmmix_simulated.ipynb
  mutate(
    x = (average_pointwise_posterior_probability - min_x) / (max_x - min_x),
    y = (loglikelihood - min_y) / (max_y - min_y)
  ) %>%
  mutate(dist = abs(x - y) / sqrt(2)) %>%
  mutate(best_alpha = dist == min(dist, na.rm = TRUE))

out_alpha <- in_artemis1$alpha[in_artemis1$best_alpha][1]

# output this & pass to hybrid decode
write.table(out_alpha, file = outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)
#---------------------------------------------------------------------------------------------------------------------