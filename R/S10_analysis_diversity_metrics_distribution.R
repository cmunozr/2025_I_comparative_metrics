library(tidyr)
library(ggplot2)


#--------------------
a <- a |>
  as.data.frame() |> 
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "alpha_richness"
    ) |> 
  mutate(state = "baseline")


a2 <- a2 |>
  as.data.frame() |> 
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "alpha_richness"
  ) |> 
  mutate(state = "scenario")

bind_rows(a, a2) |> 
  filter(sample %in% c("V1", "V10", "V45", "V74", "V86")) |> 
  ggplot() + 
  geom_density(mapping = aes(alpha_richness, fill = state), alpha = 0.4) + 
  facet_wrap(~sample, scales = "free", nrow = 5, ncol = 1)


bind_rows(a, a2) |> 