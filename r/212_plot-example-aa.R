source("r/header.R")
# library(gridExtra)
# library(grid)
# library(gridtext)

scenarios = c("No amphistomy\nadvantage",
              "Significant amphistomy\nadvantage")
treatments = c("pseudohypo\nlower surface only", "amphi\nboth surfaces")

df_segment = crossing(
  scenario = factor(scenarios, levels = scenarios),
  nesting(
    `Gas exchange\nthrough:` = treatments,
    x = c(-1, -0.5),
    xend = c(0.5, 1)
  )
) |>
  mutate(
    y = x + 0.2 * (`Gas exchange\nthrough:` == treatments[2]) +
      0.5 * (scenario == scenarios[2] &
               `Gas exchange\nthrough:` == treatments[2]),
    yend = xend + 0.2 * (`Gas exchange\nthrough:` == treatments[2]) +
      0.5 * (scenario == scenarios[2] &
               `Gas exchange\nthrough:` == treatments[2])
  )

df_polygon = df_segment |>
  reframe(
    x = case_when(
      `Gas exchange\nthrough:` == treatments[2] ~ x,
      `Gas exchange\nthrough:` == treatments[1] ~ xend
    )
  ) |>
  crossing(`Gas exchange\nthrough:` = treatments, scenario = scenarios) |>
  mutate(
    y = x + 0.2 * (`Gas exchange\nthrough:` == treatments[2]) +
      0.5 * (scenario == scenarios[2] &
               `Gas exchange\nthrough:` == treatments[2]),
    vertex = case_when(
      x == min(x) & `Gas exchange\nthrough:` == treatments[1] ~ 1,
      x == min(x) &
        `Gas exchange\nthrough:` == treatments[2] ~ 2,
      x == max(x) &
        `Gas exchange\nthrough:` == treatments[2] ~ 3,
      x == max(x) &
        `Gas exchange\nthrough:` == treatments[1] ~ 4,
    ),
    xend = NA_real_,
    yend = NA_real_
  ) |>
  arrange(scenario, vertex)

df_path = df_polygon |>
  filter(vertex %in% c(3,4)) |>
  mutate(
    x = x + 0.5,
    `Gas exchange\nthrough:` = NA,
    label = case_when(
      scenario == scenarios[1] ~ "no change",
      scenario == scenarios[2] ~ "significant\ndifference"
    )
  )

df_text = df_path |>
  summarize(x = first(x),
            y = mean(y),
            label = first(label),
            .by = scenario) |>
  mutate(xend = NA_real_,
         yend = NA_real_)

gp_aa = ggplot(df_segment,
               aes(
                 x,
                 y,
                 xend = xend,
                 yend = yend,
                 color = `Gas exchange\nthrough:`
               )) +
  facet_grid(. ~ scenario) +
  xlim(-1, 2.5) +
  geom_polygon(
    data = df_polygon,
    mapping = aes(x = x, y = y),
    alpha = 0.2,
    fill = "black",
    color = "white"
  ) +
  geom_segment(linewidth = 2, lineend = "round") +
  geom_path(data = df_path,
               lineend = "butt",
               linewidth = 1.5) +
  geom_text(
    data = df_text,
    mapping = aes(label = label),
    hjust = -0.1,
    vjust =0.5, 
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "amphi\nboth surfaces" = "steelblue",
      "pseudohypo\nlower surface only" = "tomato"
    ),
    breaks = c("amphi\nboth surfaces", "pseudohypo\nlower surface only")
  ) +
  xlab("Stomatal conductance") +
  ylab("Photosynthesis") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 12),
    title = element_text(size = 12)
  )

gp_aa

ggsave(
  "figures/example-aa.pdf",
  gp_aa,
  width = 6.5,
  height = 3.5,
  device = cairo_pdf
)
