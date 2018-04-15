##### FIGURE 3

library(cowplot)

bottom_row_fig3 = plot_grid(total_just_plot_full_bn+theme(plot.margin = unit(c(5.5, 35, 5.5, 5.5), 'pt')), 
                            evo_phage_sum_plot_t3,
                            labels = c('B', 'C'))
quartz()
bottom_row_fig3

fig_3 = plot_grid(pheno_plot_2+theme(plot.margin = unit(c(.5, .5, .8, .5), 'cm')), bottom_row_fig3,
                  ncol = 1, labels = c('A', ''),
                  rel_widths = c(1,1.3))

quartz()
fig_3

ggsave('Fig 3.png', fig_3, device='png',
       path = '../paper figures/', width = 25, height = 20, unit=c('cm'), dpi=300)

##### FIGURE 5

bottom_row_fig5 = plot_grid(total_just_plot_phage_bn+theme(plot.margin = unit(c(5.5, 35, 5.5, 5.5), 'pt')),
                            t3_facet_plot,
                            labels = c('B', 'C'))
fig_5 = plot_grid(pheno_plot_2+theme(plot.margin = unit(c(.5, .5, .8, .5), 'cm')),
          bottom_row_fig5,
          ncol = 1, labels = c('A', ''),
          rel_widths = c(1,1.3))
quartz()
fig_5

ggsave('Fig 5.png', fig_5, device='png',
       path = '../paper figures/', width = 25, height = 20, unit=c('cm'), dpi=300)

##### FIGURE 7

bottom_row_fig7 = plot_grid(total_just_plot_phage_bn+theme(plot.margin = unit(c(5.5, 35, 5.5, 5.5), 'pt')),
                            #t3_facet_plot,
                            labels = c('B', 'C'))
fig_7 = plot_grid(pheno_plot_2+theme(plot.margin = unit(c(.5, .5, .8, .5), 'cm')),
                  bottom_row_fig7,
                  ncol = 2, labels = c('A', ''),
                  rel_widths = c(1,1))
quartz()
fig_7

ggsave('Fig 7.png', fig_7, device='png',
       path = '../paper figures/', width = 15, height = 20, unit=c('cm'), dpi=300)
  