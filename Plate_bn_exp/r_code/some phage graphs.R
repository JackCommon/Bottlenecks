OG_monoculture_IDs = c('M.1', 'M.2', 'M.3', 'M.4', 'M.5', 'M.6')
OG_5clone_IDs = c('5.1', '5.2', '5.3', '5.4', '5.5', '5.6')
OG_50clone_IDs = c('50.1', '50.2', '50.3', '50.4', '50.5', '50.6')

replicate_names_legend = c('1', '2', '3', '4', '5', '6')

mono_phage_plot = ggplot(aes(y=pfu, x=timepoint, group=ID), 
                         data=subset(extreme, bottleneck == 'monoculture'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("P.f.u ml"*{}^{-1}*"")))+
  ggtitle('Monoculture')+
  
  #facet_wrap(~bottleneck, labeller = bottleneck_labeller)+
  
  scale_y_continuous(trans = 'log10',
                    breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x)),
                    limits = c(NA,1e+12))+
  
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = OG_monoculture_IDs,
                        labels = replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5'),
                   labels=c('0', '1', '2', '3', '4', '5'))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

mono_phage_plot

fiveclone_phage_plot = ggplot(aes(y=pfu, x=timepoint, group=ID), 
                         data=subset(extreme, bottleneck == '5-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("P.f.u ml"*{}^{-1}*"")))+
  ggtitle('5-clone')+
  
  #facet_wrap(~bottleneck, labeller = bottleneck_labeller)+
  
  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = OG_5clone_IDs,
                        labels = replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5'),
                   labels=c('0', '1', '2', '3', '4', '5'))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

fiveclone_phage_plot

fiftyclone_phage_plot = ggplot(aes(y=pfu, x=timepoint, group=ID), 
                         data=subset(extreme, bottleneck == '50-clone'))+
  
  geom_point(stat='identity', position=pd)+
  geom_path(stat='identity', position=pd)+
  
  labs(x='Days post-infection (d.p.i.)', y=expression(bold("P.f.u ml"*{}^{-1}*"")))+
  ggtitle('50-clone')+
  
  
  #facet_wrap(~bottleneck, labeller = bottleneck_labeller)+
  

  scale_y_continuous(trans = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     limit = c(NA, 1e+12))+
  
  theme_bw()+
  scale_colour_discrete(name='Replicate',
                        breaks = OG_50clone_IDs,
                        labels = replicate_names_legend)+
  theme(plot.title = element_text(face="bold", hjust=0, size = 16))+
  theme(axis.title = element_text(face="bold", size=16))+
  theme(legend.title = element_text(face='bold', size=14))+
  theme(legend.title.align = 0.5)+
  theme(legend.position = 'right')+
  theme(legend.key.width = unit(2, 'cm'))+
  theme(legend.key.height = unit(0.5, 'cm'))+
  theme(strip.text = element_text(face='bold', size=14))+
  
  scale_x_discrete(breaks=c('t0', 't1', 't2', 't3', 't4', 't5'),
                   labels=c('0', '1', '2', '3', '4', '5'))+
  
  theme(axis.text = element_text(size=12))+
  theme(legend.text = element_text(size=12))

fiftyclone_phage_plot

phage_lines_all = ggarrange(mono_phage_plot+labs(y='', x='')+
                              theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), 'pt')), 
                            fiveclone_phage_plot+labs(x='', y='')+
                              theme(plot.margin = unit(c(0, 0, 5.5, 5.5), 'pt')),
                            fiftyclone_phage_plot+labs(y=''),

                            ncol=1, nrow=3,
                            align='h',
                            common.legend = T, legend='right')
quartz()
phage_lines_all

library(cowplot)
mono_phage_plot = mono_phage_plot + theme(plot.margin = unit(c(2,2,0,1), 'pt'))
fiveclone_phage_plot = fiveclone_phage_plot + theme(plot.margin = unit(c(2,2,0,1), 'pt'))
fiftyclone_phage_plot = fiftyclone_phage_plot + theme(plot.margin = unit(c(0,2,1,1), 'pt'))

quartz()

legend = get_legend(mono_phage_plot)
phage_plate_all = plot_grid(mono_phage_plot+labs(x='')+theme(legend.position = 'none'), 
          fiveclone_phage_plot+labs(y='', x='')+theme(legend.position = 'none'), 
          fiftyclone_phage_plot+theme(legend.position = 'none'),
          ncol = 2, nrow=2, align = 'v', axis='l', rel_widths = c(.3,.3,.3,3), rel_heights = c(.3,.3,.3,3))

#test2 = plot_grid(test, legend, rel_widths = c(1, 1))

quartz()
test

ggsave('plate_phage_all.png', test, device = 'png',
       path = '../paper figures/', width=25, height=15, unit=c('cm'), dpi=300)

