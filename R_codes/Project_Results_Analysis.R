library(ggpubr)

#########################################################################
############################ Set Directories ############################
#########################################################################

#Saving Directory
s_dir = '../results'
dir.create(s_dir, showWarnings = F)

df = read.csv(file.path(s_dir, 'All_Kriging_Scenarios_DF.csv'))

nutrients = unique(df$Nutrient)
nutrient = nutrients[3]
for (nutrient in nutrients){
  print(nutrient)

  df_s = df[df$Nutrient == nutrient,]

  df_s = df_s[df_s$Scenario %in% c("UK_Max_5", "OK_Max_10"),]
  df_s$res = df_s$Predicted - df_s$Observed

  box_plot_ = ggboxplot(df_s, x = "Scenario", y = 'res',
                        color = "Scenario",
                        ylab = 'SErr', xlab = "Scenario")

  res.aov <- aov(res ~ Scenario , data = df_s)
  summary(res.aov)
  res.aov2 = oneway.test(V1 ~ SM, data = df_s, var.equal=FALSE)
  pr = res.aov2$p.value

  out = HSD.test(res.aov, 'SM', alpha = 0.05, group=TRUE, unbalanced=T)

}
