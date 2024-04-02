#Creador: Oscar Navarro Franco
#Gráfica de dispersión 2
#Se realiza una serie de experimentos comparativos de cerebros 
#de tres animales destetados cuyas madres estuvieron expuestas a 
#una dieta restrictivas (desnutridas) y se compararon con tres animales 
#cuyas madres se expusieron a una dieta ad libitum (eunutridas). Se desea 
#saber: ¿Qué miRNAs fueron los que más cambiaron entre los animales 
#desnutridos (tratamiento) comparados con los animales eunutridos (control)?

library("pacman")

p_load("readr",
       "ggplot2",
       "dplyr")

datos <- read_csv(file="https://raw.githubusercontent.com/ManuelLaraMVZ/Transcript-mica/main/DesnutridasvsEunutridas.csv")

head(datos)

#extracción de genes controles (referencia)

Controles <- datos %>% 
  filter(Condicion=="Control")

head(Controles)

#Sacar promedios

promedio_controles <- Controles %>% 
  summarise(Mean_C1= mean(Cx1),
            Mean_C2= mean(Cx2),
            Mean_C3= mean(Cx3),
            Mean_T1= mean(T1),
            Mean_T2= mean(T2),
            Mean_T3= mean(T3)) %>% 
  mutate(Gen="Promedio_controles") %>% 
  select(7,1,2,3,4,5,6)

promedio_controles

#####################################################
#Extraer los genes de la tabla "datos"#

genes <- datos %>% 
  filter(Condicion=="Target") %>% 
  select(-2)
head(genes)
##############

#Sacar el 2^DCT

DCT <- genes %>% 
  mutate(DCT_C1=2^-(Cx1-promedio_controles$Mean_C1),
         DCT_C2=2^-(Cx2-promedio_controles$Mean_C2),
         DCT_C3=2^-(Cx3-promedio_controles$Mean_C3),
         DCT_T1=2^-(T1-promedio_controles$Mean_T1),
         DCT_T2=2^-(T2-promedio_controles$Mean_T2),
         DCT_T3=2^-(T3-promedio_controles$Mean_T3)) %>% 
  select(-2,-3,-4,-5,-6,-7)

DCT

promedio_genes <- DCT %>% 
  mutate(Mean_DCT_Cx=(DCT_C1+DCT_C2+DCT_C3)/3,
         Mean_DCT_Tx=(DCT_T1+DCT_T2+DCT_T3)/3)

promedio_genes

grafica_dispersion <- ggplot(promedio_genes,
                             mapping = aes(x= Mean_DCT_Cx,
                                           y= promedio_genes$Mean_DCT_Tx))+
  geom_point()

grafica_dispersion

grafica_dispersion2 <- grafica_dispersion +
  labs(
    title = "Condición control vs tratamiento",
    caption = "Creador: Oscar Navarro",
    x = expression("Control " * 2^-DCT),
    y = expression("Tratamiento " * 2^-DCT)
  ) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", color = "black") +
  theme_minimal() +
  theme(
    axis.title = element_text(family = "Times New Roman", size = 12),
    axis.text = element_text(family = "Times New Roman", size = 12, face= "bold"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    legend.text = element_text(family = "Times New Roman", size = 14),
    legend.title = element_text(family = "Times New Roman", size = 14)
  )
grafica_dispersion2

###################################################################
#Identificación de los genes
head(promedio_genes)

top_10 <- promedio_genes %>% 
  select(1,8,9) %>% 
  top_n(10, Mean_DCT_Cx) %>% 
  arrange(desc(Mean_DCT_Cx))

head(top_10)

grafica_dispersion3 <- grafica_dispersion2 +
  geom_label_repel(data=top_10,mapping = 
                     aes(x = Mean_DCT_Cx, y = Mean_DCT_Tx, label = Gen),
                   label.padding = unit(0.2, "lines")
  )
grafica_dispersion3

ggsave("Dispersión2.jpeg", plot = grafica_dispersion3, height = 5, width = 7, dpi=300)
