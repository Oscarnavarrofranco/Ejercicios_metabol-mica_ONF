#Creador: Oscar Navarro
#Llamado de bases de datos

#Llamar la base de datos desde la computadora

install.packages("readr")
library(readr)

titulación <- read_csv(file="datos_titulacion (2).csv")

#####################################################
#Llamar desde un respositorio (internet)

repositorio <- read_csv(file="https://raw.githubusercontent.com/ManuelLaraMVZ/titulacion_amino_acidos/main/datos_titulacion%20(2).csv")
head(repositorio)#para ver el encabezado
View(repositorio)#para ver la tabla

###############################################
#gráfica

install.packages("ggplot2")
library(ggplot2)

grafica <- ggplot(repositorio, 
                  aes(x=Volumen,
                  y=pH))+
  geom_line()+
  labs(title="Titulación cisteína",
       x="Volumen ácido (uL)",
       y="Valor pH")+
  theme_dark()

grafica

ggsave("titulacion_repertorio.jpeg", plot = grafica, width = 6, height = 4, dpi = 500)
