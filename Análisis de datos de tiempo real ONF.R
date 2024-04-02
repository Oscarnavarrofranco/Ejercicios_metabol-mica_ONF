
# Creador: Oscar Navarro Franco
#Análisis de datos de tiempo real

install.packages("pacman") #pacman sirve para llamar e instalar a otros paquetes

library(pacman)

p_load("readr",#sirve para llamar a las bases de datos
       "dplyr")# sirve para facilitar el manejo de bases de datos

##############################################

#Llamar la base de datos

datos_pcr<- read_csv(file="https://raw.githubusercontent.com/ManuelLaraMVZ/Transcript-mica/main/Genes.csv")

head(datos_pcr)
#################################################

#Obtención de genes de referencia y de interés

actina <-  datos_pcr %>% 
  slice(1)

actina

genes_interes <- datos_pcr %>% 
  slice(2,3,4,5,6,7)

genes_interes

##################################################
#Promediar los controles y los tratamientos

promedio_actina <-actina %>% 
  mutate(Mean_Cx=(C1+C2+C3)/3) %>% 
  mutate(Mean_Tx=(T1+T2+T3)/3) %>% 
  select(Gen, Mean_Cx, Mean_Tx)

promedio_actina

promedio_GI <- genes_interes %>% 
  mutate(Mean_Cx=(C1+C2+C3)/3) %>% 
  mutate(Mean_Tx=(T1+T2+T3)/3) %>% 
  select(Gen, Mean_Cx, Mean_Tx)

promedio_GI

###############################################

#Análisis DCT

DCT <- promedio_GI %>% 
  mutate(DCT_Cx=(Mean_Cx-promedio_actina$Mean_Cx)) %>% 
  mutate(DCT_Tx=(Mean_Tx-promedio_actina$Mean_Tx)) %>% 
  select(Gen, DCT_Cx, DCT_Tx)

DCT
####################################################

#Análsisi DDCT

DDCT <- DCT %>% 
  mutate(DDCT=(DCT_Tx-DCT_Cx)) %>% 
  mutate("2^DDCT"=(2^(-DDCT)))

DDCT

################################

#Guardar Tabla

write.csv(DDCT, "2DDCT.csv", row.names=FALSE)
