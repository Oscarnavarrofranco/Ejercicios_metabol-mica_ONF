2+2
2^(-4)
Suma<-2+2
Suma
2+Suma
CHO <- 535000
VERO <- 1350000
CHO+VERO
TOTAL <- CHO+VERO
TOTAL #No se va a leer esto
y <- -22
x <- 35
w <- "34"
z <- "25"
x+w
V1 <- c(1,2,3)
V2 <- c(4,5,6)
V3 <- c(7,8,9)
V4 <- c("M","D","z")
DF_v <- data.frame(V1,V2,V3,V4)
View(DF_v)
install.packages("readr")
library(readr)
library(readr)
datos_titulacion_2_ <- read_csv("datos_titulacion (2).csv")
View(datos_titulacion_2_)
Titulacion <- read_csv(file="https://raw.githubusercontent.com/ManuelLaraMVZ/titulacion_amino_acidos/main/datos_titulacion%20(2).csv")
head(Titulacion)
install.packages("ggplot2")
library(ggplot2)
grafica <- ggplot(Titulacion,aes(x=Volumen,y=pH ))+geom_point()
grafica
