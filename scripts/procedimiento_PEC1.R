
#Para construir un objeto `SummarizedExperiment` usando los datasets seleccionados:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")


################################################################################

#1) Cargar las bibliotecas necesarias
library(BiocManager)
library(SummarizedExperiment)
library(tidyverse)
library(factoextra)
library(DMwR2)
library(naniar)

################################################################################

#CREACIÓN DEL OBJETO SUMMARYZEDEXPERIMENT

#2) Carga de datos
#Hacemos la lectura del archivo general
assay_1 <- read.csv("data/DataValues_S013.csv", row.names = 1) 

#Podemos ver que las primeras 5 columnas son, realmente, metadatos de las muestras así que los guardaremos para colData
coldata <- assay_1[,1:5]

#Eliminamos dichas columnas de la matriz de expresión.
assay_1[,1:5] <- NULL

#La matriz de expresión debe tener las muestras en columnas, por lo que transponemos y dejamos en formato matriz
assay_1 <- t(assay_1)

#Leemos el siguiente archivo
rowdata <- read.csv("data/DataInfo_S013.csv", row.names = 1) 

#Por el aspecto del anterior, parecen ser los metadatos de las features, así que eliminamos las definiciones de los
#metadatos para que las dimensiones coincidan.
rowdata <- rowdata[-(1:5), ]


#3) Creamos el objeto SummarizedExperiment:
se <- SummarizedExperiment(assays = list(counts = as.matrix(assay_1)),
                           rowData = rowdata,
                           colData = coldata)

#Vemos características del objeto creado
se

#3.1) Para crear y ver los metadatos
metadata(se)$description <- "Este es un experimento de metabolómica"
metadata(se)$fecha <- Sys.Date()  # Fecha actual
metadata(se)$title <- "Metabotypes of response to bariatric surgery independent of the magnitude of weight loss"

#OTras características que se pueden observar en el objeto creado
assay(se) #Ver la matriz de expresión
rowData(se) #Ver los metadatos de las features
colData(se) #Ver los metadatos de las muestras
assays(se)$counts #Otra forma de ver la matriz de expresión en el caso de que hubieran varios experimentos/assays
metadata(se) #Ver los metadatos asociados al experimento

#Para el informe el repositorio, guardamos el objeto .Rda
save(se, file="se.Rda")

################################################################################

#ANÁLISIS EXPLORATORIO DE LOS DATOS

#Para mostrar los nombres de los ensayos en el objeto SummarizedExperiment podemos usar:
assayNames(se)

#Colocamos el nombre del ensayo: "counts" para extraer la matriz de expresión
data_matrix <- assay(se, "counts")  

#ESTADÍSTICA DESCRIPTIVA
#Pero usaremos kableExtra en el pdf por los problemas de visualización de este comando.
summary(data_matrix)

#ANÁLISIS GRÁFICO
#Para generar un estudio gráfico de las variables de interés: necesitamos convertir la matriz a numérica
data_matrix_clean <- apply(data_matrix, 2, function(x) as.numeric(trimws(x)))

#Esta función está tomada de los ejemplos del Curso:
f <- function(x) {
  x <- na.omit(x)
  if (is.numeric(x)) {
    hist(x, breaks = 10, main = "Histograma de datos numéricos")
  } else {
    barplot(table(x), main = "Gráfico de barras de datos categóricos")
  }
}

# Configurar par() para mostrar múltiples gráficos
par(mfrow = c(3, 2))

# Aplicar la función a cada columna
apply(data_matrix_clean, 2, f)
dev.off()


#BOXPLOTS Y DETECCIÓN DE OUTLIERS
boxplot(data_matrix_clean, outline = TRUE, main = "Distribución de las variables", las = 2)

#ALGUNAS CORRELACIONES
#Matriz de correlación y visualización
corr_matrix <- cor(data_matrix_clean, use = "complete.obs")
heatmap(corr_matrix, main = "Mapa de calor de las correlaciones")

################################################################################

#PREPROCESADO DE LOS DATOS: MANEJO DE NAs
#Convertir a dataframe para ver la distribución de NAs y posterior imputación
data_matrix_clean <- as.data.frame(data_matrix_clean)
gg_miss_var(data_matrix_clean, show_pct = TRUE)

#Aplicar la imputación kNN siguiendo lo que indican en el paper original
data_matrix_clean <- knnImputation(data_matrix_clean, k = 10)

################################################################################

#REDUCCIÓN DE LA DIMENSIONALIDAD: PCA

#Realizar PCA centrando y escalando los datos
pca_result <- prcomp(t(data_matrix_clean), center = TRUE, scale. = TRUE)

# Resumen de los resultados del PCA
summary(pca_result)
plot(pca_result) #Efectivamente el primer componente es el que más varianza explica

################################################################################

#Para poder saber qué variables "forman parte" del primer componente y ver si se puede 
#generar un "patrón" para explicar el PC1.

#A) Crear el data frame de cargas de PC1 en orden descendente
loadings_pc1 <- as.data.frame(sort(pca_result$rotation[, "PC1"], decreasing = TRUE))

#B) Obtener los nombres de las features desde rowData(se)
feature_names <- rownames(rowData(se))

#C) Ordenar los nombres de las features según el orden de PC1 y asignarlos a loadings_pc1
orden_pc1 <- order(pca_result$rotation[, "PC1"], decreasing = TRUE)
rownames(loadings_pc1) <- feature_names[orden_pc1]

#Conclusión: hay demasiadas variables. Podríamos quedarnos con algunas de cara al informe.

################################################################################

#Ver las cargas de cada variable en cada componente principal
pca_result$rotation

#Convertir las cargas a un data frame
loadings_df <- as.data.frame(pca_result$rotation)

#Graficar el biplot del PCA
biplot(pca_result, main = "PCA Biplot")
#Conclusión: demasiadas variables para poder observar

#Extraer los puntajes de los componentes principales
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- rownames(pca_df)

#Graficar PC1 vs PC2
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point() +
  geom_text(nudge_y = 0.2) +
  labs(title = "PCA", x = "PC1", y = "PC2") +
  theme_minimal()

################################################################################

#HEATMAP DE DISTANCIAS
manDist <- dist(t(data_matrix_clean), method = "manhattan")

#Convertir la matriz de distancias a una matriz completa
manDist_matrix <- as.matrix(manDist)

#Generar el heatmap
heatmap(manDist_matrix, col = heat.colors(16), scale = "none")

################################################################################
#CLUSTER JERÁRQUICO
dist_matrix <- dist(t(data_matrix), method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")
plot(hc, main = "Dendrograma de clustering jerárquico")