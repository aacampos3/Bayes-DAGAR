# Bayes-DAGAR

En el repositorio contiene el código para realizar un estudio de simulación en modelos espaciales para datos de áreas con el objetivo de analizar la estimación de los parámetros. Particularmente se ajustan el modelo DAGAR y se comparan los resultados con el modelo CAR. Éstos usualmente se especifican de manera jerárquica, por lo que se utiliza un enfoque Bayesiano. Para ajustar modelos bayesianos se utiliza el lenguaje STAN.


## Ruta :file_folder:

El repositorio contiene las siguiente ruta:

* `funciones_dagar.R`: Contiene funciones auxiliares para la implementación de los modelos CAR y DAGAR.

* Normal_CAR: Contiene los archivos para simular datos utilizando un esquema CAR.
    * `car_normal.R`: Código para realizar el esquema de simulación utilizando esquema CAR. 

    * `normal_car.stan`: Código en STAN para ajustar un modelo bayesiano utilizando un modelo CAR. 

* Normal_DAGAR: Contiene archivos para simular datos utilizando un esquema DAGAR.

    * `dagar_normal.R`

    * `dagar_normal_final.R`: Código para ajustar un modelo de regresión con efecto aleatorio con enfoque bayesiano utilizando un modelo DAGAR. 

    * `normal_dagar.stan`: Código en STAN para ajustar un modelo DAGAR.

    * `random_dagar_normal.R`: Código para ajustar un modelo de regresión con efecto aleatorio con enfoque bayesiano utilizando un modelo DAGAR con matriz de adyacencia aleatoria. 

    * `random_norma_dagar.stan`: Código en STAN para ajustar un modelo DAGAR con matriz de adyacencia aleatoria.


## Referencias 📚

Peaper de referencia para el estudio

* Datta, A., Banerjee, S., Hodges, J. S., & Gao, L. (2019). Spatial disease mapping using directed acyclic graph auto-regressive (DAGAR) models. *Bayesian analysis*, 14(4), 1221.

Código adaptado desde

* Datta, A. (s/f). *advanced-spatial-statistics-2021*. Github. Recuperado el 15 de marzo de 2024, de https://github.com/abhirupdatta/advanced-spatial-statistics-2021




