# parallelOfeli

Implementación paralela de la herramienta Ofeli (Open Fast and Eficient Level set Implementation) basada en el trabajo de una aproximación al Level Set original sin resolver ecuaciones diferenciales[1].

Versiones:

************ 1.0 ************

- Añadida: División de la lista Lin y Lout en tiempo lineal cada una: O(n) + O(m). Siendo n y m las longitudes de las dos listas respectivamente.
- Mejora: Paralelización de los bluces de recorrido del las listas.
- Añadida: Recronstrucción de las listas = O(n) + O(m)

******************************


************ 2.0 ************
-Mejora: función size() orden constante.
-


******************************

|1| Y. Shi, W. C. Karl - A real-time algorithm for the approximation of level-set based curve evolution - IEEE Trans. Image Processing, vol. 17, no. 5, May 2008.


