# parallelOfeli

Implementación paralela de la herramienta Ofeli (Open Fast and Eficient Level set Implementation) basada en el trabajo de una aproximación al Level Set original sin resolver ecuaciones diferenciales[1].

Versiones:

************ 1.0 ************

- Añadida: División de la lista Lin y Lout en tiempo lineal cada una: O(n) + O(m). Siendo n y m las longitudes de las dos listas respectivamente.
- Mejora: Paralelización de los bluces de recorrido del las listas.
- Añadida: Recronstrucción de las listas = O(n) + O(m)

******************************


************ 2.0 ************
- Mejora: Las sublistas ahora son listas ligadas de la clase linked_list del propio proyecto, lo que da pie a modificaciones internas.
- Cambio drástico de la implementación de linked_list:
	- Eliminación del memory_pool ya que daba problemas en la implementación
	- Se ha insertado en cada nodo un puntero hacia el elemento anterior, es decir, se ha creado una lista doblemente ligada. 
	- Mejora: La lista doblemente ligada ha permitido crear un apuntador al final de la lista (y que las funciones de esta cambien fácilmente ese apuntador).
	- Se ha eliminado la función splice front por dar problemas a la hora de llamarse a esta mezclando dos listas.
	- Mejora: Con la anterior eliminación se ha podido cambiar la función size() de manera que sea de orden constante
	- Eliminadas varias funciones redundantes.
- Mejora: Con el apuntador al final de la lista, la función collectList ahora pasa a ser de orden constante.
- Mejora: Ahora se aprovechan los nodos creados, en vez de crear nuevos a la hora de construir las listas.
******************************


************ 3.0 ************
Cambios posibles:
	- Convertir función splitList en orden constante

Nuevo:	
	- Implementación de calculo de recubrimiento.
	- ...
******************************

|1| Y. Shi, W. C. Karl - A real-time algorithm for the approximation of level-set based curve evolution - IEEE Trans. Image Processing, vol. 17, no. 5, May 2008.


