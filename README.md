# Modelo--Momento-Curvatura
Aplicación en Python para obtener curvas Momento-Curvatura (M-ϕ) en columnas de concreto reforzado utilizando el modelo de Mander para concreto confinado y no confinado. El programa considera de manera explícita la diferencia entre el concreto confinado del núcleo estructural y el concreto no confinado del recubrimiento, empleando para ello el modelo constitutivo propuesto por Mander, Priestley y Park (1988).

La aplicación permite analizar tres configuraciones de columnas comúnmente utilizadas en la práctica estructural:
Columnas circulares con espiral.
Columnas circulares con estribos circulares.
Columnas rectangulares con estribos rectangulares.

A partir de las propiedades geométricas de la sección, las características del acero transversal y longitudinal, las propiedades mecánicas de los materiales y la carga axial aplicada, el programa determina las propiedades del concreto confinado y genera la correspondiente relación Momento-Curvatura. Asimismo, se incorpora el efecto de la tasa de deformación para evaluar la respuesta dinámica del material bajo acciones sísmicas.

Este proyecto representa una extensión del repositorio previamente desarrollado para la implementación del Modelo de Mander, integrando dicho modelo dentro de un procedimiento completo de análisis Momento-Curvatura orientado a la evaluación estructural de columnas de concreto reforzado.
