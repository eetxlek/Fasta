🧬 FASTA Analizador de secuencia genética  
Herramienta Python para análisis básico de secuencias genéticas en formato FASTA.

🚀 Uso Rápido
bash

# Instalar dependencias
pip install biopython matplotlib

# Ejecutar análisis (archivo en secuencias_ej/)
python analisis_genes.py mi_secuencia.fasta

📊 Resultados
CSV con métricas: resultados/analisis.csv

Gráfico GC: imagenes/grafica_gc.png

🔧 Funciones
Cálculo de contenido GC
Detección de codones inicio/paro
Traducción a proteína
Visualización automática

📁 Estructura

proyecto/
├── analisis_genes.py       # Script principal  
├── secuencias_ej/          # Archivos FASTA de entrada  
├── resultados/             # CSV con métricas generadas  
├── imagenes/               # Gráficos de análisis  
└── test/                   # Pruebas unitarias con Pytest  

⚡ Comandos Útiles
bash
# Análisis con rutas personalizadas
python analisis_genes.py secuencia.fasta --output-csv resultados/mi_analisis.csv --output-image imagenes/mi_grafico.png

✅ Testing

El proyecto incluye pruebas automatizadas con pytest para asegurar la calidad del análisis.

📂 Estructura de tests
test/  
└── test_analisis.py       # Tests unitarios y de integración  

🧪 Cobertura de tests

test_calcular_gc: Verifica el cálculo correcto del contenido GC.
test_traducir_proteina: Traduce una secuencia de nucleótidos a aminoácidos y verifica resultados esperados.
test_contar_codones: Cuenta codones de inicio/parada correctamente en distintas secuencias.
test_analisis_completo: Test de integración que verifica la ejecución completa del análisis y la generación de archivos.
test_secuencia_invalida: Placeholder para manejar casos de secuencias con caracteres inválidos. Puedes completarlo según la lógica de validación que hayas implementado.

▶️ Ejecutar tests

Desde la raíz del proyecto:
pytest test/test_analisis.py -v

Asegúrate de tener instalado pytest:
pip install pytest

¡Listo para usar! Coloca tus archivos FASTA en secuencias_ej/ y ejecuta.
