import pytest
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) #busca en directorio padre, no solo test
from analisis_genes import calcular_gc, traducir_proteina, contar_codones

# Añade el directorio src al path para poder importar el módulo
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

def test_calcular_gc():
    """Test para la función calcular_gc con secuencias conocidas"""
    # Secuencia con 50% GC (2 de 4 bases)
    secuencia = "ATGC"
    assert calcular_gc(secuencia) == 50.0
    
    # Secuencia con 100% GC
    secuencia = "GGCC"
    assert calcular_gc(secuencia) == 100.0
    
    # Secuencia con 0% GC
    secuencia = "ATTA"
    assert calcular_gc(secuencia) == 0.0
    
    # Secuencia vacía
    secuencia = ""
    assert calcular_gc(secuencia) == 0.0

def test_traducir_proteina():
    """Test para la función traducir_proteina con secuencia conocida"""
    # Secuencia que codifica para la proteína "MAT"
    secuencia = "ATGGCTACT"  # ATG -> M, GCT -> A, ACT -> T
    assert traducir_proteina(secuencia) == "MAT"
    
    # Secuencia con codón de paro temprano
    secuencia = "ATGTAA"  # ATG -> M, TAA -> (stop)
    assert traducir_proteina(secuencia) == "M"
    
    # Secuencia vacía
    secuencia = ""
    assert traducir_proteina(secuencia) == ""

def test_contar_codones():
    """Test para la función contar_codones"""
    secuencia = "ATGTAGATGTGAATG"
    # Contar codones de inicio ATG
    assert contar_codones(secuencia, {"ATG"}) == 2
    
    # Contar codones de paro
    codones_paro = {"TAA", "TAG", "TGA"}
    assert contar_codones(secuencia, codones_paro) == 2
    
    # Secuencia más corta que un codón
    secuencia = "AT"
    assert contar_codones(secuencia, {"ATG"}) == 0

def test_analisis_completo():
    """Test de integración: análisis completo de una secuencia pequeña"""
    # Crear un archivo FASTA de prueba temporal
    test_fasta = "secuencias_test/test_sequence.fasta"
    test_csv = "secuencias_test/resultados_test.csv"
    
    # Ejecutar la función principal
    from analisis_genes import analizar_fasta
    analizar_fasta(test_fasta, test_csv)
    
    # Verificar que el archivo CSV se creó y tiene contenido
    assert os.path.exists(test_csv)
    
    # Leer y verificar resultados
    with open(test_csv, 'r') as f:
        lines = f.readlines()
        assert len(lines) == 2  # Header + 1 línea de datos
        
        # Verificar contenido
        data = lines[1].split(',')
        assert data[0] == 'gen_test'
        assert float(data[2]) == 50.0  # %GC
        assert data[3] == '1'  # Codones inicio
        assert data[4] == '1'  # Codones paro
        assert data[5].strip() == 'MA'  # Proteína
    
    # Limpiar: eliminar archivo de resultados de test
    if os.path.exists(test_csv):
        os.remove(test_csv)

# Para manejar secuencias inválidas
def test_secuencia_invalida():
    """Test que verifica el manejo de secuencias inválidas"""
    # Secuencia con caracteres no válidos
    secuencia = "ATGNXYZ"
    # La función calcular_gc debería manejarlo
    # (depende de tu implementación actual)
    # Puedes modificar el test según tu manejo de errores
    pass

if __name__ == "__main__":
    pytest.main([__file__, "-v"])