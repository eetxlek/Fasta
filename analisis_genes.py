from Bio import SeqIO
import csv
import argparse
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import os 
from pathlib import Path  # para manejar rutas


def graficar_gc(resultados, output_path):
    #grafico de barras (cada elemento dic): id + contenido GC
    nombres = [r["ID"] for r in resultados]  # lista con los nombres de los genes (como id)
    gc_vals = [r["%GC"] for r in resultados] # porcentaje de gc

    plt.figure(figsize=(5, 3)) #inicia figura con tamaño: ancho x alto
    plt.bar(nombres, gc_vals, color="green") # x , y, color
    plt.title("Contenido GC por gen") 
    plt.xlabel("Gen") #nombre eje x
    plt.ylabel("% GC") #nombre eje y
    plt.tight_layout() #ajusta margen para que no se corten nombres
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path)
    print(f"📊 Gráfico guardado como {output_path}")

def traducir_proteina(secuencia):
    secuencia_obj = Seq(secuencia)
    return str(secuencia_obj.translate(to_stop=True))

def calcular_gc(secuencia):
    gc = secuencia.count("G") + secuencia.count("C")
    return round(gc / len(secuencia) * 100, 2)

def contar_codones(secuencia, codones):
    return sum(secuencia[i:i+3] in codones for i in range(0, len(secuencia)-2, 3))

def get_project_root():
    """Devuelve la ruta a la raíz del proyecto (una carpeta arriba del script)"""
    return Path(__file__).parent

def analizar_fasta(fasta_path, output_csv_path, output_image_path):
    resultados = []

    # Crear directorio para imágenes si no existe
    imagenes_dir = "imagenes"
    os.makedirs(imagenes_dir, exist_ok=True)

    for registro in SeqIO.parse(fasta_path, "fasta"):
        seq = str(registro.seq).upper()
        #validacion de secuencia . not seq evalua si esta vacia (true o false), not all (no contiene ningun caracter no valido, distinto a ACGT)
        # si se cumple alguna de las dos (vacio o caracter no valido), será invalidado.
        if not seq or not all(base in "ACGT" for base in seq):
            print(f"⚠️ Secuencia inválida o vacía en {registro.id}, se omite.")
            continue
        # Análisis
        gc_content = calcular_gc(seq)
        inicio_count = contar_codones(seq, {"ATG"})
        paro_count = contar_codones(seq, {"TAA", "TAG", "TGA"})
        proteina = traducir_proteina(seq)
        resultados.append({
            "ID": registro.id,
            "Longitud": len(seq),
            "%GC": gc_content,
            "Codones_inicio_ATG": inicio_count,
            "Codones_paro": paro_count,
            "Proteina_traducida": proteina
        })
        #generar  grafico
        if resultados:
            graficar_gc(resultados, output_image_path)
        else:
            print("⚠️ No se encontraron secuencias válidas para analizar")

    output_csv_path.parent.mkdir(parents=True, exist_ok=True)
    # Guardar resultados en CSV
    with open(output_csv_path, mode='w', newline='', encoding='utf-8') as csv_file:
        campo = ["ID", "Longitud", "%GC", "Codones_inicio_ATG", "Codones_paro", "Proteina_traducida"]
        writer = csv.DictWriter(csv_file, fieldnames=campo)
        writer.writeheader()
        writer.writerows(resultados)

    print(f"Análisis completado. Resultados guardados en {output_csv_path}")

# Ejecutar metodo genera un archivo csv con los resultados como output. Genera un gráfico del contenido GC.
#analizar_fasta("secuencia.fasta", "resultados.csv")  
# #input archivo del proyecto y output resutados.csv FIJOS, no coge args de terminal. 
# Hay que cambiar manualmente code en cada analisis o archivo.


# Añade parámetros para usar cualquier archivo de fuera del proyecto. Para ejecuta por línea de comandos: python analisis_genes.py secuencia.fasta resultados.csv
if __name__ == "__main__":

    project_root = get_project_root()
    parser = argparse.ArgumentParser(description="Análisis básico de secuencias FASTA")
    parser.add_argument("input_fasta", help="Ruta al archivo FASTA de entrada")
    parser.add_argument("--output-csv", help="Ruta al archivo CSV de salida",default="Resultados/analisis.csv")
    parser.add_argument("--output-image", help="Ruta para la imagen del gráfico", default="imagenes/grafica_gc.png")
    # al usar argumento posicionar para la salida del csv, debes ejecutar el script así: 
    # python analisis_genes.py secuencia.fasta --output-csv resultados.csv

    args = parser.parse_args()

    # Convertir rutas a absolutas para mejor manejo
    fasta_path = project_root / "secuencias_ej" / args.input_fasta
    # Construir ruta del CSV - verifica si es rutacon carpeta o solo fichero
    if "/" not in args.output_csv and "\\" not in args.output_csv:
        output_csv_path = project_root / "Resultados" / args.output_csv
    else:
        output_csv_path = project_root / args.output_csv

    output_image_path = project_root / args.output_image
    print(f"🔍 Buscando archivo FASTA en: {fasta_path}")
    try:
        analizar_fasta(fasta_path, output_csv_path, output_image_path)
    except FileNotFoundError as e:
        print(f"❌ Error: {e}")
        print("💡 Asegúrate de que:")
        print("   - El archivo existe en secuencias_ej/")
        print("   - El nombre del archivo es correcto (incluyendo extensión .fasta)")
        print("   - La estructura de carpetas es correcta")
        exit(1)
    except Exception as e:
        print(f"❌ Error inesperado: {e}")
        print(f"💡 Debug info:")
        print(f"   - fasta_path: {fasta_path}")
        print(f"   - output_csv_path: {output_csv_path}")
        print(f"   - output_image_path: {output_image_path}")
        print(f"   - ¿Existe fasta_path? {os.path.exists(fasta_path)}")
        exit(1)
#El argumento puede ser una ruta relativa desde el script o una absoluta del pc.