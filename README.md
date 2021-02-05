# TNBC_research

## El código principal 
Se encuentra en el archivo con extensión .Rproj en formato Rmarkdown

## Data frames
Se encuentran en carpetas zip, excepto la que contiene los datos RNASeq de todos los pacientes, ya que excede la capacidad de 25MB de github, pero se encuentra en el drive del projecto: [https://drive.google.com/drive/folders/1jslFtYwAlzpFnOmgiztPSmCKdHo6nvm5?usp=sharing]

**Actualización:** 
- Ya están filtrados los id's de los pacientes con TNBC con base en los datos clínicos
- Ya están transpuestos los data frames dónde *row = datos, col = variable o campo*
- **Revisar:** Parece que hay un error al momento de aplicar el filtro para seleccionar solamente las filas del data frame RNASeq que corresponden con el patrón de su id de paciente, al parecer selecciona todas las filas como si cumplieran la condición. **NOTA: El id del paciente solamente cuenta con los primeros 9 caracteres, mientras que el id dentro del data frame de RNASeq cuenta con más (ambos comienzan con la misma secuencia de caracteres).**
