# ProyectoTT1

Este proyecto corresponde a la asignatura de  Taller Transversal 1 de la carrera de Ingeniería Informática.
Se trata de la traducción a **C++** del código para la determinación de órbita desde **MatLab**.

---

## 🔧 Instrucciones de compilación

### Compilar la aplicación principal:

```bash
g++ tests/EKF_GEOS3.cpp src/*.cpp -lm -std=c++23 -o bin/EKF_GEOS3.exe
```
### Compilar test unitarios:
```bash
g++ tests/tests.cpp src/*.cpp -lm -std=c++23 -o bin/tests.exe
```