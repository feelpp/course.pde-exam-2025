# Examen FEM – Template Git

Ce dépôt contient le code et les tests pour l’examen de Méthode des Éléments Finis.

## Prérequis
- CMake ≥ 3.21
- Compilateur C++20 avec support Eigen

## Compilation
```bash
cmake --preset default
cmake --build --preset default
```

## Exécution des tests

Les tests unitaires utilisent Boost.UT et attendent un dossier de données `src/data/` :
```bash
build/default/src/fem_tests src/data/
```

## Notebook de l’examen
Le notebook Jupyter **`Examen.ipynb`** contient l’énoncé interactif et des cellules à compléter.


