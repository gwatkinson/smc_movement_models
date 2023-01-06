# Estimating behavioral parameters in animal movement models using a state-augmented particle filter

Projet pour le cours **Hidden Markov models and Sequential Monte-Carlo methods** de l'ENSAE.

* Gabriel Watkinson
* Gabriel Guaquiere
* Jérémie Stym-Popper
* Benjamin Maurel

## Introduction

Dans ce projet, nous implementons le modèle décrit dans le papier [Estimating behavioral parameters in animal movement models using a state-augmented particle filter](https://dalspace.library.dal.ca/bitstream/handle/10222/33464/Dowd_et_al-2011-Ecology.pdf), et l'approfondissons avec des méthodes plus récentes.

Notre rapport est disponible [ici](https://github.com/gwatkinson/smc_movement_models/blob/main/SMC_Movement_Model_Ecology.pdf) à la racine du projet.

### Le package `smc_movement_models`

Les fonctions sont définies dans le package `smc_movement_models`.

Les modules `models.py` et `models_SMC2.py` contiennent les modèles de mouvement et l'algorithme du papier et une version SMC2 respectivement.

Le module `plot.py` contient des fonctions utilitaires pour génerer des graphes.

Le module `process_data.py` contient des fonctions utilitaires pour traiter les données brutes et une CLI.

### Les notebooks

Les notebooks sont dans le dossier `notebooks` et permettent de lancer les expériences et de visualiser les résultats.

Le notebook `notebooks/implementation_papier.ipynb` contient les expériences avec l'algorithme du papier.

Le notebook `notebooks/implementation_smc2.ipynb` contient les expériences avec la version SMC2.

### Le reste

Le dossier `data` contient les données brutes et les données traitées.

Le dossier `images` contient les images utilisées dans le rapport.

Le dossier `paper` contient les papiers liés au projet.

Le dossier `r_code` contient le code R original du papier (pas utilisé par notre groupe).

Les fichiers `poetry.lock` et `pyproject.toml` sont utilisés par [poetry](https://python-poetry.org/) pour gérer les dépendances.

Le fichier `requirements.txt` est utilisé par `pip` pour gérer les dépendances et est généré par poetry.

Le dossier `latex` contient le LaTeX utilisé pour le rapport.

## Installation

Nous utilisons Python pour simuler les données et mettre en place les modèles, notamment la librairie [particles](https://github.com/nchopin/particles).

Avant toutes choses, il faut se déplacer dans le dossier `smc_movement_models`:

```bash
cd /path/to/smc_movement_models
```

### Avec [`poetry`](https://python-poetry.org/)

```bash
# Creation d'un environement virtuel et installation des packages
poetry install

# Activation de l'environement
poetry shell  # sub shell
# OU source $(poetry env info --path)/bin/activate  # activer l'environement dans le shell actuel
```

### Avec `pip`

```bash
# Creation d'un environement virtuel
python -m venv .venv

# Activation de l'environement
.venv/Script/activate  # pour Windows
# OU source .venv/bin/activate  # pour Linux / MacOS

# Installation des packages
pip install -r requirements.txt
```

### Processing des données

Pour générer les données à partir du fichier brut, il faut lancer le script `process_data.py` :

```python
python smc_movement_models/process_data.py
# Pour plus d'informations
# python smc_movement_models/process_data.py --help
```

### `pre-commit`

Pour activer les pre-commit qui formattent le code avant chaque commit :

```bash
pre-commit install
pre-commit run --all-files  # Pour installer les modules et lancer les tests
```

![Exemple de pre-commit](images/pre-commit-exemple.png)

Pour forcer un commit sans vérifier :

```bash
git commit -m "..." --no-verify
```
