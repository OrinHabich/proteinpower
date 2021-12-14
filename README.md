# Protein Pow(d)er
This project is a solution to one of the cases in the course `Heuristieken` on the University of
Amsterdam (UvA). See <http://heuristieken.nl/wiki/index.php?title=Heuristieken>. In 2017 I took this course and worked on the project `Amstelhaege'. This was a lot of fun, so I decided to do one of the other cases.


#### Description
This project is about folding proteins. A protein is a sequence of acids of various types. This protein can be folded with 90 degree angles. Some types of acids will form certain  `bonds' if they are close enough to each other (i.e. a distance of 1). This will contribute to the stability of the folding, measured with a score. Goal is to find the best folding for a given protein.


#### Project contents

- /images - Some images of previous results
- /src - The main code for the project
- /tests - The test suite for this project, i.e. a set of unittests using Pytest
- .git_original - The git history before I tried to push this project to Github. This gave some problems (due to 2FA). I solved this by starting a fresh git repo.
- presentation_video.mkv - A video in which I introduce the project and demonstrate the code.
- requirements.txt - The requirements needed for this project.

## Getting Started

- Download or clone this project.
- cd into the directory and install the requirements with the following line of code:
```
pip install -r requirements.txt
```
Consider using a virtual environment.

### Prerequisites

- Python 3.8.x

The following Python packages are used:
  - Math
  - Matplotlib
  - Numpy
  - Random
  - Sys
  - Itertools

### Run the code

The code can be executed with the following command:
```
python3 src/main.py
```

## Built With

- [Sublime Text] - Editor
- [Python 3.8](https://docs.python.org/3/) - Python3


## Versioning

- v1.0 (30-04-2021)


## Authors

**Name:** Orin Habich
