<!--
SPDX-FileCopyrightText: 2025 Thomas Förster
SPDX-License-Identifier: CC-BY-4.0
-->

<h1 align="center">Properties Package</h1>
<p align=center> 
Practice Repo for demonstrating OOP for IU class Python DLMDWPMP01.
</p>


<!--===============-->
<!--=== Chapter ===-->
<!--===============-->
# About the Project

The package creates templates properties for gas,liquid and solid substances.

<!--===============-->
<!--=== Chapter ===-->
<!--===============-->
# Getting Started

## Prerequisites

- To work with the package you need to have at least Python 3.10 installed.

- This project uses the Python dependency manager `poetry` [(Installation guide for poetry)][poetry-install]

    Example command line for Linux/Ubuntu:

    ``` bash
    curl -sSL https://install.python-poetry.org | python3 -
    ```
    
    or

    ``` bash
    pipx install poetry
    ```

## Installation

### Using the environment distributed with this model

**Step 1. Clone the repository.**

- via ssh (first, create [SSH key][ssh-key])

    ``` bash
    git clone --depth=1 --branch=main git@github.com:tfoerst3r/properties.git
    ``` 

Alternatively, you can download source code directly via the download option on top of the repository page.


**Step 2. Install required dependencies.**

To build the package run the following in the root folder of the project:

```
poetry install
```

### Using any other python environments

Poetry can build, so called wheels, which can be installed in any other python environment.

**Step 1. Create the wheel or download the wheel**

You can download any wheel via any desired CI job called "build_wheel_package". Or build from scratch via the following step:

```
poetry build --format wheel
```

This will create a wheel file `*.whl` in the folder `./dist`.

**Step 2. Installing the wheel**

- First you need to activate your desired python environment, e.g. `$ source /path/to/your/env/bin/activate`. And then install your the wheel.

``` bash
pip install NAME.whl
```

## Uninstall the module:

For uninstalling `properties` in your current environment use:

``` bash
pip uninstall properties
```

<!--===============-->
<!--=== Chapter ===-->
<!--===============-->
# Usage

## Prerequisites

First you need to activate the environment you want to work in. Either use

- I recommend using `poetry shell`. Beforehand you need to install the `poetry-plugin-shell`

    ``` bash
    poetry self add poetry-plugin-shell
    ```
    
    Then you can use the local `.venv` pyton environment

    ```
    poetry shell
    ```

- or any other environment where you installed the provided wheel

    ```
    source /path/to/your/env/bin/activate
    ```

## Usage 

This is a python package which is used via the `import` command.

The following **Syntax** creates default instances.

``` python
import properties

mygas = properties.GasProp()
myliquid = properties.LiquidProp()
mysolid = properties.SolidProp()
```

<!--===============-->
<!--=== Chapter ===-->
<!--===============-->
# Contributing

**We welcome your contribution!**

The repository is still under development and any feedback, suggestions, technical contributions are highly welcome.

General Options:

- open an issue, if you have improvement ideas
- fork the project and contribute via merge request against the main branch of this repository


<!--===============-->
<!--=== Chapter ===-->
<!--===============-->
# License

Please see the file [LICENSE.md](./LICENSE.md) for further information about how the content is licensed.

<!---- Links ---->
[ssh-key]: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent
[poetry-install]: https://python-poetry.org/docs/

