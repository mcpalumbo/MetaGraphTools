from setuptools import setup, find_packages

setup(
    name='MetaGraphTools',
    version='0.1.0',
    description='Python tool is designed to analyze and visualize metabolic networks',
    author='mcpalumbo',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'networkx',
        'cobra',
        'pyvis'
    ],
    entry_points={  
        'console_scripts': [
            'MetaGraphTools = MetaGraphTools.main:main', 
        ],
    },
    python_requires='>=3.11,<3.12',
    classifiers=[
        'Programming Language :: Python :: 3.11',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
