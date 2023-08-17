from setuptools import setup, find_packages

setup(
    name='cellsegmentationtracker',
    version="0.1.0",
    author='Simon Guldager Andersen',
    author_email='guldager.simon@gmail.com',
    url='https://github.com/simonguld/CellSegmentationTracker',
    description='Package for cell segmentation, tracking and subsequent (biophysical) statistical analysis',
    install_requires=["tifftools", "seaborn", "scikit-image"],
    packages=find_packages(include=['cellsegmentationtracker', 'cellsegmentationtracker.*']),

    classifiers=[
        'Development Status :: 1 - Planning',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.9',
    ]

)