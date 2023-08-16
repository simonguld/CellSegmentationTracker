from setuptools import setup, find_packages

setup(
    name='cellsegmentationtracker',
    version=__version__,
    author='Simon Guldager Andersen'
    author_email='guldager.simon@gmail.com'
    url='https://github.com/simonguld/CellSegmentationTracker'
    description='Package for cell segmentation, tracking and subsequent (biophysical) statistical analysis',
    packages=find_packages(),

    classifiers=[
        'Development Status :: 1 - Planning',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.8',
    ]

)