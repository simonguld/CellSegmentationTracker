from setuptools import setup, find_namespace_packages

setup(
    name='cellsegmentationtracker',
    version="0.1.11",
    author='Simon Guldager Andersen',
    author_email='guldager.simon@gmail.com',
    url='https://github.com/simonguld/CellSegmentationTracker',
    description='Module for cell segmentation, tracking and subsequent (biophysical) statistical analysis',
    install_requires=["tifftools", "seaborn", "scipy", "scikit-image", 'scikit-learn'],  
    packages=find_namespace_packages(include=['cellsegmentationtracker'], exclude=['cellsegmentationtracker.__pycache__.*']), #, 'cellsegmentationtracker.*']),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.8, 3.9',
    ],
    include_package_data=True,
    package_data={'': ['models/*']}

)