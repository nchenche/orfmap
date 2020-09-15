# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:08:04 2020

@author: nicolas
"""

import setuptools as st

packages=['orfmap', 'orfmap.lib']

st.setup(name='orfmap',
         python_requires='>=3',
         version='0.0',
         packages=packages,
         package_data={'orfmap': ['data/Acupr.*']
                      },
         entry_points={
             'console_scripts': [
                 'run_orfmap=orfmap.main:main']
                      }
        )


