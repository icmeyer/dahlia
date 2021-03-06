from scipy.constants import physical_constants, year


MOL = physical_constants['Avogadro constant'][0]
EV2JOULE = physical_constants['electron volt'][0]
YR2S = year
INV_YR2S = 1/YR2S
thermal_nubar = 2.5
fast_nubar = 3.0

THERMAL_DATA = {'922350': {'fis':  38.8,   # Uranium
                           'cap':  8.7,
                           'nubar': 2.4},
                '922380': {'fis':  0.103, 
                           'cap':  0.86 ,
                           'nubar': thermal_nubar},
                '932370': {'fis':  0.52 ,  # Neptunium
                           'cap':  33,
                           'nubar': thermal_nubar},
                '932380': {'fis':  134  , 
                           'cap':  13.6 ,
                           'nubar': thermal_nubar},
                '942380': {'fis':  2.4  ,  # Plutonium
                           'cap':  27.7 ,
                           'nubar': thermal_nubar},
                '942390': {'fis':  102, 
                           'cap':  58.7,
                           'nubar': 2.9},
                '942400': {'fis':  0.53 , 
                           'cap':  210.2,
                           'nubar': thermal_nubar},
                '942410': {'fis':  102.2 , 
                           'cap':  40.9 ,
                           'nubar': thermal_nubar},
                '942420': {'fis':  0.44, 
                           'cap':  28.8,
                           'nubar': thermal_nubar},
                '952410': {'fis':  1.1 ,  # Americium
                           'cap':  110 ,
                           'nubar': thermal_nubar},
                '952420': {'fis':  159 , 
                           'cap':  301 ,
                           'nubar': thermal_nubar},
                '952421': {'fis':  595 , 
                           'cap':  137 ,
                           'nubar': thermal_nubar},
                '952430': {'fis':  0.44, 
                           'cap':  49  ,
                           'nubar': thermal_nubar},
                '962420': {'fis':  1.14,   # Curium
                           'cap':  4.5 ,
                           'nubar': thermal_nubar},
                '962430': {'fis':  88  ,
                           'cap':  14  ,
                           'nubar': thermal_nubar},
                '962440': {'fis':  1.0 ,
                           'cap':  16  ,
                           'nubar': thermal_nubar},
                '962450': {'fis':  116 ,
                           'cap':  17  ,
                           'nubar': thermal_nubar},
                }

FAST_DATA = {'922350': {'fis':  1.98,  # Uranium
                        'cap':  0.57,
                        'nubar': 2.6},
             '922380': {'fis':  0.04, 
                        'cap':  0.30,
                        'nubar': 2.6},
             '932370': {'fis':  0.32,  # Neptunium
                        'cap':  1.7 ,
                        'nubar': fast_nubar},
             '932380': {'fis':  3.6 , 
                        'cap':  0.2 ,
                        'nubar': fast_nubar},
             '942380': {'fis':  1.1 ,  # Plutonium
                        'cap':  0.58,
                        'nubar': fast_nubar},
             '942390': {'fis':  1.86, 
                        'cap':  0.56,
                        'nubar': 3.1},
             '942400': {'fis':  0.36, 
                        'cap':  0.57,
                        'nubar': fast_nubar},
             '942410': {'fis':  2.49, 
                        'cap':  0.47,
                        'nubar': fast_nubar},
             '942420': {'fis':  0.24, 
                        'cap':  0.44,
                        'nubar': fast_nubar},
             '952410': {'fis':  0.27,  # Americium
                        'cap':  2.0 ,
                        'nubar': fast_nubar},
             '952420': {'fis':  3.2 , 
                        'cap':  0.6 ,
                        'nubar': fast_nubar},
             '952421': {'fis':  3.3 , 
                        'cap':  0.6 ,
                        'nubar': fast_nubar},
             '952430': {'fis':  0.21, 
                        'cap':  1.8 ,
                        'nubar': fast_nubar},
             '962420': {'fis':  0.58,   # Curium
                        'cap':  1.0 ,
                        'nubar': fast_nubar},
             '962430': {'fis':  7.2 ,
                        'cap':  1.0 ,
                        'nubar': fast_nubar},
             '962440': {'fis':  0.42,
                        'cap':  0.6 ,
                        'nubar': fast_nubar},
             '962450': {'fis':  5.1 ,
                        'cap':  0.9 ,
                        'nubar': fast_nubar},
             }
