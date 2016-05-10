from distutils.core import setup

setup(
      name             = 'assemblies',
      version          = '0.0',
      description      = '(meta)genome assembly tools',
      long_description = 'an other day',
      license          = 'MIT',
      url              = 'http://github.com/moritz/assembly',
      author           = 'Moritz Buck',
      author_email     = 'mrtz.buck@gmail.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bioinformatics'],
      packages         = ['generic'],
      requires         = ['sh', 'tqdm', 'pandas', 'Bio'],
)
