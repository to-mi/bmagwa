#!/usr/bin/env python

# Postprocessing utilities for BMAGWA software (v2.0)
# Author: Tomi Peltola <tomi.peltola@aalto.fi>

import struct, sys, os

def read_data(filename, fmt):
  try:
    itemsize = struct.calcsize(fmt)
    f = open(filename, 'rb')
    item = f.read(itemsize)
    while item:
      yield struct.unpack(fmt, item)[0]
      item = f.read(itemsize)
  finally:
    f.close()

def writeok(filename):
  doit = True
  if os.path.exists(filename):
    doit = False
    print 'File already exists: %s' % filename
    a = raw_input('Overwrite [y/n]? ')
    while a not in ['y', 'n']:
      a = raw_input('Overwrite [y/n]? ')
    if a == 'y':
      doit = True
  return doit

def convert(args):
  filename = args[0]
  outfile = filename + '.txt'
  if len(args) == 2:
    fmt = args[1]
  else:
    fmt = guess_fmt(args[0])
  if fmt == None:
    print >> sys.stderr, 'Cannot guess format for %s. Skipping.' % filename
    return
  if writeok(outfile):
    f = open(outfile, 'w')
    for x in read_data(filename, fmt):
      f.write(str(x) + '\n')
    f.close()

def output(args):
  filename = args[0]
  if len(args) == 2:
    fmt = args[1]
  else:
    fmt = guess_fmt(args[0])
  if fmt == None:
    print >> sys.stderr, 'Cannot guess format for %s. Skipping.' % filename
    return
  for x in read_data(filename, fmt):
    print x

def guess_fmt(filename):
  import re
  patterns = ((r'jumpdistance\.', 'B'),
              (r'loci\.', 'I'),
              (r'log_likelihood\.', 'd'),
              (r'alpha\.', 'd'),
              (r'log_prior\.', 'd'),
              (r'modelsize\.', 'I'),
              (r'move_type\.', 'B'),
              (r'move_size\.', 'B'),
              (r'pve\.', 'd'),
              (r'rao\.', 'd'),
              (r'rao_types\.', 'd'),
              (r'sigma2\.', 'd'),
              (r'types\.', 'B'))
  for p in patterns:
    if re.search(p[0], filename) != None:
      return p[1]
  return None

def mcmcpos(args):
  import re
  basename = args[0]
  nsnps = int(args[1])
  burnin = int(args[2])
  if len(args) >= 4:
    thin = max(1, int(args[3]))
  else:
    thin = 1
  
  directory, name = os.path.split(basename)
  if directory == '': directory = '.'
  filenames = filter(lambda x: re.match(name + r'\d+_loci.dat$', x) != None,
                     os.listdir(directory))
  print "mcmcpos: Using following chains to compute association probabilities:"
  for x in filenames: print " " + x
  p = [[0.0] * nsnps for x in filenames]
  
  for i,f in enumerate(filenames):
    msfile = f.replace('_loci.', '_modelsize.')
    loci = read_data(os.path.join(directory, f), guess_fmt(f))
    nsamples = 0
    for j,ms in enumerate(read_data(os.path.join(directory, msfile),
                          guess_fmt(msfile))):
      if j < burnin or (j - burnin) % thin != 0:
        for x in range(int(ms)): loci.next()
      else:
        nsamples += 1
        for x in range(int(ms)):
          p[i][int(loci.next())] += 1.0
    for j,val in enumerate(p[i]):
      p[i][j] = val / nsamples
  
    outfile = os.path.join(directory, f.replace('_loci.dat', '_mcmcpos.txt'))
    if writeok(outfile):
      fout = open(outfile, 'w')
      fout.write('\n'.join(map(str, p[i])))
      fout.write('\n')
      fout.close()

  outfile = os.path.join(directory, name + '_mcmcpos.txt')
  if writeok(outfile):
    fout = open(outfile, 'w')
    for i in xrange(nsnps):
      fout.write(str(sum((x[i] for x in p)) / len(p)) + '\n')
    fout.close()


if __name__ == "__main__":
  nargs = len(sys.argv)
  if nargs >= 5 and sys.argv[1] == "mcmcpos":
    mcmcpos(sys.argv[2:])
    sys.exit()
  if nargs >= 3 and sys.argv[1] != "mcmcpos" and \
     not os.path.exists(sys.argv[2]):
    print >> sys.stderr, "File or directory does not exist: %s" % sys.argv[2]
    sys.exit()
  if nargs == 3 and sys.argv[1] == "convert" and os.path.isdir(sys.argv[2]):
    for x in (f for f in os.listdir(sys.argv[2]) if f.endswith('.dat')):
      convert([os.path.join(sys.argv[2], x)])
  elif nargs >= 3 and sys.argv[1] == "convert":
    convert(sys.argv[2:])
  elif nargs >= 3 and sys.argv[1] == "output":
    output(sys.argv[2:])
  else:
    print "Usage: python %s options" % sys.argv[0]
    print "where options is one of following:"
    print "(nothing or invalid) : print help"
    print "convert directory    : converts all binary files (.dat) to text"
    print "                       files."
    print "convert filename     : converts binary file with name filename to"
    print "                       text file."
    print "convert filename fmt : converts binary file with name filename and"
    print "                       format fmt to text file."
    print "output filename      : similar to convert but prints to standard"
    print "                       output."
    print "output filename fmt  : similar to convert but prints to standard"
    print "                       output."
    print "mcmcpos bn nsnps b t : computes posterior association probabilities"
    print "                       from MCMC chains with basename bn (i.e., path"
    print "                       to loci.dat files up to just before the chain"
    print "                       number). nsnps is the number of SNPs in the"
    print "                       dataset, b is the number of burnin samples to"
    print "                       drop and t is thinning (in addition to the"
    print "                       one specified in the configuration file)."
    print
    print "Convert commands write output to 'filename' + '.txt'."
    print "Mcmcpos writes to files 'basenameX_mcmcpos.txt' and"
    print "'basename_mcmcpos.txt' (combined posterior of all chains)."
    print
    print "Examples:"
    print "python %s convert results" % sys.argv[0]
    print "python %s convert results/chain1_modelsize.dat" % sys.argv[0]
    print "python %s convert results/chain1_modelsize.dat I" % sys.argv[0]
    print "python %s output results/chain0_rao.dat > rao0.txt" % sys.argv[0]
    print "python %s mcmcpos results/chain 10000 50000 1" % sys.argv[0]

