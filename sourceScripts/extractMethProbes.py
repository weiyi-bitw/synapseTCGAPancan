import os
import sys
import getopt

def main():
  
  inFile = None
  mapFile = None
  outFile = None
# get options
  try:
    opts, args = getopt.getopt(sys.argv[1:], "i:m:o:", ["input=","map=", "output="])
  except getopt.GetoptError, err:
    print str(err)
    # usage()
    sys.exit(2)
  for opt, arg in opts:
    if opt in ("-i", "--input"):
      inFile = arg
    elif opt in ("-m", "--map"):
      mapFile = arg
    elif opt in ("-o", "--output"):
      outFile = arg

  print "Loading map files for probes...\n"
  try:
    f = open(mapFile)
  except IOError as e:
    print 'ERROR: cannot open file ' + f + '\n'
    sys.exit(2)

  f.readline() # first line header
  probes = []
  for line in f:
    tokens = line.split('\t')
    probes.append(tokens[0])

  print "Extracting probes...\n"
  try:
    f = open(inFile)
  except IOError as e:
    print 'ERROR: cannot open file ' + f + '\n'
    sys.exit(2)

  fo = open(outFile, 'w')
  line = f.readline()
  fo.write(line)
  cnt = 0
  for line in f:
    cnt += 1
    tokens = line.split('\t')
    if tokens[0] in probes:
      fo.write(line)
    if cnt % 10000 == 0:
      print str(cnt) + " lines read.\n"

  f.close()
  fo.close()
  
  print "Have a nice day :)\n"

if __name__ == "__main__":
    sys.exit(main())
