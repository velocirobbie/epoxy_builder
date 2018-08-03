import json
import yaml
import sys

a = yaml.load(open(sys.argv[1]+'.yaml','r'))
json.dump(a,open(sys.argv[1]+'.json','w'))
