#!/usr/bin/env python
import numpy as np
from lxml import etree
import sys

if __name__ == '__main__':

    inp = sys.argv[1]
    parser = etree.XMLParser(remove_blank_text=True)
    xml = etree.parse(inp,parser)

    detset = xml.find('.//determinantset')
    
    # move out basis and sposet builder
    spobuilder = etree.Element('sposet_builder',detset.attrib)
    detset.attrib.clear()

    for node_name in ['basisset','sposet']:
        node_list = detset.findall('./'+node_name)
        for node in node_list:
            detset.remove(node)
            spobuilder.append(node)
        # end for
    # end for
    
    wf = xml.find('.//wavefunction')
    wf.insert(0,spobuilder)
    xml.write('spo-'+inp,pretty_print=True)

# end __main__
