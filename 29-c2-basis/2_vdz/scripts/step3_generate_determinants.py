def step3_generate_determinants(ndet,det_dir,phfmol_inp_template='./templates/phfrun.inp',submit_file='./templates/dets.msub'):
  import subprocess as sp
  import os
  # parse fcidump.dat for neri
  def line_num(expression,fname):
      proc = sp.Popen(['grep','-n',expression,fname],stdout=sp.PIPE,stderr=sp.PIPE)
      out,err = proc.communicate()
      first_line = out.split('\n')[0]
      idx = int(first_line.split(':')[0])
      return idx
  # end def line_num
  start = line_num('&END','fcidump.dat')
  end   = line_num('0  0','fcidump.dat')
  neri = end-start-1 # number of 2-electron integrals (i.e. electron repulsion integrals)

  with open('fcidump.dat','r') as f:
    header = f.readline() # e.g. &FCI NORB=   4,NELEC= 2,MS2=0,
  # end with
  nbas  = int( header.split('NORB')[1].split('=')[1].split(',')[0] )
  nelec = int( header.split('NELEC')[1].split('=')[1].split(',')[0] )
  nup = nelec/2
  ndn = nelec/2
  if nup+ndn!=nelec:
    raise RuntimeError('cannot setup RHF')
  # end if

  # setup phfmol run
  if not os.path.isdir(det_dir):
    os.system('mkdir '+det_dir)
  else:
    raise RuntimeError('determinants already generated? remove %s to rerun' % det_dir)
  # end if
  #  first input
  inp1 = os.path.join(det_dir,'phfrun0.inp')
  os.system( ' '.join(['cp ',phfmol_inp_template,inp1]) )
  os.system( 'sed -i "s/mynbas/%d/" %s' % (nbas,inp1) )
  os.system( 'sed -i "s/mynup/%d/" %s' % (nup,inp1) )
  os.system( 'sed -i "s/myndn/%d/" %s' % (ndn,inp1) )
  os.system( 'sed -i "s/myni2s/%d/" %s' % (neri,inp1) )
  os.system( 'sed -i "s/mylread/.false./" %s' % inp1 )
  os.system( 'sed -i "s/myigsmix/0/" %s' % inp1 ) # do not mix orbitals to keep HF determinant
  #  next inputs
  inp2 = os.path.join(det_dir,'phfrun1.inp')
  os.system( ' '.join(['cp ',phfmol_inp_template,inp2]) )
  os.system( 'sed -i "s/mynbas/%d/" %s' % (nbas,inp2) )
  os.system( 'sed -i "s/mynup/%d/" %s' % (nup,inp2) )
  os.system( 'sed -i "s/myndn/%d/" %s' % (ndn,inp2) )
  os.system( 'sed -i "s/myni2s/%d/" %s' % (neri,inp2) )
  os.system( 'sed -i "s/mylread/.true./" %s' % inp2 )
  os.system( 'sed -i "s/myigsmix/3/" %s' % inp2 ) # medium mixing
  #  submit job
  fsub = os.path.join(det_dir,os.path.basename(submit_file))
  os.system( 'cp %s %s' % (submit_file,det_dir) )
  os.system( 'sed -i "s/myndet/%d/" %s' % (ndet,fsub) )
  os.system( 'cd %s; sbatch %s' % (det_dir,os.path.basename(fsub)) )
