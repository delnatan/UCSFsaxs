# -*- mode: python -*-
a = Analysis(['main_saxsgui.py'],
             pathex=['/Users/delnatan/Codes/UCSFsaxs'],
             hiddenimports=['scipy.special._ufuncs_cxx'],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='UCSFsaxs',
          debug=True,
          strip=None,
          upx=True,
          console=True , icon='Resources/ucsfsaxsicon.icns')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='UCSFsaxs')

for i in a.binaries:
	print i
