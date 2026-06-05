# run_m.py — Claude CodeがMATLABスクリプトを実行するためのブリッジ
import matlab.engine
import sys
import os

script_path = os.path.abspath(sys.argv[1])
script_dir = os.path.dirname(script_path)
script_name = os.path.splitext(os.path.basename(script_path))[0]

eng = matlab.engine.start_matlab()
eng.addpath(script_dir, nargout=0)
eng.eval(script_name, nargout=0)
eng.quit()