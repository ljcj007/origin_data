"""在新目录内编译，避免读取工作区旧模块接口。"""
import argparse
from pathlib import Path
import subprocess
from run_thermal_scan import save,sha


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--build-dir',type=Path,required=True)
    args=p.parse_args()
    source=args.source.resolve();folder=args.build_dir.resolve()
    folder.mkdir(parents=True,exist_ok=False)
    executable=folder/('sse_thermal.exe' if __import__('os').name=='nt' else 'sse_thermal')
    command=['gfortran','-O2','-Wall','-Wextra','-fcheck=all',str(source),'-o',str(executable)]
    result=subprocess.run(command,cwd=folder,capture_output=True,text=True,timeout=120)
    for name,content in [('compile.stdout',result.stdout),('compile.stderr',result.stderr)]:
        with (folder/name).open('x',encoding='utf-8') as handle: handle.write(content)
    save(folder/'build_manifest.json',dict(command=command,returncode=result.returncode,source_sha256=sha(source),
         executable_sha256=sha(executable) if result.returncode==0 else None,isolated_working_directory=str(folder)))
    raise SystemExit(result.returncode)
