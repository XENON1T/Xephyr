from pathlib import Path
import json

p = Path('.')
l = list(p.glob('**/info.json'))

cmake_list_file = open('CMakeLists.txt','a')

print("Xephyr package manager found %i packages:" %len(l))
for i in l:
    with open (str(i)) as pkg_info_file:
        info = json.load(pkg_info_file)
    print("\t",info['pkg_name'])
    print("\t ---- Version - ",info['pkg_version'])
    print("\t ---- Require - ",info['dependencies'])
    if info['pkg_name'] == 'xephyr':
        pass
    else:
        src = i.parent / 'src'
        print(src)
        src_list = list(src.glob('*.cxx'))
        if(src.exists() and len(src_list) > 0 ):
            source = info['pkg_name']+ '_source'
            lib    = info['pkg_name']+ '_lib'
            cmake_list_file.write('file(GLOB ' + source + ' "' + str(src) + '/*.cxx")\n')
            cmake_list_file.write('add_library('+ lib + ' SHARED ${'+ source +'})\n')
            cmake_list_file.write('target_link_libraries('+ lib + ' xelib)\n')

        exe_list = list((i.parent).glob('*main.cxx'))
        for exe in exe_list :
            cmake_list_file.write('add_executable(' + info['pkg_name']  + ' ' + str(exe) +')\n')
            if(len(src_list) > 0):
                cmake_list_file.write('target_link_libraries(' + info['pkg_name']  +  ' ' + lib +')\n')
            else:
                cmake_list_file.write('target_link_libraries(' + info['pkg_name']  +  ' xelib)\n'  )

