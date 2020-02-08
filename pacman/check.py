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
        print('.... Scanning source files in: ' + str(src) )
        src_list = list(src.glob('*.cxx'))
        if(src.exists() and len(src_list) > 0 ):
            print('.... Building shared library for ' + info['pkg_name'])
            source = info['pkg_name']+ '_source'
            lib    = info['pkg_name']+ '_lib'
            #cmake_list_file.write('include_directories('+ str(src) + ')\n')  # maybe this will be needed when compiling including also other lib fro other package
            cmake_list_file.write('file(GLOB ' + source + ' "' + str(src) + '/*.cxx")\n')
            cmake_list_file.write('add_library('+ lib + ' SHARED ${'+ source +'})\n')
            cmake_list_file.write('target_link_libraries('+ lib + ' xelib)\n')

        exe_list = list((i.parent).glob('*_main.cxx'))
        for exe in exe_list :
            #remove the main.cxx suffix and split by / 
            exe_suffix = (str(exe).split('/')[-1])[:-9] 
            print('.... Building executable: ' + exe_suffix)
            cmake_list_file.write('add_executable(' + info['pkg_name'] + '_' + exe_suffix  + ' ' + str(exe) +')\n')
            if(len(src_list) > 0):
                cmake_list_file.write('target_link_libraries(' + info['pkg_name']  + '_' + exe_suffix  +  ' ' + lib +')\n')
                cmake_list_file.write('target_include_directories(' + info['pkg_name'] + '_' + exe_suffix  + ' PUBLIC ' + str(src) +')\n')
            else:
                cmake_list_file.write('target_link_libraries(' + info['pkg_name']  + '_' + exe_suffix  +  ' xelib)\n'  )

