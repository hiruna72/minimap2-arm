#!/bin/sh

toolchain_file="/home/sanoj/Android/Sdk/ndk-bundle/build/cmake/android.toolchain.cmake"

# # to create a minimap2 library
cp main.c tempmain
# sed -i 's/int main(int argc/int init_samtools(int argc/g' bamtk.c
sed -i ':a;N;$!ba;s/int main(int argc, char \*argv\[\])\n{/int init_minimap2(int argc, char *argv[])\n{/g' main.c

touch interface.h
echo "int init_minimap2(int argc, char *argv[]);" > interface.h

mkdir -p build
rm -rf build
mkdir build
cd build

# for architecture x86
 # cmake .. -DDEPLOY_PLATFORM=x86
 # make -j 8

# # for architecture armeabi-V7a
cmake .. -G Ninja -DCMAKE_TOOLCHAIN_FILE:STRING=$toolchain_file -DANDROID_PLATFORM=android-21 -DDEPLOY_PLATFORM:STRING="armeabi-v7a" -DANDROID_ABI="armeabi-v7a" -DANDROID_ARM_NEON=ON

# # for architecture arm64-v8a
# cmake .. -G Ninja -DCMAKE_TOOLCHAIN_FILE:STRING=$toolchain_file -DANDROID_PLATFORM=android-23 -DDEPLOY_PLATFORM:STRING="arm64-v8a" -DANDROID_ABI="arm64-v8a" -DANDROID_ARM_NEON=ON

ninja
cd -
mv tempmain main.c