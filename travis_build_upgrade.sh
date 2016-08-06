#!/bin/sh

apt-get install -y python-software-properties
add-apt-repository -y ppa:george-edison55/cmake-3.x
add-apt-repository -y ppa:ubuntu-toolchain-r/test
apt-get update -qq

#echo "Package: *\nPin: release o=LP-PPA-george-edison55-cmake-3.x\nPin-Priority: 1000" > /etc/apt/preferences.d/cmake_prefs

#cat /etc/apt/preferences.d/cmake_prefs
sudo apt-cache policy cmake

apt-get install -y --reinstall cmake
apt-get install -y gcc-5 g++-5


