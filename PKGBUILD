# Maintainer: Dictor <dictor@terminator>
pkgname=orthrus-git
pkgver=0.3.0
pkgrel=1
pkgdesc="High-Assurance KLPT Quaternion Navigation Core for SQISign (C++20/C-API)"
arch=('x86_64')
url="https://github.com/Dictor457/orthrus"
license=('MIT')
depends=('gmp')
makedepends=('cmake' 'gcc' 'git' 'pkgconf')
provides=('orthrus' 'liborthrus.so')
conflicts=('orthrus')

build() {
    cmake -B build -S "$startdir" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/usr
    cmake --build build
}

package() {
    DESTDIR="$pkgdir" cmake --install build
}
