fn main() {
    println!("cargo:rustc-link-search=/usr/local/lib");
    println!("cargo:rustc-link-lib=static=GKlib");
    println!("cargo:rustc-link-lib=static=metis");
}