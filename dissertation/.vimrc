set number
syntax on
set autoindent


filetype plugin on
filetype indent on
let g:tex_flavor='latex'
let g:TTarget='pdf'

let g:Tex_BibtexFlavor = 'biber'
let g:Tex_DefaultTargetFormat = 'pdf'
let g:Tex_CompileRule_pdf = 'arara -v $*'
let g:Tex_MultipleCompileFormats='pdf'

set tabstop=4
set shiftwidth=4
set expandtab


set grepprg=grep\ -nH\ $*
set runtimepath=~/.vim,$VIM/vimfiles,$VIMRUNTIME,$VIM/vimfiles/after,~/.vim/after
