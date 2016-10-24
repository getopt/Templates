" set tags=tags
set tags=tags;
"  search for tags files in directories higher up.

filetype on

"   syntax on
set et
set expandtab
set nowrap
set tabstop=4
set shiftwidth=4
set sidescrolloff=3
set scrolloff=3
set backspace=2
set background=dark
set backspace=indent,eol,start " this makes bakspace to work normally
hi Matchparen none

set smarttab
set ai
set ic
set nocindent
" set autochdir

" hi Search NONE
" set hls
" :nnoremap <CR> :nohlsearch<CR>/<BS><CR>

set writebackup
set nobackup

set laststatus=2
set statusline=%F%m%r%h%w\ [POS=%04l,%04v][%p%%]\ [LEN=%L]
hi StatusLine ctermfg=black ctermbg=green

command NS :%s/  / /g 
command TU :TlistUpdate
command TT :TlistToggle
command TC :TlistClose
command TO :TlistOpen

setlocal spell spelllang=en
set spellcapcheck=
set nospell
set spf=$HOME/.spell.en.add

let Tlist_Exit_OnlyWindow = 1
let Tlist_WinWidth = 80
let Tlist_Use_Right_Window = 1
let Tlist_Display_Prototype = 1

" SecureCRT versions prior to 6.1.x do not support 4-digit DECSET
"     let &t_ti = "\<Esc>[?1049h"
"     let &t_te = "\<Esc>[?1049l"
" Use 2-digit DECSET instead
" let &t_ti = "\<Esc>[?47h"
" let &t_te = "\<Esc>[?47l"
