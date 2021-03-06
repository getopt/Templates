set nocompatible              " required
filetype off                  " required

" set the runtime path to include Vundle and initialize
" set rtp+=~/.vim/bundle/Vundle.vim
" call vundle#begin()

set runtimepath^=/Users/smanakov/.vim/majutsushi-tagbar-c004652/


" alternatively, pass a path where Vundle should install plugins
"call vundle#begin('~/some/path/here')

" let Vundle manage Vundle, required
" Plugin 'gmarik/Vundle.vim'

" Add all your plugins here (note older versions of Vundle used Bundle instead of Plugin)


" All of your Plugins must be added before the following line
" call vundle#end()            " required
" filetype plugin indent on    " required


set hlsearch

" set tags=tags
set tags=tags;
"  search for tags files in directories higher up.

filetype on

syntax on
let Tlist_Ctags_Cmd = '/usr/bin/ctags'
set et
set expandtab
set nowrap
set tabstop=4
set shiftwidth=4
set sidescrolloff=3
set scrolloff=3
set backspace=2
set background=dark
set noerrorbells visualbell t_vb=
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

" " taglist specific options
" command NS :%s/  / /g 
" command TU :TlistUpdate
" command TT :TlistToggle
" command TC :TlistClose
" command TO :TlistOpen
" let Tlist_Exit_OnlyWindow = 1
" let Tlist_WinWidth = 80
" let Tlist_Use_Right_Window = 1
" let Tlist_Display_Prototype = 1

" tagbar specific options
command TO :TagbarOpen
command TC :TagbarClose
let tagbar_width=80

setlocal spell spelllang=en
set spellcapcheck=
set nospell
set spf=$HOME/.spell.en.add

" Can cause issues with clearing the screen after closing
" SecureCRT versions prior to 6.1.x do not support 4-digit DECSET
"     let &t_ti = "\<Esc>[?1049h"
"     let &t_te = "\<Esc>[?1049l"
" Use 2-digit DECSET instead
" let &t_ti = "\<Esc>[?47h"
" let &t_te = "\<Esc>[?47l"


" ####
" ##  Cygwin specific: to paste into windows clipboard
" ####

function! Putclip(type, ...) range
  let sel_save = &selection
  let &selection = "inclusive"
  let reg_save = @@
  if a:type == 'n'
    silent exe a:firstline . "," . a:lastline . "y"
  elseif a:type == 'c'
    silent exe a:1 . "," . a:2 . "y"
  else
    silent exe "normal! `<" . a:type . "`>y"
  endif
  "call system('putclip', @@)
  "As of Cygwin 1.7.13, the /dev/clipboard device was added to provide
  "access to the native Windows clipboard. It provides the added benefit
  "of supporting utf-8 characters which putclip currently does not. Based
  "on a tip from John Beckett, use the following:
  call writefile(split(@@,"\n"), '/dev/clipboard')
  let &selection = sel_save
  let @@ = reg_save
endfunction

vnoremap <silent> <leader>y :call Putclip(visualmode(), 1)<CR>
nnoremap <silent> <leader>y :call Putclip('n', 1)<CR>

com! -nargs=0 -range=% Putclip call Putclip('c', <line1>, <line2>)

" ####
" ##  Cygwin specific: to paste from windows clipboard
" ####

" # use \p to paste in command mode

function! Getclip()
  let reg_save = @@
  "let @@ = system('getclip')
  "Much like Putclip(), using the /dev/clipboard device to access to the
  "native Windows clipboard for Cygwin 1.7.13 and above. It provides the
  "added benefit of supporting utf-8 characters which getclip currently does
  "not. Based again on a tip from John Beckett, use the following:
  let @@ = join(readfile('/dev/clipboard'), "\n")
  setlocal paste
  exe 'normal p'
  setlocal nopaste
  let @@ = reg_save
endfunction

nnoremap <silent> <leader>p :call Getclip()<CR>

highlight DiffChange cterm=bold ctermbg=8 ctermfg=2
highlight DiffText   cterm=reverse ctermbg=2 ctermfg=17
highlight DiffAdd    cterm=bold ctermbg=22 ctermfg=255
highlight DiffDelete cterm=bold ctermfg=1 ctermbg=1
