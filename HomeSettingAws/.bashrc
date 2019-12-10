
alias rm='rm -i'
alias cp='cp -i'
alias mv='mv -i'
alias grep='grep --colour=auto'
alias via='vi ~/.apparixrc'
alias c='clear'

source $HOME/.bash_path
source $HOME/.bash_generic
#source $HOME/.bourne-apparish
source $HOME/.bash_apparix


export CLICOLOR=1
export LSCOLORS=gxBxhxDxfxhxhxhxhxcxcx
export INPUTRC=~/.inputrc

bind '"\t":menu-complete'

## bind '^F:reverse-search-history'
## bind '^G:forward-search-history'


# enable color support of ls and also add handy aliases
if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
   alias dir='dir --color=auto'
   alias vdir='vdir --color=auto'

    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
fi


# test
