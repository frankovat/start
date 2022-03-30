echo ----Running .bash_profile

conda -V
python -V
echo $PATH
echo $PATH | tr ":" "\n"

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/storage/plzen1/home/frankovat/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/storage/plzen1/home/frankovat/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/storage/plzen1/home/frankovat/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/storage/plzen1/home/frankovat/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

echo ----Completed running .bash_profile

if [ -f ~/.bashrc ]; then
   source ~/.bashrc
fi
