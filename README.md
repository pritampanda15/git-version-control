# git-version-control
Tutorial on git version control with basic and advanced commands

Configure globally on your 

'''
git –version
git config --global user.name "pritam”
git config --global user.email pritampkp15@gmail.com
git config --global core.editor "code --wait” (in vscode , open the command pallete and type >shell.commands and add vscode to path. in windows its default)
git config --global –e (open the editor and see the config file)
git config --global core.autocrlf input (setting carriage characters ; for windows set this as true.
git config --global diff.tool vscode (set vscode as difftool)
git config --global difftool.vscode.cmd "code --wait --diff" $LOCAL $REMOTE (open vscode and add $LOCAL $REMOTE)
'''

