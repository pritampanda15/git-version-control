# Git Version Control Tutorial

## Introduction to Git
Git is a distributed version control system used for tracking changes in source code during software development. It offers coordination, scalability, and speed for small to very large projects.

## Getting Started with Git

### Installation
Before using Git, ensure it is installed on your system. To check if Git is installed, run:
```
git --version
```

### Global Configuration
After installing Git, configure your username and email, as these will be associated with your commits.

1. **Set Username**
   ```
   git config --global user.name "Your Name"
   ```
   Replace `"Your Name"` with your actual name.

2. **Set Email**
   ```
   git config --global user.email "your.email@example.com"
   ```
   Replace `"your.email@example.com"` with your actual email.

3. **Set Default Editor**
   To set Visual Studio Code as your default editor:
   ```
   git config --global core.editor "code --wait"
   ```
   - Note: In VS Code, open the command palette and type `>Shell Command: Install 'code' command in PATH`.

4. **View Configurations**
   To view your configurations:
   ```
   git config --global --edit
   ```

5. **Set Carriage Return and Line Feed (CRLF) Characters**
   - For Unix/Linux/macOS:
     ```
     git config --global core.autocrlf input
     ```
   - For Windows:
     ```
     git config --global core.autocrlf true
     ```

6. **Set Visual Studio Code as the Diff Tool**
   ```
   git config --global diff.tool vscode
   git config --global difftool.vscode.cmd "code --wait --diff $LOCAL $REMOTE"
   ```

### Basic Git Commands

1. **Initialize a New Repository**
   ```
   git init
   ```

2. **Clone an Existing Repository**
   ```
   git clone <repository_url>
   ```

3. **Add Changes to Staging Area**
   ```
   git add <file_name>
   ```

4. **Commit Changes**
   ```
   git commit -m "Commit message"
   ```

5. **Push Changes to Remote Repository**
   ```
   git push origin <branch_name>
   ```

6. **Pull Changes from Remote Repository**
   ```
   git pull
   ```

7. **Check Status**
   ```
   git status
   ```

### Advanced Git Commands

1. **Branching**
   - Create a new branch:
     ```
     git branch <branch_name>
     ```
   - Switch to a branch:
     ```
     git checkout <branch_name>
     ```
   - Merge a branch:
     ```
     git merge <branch_name>
     ```

2. **Tagging**
   - Create a tag:
     ```
     git tag <tag_name>
     ```

3. **Stashing**
   - Stash changes:
     ```
     git stash
     ```
   - Apply stashed changes:
     ```
     git stash pop
     ```

4. **Viewing Commit History**
   ```
   git log
   ```

5. **Using Diff**
   - Compare changes:
     ```
     git diff
     ```

6. **Rebasing**
   ```
   git rebase <branch_name>
   ```

## Conclusion
This tutorial covers the basic and advanced usage of Git for version control. Consistent practice and exploration of Git's features will enhance your proficiency in managing software development projects effectively.
