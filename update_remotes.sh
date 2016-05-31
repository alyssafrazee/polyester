#!/usr/bin/env bash

# This script should be run after cloning, it will update the remote
# information so that git svn commands work properly for both the master and
# release branches.

set -eou pipefail
IFS=$'\n\t'

set +u
package=$1

# otherwise use the default repo name
if [ -z "$package" ];then
  package=$(git remote -v | perl -ne 'if (m!/([^/]+?)(?:.git)?\s!) { print $1; exit}')
fi
set -u

base_url="https://hedgehog.fhcrc.org/bioconductor/"

is_mirror_clone () {
  git remote show origin | grep -i -q bioconductor-mirror
  return $?
}

add_release_tracking () {
  local branch=$1
  shift
  local remote=$1
  shift
  for release_branch in $@; do
    svn_branch=$(echo $release_branch | perl -ne 'if (/release-(\d+)\.(\d+)/) { print "RELEASE_$1_$2"; }')
    svn_url="$base_url/branches/$svn_branch/madman/Rpacks/$package"
    git config --add svn-remote.$release_branch.url $svn_url
    git config --add svn-remote.$release_branch.fetch :refs/remotes/git-svn-$release_branch
    git update-ref refs/remotes/git-svn-$release_branch refs/$remote/$release_branch
  done
}

add_branch () {
  set +eu
  local local_branch=$1
  local remote_branch=${2-$local_branch}
  if ! git branch --track $local_branch bioc/$remote_branch  2>/dev/null 1>&2; then
    1>&2 cat <<END
$local_branch already exists, create a custom branch to track bioc/$remote_branch with
  \`git branch --track NEW_NAME bioc/$remote_branch\`
END
  fi
}

if is_mirror_clone; then
    git checkout master
    git svn init "$base_url/trunk/madman/Rpacks/$package"
    git update-ref refs/remotes/git-svn refs/remotes/origin/master
    git svn rebase
    git remote add bioc https://github.com/Bioconductor-mirror/${package}.git
    git fetch bioc 2>/dev/null 1>&2

    release_branches=$(git branch -r | perl -ne 'if (m!origin/(release-.*)!) { print $1, "\n" }')
    for release_branch in ${release_branches[@]}; do
      add_branch $release_branch
    done

    add_release_tracking origin heads $release_branches

    cat <<\END
Commit to git as normal, when you want to push your commits to svn
  1. `git svn rebase` to get the latest SVN changes.
  2. `git svn dcommit --add-author-from` to commit your changes to SVN.
END

else
    git remote add bioc "https://github.com/Bioconductor-mirror/${package}.git"
    git fetch bioc 2>/dev/null 1>&2
    git config --add svn-remote.devel.url "$base_url/trunk/madman/Rpacks/$package"
    git config --add svn-remote.devel.fetch :refs/remotes/git-svn-devel
    git update-ref refs/remotes/git-svn-devel refs/remotes/bioc/master

    release_branches=$(git branch -r | perl -ne 'if (m!bioc/(release-.*)!) { print $1, "\n" }')
    add_release_tracking bioc remotes/bioc $release_branches

    for release_branch in ${release_branches[@]}; do
      add_branch $release_branch
    done
    add_branch devel master
    cat <<\END
Commit to git as normal, when you want to push your commits to svn
  1. `git checkout devel` to switch to the devel branch. (use release-X.X for
        release branches)
  2. `git svn rebase` to get the latest SVN changes.
  3. `git merge master --log` to merge your changes from the master branch
        or skip this step and work directly on the current branch.
  4. `git svn dcommit --add-author-from` to sync and commit
        your changes to svn.
END
fi
