#!/bin/sh

## This might be a bit flakey, but should do the job.

WIKI_DIR=../wiki
GIT_WIKI="--git-dir=$WIKI_DIR/.git --work-tree=$WIKI_DIR"

if test -n "$(git $GIT_WIKI status -s 2> /dev/null)"
then
    echo "Error: wiki git is dirty"
    echo "(try git reset --hard HEAD, *in the wiki dir*, perhaps)"
    exit 1
fi

git $GIT_WIKI pull --rebase --quiet origin
git push
