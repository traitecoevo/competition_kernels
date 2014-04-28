#!/bin/sh

WIKI_DIR=../wiki
FIGURE_DIR="figure"
GIT_WIKI="--git-dir=$WIKI_DIR/.git --work-tree=$WIKI_DIR"
SHORT_SHA=$(git rev-parse --short HEAD)

# Check that we don't have any residual "unnamed chunk" figures that I
# never want committed.
if test -n "$(find $FIGURE_DIR -maxdepth 1 -name '*unnamed*' -print -quit)"
then
    echo "Error: Unnamed figure chunks found"
    exit 1
fi

if test -n "$(git $GIT_WIKI status -s 2> /dev/null)"
then
    echo "Error: wiki git is dirty"
    echo "(try git reset --hard HEAD, *in the wiki dir*, perhaps)"
    exit 1
fi

# It will be good to work out where images are deleted (they'll remain
# in the wiki even when unused).  Better would be to use rsync, but
# delete and recreate also works.

# Once there are multiple files to track here, this will be looped
# over probably.

cp simple_models.md $WIKI_DIR
git $GIT_WIKI rm -f --quiet --ignore-unmatch -- "$FIGURE_DIR/simple_models*"
mkdir -p $WIKI_DIR/$FIGURE_DIR
cp $FIGURE_DIR/simple_models* $WIKI_DIR/$FIGURE_DIR
git $GIT_WIKI add simple_models.md "$FIGURE_DIR/simple_models*"

if git $GIT_WIKI status --porcelain --untracked-files=no | grep --quiet '^[A-Z]'
then
    git $GIT_WIKI commit -q -m "Updated wiki [at $SHORT_SHA]"
else
    echo "Wiki already up to date"
fi
