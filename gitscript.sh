
SRC=$(echo "main.py \
  pibuck/main.py \
  pibuck/tools.py \
  pibuck/pibuck.py \
  pibuck/__init__.py")


untFiles=($(git ls-files --others --exclude-standard))


for f in ${SRC[@]}; do
  if [[ "${untFiles[@]}" =~ ( |^)"$f"( |$) ]] ; then
    echo "$f have not been added to track list"
    git add "$f"
  else
    echo "$f is in the track list"
  fi
done

