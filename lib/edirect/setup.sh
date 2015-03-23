#!/bin/bash -norc

DIR="$( cd "$( dirname "$0" )" && pwd )"

cat <<EOF

Trying to establish local installations of any missing Perl modules
(as logged in $DIR/setup-deps.log).
Please be patient, as this step may take a little while.
EOF
mkdir -p "$DIR/_cpan/CPAN"
echo '1;' >> "$DIR/_cpan/CPAN/MyConfig.pm"
perl -I"$DIR/_cpan" "$DIR/setup-deps.pl" </dev/null >"$DIR/setup-deps.log" 2>&1
rm -rf "$DIR/_cpan"

target=bash_profile
if ! grep "$target" "$HOME/.bashrc" >/dev/null 2>&1
then
  if grep 'bashrc' "$HOME/.target" >/dev/null 2>&1
  then
    target=bashrc
  else
    echo 'source ~/.bash_profile' >> $HOME/.bashrc
  fi
fi
if ! grep "PATH.*edirect" "$HOME/.$target" >/dev/null 2>&1
then
  echo "export PATH=\$PATH:$DIR" >> $HOME/.$target
fi

echo ""
echo "ENTREZ DIRECT HAS BEEN SUCCESSFULLY INSTALLED AND CONFIGURED"
echo ""
