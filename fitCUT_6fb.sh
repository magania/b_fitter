./WsBs
./DataBs -i ws.root -o ws_CUT.root -r /work/elsanto-clued0/Z/Bs/sample/PRL.root -c "run<254000"
./FitBs -i ws_CUT.root -o ws_CUT_fit.root -n 8 
