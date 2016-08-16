function Conn = findConn(NURBS, Knt, dir)

mcp = NURBS.NCtrlPts(1);
ncp = NURBS.NCtrlPts(2);
k = FindSpan(NURBS.NCtrlPts(1), NURBS.Order(dir), Knt, NURBS.KntVect{dir});
FDofs = sub2ind([mcp, ncp], (k - NURBS.Order(dir) - 1) * ones(1, ncp), 1 : ncp)';
SDofs = sub2ind([mcp, ncp], (k - NURBS.Order(dir)) * ones(1, ncp), 1 : ncp)';
TDofs = sub2ind([mcp, ncp], (k - NURBS.Order(dir) + 1) * ones(1, ncp), 1 : ncp)';

Conn = [FDofs SDofs; SDofs TDofs];
end