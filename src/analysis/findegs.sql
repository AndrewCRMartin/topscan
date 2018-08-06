select d.domid, p.resol
from domains d, pdb p
where d.cat = '2.60.120'
and d.pdbcode = p.pdbcode
order by resol;
