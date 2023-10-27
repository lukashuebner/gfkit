# coding: utf-8
import io
import tskit

num_sites = 10**10

nodes = io.StringIO(
    "is_sample   time\n"
    "1  0 \n"
    "1  0 \n"
    "1  0 \n"
    "1  0 \n"
    "0  1 \n"
    "0  2 \n"
    "0  3 \n"
)

edges = io.StringIO(
    "left    right   parent  child\n"
    "0  1   4   0,1\n"
    "0  1   5   2,3\n"
    "0  1   6   4,5\n"
)

sites = io.StringIO(
    "position ancestral_state\n" + "\n".join([f"{1/num_sites * site} 0" for site in range(0, num_sites)])
)

mutations = io.StringIO(
    "site node derived_state\n" + "\n".join(f"{site} 1 1\n{site} 5 2" for site in range(0, num_sites))
)

#ts = tskit.load_text(nodes, edges, strict=False, sequence_length=1)
#for site in range(0, num_sites):
#    ts.tables.sites.add_row(position = site, ancestral_state='0')
#    ts.tables.mutations.add_row(site = site, node=1, derived_state="1")
#    ts.tables.mutations.add_row(site = site, node=5, derived_state="2")

#ts.tables.compute_mutation_parens()

print(f"Num sites:      {num_sites}")
print(f"Num seg. sites: {ts.segregating_sites(span_normalise=False)}")
