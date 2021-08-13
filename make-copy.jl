#!/home/marcelovmaciel/.local/bin/julia
using Base.Filesystem

roam_path = "~/Drive/Org/org-roam-mvm/"
bib_path = "~/Drive/Org/bib/refs.bib"

files = [
    "acemoglu10_opinion_dynam_learn_social_networ.org",
    "grofman04_downs_two_party_conver.org",
    "falmagne2019stochastic.org",
    "stasser1981group.org",
    "grofman1981mathematical.org",
    "degroot1974reaching.org",
    "becker19_wisdom_partis_crowd.org",
    "guilbeault18_social_learn_partis_bias_inter_climat_trend.org",
    "centola2018truth.org",
    "holland1996hidden.org",
    "OSTROM2013155.org",
    "kollman1998political.org",
    "laver2011party.org"
]


notespaths = map(x->joinpath(expanduser(roam_path),x), files)

function copytonotes(file)
        cp(file, string(@__DIR__,
                      "/notes/zettelkasten-copies/",
                      last(splitpath(file))), force = true)
end

function copybib(file)
            cp(file, string(@__DIR__,
                      "/notes/",
                      last(splitpath(file))), force = true)
end


foreach(copytonotes, notespaths)
copybib(expanduser(bib_path))
