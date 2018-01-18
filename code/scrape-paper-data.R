library("pdftools")
library("tidyverse")
library("poppr")

txt <- pdf_text("~/Documents/qual_papers/hambelton2002clonal.pdf")
tab <- unlist(strsplit(txt[5:6], "\n"))
prints <- tab %>% 
  str_extract("[01. ]{20,}") %>% 
  trimws() %>% 
  str_replace_all("(?<=[.01])[ ]{1}([.01])", "\\1") %>% # Zhian's first look-behind!
  str_replace_all("(?<=[ ])[1 ]{10,}", "") %>%
  trimws()
haps <- tab %>% 
  str_extract("[A-Z](?=[[:blank:]]{2,})")
nind <- tab %>%
  str_extract("(?<= )[0-9]{1,2}$") %>%
  as.integer()

replace_dots <- function(x){
  x <- strsplit(x, "")
  if (length(x) > 1){
    for (i in seq(from = 2, to = length(x))){
      to_replace <- x[[i]] == "."
      x[[i]][to_replace] <- x[[1]][to_replace]
    }
  } 
  vapply(x, paste, FUN.VALUE = character(1), collapse = ".")
}

df <- data_frame(prints, haps, nind) %>%
  mutate(n = nchar(prints)) %>%
  filter(n > 10) %>%
  mutate(prints = gsub(" ", "", prints)) %>% 
  mutate(prints = substr(prints, 1, floor(mean(n)))) %>%
  mutate(n = nchar(prints)) %>%
  group_by(haps) %>%
  mutate(prints = replace_dots(prints)) %>%
  separate(prints, into = paste0("L", seq(min(.$n))), sep = "\\.") %>%
  select(-n) %>%
  ungroup()

df$haps[is.na(df$haps)] <- c("R", "S", "T")
df$nind[is.na(df$nind)] <- 1

gid <- df2genind(select(df, starts_with("L")), pop = df$haps, type = "PA", ploidy = 1) %>% as.genclone()
strata(gid) <- data.frame(pop = pop(gid))
poppr.msn(gid, dist = diss.dist(gid), palette = "funky", vertex.label = "inds")
glue::glue("
Number of samples   : {sum(df$nind) + 1}
Number of haplotypes: {nmll(gid) + 1}
Number of MCG       : {nPop(gid) + 1}
Concordant MCGs     : {group_by(df, haps) %>% summarize(N = sum(nind > 0) == 1) %>% pull(N) %>% sum() + 1}")


txt <- pdf_text("~/Documents/Everhart/differentiation-papers/phillips2002phylogeography.pdf")
txt2 <- paste(txt[4:5], collapse = "\n") %>%
  strsplit("\n") %>% 
  `[[`(1) %>%
  str_extract("^([0-9-]{3}|[ ]{3})[ ]+[0-9]{3} .+$") %>%
  na.omit()
phil <- txt2 %>% 
  strsplit("[[:blank:]]{2,}", perl = TRUE) %>%
  vapply(function(i) paste(i[c(2, 8, 10:11)], collapse = ";"), character(1)) 
pop <- substring(txt2, 1, 3)
for (i in seq(from = 2, to = length(pop))){
  if (trimws(pop[i]) == ""){
    pop[i] <- pop[i - 1]
  }
}

pdf <- paste(pop, phil, sep = ";") %>% 
  paste(collapse = "\n") %>%
  read.table(text = ., header = FALSE, sep = ";") %>%
  setNames(c("Pop", "Strain", "MLH", "MCG", "Fingerprint"))

mlhmcg <- table(MLH = paste(pdf$MLH, pdf$Fingerprint), MCG = pdf$MCG) > 0
rowSums(mlhmcg)
conco <- colSums(mlhmcg)
length(conco)
sum(conco > 1)

# How are the MLHs connected to each other?
library("igraph")
mlhmcg %*% t(mlhmcg) %>% 
  graph_from_adjacency_matrix(diag = FALSE, mode = "undirected", weighted = TRUE) %>% 
  plot.igraph(layout = layout_with_fr, arrow.mode = "-", vertex.size = 3, edge.width = E(.)$weight)

# How are the MCGs connected to each other?
t(mlhmcg) %*% mlhmcg %>% 
  graph_from_adjacency_matrix(diag = TRUE, mode = "undirected", weighted = TRUE) %>% 
  plot.igraph(layout = layout_with_fr, arrow.mode = "-", vertex.size = 3, edge.width = E(.)$weight)

# What about the bipartite graph?
mlhmcg %>% 
  graph_from_incidence_matrix(weighted = TRUE, directed = TRUE, mode = "out") %>% 
  plot.igraph(layout = layout_with_fr, vertex.size = 3, edge.width = E(.)$weight, edge.arrow.size = 0.15)


# Ford 1995
tribble(
  ~strain, ~VCG, ~MCG,
"84.18 ",  "I  ",  "A",
"KAZAZ ",  "I  ",  "B",
"C-11  ",  "I  ",  "B",
"PA-2  ",  "II ",  "B",
"CM-813",  "II ",  "C",
"160   ",  "III",  "D"
) %>% 
  mutate_all(trimws) %>%
  select(-strain) %>%
  table() %T>%
  print() %>%
  graph_from_incidence_matrix(weighted = TRUE, directed = TRUE, mode = "in") %>%
  plot.igraph(edge.arrow.size = 0.5)
