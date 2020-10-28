library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("gghighlight")
library(data.table)
world <- ne_countries(scale = "medium", returnclass = "sf")

lakes = fread("sra_table.csv")
lakes = lakes[lon != "NA"]
lakes = lakes[LIBRARY_SELECTION != "PCR"]
lakes = lakes[abs(lakes$lat)<90,]

lakes_derep = lakes[,.(lon = lon[1], lat = lat[1], nb_bases = as.numeric(sum(nb_bases)), taxon = taxon[1], types = length(levels(factor(taxon))), nb_samples = length(taxon)),  by=coord]
lakes_derep[, sediment := grepl("sediment", taxon)]


lakes_sf = st_as_sf(x = lakes_derep, coords = c("lon", "lat"), crs=  "+proj=longlat")

ggplot(data = world) + geom_sf(fill="floralwhite", alpha=0.5, col="lightgray")+geom_sf(data=lakes_sf, mapping = aes(col=taxon, size = nb_bases), shape=18)+coord_sf(xlim=c(-180, 180), ylim=c(-90,90))+scale_size_continuous(range = c(2, 10))
#  geom_sf_text_repel(data=lakes_sf, aes(label=Lake), size=3, force = 20, nudge_x = -4, seed = 17)+xlab("")+ylab("")
ggsave("lat_map.pdf")

proj="+proj=laea +lat_0=45 +lon_0=-30 +ellps=WGS84 +units=m"
world_pretty = st_transform(world, crs = proj)
lakes_pretty = st_transform(lakes_sf, crs = proj)

zoom = make_box(-100,-10,-170, 20, 50,-15, 120,20)
xwind = c(min(st_coordinates(lakes_pretty)[,"X"]), max(st_coordinates(lakes_pretty)[,"X"]))*1.2
ywind = c(min(st_coordinates(lakes_pretty)[,"Y"]), max(st_coordinates(lakes_pretty)[,"Y"]))*1.2

ggplot(data = world_pretty) + geom_sf(fill="floralwhite", alpha=0.5, col="lightgray")+
  geom_sf(data=lakes_pretty, mapping = aes(col=taxon, size = nb_bases), shape=18)+
  coord_sf(xlim=c(-180, 180), ylim=c(-90,90))+scale_size_continuous(range = c(2, 10))+
  coord_sf(xlim=xwind, ylim=ywind)
ggsave("north_map.pdf")


globe = function(lat, lon, zoom = 1){

  proj=paste0("+proj=laea +lat_0=",lat," +lon_0=",lon, " +ellps=WGS84 +units=m")

  world_pretty = st_transform(world, crs = proj)
  lakes_pretty = st_transform(lakes_sf, crs = proj)

  xwind = c(min(st_coordinates(lakes_pretty)[,"X"]), max(st_coordinates(lakes_pretty)[,"X"]))*zoom
  ywind = c(min(st_coordinates(lakes_pretty)[,"Y"]), max(st_coordinates(lakes_pretty)[,"Y"]))*zoom

  ggplot(data = world_pretty) + geom_sf(fill="floralwhite", alpha=0.5, col="lightgray")+
    geom_sf(data=lakes_pretty, mapping = aes(col=taxon, size = nb_bases), shape=18)+
    scale_size_continuous(range = c(2, 10))+
    coord_sf(xlim=xwind, ylim=ywind, expand=TRUE)

}
