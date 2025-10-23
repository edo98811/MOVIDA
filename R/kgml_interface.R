#' Parse KEGG KGML files to extract relations and edges data frames.
#' 
#' @param file Path to the KGML XML file.
#' @return A tibble with columns: from, to, type, name, value.
#' @importFrom xml2 read_xml xml_find_all xml_attr
#' @importFrom tibble tibble
#' @importFrom purrr map_dfr
parse_kgml_relations <- function(file) {
  doc <- read_xml(file)

  rels <- xml_find_all(doc, ".//relation")

  purrr::map_dfr(rels, function(rel) {
    entry1 <- xml_attr(rel, "entry1")
    entry2 <- xml_attr(rel, "entry2")
    type <- xml_attr(rel, "type")
    subnodes <- xml_find_all(rel, ".//subtype")

    if (length(subnodes) == 0) {
      tibble(
        from = entry1,
        to = entry2,
        type = type,
        subtype = NA_character_,
        value = NA_character_
      )
    } else {
      tibble(
        from = entry1,
        to = entry2,
        type = type,
        subtype = xml_attr(subnodes, "name"),
        value = xml_attr(subnodes, "value")
      )
    }
  })
}


#' Parse KEGG KGML files to extract nodes and edges data frames.
#' 
#' @param file Path to the KGML XML file.
#' @return A tibble with columns: id, name, type, link, reaction, graphics_name, fgcolor, bgcolor, graphics_type, x, y, width, height.
#' @importFrom xml2 read_xml xml_find_all xml_attr
#' @importFrom tibble tibble
#' @importFrom purrr map_dfr
parse_kgml_entries <- function(file) {
  # read the KGML file
  doc <- read_xml(file)

  # find all entry nodes
  entries <- xml_find_all(doc, ".//entry")
  graphics <- xml_find_all(entries, ".//graphics")

  # extract attributes into a data.frame
  df <- tibble(
    id = xml_attr(entries, "id"),
    name = xml_attr(entries, "name"),
    type = xml_attr(entries, "type"),
    link = xml_attr(entries, "link"),
    reaction = xml_attr(entries, "reaction"),
    graphics_name = xml_attr(graphics, "name"),
    label = xml_attr(graphics, "name"), # for visNetwork
    fgcolor = xml_attr(graphics, "fgcolor"),
    bgcolor = xml_attr(graphics, "bgcolor"),
    graphics_type = xml_attr(graphics, "type"),
    x = xml_attr(graphics, "x"),
    y = xml_attr(graphics, "y"),
    width = xml_attr(graphics, "width"),
    height = xml_attr(graphics, "height")
  )
  # Example entry:
  #     <entry id="184" name="ko:K15359 ko:K18276" type="ortholog" reaction="rn:R09472"
  # link="https://www.kegg.jp/dbget-bin/www_bget?K15359+K18276">
  # <graphics name="K15359..." fgcolor="#000000" bgcolor="#FFFFFF"
  #      type="rectangle" x="303" y="561" width="46" height="17"/>

  return(df)
}