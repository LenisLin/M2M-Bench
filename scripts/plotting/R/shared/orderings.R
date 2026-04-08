m2m_orderings <- list(
  representation = c("Gene", "Pathway", "FM"),
  analysis_family = c("group_concordance", "retrieval", "cross_eligibility"),
  direction = c("ALL", "C2G", "G2C"),
  scope_status = c("materialized", "available", "excluded_by_support_gate", "not_applicable_scope", "historical_only")
)

representation_class <- function(values) {
  ifelse(values %in% c("Gene", "Pathway", "FM"), values, "FM")
}

apply_representation_order <- function(data, column = "representation") {
  data[[column]] <- factor(representation_class(data[[column]]), levels = m2m_orderings$representation)
  data
}
