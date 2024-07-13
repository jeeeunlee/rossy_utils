#include "yaml/node/node.h"
#include "yaml/nodebuilder.h"
#include "yaml/nodeevents.h"

namespace YAML {
Node Clone(const Node& node) {
  NodeEvents events(node);
  NodeBuilder builder;
  events.Emit(builder);
  return builder.Root();
}
}
