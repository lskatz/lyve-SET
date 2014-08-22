package Graph::Centrality::Pagerank;

require 5.006_000;
use strict;
use warnings;
use Graph;
use Data::Dump qw(dump);

BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
    $VERSION     = '1.05';
    @ISA         = qw(Exporter);
    @EXPORT      = qw();
    @EXPORT_OK   = qw();
    %EXPORT_TAGS = ();
}

=head1 NAME

C<Graph::Centrality::Pagerank> - Computes pagerank of all nodes in a graph.

=head1 SYNOPSIS

  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $ranker = Graph::Centrality::Pagerank->new();
  my $listOfEdges = [[1,2],[3,4]];
  dump $ranker->getPagerankOfNodes (listOfEdges => $listOfEdges);
  # dumps:
  # {
  #   1 => "0.175438596989046",
  #   2 => "0.324561403010954",
  #   3 => "0.175438596989046",
  #   4 => "0.324561403010954",
  # }

=head1 DESCRIPTION

C<Graph::Centrality::Pagerank> computes the pagerank of the all nodes in a graph.
The input can be a list of edges or a L<Graph>. C<Graph::Centrality::Pagerank> is
written entirely in Perl and is not recommended for use in high performance applications.

=head1 CONSTRUCTOR

=head2 C<new>

The method C<new> creates an instance of the C<Graph::Centrality::Pagerank>
class with the following parameters:

=over

=item C<dampeningFactor>

 dampeningFactor => 0.85

C<dampeningFactor> is the dampening factor used when computing pagerank. It
must be a value ranging from zero to one; the default is 0.85. Note the
incidence matrix
generated from the graph is multiplied (scaled) by C<dampeningFactor>, I<not>
by C<1 - dampeningFactor>.

=item C<maxRelError>

 maxRelError => sqrt (machine-epsilon)

C<maxRelError> is the maximum I<average> relative error that is permitted between
successive pagerank vectors before the iterative process to approximate the pagerank vector
should be stopped. The default is the square root of the systems machine epsilon.
Usually, most pagerank values computed will have C<-log10(maxRelError)> digits of accuracy.
C<maxRelError> must be positive and less than or equal to 0.01.

=item C<minIterations>

 minIterations => 0

C<minIterations> is the minimum number of iterations that will be computed
before the pagerank iterations are stopped, even if C<maxRelError> is achieved.
The default is zero.

=item C<maxIterations>

 maxIterations => int (2 * ((maxRelError / ln (dampeningFactor) + 1))

C<maxIterations> is the maximum number of iterations that can be performed
to approximate the pagerank vector even if C<maxRelError> is not achieved.
The default is C<2 * ((maxRelError / ln (dampeningFactor) + 1)>. If
C<dampeningFactor> is zero, then C<maxIterations> is one. If C<dampeningFactor>
is one, then C<maxIterations> is equal to the total nodes in the graph.

=item C<linkSinkNodes>

 linkSinkNodes => 1

In a directed graph sink nodes are the nodes with no edges emanating out
from them. In the pagerank algorithm these nodes are automatically linked
to all other nodes in the graph. To prevent this set C<linkSinkNodes> to zero;
the default is one.

=item C<directed>

 directed => 1

If C<directed> is true, the pagerank computations are done with the graph
edges being directed. If C<directed> is false, the pageranks are computed
treating the graph as undirected; the default value of C<directed> is one.

=item C<useEdgeWeights>

 useEdgeWeights => 0

If C<useEdgeWeights> is true, then any weight associated with an edge is
used in the computation of pagerank. The default weight for any edge without an
assigned weight is one. The default value of C<useEdgeWeights> is zero,
which forces all edge weights to be one.

=back

=cut

sub new
{
  my ($Class, %Parameters) = @_;
  my $Self = bless ({}, ref ($Class) || $Class);

  # set the default parameters.
  my %parameters = $Self->_setDefaultParameters (%Parameters);
  $Self->{defaultParameters} = \%parameters;

  return $Self;
}

sub _setDefaultParameters
{
  my ($Self, %Parameters) = @_;

  # set dampening factor, default is .85;
  $Parameters{dampeningFactor} = 0.85 unless (exists ($Parameters{dampeningFactor}));
  $Parameters{dampeningFactor} = abs ($Parameters{dampeningFactor}) if (exists ($Parameters{dampeningFactor}));
  $Parameters{dampeningFactor} = 1 if ($Parameters{dampeningFactor} > 1);

  # set max relative error, default is sqrt (machine epsilon).
  $Parameters{maxRelError} = sqrt ($Self->_getMachineEpsilon ()) unless (exists ($Parameters{maxRelError}));
  $Parameters{maxRelError} = abs ($Parameters{maxRelError}) if (exists ($Parameters{maxRelError}));
  $Parameters{maxRelError} = 0.01 if ($Parameters{maxRelError} > 0.01);

  # set default min iterations.
  $Parameters{minIterations} = 1 unless (exists ($Parameters{minIterations}));
  $Parameters{minIterations} = abs ($Parameters{minIterations});

  # set default max iterations.
  unless (exists ($Parameters{maxIterations}))
  {
    $Parameters{maxIterations} = 1;
    $Parameters{maxIterations} = 1 if ($Parameters{dampeningFactor} <= 0);
    $Parameters{maxIterationsIsTotalNodes} = 1 if ($Parameters{dampeningFactor} >= 1);
    if ($Parameters{dampeningFactor} < 1)
    {
      $Parameters{maxIterations} = 2 * (int (log ($Parameters{maxRelError}) / log ($Parameters{dampeningFactor})) + 1);
    }
  }
  $Parameters{maxIterations} = abs ($Parameters{maxIterations});
  $Parameters{maxIterations} = 1 if ($Parameters{maxIterations} < 1);
  $Parameters{maxIterations} = $Parameters{minIterations} if ($Parameters{maxIterations} < $Parameters{minIterations});

  # set type of graph, default is directed.
  unless (exists ($Parameters{directed}))
  {
    if (exists ($Parameters{undirected})) { $Parameters{directed} = !$Parameters{undirected}; }
    else { $Parameters{directed} = 1; }
  }

  # set use of edge weights in graph, default is 0.
  $Parameters{useEdgeWeights} = 0 unless (exists ($Parameters{useEdgeWeights}));

  # set linking of sinks nodes in graph, default is 1. if graph is
  # undirected there are no sink nodes.
  $Parameters{linkSinkNodes} = 1 unless (exists ($Parameters{linkSinkNodes}));

  # when forcing dangling nodes to link to all other nodes, this should be
  # with a probability of (1 - dampenFactor) / totalNodes to keep the matrix
  # stocastic; however, some implementation of Pagerank do not do this. the
  # default is 1.
  $Parameters{scaleDampeningFactor} = 1 unless (exists ($Parameters{scaleDampeningFactor}));

  return %Parameters;
}

# gets all the parameters needed to compute the pagerank. uses the
# values set when the object was insantiated to set any missing values.
sub _getAllParameters
{
  my ($Self, %Parameters) = @_;

  # set dampening factor.
  $Parameters{dampeningFactor} = $Self->{defaultParameters}{dampeningFactor} unless (exists ($Parameters{dampeningFactor}));
  $Parameters{dampeningFactor} = abs ($Parameters{dampeningFactor}) if (exists ($Parameters{dampeningFactor}));
  $Parameters{dampeningFactor} = 1 if ($Parameters{dampeningFactor} > 1);

  # set max relative error.
  $Parameters{maxRelError} = $Self->{defaultParameters}{maxRelError} unless (exists ($Parameters{maxRelError}));
  $Parameters{maxRelError} = abs ($Parameters{maxRelError}) if (exists ($Parameters{maxRelError}));
  $Parameters{maxRelError} = 0.01 if ($Parameters{maxRelError} > 0.01);

  # set default min iterations.
  $Parameters{minIterations} = $Self->{defaultParameters}{minIterations} unless (exists ($Parameters{minIterations}));
  $Parameters{minIterations} = abs ($Parameters{minIterations});

  # set default max iterations.
  unless (exists ($Parameters{maxIterations}))
  {
    $Parameters{maxIterations} = $Self->{defaultParameters}{maxIterations};
    if ((0 < $Parameters{dampeningFactor}) && ($Parameters{dampeningFactor} < 1))
    {
      $Parameters{maxIterations} = 2 * (int (log ($Parameters{maxRelError}) / log ($Parameters{dampeningFactor})) + 1);
    }
  }
  $Parameters{maxIterations} = abs ($Parameters{maxIterations});
  $Parameters{maxIterations} = 1 if ($Parameters{maxIterations} < 1);
  $Parameters{maxIterations} = $Parameters{minIterations} if ($Parameters{maxIterations} < $Parameters{minIterations});

  # set the type of graph.
  my $directed;
  $directed = !$Parameters{undirected} if (exists ($Parameters{undirected}));
  $directed = $Parameters{directed} if (exists ($Parameters{directed}));
  $directed = $Parameters{graph}->is_directed () if (!defined ($directed) && exists ($Parameters{graph}));
  $directed = $Self->{defaultParameters}{directed} unless defined $directed;
  $Parameters{directed} = $directed;

  # set use of edge weights in graph, default is 0.
  $Parameters{useEdgeWeights} = $Self->{defaultParameters}{useEdgeWeights} unless (exists ($Parameters{useEdgeWeights}));

  # set linking of sinks nodes in graph, default is 1. if graph is
  # undirected there are no sink nodes.
  $Parameters{linkSinkNodes} = $Self->{defaultParameters}{linkSinkNodes} unless (exists ($Parameters{linkSinkNodes}));

  # when forcing dangling nodes to link to all other nodes, this should be
  # with a probability of (1 - dampenFactor) / totalNodes to keep the matrix
  # stocastic; however, some implementation of Pagerank do not do this. the
  # default is 1.
  $Parameters{scaleDampeningFactor} = $Self->{defaultParameters}{scaleDampeningFactor} unless (exists ($Parameters{scaleDampeningFactor}));

  return %Parameters;
}

=head1 METHOD

=head2 C<getPagerankOfNodes>

The method C<getPagerankOfNodes> computes the pagerank of each node in the graph.
The graph can be defined using the C<graph> parameter or by supplying a list of edges.
All the parameters used by the constructor C<new> can also be set here and they will override
the values used with C<new>. C<getPagerankOfNodes> returns a reference to a hash where the
keys are the graph nodes and the values are the pageranks of the node.

=over

=item C<graph>

 graph => Graph

C<graph> must be a L<Graph> object. If the C<directed> parameter was not set
with the constructor C<new> or with this method, then C<directed>
is set to the value of L<Graph>->L<is_directed|Graph/Accessors>().

=item C<listOfEdges>

 listOfEdges => [['a',10],[10,11]]

C<listOfEdges> must be a list of edges, where an edge is
a pair of strings of the form C<[from-node, to-node]> or a triple of the
form C<[from-node, to-node, numeric-edge-weight]>. Note that C<graph> and C<listOfEdges> can
both be defined, in which case the union of their list of edges is used to compute the
pageranks of the nodes.

=item C<listOfNodes>

 listOfNodes => ['a',10, 'b']

C<listOfNodes> is optional but, must be the list of nodes in the graph when provided;
it defaults to all the nodes comprising the edges in C<listOfEdges> or C<graph>.

=item C<nodeWeights>

  nodeWeights => {}

C<nodeWeights> is an optional hash reference that can provide a weight for the
nodes. If C<nodeWeights> is not defined for any node in the graph, then each
node has a weight of C<1/scale(@listOfNodes)>. If C<nodeWeights> is defined for
at least one node in the graph, then the default weight of any undefined
node is zero.

=back

=cut

sub getPagerankOfNodes
{
  my ($Self, %Parameters) = @_;

  # set any missing parameters to their default values.
  %Parameters = $Self->_getAllParameters (%Parameters);

  # get the list of edges from the graph.
  my @listOfEdges;
  my @listOfNodes;
  if (exists ($Parameters{graph}))
  {
    # get the graph.
    my $graph = $Parameters{graph};

    # get the list of graph edges.
    @listOfEdges = $graph->edges();

    # get the list of vertices.
    @listOfNodes = $graph->vertices();

    # get the graph edge weights if they are to be used and exist.
    if ($Parameters{useEdgeWeights})
    {
      foreach my $edge (@listOfEdges)
      {
        my $weight = $graph->get_edge_weight (@$edge);
        $edge->[2] = $weight if (defined ($weight));
      }
    }
  }

  # add edges from the parameter listOfEdges; assuming they are unique.
  push @listOfEdges, @{$Parameters{listOfEdges}} if exists $Parameters{listOfEdges};
  push @listOfNodes, @{$Parameters{listOfNodes}} if exists $Parameters{listOfNodes};

  return $Self->_getPageranksOfNodesFromEdgeList
    (
      %Parameters,
      listOfEdges => \@listOfEdges,
      listOfNodes => \@listOfNodes
    );
}

sub _getPageranksOfNodesFromEdgeList
{
  my ($Self, %Parameters) = @_;

  # get the list of edges which may have weights, default weight is 1.
  my $listOfEdges = $Parameters{listOfEdges};

  # get the list of edges which may have weights, default weight is 1.
  my $listOfNodes;
  $listOfNodes = $Parameters{listOfNodes} if exists $Parameters{listOfNodes};

  # store the type of graph.
  my $directed = $Parameters{directed};

  # used to build the row oriented matrix from the list of graph edges.
  my %matrixRows;
  my %columnSum;
  my $addEdgeSub = sub
  {
    my ($from, $to, $weight) = @_;

    # add the edge weight to its row.
    unless (exists ($matrixRows{$to}))
    {
      $matrixRows{$to} = {};
      $matrixRows{$to}->{$from} = $weight;
    }
    else
    {
      $matrixRows{$to}->{$from} += $weight;
    }

    # accumulate the column sums, which are used to normalize them later.
    unless (exists ($columnSum{$from}))
    {
      $columnSum{$from} = $weight;
    }
    else
    {
      $columnSum{$from} += $weight;
    }

    # set the column sum of the $to node if not already set, this is used
    # to keep track of sink nodes.
    unless (exists ($columnSum{$to}))
    {
      $columnSum{$to} = 0;
    }
  };

  # convert the list of edges into a row oriented matrix.
  foreach my $edge (@$listOfEdges)
  {
    # get the edge weight, it must be positive.
    my $weight = 1;
    $weight = abs $edge->[2] if (defined ($edge->[2]));
    next if ($weight == 0);

    # add the edge to the matrix of rows.
    &$addEdgeSub ($edge->[0], $edge->[1], $weight);
    &$addEdgeSub ($edge->[1], $edge->[0], $weight) unless ($directed);
  }

  # if $listOfNodes is defined, ensure any missing nodes are added.
  if (defined $listOfNodes)
  {
    foreach my $node (@$listOfNodes)
    {
      unless (exists ($columnSum{$node}))
      {
        $columnSum{$node} = 0;
      }
    }
  }

  # normalize the column sums of the matrix and find all node sinks.
  my @sinkNodes;
  my @rowNodes = keys %matrixRows;
  while (my ($node, $sum) = each %columnSum)
  {
    if ($sum == 0)
    {
      # if the column sum for a node is zero, it is a sink.
      push @sinkNodes, $node;
    }
    else
    {
      # normalize the columns of the node so it sums to 1.
      foreach my $rowNode (@rowNodes)
      {
        if (exists ($matrixRows{$rowNode}->{$node}))
        {
          $matrixRows{$rowNode}->{$node} /= $sum;
        }
      }
    }
  }

  # get the list of all the nodes.
  my @allNodes = keys %columnSum;

  # get the total number of nodes.
  my $totalNodes = @allNodes;

  # get or make the nodeWeights vector.
  my %nodeWeights;
  %nodeWeights = %{$Parameters{nodeWeights}} if exists $Parameters{nodeWeights};

  # ensure $nodeWeights is nonegative for all nodes.
  my $sum = 0;
  foreach my $node (@allNodes)
  {
    # if a value is not defined for a node, make it zero at first.
    $nodeWeights{$node} = 0 unless exists $nodeWeights{$node};

    # force values to be nonnegative.
    $nodeWeights{$node} = -$nodeWeights{$node} if ($nodeWeights{$node} < 0);

    # sum the positive values.
    $sum += $nodeWeights{$node};
  }

  # ensure $nodeWeights sum to one.
  if ($sum > 0)
  {
    foreach my $node (@allNodes) { $nodeWeights{$node} /= $sum; }
  }
  else
  {
    foreach my $node (@allNodes) { $nodeWeights{$node} = 1 / $totalNodes; }
  }


  # initialize the pagerank vector;
  my $pagerank = {};
  return $pagerank if ($totalNodes == 0);
  foreach my $node (@allNodes) { $pagerank->{$node} = $nodeWeights{$node}; }
  my $newPageRank = {};

  # set the maximum number of iterations.
  my $maxIterations = $Parameters{maxIterations};
  $maxIterations = $totalNodes if exists $Parameters{maxIterationsIsTotalNodes};

  for (my $iteration = 0; $iteration < $maxIterations; $iteration++)
  {
    # first set the new page rank to the average pageranks times (1 - $dampeningFactor).
    if ($Parameters{scaleDampeningFactor})
    {
      #my $sum = 0;
      #foreach my $node (@allNodes) { $sum += $pagerank->{$node}; }
      #$sum *= (1 - $Parameters{dampeningFactor}) / $totalNodes;
      foreach my $node (@allNodes) { $newPageRank->{$node} = (1 - $Parameters{dampeningFactor}) * $nodeWeights{$node}; }
    }
    else
    {
      foreach my $node (@allNodes) { $newPageRank->{$node} = 1 - $Parameters{dampeningFactor}; }
    }

    # add in the values for the sink nodes.
    if ($Parameters{linkSinkNodes})
    {
      my $sinkSum = 0;
      foreach my $sinkNode (@sinkNodes)
      {
        $sinkSum += $pagerank->{$sinkNode};
      }
      $sinkSum *= $Parameters{dampeningFactor} / $totalNodes;
      foreach my $node (@allNodes) { $newPageRank->{$node} += $sinkSum; }
    }

    # add in the rank from the graph links.
    foreach my $rowNode (@rowNodes)
    {
      my $sum = 0;
      while (my ($colNode, $value) = each %{$matrixRows{$rowNode}})
      {
        $sum += $value * $pagerank->{$colNode};
      }
      $newPageRank->{$rowNode} += $Parameters{dampeningFactor} * $sum;
    }

    # normalize (for rounding error stability -- i hope) then swap pagerank vectors.
    _normalizeByNorm1 ($newPageRank, \@allNodes);
    ($pagerank, $newPageRank) = ($newPageRank, $pagerank);

    # compute the error.
    my $error = 0;
    my $totalNonzero = 0;
    foreach my $node (@allNodes)
    {
      if ($pagerank->{$node} != 0)
      {
        $error += abs (($newPageRank->{$node} - $pagerank->{$node}) / $pagerank->{$node});
      }
      else
      {
        $error += abs ($newPageRank->{$node} - $pagerank->{$node});
      }
    }
    $error /= $totalNodes;
    #print $iteration . ' ' . $error . "\n";

    # stop iterating if the error is small enough, unless the minIterations has
    # not been reached.
    last if (($error < $Parameters{maxRelError}) && ($iteration >= $Parameters{minIterations}));
  }

  # normalize the pagerank vector (for rounding error stability -- i hope).
  _normalizeByNorm1 ($pagerank, \@allNodes);

  return $pagerank;
}

# normalize $hashVector so it sums to one (or zero).
sub _normalizeByNorm1
{
  my ($hashVector, $indices) = @_;
  my $sum = 0;
  foreach my $node (@$indices) { $sum += $hashVector->{$node}; }
  $sum = 1 if ($sum == 0);
  foreach my $node (@$indices) { $hashVector->{$node} /= $sum; }
  return 1;
}

# get the value of machine epsilon, sort of, really
# the unit roundoff value.
sub _getMachineEpsilon
{
  my $one = 1;
  my $epsilon = 2;
  my $halfOfepsilon = 1;
  my $powerOf2 = 0;
  my $sum;
  do
  {
    $epsilon = $halfOfepsilon;
    $halfOfepsilon = $epsilon / 2;
    $sum = 1 + $halfOfepsilon;
    ++$powerOf2;
  }
  until (($sum == $one) || ($powerOf2 > 2048)) ;
  return $epsilon;
}

=head1 EXAMPLES

A rather dull example with one node and no edges:

  use Graph;
  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $ranker = Graph::Centrality::Pagerank->new();
  my $listOfNodes = [1];
  dump $ranker->getPagerankOfNodes (listOfNodes => $listOfNodes);
  # dumps:
  # {
  #   1 => 1
  # }

An example of a graph with two components:

  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $ranker = Graph::Centrality::Pagerank->new();
  my $listOfEdges = [[1,2],[3,4]];
  dump $ranker->getPagerankOfNodes (listOfEdges => $listOfEdges);
  # dumps:
  # {
  #   1 => "0.175438596989046",
  #   2 => "0.324561403010954",
  #   3 => "0.175438596989046",
  #   4 => "0.324561403010954",
  # }

In this case the edges are placed in a L<Graph>:

  use Graph;
  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $listOfEdges = [[1,2],[3,4]];
  my $graph = Graph->new (edges => $listOfEdges);
  my $ranker = Graph::Centrality::Pagerank->new();
  dump $ranker->getPagerankOfNodes (graph => $graph);
  # dumps:
  # {
  #   1 => "0.175438596989046",
  #   2 => "0.324561403010954",
  #   3 => "0.175438596989046",
  #   4 => "0.324561403010954",
  # }

Below is the first example in the paper
I<How Google Finds Your Needle in the Web's Haystack> by David Austin.

  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $ranker = Graph::Centrality::Pagerank->new();
  my $listOfEdges = [[1,2],[1,3],[2,4],[3,2],[3,5],[4,2],[4,5],[4,6],[5,6],
    [5,7],[5,8],[6,8],[7,5],[7,1],[7,8],[8,6],[8,7]];
  dump $ranker->getPagerankOfNodes (listOfEdges => $listOfEdges,
    dampeningFactor => 1);
  # dumps:
  # {
  #   1 => "0.0599999994835539",
  #   2 => "0.0675000002254998",
  #   3 => "0.0300000002967361",
  #   4 => "0.0674999997408677",
  #   5 => "0.0974999994123176",
  #   6 => "0.202500001447512",
  #   7 => "0.180000001348251",
  #   8 => "0.294999998045262",
  # }

Below is the second example in the paper. Notice C<linkSinkNodes> is
set to zero.

  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $ranker = Graph::Centrality::Pagerank->new();
  my $listOfEdges = [[1,2]];
  dump $ranker->getPagerankOfNodes (listOfEdges => $listOfEdges,
    dampeningFactor => 1, linkSinkNodes => 0);
  # dumps:
  # { 1 => 0, 2 => 0 }

Below is the third example in the paper. Notice in this case
C<linkSinkNodes> is set to one, the default value.

  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $ranker = Graph::Centrality::Pagerank->new();
  my $listOfEdges = [[1,2]];
  dump $ranker->getPagerankOfNodes (listOfEdges => $listOfEdges,
    dampeningFactor => 1, linkSinkNodes => 1);
  # dumps:
  # { 1 => "0.33333333209157", 2 => "0.66666666790843" }

Below is the fourth example in the paper. The
result is different from the paper since the starting vector for
L<Graph::Centrality::Pagerank> is

  { 1 => "0.2", 2 => "0.2", 3 => "0.2", 4 => "0.2", 5 => "0.2" }

while the starting vector in the paper is

  { 1 => 1, 2 => 0, 3 => 0, 4 => 0, 5 => 0 }.

  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $ranker = Graph::Centrality::Pagerank->new();
  my $listOfEdges = [[1,2],[2,3],[3,4],[4,5],[5,1]];
  dump $ranker->getPagerankOfNodes (listOfEdges => $listOfEdges,
    dampeningFactor => 1, linkSinkNodes => 0);
  # dumps:
  # { 1 => "0.2", 2 => "0.2", 3 => "0.2", 4 => "0.2", 5 => "0.2" }

Below is the fifth example in the paper.

  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $ranker = Graph::Centrality::Pagerank->new();
  my $listOfEdges = [[1,3],[1,2],[2,4],[3,2],[3,5],[4,2],[4,5],[4,6],[5,6],
    [5,7],[5,8],[6,8],[7,5],[7,8],[8,6],[8,7]];
  dump $ranker->getPagerankOfNodes (listOfEdges => $listOfEdges,
    dampeningFactor => 1, linkSinkNodes => 0);
  # dumps:
  # {
  #   1 => 0,
  #   2 => "2.39601089109228e-54",
  #   3 => 0,
  #   4 => "5.47659632249665e-54",
  #   5 => "0.119999999997811",
  #   6 => "0.240000000003975",
  #   7 => "0.240000000003975",
  #   8 => "0.399999999994238",
  # }

An example of the effect of including edge weights:

  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $ranker = Graph::Centrality::Pagerank->new();
  my $listOfEdges = [[2,1],[2,3]];
  dump $ranker->getPagerankOfNodes (listOfEdges => $listOfEdges);
  $listOfEdges = [[2,1,2],[2,3,1]];
  dump $ranker->getPagerankOfNodes (listOfEdges => $listOfEdges,
    useEdgeWeights => 1);

  # dumps:
  # all edges have weight 1.
  # {
  #   1 => "0.370129870353883",
  #   2 => "0.259740259292235",
  #   3 => "0.370129870353883",
  # }
  # edge [2, 1] has twice the weight of edge [2,3].
  # {
  #   1 => "0.406926407374432",
  #   2 => "0.259740259292235",
  #   3 => "0.333333333333333",
  # }

An example of the effect of including node weights:

  use Graph;
  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $ranker = Graph::Centrality::Pagerank->new();
  my $listOfEdges = [[1,2],[2,3]];
  dump $ranker->getPagerankOfNodes (listOfEdges => $listOfEdges);
  dump $ranker->getPagerankOfNodes (listOfEdges => $listOfEdges,
    nodeWeights => {2 => .9, 3 => .1 });

  # dumps:
  # {
  #   1 => "0.184416783248514",
  #   2 => "0.341171047056969",
  #   3 => "0.474412169694517",
  # }
  # {
  #   1 => "0.135592438389592",
  #   2 => "0.385846009631034",
  #   3 => "0.478561551979374",
  # }

A example of the modules speed, or lack of.

  use Graph::Centrality::Pagerank;
  use Data::Dump qw(dump);
  my $ranker = Graph::Centrality::Pagerank->new();
  my @listOfEdges;
  for (my $i = 0; $i < 1000000; $i++)
    { push @listOfEdges, [int rand 10000, int rand 10000]; }
  my $startTime = time;
  my $pageranks = $ranker->getPagerankOfNodes (listOfEdges => \@listOfEdges);
  print time()-$startTime . "\n";
  # prints:
  # a non-negative integer after a long time.

=head1 INSTALLATION

To install the module run the following commands:

  perl Makefile.PL
  make
  make test
  make install

If you are on a windows box you should use 'nmake' rather than 'make'.

=head1 BUGS

Please email bugs reports or feature requests to C<bug-graph-centrality-pagerank@rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Graph-Centrality-Pagerank>.  The author
will be notified and you can be automatically notified of progress on the bug fix or feature request.

=head1 AUTHOR

 Jeff Kubina<jeff.kubina@gmail.com>

=head1 COPYRIGHT

Copyright (c) 2009 Jeff Kubina. All rights reserved.
This program is free software; you can redistribute
it and/or modify it under the same terms as Perl itself.

The full text of the license can be found in the
LICENSE file included with this module.

=head1 KEYWORDS

centrality measure, eigenvector centrality, graph, network, pagerank

=head1 SEE ALSO

=begin html

<p>A good tutorial on the <a href="http://en.wikipedia.org/wiki/PageRank">pagerank</a> algorithm is the article
<a href="http://www.ams.org/featurecolumn/archive/pagerank.html">How Google Finds Your Needle in the Web's Haystack</a> by David Austin.</p>

=end html

L<Graph>

=begin html

<a href="http://en.wikipedia.org/wiki/Centrality">centrality measure</a>,
<a href="http://en.wikipedia.org/wiki/Centrality#Eigenvector_centrality">eigenvector centrality</a>,
<a href="http://en.wikipedia.org/wiki/Graph_(mathematics)">graph</a>,
<a href="http://en.wikipedia.org/wiki/Network_theory">network</a>,
<a href="http://en.wikipedia.org/wiki/PageRank">pagerank</a>

=end html

=cut

1;
# The preceding line will help the module return a true value
