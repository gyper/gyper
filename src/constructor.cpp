#include <gyper/constructor.hpp>

namespace gyper
{

/*********
 * BASIC *
 *********/

Constructor::Constructor ()
{
	TGraph graph();
}


Constructor::Constructor (Options & CO)
{
  TGraph graph();
  Constructor::Constructor::CO = CO; // Replaces default options
}

void
Constructor::set_genomic_region(const char * region)
{
  seqan::parse(genomic_region, region);
  SEQAN_ASSERT_MSG(fasta_set, "You need to read reference genome before setting region in it.");
  bool ret = seqan::getIdByName(genomic_region.rID, fasta_index, genomic_region.seqName);
  SEQAN_ASSERT_MSG(ret, "No sequence in FASTA index with id = %s", genomic_region.seqName);

  // Fix invalid positions after parsing, this should probably be a part of seqan::readRegion(..) to do this though.
  if (static_cast<unsigned>(genomic_region.endPos) > seqan::sequenceLength(fasta_index, genomic_region.rID) ||
      genomic_region.endPos == seqan::GenomicRegion::INVALID_POS
     )
  {
    genomic_region.endPos = seqan::sequenceLength(fasta_index, genomic_region.rID);
  }

  if (genomic_region.beginPos == seqan::GenomicRegion::INVALID_POS ||
      genomic_region.beginPos < 0
     )
  {
    genomic_region.beginPos = 0;
  }
}


/*************
 * FASTA I/O *
 *************/
bool
Constructor::read_reference_genome(const char * fasta_filename)
{
  if (fasta_set)
  {
    std::cerr << "Warning: A previously read FASTA will be discarded" << std::endl;
    seqan::clear(fasta_index);
    fasta_set = false;
  }

  if (!seqan::open(fasta_index, fasta_filename))
  {
    if (!seqan::build(fasta_index, fasta_filename))
    {
      std::cerr << "ERROR: Index could not be loaded or built." << std::endl;
      return false;
    }

    if (!seqan::save(fasta_index))
    {
      std::cerr << "ERROR: Index could not be written do disk." << std::endl;
      return false;
    }
  }

  fasta_set = true;
  return true;
}

void
Constructor::extract_reference_sequence(const char * region)
{
  SEQAN_ASSERT_MSG(seqan::numSeqs(fasta_index) != 0, "Please read the reference genome before extr like the FASTA index has not been read.");

  set_genomic_region(region);
  extract_reference_sequence();
}

void
Constructor::extract_reference_sequence(void)
{
  SEQAN_ASSERT_MSG(seqan::numSeqs(fasta_index) != 0, "Please read the reference genome before extr like the FASTA index has not been read.");
  SEQAN_ASSERT_MSG(genomic_region.endPos != seqan::GenomicRegion::INVALID_POS, "Please set region before extracting reference sequence.");

  seqan::readRegion(
    reference_sequence,
    fasta_index,
    genomic_region.rID,
    genomic_region.beginPos,
    genomic_region.endPos
  );
}


/***********
 * VCF I/O *
 ***********/
void
Constructor::open_tabix(const char * tabix_filename)
{
  seqan::open(tabix_file, tabix_filename);
  seqan::String<char> header;
  seqan::getHeader(header, tabix_file);
  vcf_set = true;

  if (genomic_region.endPos != seqan::VcfRecord::INVALID_POS)
  {
    // This means the genomic region has been set
    seqan::String<char> region_str;
    genomic_region.toString(region_str);
    // std::cout << "region_str = " << region_str << std::endl;
    seqan::setRegion(tabix_file, seqan::toCString(region_str));
  }
}

void
Constructor::open_tabix(const char * tabix_filename, const char * region)
{
  SEQAN_ASSERT_MSG(genomic_region.endPos == seqan::VcfRecord::INVALID_POS, "Genomic region was already set.");
  open_tabix(tabix_filename);
  seqan::setRegion(tabix_file, region);
}

bool
Constructor::read_tabix_record()
{
  return seqan::readRecord(vcf_record, tabix_file);
}

bool
Constructor::read_first_tabix_region()
{
  bool success = seqan::readRegion(vcf_record, tabix_file);

  while(vcf_record.beginPos < genomic_region.beginPos)
  {
    success = seqan::readRegion(vcf_record, tabix_file);
  }

  return success;
}

bool
Constructor::read_tabix_region()
{
  return seqan::readRegion(vcf_record, tabix_file);
}


/**********************
 * GRAPH CONSTRUCTION *
 **********************/
TVertex
Constructor::insert_reference_vertex(unsigned const & order, seqan::String<seqan::Dna> const & value)
{
  // std::cout << "Adding order, value = " << order << ", " << value << std::endl;
  TVertex new_vertex = seqan::addVertex(graph);
  VertexLabel new_vertex_label(order, value);
  vertex_label_map[new_vertex_label] = new_vertex;
  vertex_labels.push_back(new_vertex_label);
  return new_vertex;
}

TVertex
Constructor::insert_reference_vertex(unsigned const & order, seqan::String<seqan::Dna> const & value, TVertex const & previous_vertex)
{
  TVertex new_vertex = insert_reference_vertex(order, value);
  seqan::addEdge(graph, previous_vertex, new_vertex);
  return new_vertex;
}

void
Constructor::add_first_reference_sequence()
{
  if (reference_sequence[0] != 'N')
  {
    vertex_head = insert_reference_vertex(0, "");
    return;
  }

  unsigned starting_pos = 0;

  while (reference_sequence[starting_pos] == 'N' && starting_pos < static_cast<unsigned>(genomic_region.endPos))
  {
    ++starting_pos;
  }

  vertex_head = insert_reference_vertex(starting_pos, "");
}

void
Constructor::add_last_reference_sequence()
{
  add_reference_sequence_preceding_a_point(vertex_labels[vertex_head].order);
}

void
Constructor::add_reference_sequence_preceding_a_point(unsigned const & point)
{
  SEQAN_ASSERT_MSG(point <= seqan::length(reference_sequence) + genomic_region.beginPos,
                   "Point = %d, ref_seq ends at %d",
                   point,
                   seqan::length(reference_sequence) + genomic_region.beginPos
                  );

  // add_first_reference_sequence_preceding_a_point(point);
  bool sequence_started = true;
  unsigned starting_pos = vertex_labels[vertex_head].order + seqan::length(vertex_labels[vertex_head].dna);
  seqan::String<seqan::Dna> current_dna = "";

  for (unsigned pos = starting_pos; pos < point - genomic_region.beginPos; ++pos)
  {
    if (reference_sequence[pos] == 'N')
    {
      if (sequence_started)
      {
        vertex_head = insert_reference_vertex(starting_pos, current_dna, vertex_head);
        starting_pos = pos;
        current_dna = "";
        sequence_started = false;
      }

      continue;
    }

    if (!sequence_started)
    {
      vertex_head = insert_reference_vertex(pos, "", vertex_head);
      starting_pos = pos;
      sequence_started = true;
    }
    
    appendValue(current_dna, reference_sequence[pos]);
  }

  if (sequence_started)
  {
    vertex_head = insert_reference_vertex(starting_pos, current_dna, vertex_head);
  }
  else
  {
    // The sequence ends with a 'N'
    // std::cout << "point = " << point << std::endl;
    vertex_head = insert_reference_vertex(point - genomic_region.beginPos, current_dna, vertex_head);
  }
  
}

void
Constructor::add_sequence_preceding_a_vcf_record(void)
{
  SEQAN_ASSERT_MSG(vcf_record.beginPos != seqan::VcfRecord::INVALID_POS, "No VCF record.");
  add_reference_sequence_preceding_a_point(vcf_record.beginPos);
}

void
Constructor::add_vcf_record_to_graph(void)
{
  unsigned const & order = vcf_record.beginPos;
  insert_reference_vertex(order, static_cast<seqan::String<seqan::Dna> >(vcf_record.ref), vertex_head);

  seqan::StringSet<seqan::String<char> > alternatives;
  seqan::strSplit(alternatives, vcf_record.alt, seqan::EqualsChar<','>());

  for (auto alt_it = seqan::begin(alternatives); !seqan::atEnd(alt_it); ++alt_it)
  {
    std::cout << *alt_it << std::endl;
  }
  
}

void
Constructor::construct_graph()
{
  
}

} // namespace gyper
