package sam

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"

	"github.com/exascience/elprep/internal"
	"github.com/exascience/elprep/utils"
)

/*
Parse a field in a header line in a SAM file. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
*/
func (sc *StringScanner) ParseHeaderField() (tag, value string) {
	if sc.err != nil {
		return
	}
	tag, ok := sc.readUntil(':')
	if !ok || (len(tag) != 2) {
		if sc.err == nil {
			sc.err = fmt.Errorf("Invalid field tag %v", tag)
		}
		return "", ""
	} else {
		value, _ = sc.readUntil('\t')
		return tag, value
	}
}

/*
Parse a header line in a SAM file. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.

The @ record type code must have already been scanned. ParseHeaderLine
cannot be used for @CO lines.
*/
func (sc *StringScanner) ParseHeaderLine() utils.StringMap {
	if sc.err != nil {
		return nil
	}
	record := make(utils.StringMap)
	for sc.Len() > 0 {
		tag, value := sc.ParseHeaderField()
		if !record.SetUniqueEntry(tag, value) {
			if sc.err == nil {
				sc.err = fmt.Errorf("Duplicate field tag %v in a SAM header line", tag)
			}
			break
		}
	}
	return record
}

/*
Parse a complete header in a SAM file. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.

Returns a freshly allocated header, the number of header lines, and a
non-nil error value if an error occurred during parsing.
*/
func ParseHeader(reader *bufio.Reader) (hdr *Header, lines int, err error) {
	hdr = NewHeader()
	var sc StringScanner
	for first := true; ; first = false {
		switch data, err := reader.Peek(1); {
		case err == io.EOF:
			return hdr, lines, sc.err
		case err != nil:
			return hdr, lines, err
		case data[0] != '@':
			return hdr, lines, sc.err
		}
		bytes, err := reader.ReadSlice('\n')
		length := len(bytes)
		switch {
		case err == nil:
			length--
		case err != io.EOF:
			return hdr, lines, err
		}
		lines++
		line := string(bytes[4:length])
		sc.Reset(line)
		switch string(bytes[0:4]) {
		case "@HD\t":
			if !first {
				return hdr, lines, errors.New("@HD line not in first line when parsing a SAM header")
			}
			hdr.HD = sc.ParseHeaderLine()
		case "@SQ\t":
			hdr.SQ = append(hdr.SQ, sc.ParseHeaderLine())
		case "@RG\t":
			hdr.RG = append(hdr.RG, sc.ParseHeaderLine())
		case "@PG\t":
			hdr.PG = append(hdr.PG, sc.ParseHeaderLine())
		case "@CO\t":
			hdr.CO = append(hdr.CO, line)
		default:
			switch code := string(bytes[0:3]); {
			case code == "@CO":
				hdr.CO = append(hdr.CO, string(bytes[3:]))
			case IsHeaderUserTag(code):
				if bytes[3] != '\t' {
					return hdr, lines, fmt.Errorf("Header code %v not followed by a tab when parsing a SAM header", code)
				}
				hdr.AddUserRecord(code, sc.ParseHeaderLine())
			default:
				return hdr, lines, fmt.Errorf("Unknown SAM record type code %v", code)
			}
		}
	}
}

/*
Skip the complete header in a SAM file. This is more efficient than
calling ParseHeader and ignoring its result.

Returns the number of header lines and a non-nil error value if an
error occurred.
*/
func SkipHeader(reader *bufio.Reader) (lines int, err error) {
	for {
		data, err := reader.Peek(1)
		if err != nil {
			if err == io.EOF {
				return lines, nil
			} else {
				return lines, err
			}
		}
		if data[0] != '@' {
			break
		}
		for {
			byte, err := reader.ReadByte()
			if err != nil {
				if err == io.EOF {
					return lines, nil
				} else {
					return lines, err
				}
			}
			if byte == '\n' {
				break
			}
		}
		lines++
	}
	return lines, nil
}

func splitHeaderField(field string) (tag, value string, err error) {
	if field[2] != ':' {
		return "", "", fmt.Errorf("Incorrectly formatted SAM file field %v", field)
	}
	return field[:2], field[3:], nil
}

/*
Parse a SAM header line from a string, except that entries are
separated by white space instead of tabulators. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.

The @ record type code must have already been
scanned. ParseHeaderLineFromString cannot be used for @CO lines.
*/
func ParseHeaderLineFromString(line string) (utils.StringMap, error) {
	record := make(utils.StringMap)
	fields := strings.Fields(line)
	for _, field := range fields {
		switch tag, value, err := splitHeaderField(field); {
		case err != nil:
			return record, err
		case !record.SetUniqueEntry(tag, value):
			return record, fmt.Errorf("Duplicate field tag %v in a SAM header line", tag)
		}
	}
	return record, nil
}

/*
All parsers for optional fields in read alignment lines in SAM files
have this signature.
*/
type FieldParser func(*StringScanner, utils.Symbol) (utils.Symbol, interface{})

/*
Parses a single tab-delimited character and returns it as a byte. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
*/
func (sc *StringScanner) ParseChar(tag utils.Symbol) (utils.Symbol, interface{}) {
	if sc.err != nil {
		return tag, nil
	}
	value, _ := sc.readByteUntil('\t')
	return tag, value
}

/*
Parses a single tab-delimited integer and returns it as an int32. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
*/
func (sc *StringScanner) ParseInteger(tag utils.Symbol) (utils.Symbol, interface{}) {
	if sc.err != nil {
		return tag, nil
	}
	value, _ := sc.readUntil('\t')
	val, err := strconv.ParseInt(value, 10, 32)
	if (err != nil) && (sc.err == nil) {
		sc.err = err
	}
	return tag, int32(val)
}

/*
Parses a single tab-delimited float and returns it as a float32. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
*/
func (sc *StringScanner) ParseFloat(tag utils.Symbol) (utils.Symbol, interface{}) {
	if sc.err != nil {
		return tag, nil
	}
	value, _ := sc.readUntil('\t')
	val, err := strconv.ParseFloat(value, 32)
	if (err != nil) && (sc.err == nil) {
		sc.err = err
	}
	return tag, float32(val)
}

/*
Parses a single tab-delimited string and returns it. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
*/
func (sc *StringScanner) ParseString(tag utils.Symbol) (utils.Symbol, interface{}) {
	if sc.err != nil {
		return tag, nil
	}
	value, _ := sc.readUntil('\t')
	switch tag {
	case CC, LB, PG, PU, RG:
		return tag, *utils.Intern(value)
	default:
		return tag, value
	}
}

/*
Parses a byte array in the tab-delimited Hex format and returns it as
a ByteArray. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.5.
*/
func (sc *StringScanner) ParseByteArray(tag utils.Symbol) (utils.Symbol, interface{}) {
	if sc.err != nil {
		return tag, nil
	}
	value, _ := sc.readUntil('\t')
	result := ByteArray(make([]byte, 0, len(value)>>1))
	for i := 0; i < len(value); i += 2 {
		val, err := strconv.ParseUint(value[i:i+2], 16, 8)
		if err != nil {
			if sc.err == nil {
				sc.err = err
			}
			return tag, nil
		}
		result = append(result, byte(val))
	}
	return tag, result
}

/*
Parses a typed, tab-delimited, and comma-separated integer or numeric
array and returns it as a []int8, []uint8, []int16, []uint16, []int32,
[]uint32, or []float32. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.
*/
func (sc *StringScanner) ParseNumericArray(tag utils.Symbol) (utils.Symbol, interface{}) {
	if sc.err != nil {
		return tag, nil
	}
	ntype, ok := sc.readByteUntil(',')
	if !ok {
		if sc.err == nil {
			sc.err = errors.New("Missing entry in numeric array")
		}
		return tag, nil
	}
	switch ntype {
	case 'c':
		var result []int8
		for {
			entry, sep := sc.readUntil2(',', '\t')
			val, err := strconv.ParseInt(entry, 10, 8)
			if err != nil {
				if sc.err == nil {
					sc.err = err
				}
				return tag, nil
			}
			result = append(result, int8(val))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 'C':
		var result []uint8
		for {
			entry, sep := sc.readUntil2(',', '\t')
			val, err := strconv.ParseUint(entry, 10, 8)
			if err != nil {
				if sc.err == nil {
					sc.err = err
				}
				return tag, nil
			}
			result = append(result, uint8(val))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 's':
		var result []int16
		for {
			entry, sep := sc.readUntil2(',', '\t')
			val, err := strconv.ParseInt(entry, 10, 16)
			if err != nil {
				if sc.err == nil {
					sc.err = err
				}
				return tag, nil
			}
			result = append(result, int16(val))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 'S':
		var result []uint16
		for {
			entry, sep := sc.readUntil2(',', '\t')
			val, err := strconv.ParseUint(entry, 10, 16)
			if err != nil {
				if sc.err == nil {
					sc.err = err
				}
				return tag, nil
			}
			result = append(result, uint16(val))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 'i':
		var result []int32
		for {
			entry, sep := sc.readUntil2(',', '\t')
			val, err := strconv.ParseInt(entry, 10, 32)
			if err != nil {
				if sc.err == nil {
					sc.err = err
				}
				return tag, nil
			}
			result = append(result, int32(val))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 'I':
		var result []uint32
		for {
			entry, sep := sc.readUntil2(',', '\t')
			val, err := strconv.ParseUint(entry, 10, 32)
			if err != nil {
				if sc.err == nil {
					sc.err = err
				}
				return tag, nil
			}
			result = append(result, uint32(val))
			if sep != ',' {
				break
			}
		}
		return tag, result
	case 'f':
		var result []float32
		for {
			entry, sep := sc.readUntil2(',', '\t')
			val, err := strconv.ParseFloat(entry, 32)
			if err != nil {
				if sc.err == nil {
					sc.err = err
				}
				return tag, nil
			}
			result = append(result, float32(val))
			if sep != ',' {
				break
			}
		}
		return tag, result
	default:
		if sc.err == nil {
			sc.err = fmt.Errorf("Invalid numeric array type %v", ntype)
		}
		return tag, nil
	}
}

/*
Parses a single tab-delimited mandatory field in a SAM read alignment
line and returns it as a string. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.4.
*/
func (sc *StringScanner) ParseMandatoryField() string {
	s, _ := sc.readUntil('\t')
	return s
}

var optionalFieldParseTable = map[byte]FieldParser{
	'A': (*StringScanner).ParseChar,
	'i': (*StringScanner).ParseInteger,
	'f': (*StringScanner).ParseFloat,
	'Z': (*StringScanner).ParseString,
	'H': (*StringScanner).ParseByteArray,
	'B': (*StringScanner).ParseNumericArray,
}

/*
Parses a single tab-delimited optional field in a SAM read alignment
line and returns it as a tag/value pair. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.5.

The second return value is one of byte (representing an ASCII
character), int32, float32, string, ByteArray, []int8, []uint8,
[]int16, []uint16, []int32, []uint32, or []float32.
*/
func (sc *StringScanner) ParseOptionalField() (tag utils.Symbol, value interface{}) {
	if sc.err != nil {
		return nil, nil
	}
	tagname, ok := sc.readUntil(':')
	if !ok || (len(tagname) != 2) {
		if sc.err == nil {
			sc.err = fmt.Errorf("Invalid field tag %v in SAM alignment line", tagname)
		}
		return nil, nil
	}
	tag = utils.Intern(tagname)
	typebyte, ok := sc.readByteUntil(':')
	if !ok {
		if sc.err == nil {
			sc.err = fmt.Errorf("Invalid field type %v in SAM alignment line", typebyte)
		}
		return nil, nil
	}
	return optionalFieldParseTable[typebyte](sc, tag)
}

func (sc *StringScanner) doString() string {
	if sc.err != nil {
		return ""
	}
	value, ok := sc.readUntil('\t')
	if !ok {
		if sc.err == nil {
			sc.err = errors.New("Missing tabulator in SAM alignment line")
		}
		return ""
	}
	return value
}

func (sc *StringScanner) doInt32() int32 {
	if sc.err != nil {
		return 0
	}
	value, err := strconv.ParseInt(sc.doString(), 10, 32)
	if (err != nil) && (sc.err == nil) {
		sc.err = err
	}
	return int32(value)
}

func (sc *StringScanner) doUint(bitSize int) uint64 {
	if sc.err != nil {
		return 0
	}
	value, err := strconv.ParseUint(sc.doString(), 10, bitSize)
	if (err != nil) && (sc.err == nil) {
		sc.err = err
	}
	return value
}

/*
Parses a read alignment line in a SAM file and returns a freshly allocated alignment. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Sections 1.4 and 1.5.
*/
func (sc *StringScanner) ParseAlignment() *Alignment {
	aln := NewAlignment()

	aln.QNAME = sc.doString()
	aln.FLAG = uint16(sc.doUint(16))
	aln.RNAME = *utils.Intern(sc.doString())
	aln.POS = sc.doInt32()
	aln.MAPQ = byte(sc.doUint(8))
	aln.CIGAR = sc.doString()
	aln.RNEXT = *utils.Intern(sc.doString())
	aln.PNEXT = sc.doInt32()
	aln.TLEN = sc.doInt32()
	aln.SEQ = sc.doString()
	aln.QUAL, _ = sc.readUntil('\t')

	for sc.Len() > 0 {
		aln.TAGS.Set(sc.ParseOptionalField())
	}

	return aln
}

/*
Write a SAM file TAG of type string.
*/
func FormatString(out *bufio.Writer, tag, value string) {
	out.WriteByte('\t')
	out.WriteString(tag)
	out.WriteByte(':')
	out.WriteString(value)
}

/*
Write a header line in a SAM file header section. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
*/
func FormatHeaderLine(out *bufio.Writer, code string, record utils.StringMap) {
	out.WriteString(code)
	for key, value := range record {
		FormatString(out, key, value)
	}
	out.WriteByte('\n')
}

/*
Write a header comment line in a SAM file header section. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
*/
func FormatComment(out *bufio.Writer, code, comment string) {
	out.WriteString(code)
	out.WriteByte('\t')
	out.WriteString(comment)
	out.WriteByte('\n')
}

/*
Write the header section of a SAM file. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.3.
*/
func (hdr *Header) Format(out *bufio.Writer) {
	if hdr.HD != nil {
		FormatHeaderLine(out, "@HD", hdr.HD)
	}
	for _, record := range hdr.SQ {
		FormatHeaderLine(out, "@SQ", record)
	}
	for _, record := range hdr.RG {
		FormatHeaderLine(out, "@RG", record)
	}
	for _, record := range hdr.PG {
		FormatHeaderLine(out, "@PG", record)
	}
	for _, comment := range hdr.CO {
		FormatComment(out, "@CO", comment)
	}
	for code, records := range hdr.UserRecords {
		for _, record := range records {
			FormatHeaderLine(out, code, record)
		}
	}
}

/*
Write a SAM file TAG by appending its ASCII-string representation to
out and returning the result, dispatching on the actual type of the
given value. See http://samtools.github.io/hts-specs/SAMv1.pdf -
Section 1.5.

The following types are accepted: byte (A), int32 (i), float32 (f),
string (Z), ByteArray (H), []int8 (B:c), []uint8 (B:C), []int16 (B:s),
[]uint16 (B:S), []int32 (B:i), []uint32 (B:I), and []float32 (B:f).
*/
func FormatTag(out []byte, tag utils.Symbol, value interface{}) ([]byte, error) {
	out = append(out, '\t')
	out = append(out, *tag...)

	switch val := value.(type) {
	case byte:
		out = append(append(out, ":A:"...), val)
	case int32:
		out = strconv.AppendInt(append(out, ":i:"...), int64(val), 10)
	case float32:
		out = strconv.AppendFloat(append(out, ":f:"...), float64(val), 'g', -1, 32)
	case string:
		out = append(append(out, ":Z:"...), val...)
	case ByteArray:
		out = append(out, ":H:"...)
		for _, b := range val {
			if b < 16 {
				out = append(out, '0')
			}
			out = strconv.AppendUint(out, uint64(b), 16)
		}
	case []int8:
		out = append(out, ":B:c"...)
		for _, v := range val {
			out = strconv.AppendInt(append(out, ','), int64(v), 10)
		}
	case []uint8:
		out = append(out, ":B:C"...)
		for _, v := range val {
			out = strconv.AppendUint(append(out, ','), uint64(v), 10)
		}
	case []int16:
		out = append(out, ":B:s"...)
		for _, v := range val {
			out = strconv.AppendInt(append(out, ','), int64(v), 10)
		}
	case []uint16:
		out = append(out, ":B:S"...)
		for _, v := range val {
			out = strconv.AppendUint(append(out, ','), uint64(v), 10)
		}
	case []int32:
		out = append(out, ":B:i"...)
		for _, v := range val {
			out = strconv.AppendInt(append(out, ','), int64(v), 10)
		}
	case []uint32:
		out = append(out, ":B:I"...)
		for _, v := range val {
			out = strconv.AppendUint(append(out, ','), uint64(v), 10)
		}
	case []float32:
		out = append(out, ":B:f"...)
		for _, v := range val {
			out = strconv.AppendFloat(append(out, ','), float64(v), 'g', -1, 32)
		}
	default:
		return nil, fmt.Errorf("Unknown SAM alignment TAG type %v", value)
	}

	return out, nil
}

/*
Write a SAM file read alignment line by appending its ASCII-string
representation to out and return the result. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Sections 1.4 and 1.5.
*/
func (aln *Alignment) Format(out []byte) ([]byte, error) {
	out = append(append(out, aln.QNAME...), '\t')
	out = append(strconv.AppendUint(out, uint64(aln.FLAG), 10), '\t')
	out = append(append(out, aln.RNAME...), '\t')
	out = append(strconv.AppendInt(out, int64(aln.POS), 10), '\t')
	out = append(strconv.AppendUint(out, uint64(aln.MAPQ), 10), '\t')
	out = append(append(out, aln.CIGAR...), '\t')
	out = append(append(out, aln.RNEXT...), '\t')
	out = append(strconv.AppendInt(out, int64(aln.PNEXT), 10), '\t')
	out = append(strconv.AppendInt(out, int64(aln.TLEN), 10), '\t')
	out = append(append(out, aln.SEQ...), '\t')
	out = append(out, aln.QUAL...)

	var err error
	for _, entry := range aln.TAGS {
		if out, err = FormatTag(out, entry.Key, entry.Value); err != nil {
			return nil, err
		}
	}

	return append(out, '\n'), nil
}

/*
Write a complete SAM file. See
http://samtools.github.io/hts-specs/SAMv1.pdf - Section 1.
*/
func (sam *Sam) Format(out *bufio.Writer) error {
	sam.Header.Format(out)
	buf := internal.ReserveByteBuffer()
	defer internal.ReleaseByteBuffer(buf)
	for _, aln := range sam.Alignments {
		var err error
		buf, err = aln.Format(buf)
		if err != nil {
			return err
		}
		out.Write(buf)
		buf = buf[:0]
	}
	return nil
}

type (
	InputFile struct {
		rc io.ReadCloser
		*bufio.Reader
		*exec.Cmd
	}

	OutputFile struct {
		wc io.WriteCloser
		*bufio.Writer
		*exec.Cmd
	}

	Reader bufio.Reader
	Writer bufio.Writer
)

func (input *InputFile) SamReader() *Reader {
	return (*Reader)(input.Reader)
}

func (output *OutputFile) SamWriter() *Writer {
	return (*Writer)(output.Writer)
}

/*
Open a SAM file for input.

If the filename extension is .bam or .cram,
use samtools view for input. Tell samtools view to only return the
header section for input when headerOnly is true.

samtools must be visible in the directories named by the PATH
environment variable for .bam or .cram input.

If the filename extension is not .bam or .cram, then .sam is always
assumed.

If the name is "/dev/stdin", then the input is read from os.Stdin
*/
func Open(name string, headerOnly bool) (*InputFile, error) {
	switch filepath.Ext(name) {
	case ".bam", ".cram":
		if _, err := os.Stat(name); os.IsNotExist(err) {
			return nil, err
		} else {
			args := []string{"view"}
			if headerOnly {
				args = append(args, "-H")
			} else {
				args = append(args, "-h")
			}
			args = append(args, []string{"-@", strconv.FormatInt(int64(runtime.GOMAXPROCS(0)), 10)}...)
			args = append(args, name)
			cmd := exec.Command("samtools", args...)
			outPipe, err := cmd.StdoutPipe()
			if err != nil {
				return nil, err
			}
			err = cmd.Start()
			if err != nil {
				return nil, err
			}
			return &InputFile{outPipe, bufio.NewReader(outPipe), cmd}, nil
		}
	default:
		if name == "/dev/stdin" {
			return &InputFile{os.Stdin, bufio.NewReader(os.Stdin), nil}, nil
		} else {
			file, err := os.Open(name)
			if err != nil {
				return nil, err
			}
			return &InputFile{file, bufio.NewReader(file), nil}, nil
		}
	}
}

/*
Create a SAM file for output.

If the filename extension is .bam or .cram, use samtools view for
output. If the filename extension is .cram, then either fai or fasta
must be a filename, and the other must be "". If fai is a filename, it
is passed as the -t option to samtools view. If fasta is a filename,
it is passed as the -T option to samtools view.

samtools must be visible in the directories named by the PATH
environment variable for .bam or .cram output.

If the filename extension is not .bam or .cram, then .sam is always
assumed.

If the name is "/dev/stdout", then the output is written to os.Stdout.
*/
func Create(name, fai, fasta string) (*OutputFile, error) {
	switch filepath.Ext(name) {
	case ".bam":
		args := append([]string{"view", "-Sb", "-@"}, strconv.FormatInt(int64(runtime.GOMAXPROCS(0)), 10))
		args = append(args, []string{"-o", name, "-"}...)
		cmd := exec.Command("samtools", args...)
		inPipe, err := cmd.StdinPipe()
		if err != nil {
			return nil, err
		}
		err = cmd.Start()
		if err != nil {
			return nil, err
		}
		return &OutputFile{inPipe, bufio.NewWriter(inPipe), cmd}, nil
	case ".cram":
		if (fai == "") && (fasta == "") {
			log.Fatal("When creating CRAM output, either a reference-fai or a reference-fai must be provided.")
		} else if (fai != "") && (fasta != "") {
			log.Fatal("When creating CRAM output, only either a reference-fasta or a reference-fai must be provided, but not both.")
		}
		args := append([]string{"view", "-C", "-@"}, strconv.FormatInt(int64(runtime.GOMAXPROCS(0)), 10))
		if fai != "" {
			args = append(args, []string{"-t", fai}...)
		} else {
			args = append(args, []string{"-T", fasta}...)
		}
		args = append(args, []string{"-o", name, "-"}...)
		cmd := exec.Command("samtools", args...)
		inPipe, err := cmd.StdinPipe()
		if err != nil {
			return nil, err
		}
		return &OutputFile{inPipe, bufio.NewWriter(inPipe), cmd}, nil
	default:
		if name == "/dev/stdout" {
			return &OutputFile{os.Stdout, bufio.NewWriter(os.Stdout), nil}, nil
		}
		file, err := os.Create(name)
		if err != nil {
			return nil, err
		}
		return &OutputFile{file, bufio.NewWriter(file), nil}, nil
	}
}

/*
Close the SAM input file. If samtools view is used for input, wait for
its process to finish.
*/
func (sam *InputFile) Close() error {
	if sam.rc != os.Stdin {
		err := sam.rc.Close()
		if err != nil {
			return err
		}
	}
	if sam.Cmd != nil {
		err := sam.Wait()
		if err != nil {
			return err
		}
	}
	return nil
}

/*
Close the SAM output file. If samtools view is used for output, wait
for its process to finish.
*/
func (sam *OutputFile) Close() error {
	err := sam.Flush()
	if err != nil {
		return err
	}
	if sam.wc != os.Stdout {
		err := sam.wc.Close()
		if err != nil {
			return err
		}
	}
	if sam.Cmd != nil {
		err := sam.Wait()
		if err != nil {
			return err
		}
	}
	return nil
}
