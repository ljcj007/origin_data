"""把导师汇报与后续计划排成中文 PDF。"""
from pathlib import Path
from html import escape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak

ROOT = Path(__file__).resolve().parents[1]


def pdf_text(text):
    """用排版标签替代字体缺失的下标。"""
    text = text.replace('–', '-').replace('—', '-').replace('≪', ' << ')
    return escape(text).replace('₀', '<sub>0</sub>').replace('ₛ', '<sub>s</sub>')


def main():
    font_path = Path('C:/Windows/Fonts/msyh.ttc')
    if font_path.exists():
        pdfmetrics.registerFont(TTFont('Chinese', str(font_path), subfontIndex=0))
        font = 'Chinese'
    else:
        pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
        font = 'STSong-Light'
    body = ParagraphStyle('Body', fontName=font, fontSize=10, leading=16.5,
                          spaceAfter=9, wordWrap='CJK', alignment=TA_LEFT,
                          textColor=colors.HexColor('#263444'))
    heading = ParagraphStyle('Heading', parent=body, fontSize=12, leading=18,
                             spaceBefore=10, spaceAfter=7, textColor=colors.HexColor('#126D78'),
                             keepWithNext=True)
    title = ParagraphStyle('Title', parent=body, fontSize=19, leading=27,
                           spaceAfter=18, textColor=colors.HexColor('#162C46'), keepWithNext=True)
    story = []
    for index, source in enumerate(['teacher_review_20260906.md', 'next_steps_20260906.md']):
        if index:
            story.append(PageBreak())
        for block in (ROOT/'docs'/source).read_text(encoding='utf-8').split('\n\n'):
            block = block.strip()
            if not block:
                continue
            if block == '<!-- pagebreak -->':
                story.append(PageBreak())
            elif block.startswith('# '):
                story.append(Paragraph(pdf_text(block[2:]), title))
            elif block.startswith('## '):
                story.append(Paragraph(pdf_text(block[3:]), heading))
            else:
                story.append(Paragraph(pdf_text(block.replace('\n', ' ')), body))
    target = ROOT/'output/pdf/teacher_brief_20260906.pdf'
    target.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(str(target), pagesize=(595.28, 841.89),
                            leftMargin=46, rightMargin=46, topMargin=46, bottomMargin=44,
                            title='Super-honeycomb 自旋模型：阶段汇报与后续计划',
                            author='研究讨论稿')

    def footer(canvas, document):
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor('#D6DEE5'))
        canvas.line(46, 34, 549, 34)
        canvas.setFont(font, 8)
        canvas.setFillColor(colors.HexColor('#667787'))
        canvas.drawString(46, 22, 'Super-honeycomb | 2026-09-06 | 导师讨论稿')
        canvas.drawRightString(549, 22, str(document.page))
        canvas.restoreState()

    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(target)


if __name__ == '__main__':
    main()
