<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl  ="http://www.w3.org/1999/XSL/Transform"   version="2.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="config_horiz_grid">
  <head>
    <title>CESM Grids</title>
  </head>
  <body>
    <h2>CESM Grids</h2>

    <p></p>
    <table BORDER="2" CELLPADDING="5">
    <thead>
      <tr>
      <th>Name (Shortname)</th>
      <th>Description</th>
      <th>Atm</th>
      <th>Ice</th>
      <th>Lnd</th>
      <th>Ocn</th>
      </tr>
    </thead>
    <tbody>
      <xsl:apply-templates select="horiz_grid"/>
    </tbody>
    </table>

  </body>
</xsl:template>

<xsl:template match="horiz_grid">
<xsl:if test="string-length(@GRID)>0">
  <tr>
    <td><font color="#ff0000"><xsl:value-of select="@GRID"/></font>
    (<xsl:value-of select="@SHORTNAME"/>)
    </td>
    <td>
    <xsl:value-of select="@DESC"/>
    <xsl:if test="string-length(@VALID_COMPSET_MATCH)>0">
    <p>
    (Only valid for the following type of compsets: <xsl:value-of select="@VALID_COMPSET_MATCH"/>)
    </p>
    </xsl:if>
    </td>
    <td>
    <xsl:value-of select="@ATM_GRID"/>
    </td>
    <td>
    <xsl:value-of select="@ICE_GRID"/>
    </td>
    <td>
    <xsl:value-of select="@LND_GRID"/>
    </td>
    <td>
    <xsl:value-of select="@OCN_GRID"/>
    </td>
  </tr>
</xsl:if>
</xsl:template>


</xsl:stylesheet>
