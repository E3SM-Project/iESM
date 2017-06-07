<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="config_definition">
  <head>
    <title>Configuration Variable Definition</title>
  </head>
  <body>

    <h2>Case Definition Variables [can NOT be changed after case defined]</h2>

    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name (type-of-variable)
        <p>[default=default_value_if_exists]</p></th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid-Values</th>
      </tr>
      <xsl:apply-templates select="entry[contains(@group,'case_')]"/>
    </table>

    <h2>Build Variables (env_build.xml)</h2>

    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name (type-of-variable)
        <p>[default=default_value_if_exists]</p></th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid-Values</th>
      </tr>
      <xsl:apply-templates select="entry[contains(@group,'build_')]"/>
    </table>

    <h2>Machine Variables (env_machine*)</h2>

    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name (type-of-variable)
        <p>[default=default_value_if_exists]</p></th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid-Values</th>
      </tr>
      <xsl:apply-templates select="entry[contains(@group,'mach_')]"/>
    </table>

    <h2>Run Variables (env_run.xml) [can be changed at run-time]</h2>

    <table BORDER="1" CELLPADDING="10">
      <tr>
      <th rowspan="2">Name (type-of-variable)
        <p>[default=default_value_if_exists]</p></th>
      <th>Description</th>
      </tr>
      <tr>
      <th>Valid-Values</th>
      </tr>
      <xsl:apply-templates select="entry[contains(@group,'run_')]"/>
    </table>

  </body>
</xsl:template>

<xsl:template match="entry">
  <tr>
    <td rowspan="2"><font color="#ff0000"><xsl:value-of select="@id"/> 
        (<xsl:value-of select="@type"/>)</font>
        <xsl:if test="string-length(@value)>0 and @value != 'UNSET' and @value !='/UNSET'
and @value != 'null' and @value != 'NULL'">
        <p>[default=<xsl:value-of select="@value"/>]</p>
        </xsl:if>
    </td>
    <td><xsl:value-of select="@ldesc"/></td>
  </tr>
  <tr>
    <td><xsl:value-of select="@valid_values"/></td>
  </tr>
</xsl:template>


</xsl:stylesheet>
