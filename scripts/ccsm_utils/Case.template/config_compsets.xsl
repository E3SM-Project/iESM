<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl  ="http://www.w3.org/1999/XSL/Transform"   version="2.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="config_compset">
  <head>
    <title>Configuration Component Sets</title>
  </head>
  <body>
    <h2>Configuration Component Sets</h2>

    <table BORDER="1" CELLPADDING="10">
      <th>Name (Shortname)</th>
      <th>Description, Including components, CCSM variables</th>
      <tr>
      <th colspan="3">A (All Data Models)</th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'A_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PA_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RA_')]"/>
      <tr>
      <th colspan="3">B (All Active Models)</th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'B_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PB_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RB_')]"/>
      <tr>
      <th colspan="3">C (Standalone POP)</th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'C_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PC_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RC_')]"/>
      <tr>
      <th colspan="3">D (Active sea-ice and ocean, data atmosphere and stub land)</th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'D_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PD_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RD_')]"/>
      <tr>
      <th colspan="3">E (Active land and atmosphere with slab ocean) </th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'E_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PE_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RE_')]"/>
      <tr>
      <th colspan="3">F (Active land and atmosphere with data ocean) </th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'F_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PF_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RF_')]"/>
      <tr>
      <th colspan="3">G (Active sea-ice and ocean, data atmosphere and land) </th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'G_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PG_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RG_')]"/>
      <tr>
      <th colspan="3">H (Active sea-ice and ocean, data atmosphere and stub land) </th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'H_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PH_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RH_')]"/>
      <tr>
      <th colspan="3">I (Standalone CLM) </th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'I_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PI_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RI_')]"/>
      <tr>
      <th colspan="3">S (All stub models with dead atmosphere model) </th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'S_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PS_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RS_')]"/>
      <tr>
      <th colspan="3">T (All stub models with data land model) </th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'T_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PT_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RT_')]"/>
      <tr>
      <th colspan="3">X (All dead models) </th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'X')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'PX_')]"/>
      <xsl:apply-templates select="compset[starts-with(@NAME,'RX_')]"/>
      <tr>
      <th colspan="3">GLC (All compsets from above that have the active land-ice component) </th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'_GLC')]"/>
      <tr>
      <th colspan="3">P (All compsets using WRF regional atmosphere model with CLM) </th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'P')]"/>
      <tr>
      <th colspan="3">R (All compsets using WRF regional atmosphere model with VIC 
                      regional land model) </th>
      </tr>
      <xsl:apply-templates select="compset[starts-with(@NAME,'R')]"/>
    </table>

  </body>
</xsl:template>

<xsl:template match="compset">
  <tr>
    <td><font color="#ff0000"><xsl:value-of select="@NAME"/></font>
    (<xsl:value-of select="@SHORTNAME"/>)
        <xsl:if test="string-length(@GRID_MATCH)>0">
        [Grid=<xsl:value-of select="@GRID_MATCH"/>]
        </xsl:if>
    </td>
    <td><xsl:value-of select="@DESC"/>
    <xsl:if test="string-length(@COMP_ATM)>0">
    <p>
    (atm=<xsl:value-of select="@COMP_ATM"/>
     lnd=<xsl:value-of select="@COMP_LND"/>
     rof=<xsl:value-of select="@COMP_ROF"/>
     glc=<xsl:value-of select="@COMP_GLC"/>
     ice=<xsl:value-of select="@COMP_ICE"/>
     ocn=<xsl:value-of select="@COMP_OCN"/>)
    </p>
    </xsl:if>
        <xsl:if test="string-length(@CCSM_CO2_PPMV)>0">
        CCSM_CO2_PPMV=<xsl:value-of select="@CCSM_CO2_PPMV"/>
        </xsl:if>
        <xsl:if test="string-length(@CCSM_BGC)>0">
        CCSM_BGC=<xsl:value-of select="@CCSM_BGC"/>
        </xsl:if>
        <xsl:if test="string-length(@RUN_REFCASE)>0">
        ref_case=<xsl:value-of select="@RUN_REFCASE"/>
        </xsl:if>
        <xsl:if test="string-length(@RUN_REFCASE)>0">
        ref_date=<xsl:value-of select="@RUN_REFDATE"/>
        </xsl:if>

    </td>
  </tr>
</xsl:template>


</xsl:stylesheet>
