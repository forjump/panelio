$p = "C:\Users\Zhang\Desktop\clsc_r\paper\panelio_paper.tex"
$t = [IO.File]::ReadAllText($p)

# 1) remove the malformed duplicate line 7 (literal "\n" artifact)
$bad = '\newcommand{\doi}[1]{\href{https://doi.org/#1}{doi:#1}}\n\newcommand{\doi}[1]{\href{https://doi.org/#1}{doi:#1}}'
$t = $t.Replace($bad, '')
# ensure exactly one doi macro remains
if (([regex]::Matches($t, '\\newcommand\{\\doi\}')).Count -ne 1) {
  $t = $t.Replace('\usepackage{graphicx}', '\usepackage{graphicx}' + [Environment]::NewLine + '\newcommand{\doi}[1]{\href{https://doi.org/#1}{doi:#1}}')
}

# 2) replace the sign-recovery sentence by literal IndexOf slicing
$start_marker = 'but it is not near one:'
$end_marker = 'This is a'
$si = $t.IndexOf($start_marker)
$ei = $t.IndexOf($end_marker, $si)
if ($si -ge 0 -and $ei -gt $si) {
  $repl = 'but far from one: the softmax closure redistributes the six simultaneous log-ratio increments among taxa, and the AR(1) innovation noise at the intervention time point is large because the inter-visit intervals at the end of each panel are wide (up to 28 days), so on every design a sizeable fraction of replicates produce observed CLR post-minus-pre differences whose sign opposes the latent truth (lowest sign hit on A4/A7 at 0.42--0.44). '
  $t = $t.Substring(0, $si) + $repl + $t.Substring($ei)
}

[IO.File]::WriteAllText($p, $t, (New-Object Text.UTF8Encoding $false))
$chk = [IO.File]::ReadAllText($p)
Write-Output ("doi macros: " + ([regex]::Matches($chk, '\\newcommand\{\\doi\}')).Count)
Write-Output ("far from one: " + ([regex]::Matches($chk, 'far from one')).Count)
