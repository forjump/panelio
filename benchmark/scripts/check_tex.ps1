$p = "C:\Users\Zhang\Desktop\clsc_r\paper\panelio_paper.tex"
$t = [IO.File]::ReadAllText($p)

Write-Output ("total length: " + $t.Length)
Write-Output ("ends with end doc: " + $t.TrimEnd().EndsWith('\end{document}'))
Write-Output ("CONT markers left: " + ([regex]::Matches($t, '%%%CONT')).Count)

# collect cited keys and defined bibitem keys
$cites = [regex]::Matches($t, '\\cite\{([^}]+)\}') |
  ForEach-Object { $_.Groups[1].Value -split ',' } |
  ForEach-Object { $_.Trim() } | Select-Object -Unique
$bibs = [regex]::Matches($t, '\\bibitem\{([^}]+)\}') | ForEach-Object { $_.Groups[1].Value }

Write-Output "cited keys: $($cites -join ', ')"
$missing = $cites | Where-Object { $_ -notin $bibs }
Write-Output ("missing bibitem for cited: " + $(if ($missing) { $missing -join ', ' } else { 'none' }))

# rough brace balance
$opens = ([regex]::Matches($t, '\{')).Count
$closes = ([regex]::Matches($t, '\}')).Count
Write-Output ("braces open/close: $opens / $closes")

# tables count
Write-Output ("tables: " + ([regex]::Matches($t, '\\begin\{table\}')).Count + "  figures: " + ([regex]::Matches($t, '\\begin\{figure\}')).Count)
