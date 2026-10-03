$csv = "C:\Users\Zhang\Desktop\clsc_r\benchmark\results"

# read summary CSVs as data
function Get-Rows($f) {
  Import-Csv (Join-Path $csv $f)
}
$sim = Get-Rows "summary_simulation.csv"
$abl = Get-Rows "summary_ablation.csv"
$senn = Get-Rows "summary_sensitivity_normalize.csv"
$senb = Get-Rows "summary_sensitivity_band.csv"

$checks = @()

# Table 2: A8 pooled_ols power 0.267, A12 pooled 0.272, A9 scm 0.172
$r = $sim | Where-Object { $_.dataset -eq "A8" -and $_.method -eq "pooled_ols" }
$checks += "A8 pooled power = $($r.power_fdr) (tex: 0.267)"
$r = $sim | Where-Object { $_.dataset -eq "A12" -and $_.method -eq "pooled_ols" }
$checks += "A12 pooled power = $($r.power_fdr) (tex: 0.272)"
$r = $sim | Where-Object { $_.dataset -eq "A9" -and $_.method -eq "scm" }
$checks += "A9 scm power = $($r.power_fdr) (tex: 0.172)"

# Table 3: A11 typeI
foreach ($m in @("did","pooled_ols","event_study","scm")) {
  $r = $sim | Where-Object { $_.dataset -eq "A11" -and $_.method -eq $m }
  $checks += "A11 $m typeI = $($r.typeI_fdr)"
}

# Table 5: A6 no_offset power 0.028, A11 nnls typeI 0.186, A12 nnls weight_sum
$r = $abl | Where-Object { $_.dataset -eq "A6" -and $_.variant -eq "ab_no_offset" }
$checks += "A6 no_offset power = $($r.power_fdr) (tex: 0.028)"
$r = $abl | Where-Object { $_.dataset -eq "A11" -and $_.variant -eq "ab_nnls" }
$checks += "A11 nnls typeI = $($r.typeI_fdr) (tex: 0.186)"
$r = $abl | Where-Object { $_.dataset -eq "A12" -and $_.variant -eq "ab_nnls" }
$checks += "A12 nnls weight_sum = $($r.weight_sum) (tex: 3.01)"
$r = $abl | Where-Object { $_.dataset -eq "A1" -and $_.variant -eq "ab_nb_hier" }
$checks += "A1 nb_hier typeI = $($r.typeI_fdr) (tex: 0.924)"

# Table 7: A11 min_max typeI 0.541
$r = $senn | Where-Object { $_.dataset -eq "A11" -and $_.config -eq "norm_min_max" }
$checks += "A11 min_max typeI = $($r.typeI_fdr) (tex: 0.541)"
$r = $senn | Where-Object { $_.dataset -eq "A1" -and $_.config -eq "norm_min_max" }
$checks += "A1 min_max power = $($r.power_fdr) (tex: 0.450)"

# Table 9: A1 spline 2.0 power 0.250
$r = $senb | Where-Object { $_.dataset -eq "A1" -and $_.config -eq "spline_l2.0" }
$checks += "A1 spl2.0 power = $($r.power_fdr) (tex: 0.250)"

$checks | ForEach-Object { Write-Output $_ }
