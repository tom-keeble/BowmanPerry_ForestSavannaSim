### Description
Recreation of the Bowman and Perry (2017) model of forest-savanna transitions in R.

### Results (example)

With the model condition of `Fire-Soil feedback` toggled on, while the `Edaphic boundary` condition is off, the model produces the following state transition dynamics:

![Fire-soil feedback | NO Edaphic boundary](www/Example_FSfb_NOEdaphBound.gif)

*Fire-soil feedback | NO Edaphic boundary. Forest is represented by green cells, savanna by brown, savanna colonised by forest propagule in yellow, and fire in red.*

<br>
further examples can be viewed in the below sections...

<br/>
<br/>
<br/>

<details><summary>Fire-soil feedback | Edaphic boundary</summary>

![Fire-soil feedback | Edaphic boundary](www/Example_FSfb_EdaphBound.gif)

</details>

<details><summary>NO Fire-soil feedback | Edaphic boundary</summary>

![Fire-soil feedback | Edaphic boundary](www/Example_NOFSfb_EdaphBound.gif)

</details>

<details><summary>NO Fire-soil feedback | NO Edaphic boundary</summary>

![Fire-soil feedback | Edaphic boundary](www/Example_NOFSfb_NOEdaphBound.gif)

</details>

### Usage
The `run_simulation` function within the script allows for "turning on/off" the *Edaphic boundary* and *Fire-Soil feedback* conditions, along with on-the-fly plotting of the state transition dynamics.

### Reference:
- Bowman, D.M.J.S., Perry, G.L.W. (2017). Soil or fire: what causes treeless sedgelands in Tasmanian wet forests?. *Plant Soil*, *420*, 1–18. https://doi.org/10.1007/s11104-017-3386-7
