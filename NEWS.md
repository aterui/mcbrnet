Changelog
================

# mcbrnet 1.4.2 (10/28/24)

- `ggbrnet()` has been updated to increase flexibility (some
  incompatibility with earlier versions). See **Articles** for details.
- disabled functions `findr()`, `ppm()`, `extra_prey()`, `to_alpha()`,
  `foodweb()`, `stability()`; these functions were moved to `ecotools`
  package at <https://github.com/aterui/ecotools>

# mcbrnet 1.4.1 (9/10/24)

- disabled functions `foodchain()`, `max_tp()`
- disabled `disturb` argument in `sglv()`

# mcbrnet 1.4.0 (10/24/23)

- add new major functions `sglv()`, `foodweb()`, `foodchain()`,
  `findr()`; full documentation to follow

- full description for `igpsim()` is now available; see **Articles**

- moved utility functions to `utils.R`

# mcbrnet 1.3.1 (07/10/23)

- add a new function `frgm()`

- update `igpsim()`: full description coming soon

# mcbrnet 1.3.0 (06/15/22)

- add a new function `ptsource()`

- add new arguments to `mcsim()`

- update `ggbrnet()` to be compatible with piping

# mcbrnet 1.2.3 (04/13/22)

- fix a bug in `fun_disp_mat()`

- add `ggbrnet()`

- add disturbance arguments to `mcsim()` (`p_disturb` & `m_disturb`)

# mcbrnet 1.2.2 (03/24/22)

- fix a bug in `fun_igp()`

# mcbrnet 1.2.1 (03/09/22)

- implement internal functions to `mcsim()` and `brnet()`

- remove argument weighted_distance_matrix from `mcsim()` and `brnet()`

- add argument dispersal_matrix to `mcsim()`

# mcbrnet 1.2.0 (03/09/22)

- add a major function `igpsim()`

- simplified `brnet()` and `mcsim()` by introducing internal
  sub-functions

# mcbrnet 1.1.1 (12/07/21)

- add a local noise parameter for disturbance values to `brnet()`
  (argument `sd_disturb_lon`)

# mcbrnet 1.1.0 (08/02/21)

- add disturbance arguments to `brnet()` added function `adjtodist()`

# mcbrnet 1.0.0 (05/03/21)

- initial release
