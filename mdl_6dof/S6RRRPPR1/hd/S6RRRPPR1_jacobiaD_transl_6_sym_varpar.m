% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:19
% EndTime: 2019-02-26 22:03:19
% DurationCPUTime: 0.34s
% Computational Cost: add. (557->75), mult. (469->99), div. (0->0), fcn. (360->12), ass. (0->65)
t279 = pkin(11) + qJ(6);
t273 = sin(t279);
t274 = cos(t279);
t281 = qJ(2) + qJ(3);
t275 = pkin(10) + t281;
t270 = sin(t275);
t318 = qJD(6) * t270;
t271 = cos(t275);
t280 = qJD(2) + qJD(3);
t325 = t271 * t280;
t342 = t273 * t325 + t274 * t318;
t272 = cos(pkin(11)) * pkin(5) + pkin(4);
t341 = r_i_i_C(1) * t274 + t272;
t269 = t270 * qJD(5);
t276 = sin(t281);
t284 = sin(qJ(2));
t327 = pkin(2) * qJD(2);
t315 = t284 * t327;
t283 = -pkin(9) - qJ(5);
t328 = r_i_i_C(3) - t283;
t332 = pkin(3) * t280;
t340 = (-t270 * t272 + t328 * t271) * t280 + (pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + qJ(4)) * qJD(1) - t276 * t332 + t269 - t315;
t309 = t273 * t318;
t339 = r_i_i_C(1) * t309 + t342 * r_i_i_C(2) + qJD(5) * t271;
t287 = cos(qJ(1));
t317 = qJD(6) * t271;
t301 = -qJD(1) + t317;
t338 = t287 * t301;
t300 = qJD(1) * t271 - qJD(6);
t285 = sin(qJ(1));
t323 = t280 * t285;
t314 = t270 * t323;
t336 = t300 * t287 - t314;
t334 = pkin(3) * t276;
t277 = cos(t281);
t333 = pkin(3) * t277;
t330 = r_i_i_C(2) * t273;
t329 = r_i_i_C(3) * t271;
t324 = t271 * t283;
t322 = t280 * t287;
t321 = qJD(1) * t285;
t320 = qJD(1) * t287;
t316 = t270 * t330;
t313 = t270 * t322;
t311 = t270 * t320;
t310 = t270 * t321;
t299 = t301 * t285;
t298 = -t316 - t329;
t297 = t283 * t314 + t339 * t285 + t311 * t330 + t320 * t329;
t296 = -r_i_i_C(3) * t270 - t271 * t341;
t295 = -t270 * t341 - t324;
t294 = t283 * t313 + t339 * t287 + t341 * t310 + t321 * t324;
t293 = t300 * t285 + t313;
t292 = t295 - t334;
t286 = cos(qJ(2));
t291 = qJD(4) + (-pkin(2) * t286 - t328 * t270 - t271 * t272 - pkin(1) - t333) * qJD(1);
t290 = -t277 * t332 + t296 * t280 - t286 * t327;
t289 = t280 * (t296 - t333);
t288 = r_i_i_C(3) * t325 + t269 + (-r_i_i_C(1) * t273 - r_i_i_C(2) * t274) * t317 + (t316 + t292) * t280;
t266 = -pkin(2) * t284 - t334;
t247 = t273 * t299 - t336 * t274;
t246 = t336 * t273 + t274 * t299;
t245 = t273 * t338 + t293 * t274;
t244 = t293 * t273 - t274 * t338;
t1 = [t247 * r_i_i_C(1) + t246 * r_i_i_C(2) - t340 * t285 + t291 * t287 (-t266 + t298) * t321 + t290 * t287 + t294 (t298 + t334) * t321 + t287 * t289 + t294, t320, t271 * t322 - t310, r_i_i_C(1) * t244 + r_i_i_C(2) * t245; -t245 * r_i_i_C(1) + t244 * r_i_i_C(2) + t291 * t285 + t340 * t287, t290 * t285 + (t266 + t295) * t320 + t297, t285 * t289 + t292 * t320 + t297, t321, t271 * t323 + t311, -r_i_i_C(1) * t246 + r_i_i_C(2) * t247; 0, t288 - t315, t288, 0, t280 * t270 (-t274 * t325 + t309) * r_i_i_C(2) - t342 * r_i_i_C(1);];
JaD_transl  = t1;
