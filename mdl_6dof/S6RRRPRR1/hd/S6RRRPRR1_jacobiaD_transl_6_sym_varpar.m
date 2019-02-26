% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:04
% EndTime: 2019-02-26 22:16:04
% DurationCPUTime: 0.32s
% Computational Cost: add. (763->76), mult. (550->107), div. (0->0), fcn. (404->12), ass. (0->63)
t281 = qJ(2) + qJ(3);
t276 = pkin(11) + t281;
t274 = qJ(5) + t276;
t270 = sin(t274);
t285 = cos(qJ(6));
t338 = r_i_i_C(1) * t285 + pkin(5);
t303 = t338 * t270;
t282 = sin(qJ(6));
t322 = qJD(6) * t285;
t271 = cos(t274);
t280 = qJD(2) + qJD(3);
t275 = qJD(5) + t280;
t330 = t271 * t275;
t340 = t270 * t322 + t282 * t330;
t335 = pkin(10) + r_i_i_C(3);
t317 = t335 * t271;
t264 = -pkin(3) * sin(t281) - pkin(4) * sin(t276);
t296 = t264 * t280;
t283 = sin(qJ(2));
t331 = pkin(2) * qJD(2);
t319 = t283 * t331;
t334 = pkin(5) * t270;
t339 = (t317 - t334) * t275 + t296 - t319;
t323 = qJD(6) * t282;
t310 = t270 * t323;
t336 = r_i_i_C(1) * t310 + t340 * r_i_i_C(2);
t300 = -pkin(3) * cos(t281) - pkin(4) * cos(t276);
t302 = t338 * t271;
t318 = t335 * t270;
t290 = t300 * t280 + (-t302 - t318) * t275;
t332 = r_i_i_C(2) * t282;
t284 = sin(qJ(1));
t329 = t275 * t284;
t328 = t275 * t285;
t287 = cos(qJ(1));
t327 = t275 * t287;
t326 = t285 * t287;
t325 = qJD(1) * t284;
t324 = qJD(1) * t287;
t321 = t270 * t332;
t320 = qJD(1) * t332;
t316 = t335 * t284;
t315 = t270 * t328;
t305 = qJD(6) * t271 - qJD(1);
t304 = qJD(1) * t271 - qJD(6);
t301 = t338 * t287;
t299 = t336 * t287 + t325 * t303;
t298 = t305 * t282;
t297 = t287 * t270 * t320 + t336 * t284 + t324 * t317;
t295 = -t317 - t321;
t286 = cos(qJ(2));
t294 = -t286 * pkin(2) - pkin(5) * t271 - pkin(1) + t300 - t318;
t293 = t270 * t327 + t304 * t284;
t291 = -t286 * t331 + t290;
t289 = -t271 * r_i_i_C(2) * t322 + (-t271 * t323 - t315) * r_i_i_C(1) + t335 * t330 + (-t334 + t321) * t275;
t288 = t296 + t289;
t279 = -pkin(9) - qJ(4) - pkin(8) - pkin(7);
t256 = -t283 * pkin(2) + t264;
t249 = -t304 * t326 + (t298 + t315) * t284;
t248 = t305 * t285 * t284 + (-t270 * t329 + t304 * t287) * t282;
t247 = t293 * t285 + t287 * t298;
t246 = t293 * t282 - t305 * t326;
t1 = [t249 * r_i_i_C(1) + t248 * r_i_i_C(2) + t287 * qJD(4) - t339 * t284 + (t279 * t284 + t294 * t287) * qJD(1) (-t256 + t295) * t325 + t291 * t287 + t299 (-t264 + t295) * t325 + t290 * t287 + t299, t324 (-t284 * t320 - t335 * t327) * t270 + (-qJD(1) * t316 - t275 * t301) * t271 + t299, t246 * r_i_i_C(1) + t247 * r_i_i_C(2); -t247 * r_i_i_C(1) + t246 * r_i_i_C(2) + t284 * qJD(4) + t339 * t287 + (-t279 * t287 + t294 * t284) * qJD(1) (t256 - t303) * t324 + t291 * t284 + t297 (t264 - t303) * t324 + t290 * t284 + t297, t325, -t302 * t329 + (-qJD(1) * t301 - t275 * t316) * t270 + t297, -t248 * r_i_i_C(1) + t249 * r_i_i_C(2); 0, t288 - t319, t288, 0, t289 (-t271 * t328 + t310) * r_i_i_C(2) - t340 * r_i_i_C(1);];
JaD_transl  = t1;
