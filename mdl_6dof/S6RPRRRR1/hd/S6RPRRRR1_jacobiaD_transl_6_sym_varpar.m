% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:52
% EndTime: 2019-02-26 21:14:52
% DurationCPUTime: 0.36s
% Computational Cost: add. (716->71), mult. (528->101), div. (0->0), fcn. (388->12), ass. (0->63)
t282 = qJ(3) + qJ(4);
t278 = qJ(5) + t282;
t271 = sin(t278);
t285 = cos(qJ(6));
t339 = r_i_i_C(1) * t285 + pkin(5);
t299 = t339 * t271;
t272 = cos(t278);
t320 = qJD(6) * t285;
t279 = qJD(3) + qJD(4);
t275 = qJD(5) + t279;
t283 = sin(qJ(6));
t325 = t275 * t283;
t343 = t271 * t320 + t272 * t325;
t333 = pkin(10) + r_i_i_C(3);
t341 = t333 * t272;
t342 = (-pkin(5) * t271 + t341) * t275;
t284 = sin(qJ(3));
t327 = pkin(3) * qJD(3);
t316 = t284 * t327;
t276 = sin(t282);
t331 = pkin(4) * t279;
t319 = t276 * t331;
t340 = t342 - t316 - t319;
t321 = qJD(6) * t283;
t306 = t271 * t321;
t336 = r_i_i_C(1) * t306 + t343 * r_i_i_C(2);
t300 = qJD(1) * t272 - qJD(6);
t335 = t285 * t300;
t277 = cos(t282);
t315 = t333 * t271;
t289 = (-t272 * t339 - t315) * t275 - t277 * t331;
t301 = qJD(6) * t272 - qJD(1);
t312 = t271 * t325;
t334 = t301 * t285 - t312;
t332 = pkin(4) * t276;
t328 = r_i_i_C(2) * t283;
t324 = t275 * t285;
t280 = qJ(1) + pkin(11);
t273 = sin(t280);
t323 = qJD(1) * t273;
t274 = cos(t280);
t322 = qJD(1) * t274;
t317 = qJD(1) * t328;
t313 = t333 * t275;
t311 = t271 * t324;
t298 = t339 * t275;
t297 = t336 * t274 + t323 * t299;
t296 = t300 * t283;
t295 = t274 * t271 * t317 + t336 * t273 + t322 * t341;
t294 = -t271 * t328 - t341;
t286 = cos(qJ(3));
t293 = -t286 * pkin(3) - pkin(4) * t277 - pkin(5) * t272 - pkin(2) - t315;
t292 = t301 * t283 + t311;
t290 = -t286 * t327 + t289;
t288 = (-t272 * t321 - t311) * r_i_i_C(1) + (-t272 * t320 + t312) * r_i_i_C(2) + t342;
t287 = t288 - t319;
t281 = -pkin(9) - pkin(8) - pkin(7);
t270 = -t284 * pkin(3) - t332;
t252 = t292 * t273 - t274 * t335;
t251 = t334 * t273 + t274 * t296;
t250 = t273 * t335 + t292 * t274;
t249 = t273 * t296 - t334 * t274;
t1 = [t252 * r_i_i_C(1) + t251 * r_i_i_C(2) - t340 * t273 + (-cos(qJ(1)) * pkin(1) + t273 * t281 + t293 * t274) * qJD(1), 0 (-t270 + t294) * t323 + t290 * t274 + t297 (t294 + t332) * t323 + t289 * t274 + t297 (-t273 * t317 - t274 * t313) * t271 + (-t274 * t298 - t333 * t323) * t272 + t297, t249 * r_i_i_C(1) + t250 * r_i_i_C(2); -t250 * r_i_i_C(1) + t249 * r_i_i_C(2) + t340 * t274 + (-sin(qJ(1)) * pkin(1) - t274 * t281 + t293 * t273) * qJD(1), 0 (t270 - t299) * t322 + t290 * t273 + t295 (-t299 - t332) * t322 + t289 * t273 + t295, -t273 * t272 * t298 + (-t273 * t313 - t322 * t339) * t271 + t295, -t251 * r_i_i_C(1) + t252 * r_i_i_C(2); 0, 0, t287 - t316, t287, t288 (-t272 * t324 + t306) * r_i_i_C(2) - t343 * r_i_i_C(1);];
JaD_transl  = t1;
