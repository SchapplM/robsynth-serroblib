% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP12
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP12_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP12_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:14
% EndTime: 2019-02-26 21:52:15
% DurationCPUTime: 0.34s
% Computational Cost: add. (541->70), mult. (808->100), div. (0->0), fcn. (672->8), ass. (0->55)
t333 = pkin(5) + r_i_i_C(1);
t284 = qJ(4) + qJ(5);
t281 = sin(t284);
t330 = r_i_i_C(3) + qJ(6);
t340 = t330 * t281;
t286 = sin(qJ(2));
t285 = sin(qJ(4));
t308 = pkin(4) * t285 + qJ(3);
t289 = cos(qJ(2));
t320 = qJD(2) * t289;
t288 = cos(qJ(4));
t332 = pkin(4) * t288;
t282 = cos(t284);
t318 = pkin(2) + r_i_i_C(2) + pkin(9) + pkin(8);
t329 = pkin(4) * qJD(4);
t335 = t318 * qJD(2) + qJD(6) * t282 - t288 * t329 - qJD(3);
t339 = t335 * t286 - (pkin(7) + pkin(3) + t332) * qJD(1) - t308 * t320;
t321 = qJD(2) * t286;
t283 = qJD(4) + qJD(5);
t328 = t281 * t283;
t338 = t282 * t321 + t289 * t328;
t309 = t330 * t282;
t297 = -t333 * t281 - t308 + t309;
t287 = sin(qJ(1));
t311 = t287 * t320;
t290 = cos(qJ(1));
t322 = qJD(1) * t290;
t337 = t286 * t322 + t311;
t327 = t281 * t290;
t326 = t282 * t290;
t325 = t286 * t287;
t324 = qJD(1) * t286;
t323 = qJD(1) * t287;
t319 = t281 * qJD(6);
t317 = t285 * t329;
t315 = t286 * t323;
t310 = t290 * t320;
t307 = t283 * t286 + qJD(1);
t306 = qJD(4) + t324;
t305 = t307 * t287;
t304 = (-qJD(4) * t286 - qJD(1)) * t285;
t303 = t321 * t340 + t333 * t338;
t302 = t282 * t287 + t286 * t327;
t300 = -t283 * t309 - t319;
t262 = t281 * t305 - t337 * t282 - t283 * t326;
t263 = t282 * t305 + (t311 + (t283 + t324) * t290) * t281;
t299 = -(-t281 * t325 + t326) * qJD(6) + t330 * t263 - t333 * t262;
t264 = -t282 * t310 + t302 * t283 + (t282 * t325 + t327) * qJD(1);
t265 = -t281 * t315 - t287 * t328 + (t281 * t320 + t307 * t282) * t290;
t298 = t302 * qJD(6) - t333 * t264 + t330 * t265;
t295 = t297 * t289;
t294 = qJD(2) * t297;
t293 = -t317 + t319 + (-t308 * t286 - t318 * t289 - pkin(1)) * qJD(1);
t292 = (t333 * t282 + t340) * t283 - t335;
t1 = [-t330 * t262 - t333 * t263 + t339 * t287 + t293 * t290 (t290 * t294 + t318 * t323) * t286 + (t292 * t290 + t297 * t323) * t289, t310 - t315 (t290 * t304 + (-t306 * t287 + t310) * t288) * pkin(4) + t298, t298, t264; t330 * t264 + t333 * t265 + t293 * t287 - t339 * t290 (-t318 * t286 - t295) * t322 + (t286 * t294 + t292 * t289) * t287, t337 (t306 * t290 * t288 + (t288 * t320 + t304) * t287) * pkin(4) + t299, t299, t262; 0, -qJD(2) * t295 + t292 * t286, t321, t321 * t332 + (t300 + t317) * t289 + t303, t300 * t289 + t303, -t338;];
JaD_transl  = t1;
