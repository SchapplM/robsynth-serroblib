% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:32:35
% EndTime: 2019-02-26 21:32:36
% DurationCPUTime: 0.51s
% Computational Cost: add. (453->88), mult. (984->141), div. (0->0), fcn. (906->10), ass. (0->58)
t286 = sin(qJ(2));
t283 = sin(pkin(10));
t284 = cos(pkin(10));
t285 = sin(qJ(5));
t282 = qJ(5) + qJ(6);
t279 = sin(t282);
t280 = cos(t282);
t288 = cos(qJ(5));
t328 = t288 * pkin(5) + pkin(3) + pkin(4);
t297 = t280 * r_i_i_C(1) - t279 * r_i_i_C(2) + t328;
t298 = -t279 * r_i_i_C(1) - t280 * r_i_i_C(2) - qJ(4);
t330 = (t285 * pkin(5) - t298) * t283 + t297 * t284 + pkin(2);
t289 = cos(qJ(2));
t312 = r_i_i_C(3) - qJ(3) + pkin(9) + pkin(8);
t333 = t312 * t289;
t337 = t330 * t286 + t333;
t313 = t286 * qJD(3);
t336 = (t286 * pkin(2) + t333) * qJD(2) - t313;
t281 = qJD(5) + qJD(6);
t299 = t283 * t288 - t284 * t285;
t300 = t279 * t283 + t280 * t284;
t301 = t279 * t284 - t280 * t283;
t334 = -t283 * qJD(4) - t299 * qJD(5) * pkin(5) + (t301 * r_i_i_C(1) + t300 * r_i_i_C(2)) * t281;
t326 = t281 * t286;
t287 = sin(qJ(1));
t325 = t287 * t289;
t290 = cos(qJ(1));
t324 = t290 * t283;
t323 = t290 * t284;
t271 = t284 * t325 - t324;
t315 = qJD(2) * t290;
t309 = t286 * t315;
t267 = -t271 * qJD(1) - t284 * t309;
t311 = t289 * t324;
t272 = -t287 * t284 + t311;
t307 = t272 * t281 + t267;
t270 = t283 * t325 + t323;
t266 = t270 * qJD(1) + t283 * t309;
t273 = t287 * t283 + t289 * t323;
t308 = -t273 * t281 - t266;
t262 = -t307 * t279 + t308 * t280;
t263 = t308 * t279 + t307 * t280;
t322 = t262 * r_i_i_C(1) - t263 * r_i_i_C(2);
t317 = qJD(2) * t286;
t310 = t287 * t317;
t269 = t273 * qJD(1) - t284 * t310;
t305 = -t270 * t281 - t269;
t319 = qJD(1) * t287;
t268 = qJD(1) * t311 - t283 * t310 - t284 * t319;
t306 = t271 * t281 - t268;
t321 = (t305 * t279 - t306 * t280) * r_i_i_C(1) + (t306 * t279 + t305 * t280) * r_i_i_C(2);
t316 = qJD(2) * t289;
t320 = (-t300 * t326 - t301 * t316) * r_i_i_C(1) + (-t300 * t316 + t301 * t326) * r_i_i_C(2);
t318 = qJD(1) * t290;
t303 = t312 * t286;
t296 = -pkin(2) * t289 - pkin(1) + t303;
t292 = t289 * qJD(3) + t334 * t286 + (-t289 * t330 + t303) * qJD(2);
t1 = [-t270 * qJD(4) + t298 * t268 + ((-t270 * t280 + t271 * t279) * r_i_i_C(1) + (t270 * t279 + t271 * t280) * r_i_i_C(2)) * t281 - t297 * t269 + (-t268 * t285 + (-t270 * t288 + t271 * t285) * qJD(5)) * pkin(5) + t336 * t287 + (-t287 * pkin(7) + t296 * t290) * qJD(1), t292 * t290 + t319 * t337, -t286 * t319 + t289 * t315, -t266 (-t266 * t288 - t267 * t285 + (-t272 * t285 - t273 * t288) * qJD(5)) * pkin(5) + t322, t322; t263 * r_i_i_C(1) + t262 * r_i_i_C(2) - t266 * qJ(4) + t272 * qJD(4) + t328 * t267 + (-t266 * t285 + (t272 * t288 - t273 * t285) * qJD(5)) * pkin(5) - t336 * t290 + (pkin(7) * t290 + t296 * t287) * qJD(1), t292 * t287 - t318 * t337, t286 * t318 + t287 * t316, t268 (t268 * t288 - t269 * t285 + (-t270 * t285 - t271 * t288) * qJD(5)) * pkin(5) + t321, t321; 0, -qJD(2) * t337 - t334 * t289 + t313, t317, t283 * t316 ((-t283 * t285 - t284 * t288) * t286 * qJD(5) + t299 * t316) * pkin(5) + t320, t320;];
JaD_transl  = t1;
