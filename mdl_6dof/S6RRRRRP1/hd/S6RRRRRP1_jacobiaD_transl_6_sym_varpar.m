% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:39:34
% EndTime: 2019-02-26 22:39:35
% DurationCPUTime: 0.42s
% Computational Cost: add. (736->84), mult. (666->113), div. (0->0), fcn. (501->10), ass. (0->66)
t341 = pkin(5) + r_i_i_C(1);
t275 = qJD(2) + qJD(3);
t271 = qJD(4) + t275;
t279 = sin(qJ(5));
t333 = r_i_i_C(2) * t279;
t348 = t271 * t333 + qJD(6);
t277 = qJ(2) + qJ(3);
t274 = qJ(4) + t277;
t268 = sin(t274);
t269 = cos(t274);
t282 = cos(qJ(5));
t321 = qJD(5) * t282;
t322 = qJD(5) * t279;
t347 = t348 * t269 + (r_i_i_C(2) * t321 + t322 * t341) * t268;
t281 = sin(qJ(1));
t284 = cos(qJ(1));
t303 = qJD(1) * t269 - qJD(5);
t297 = t303 * t284;
t304 = qJD(5) * t269 - qJD(1);
t299 = t304 * t282;
t326 = t271 * t281;
t314 = t268 * t326;
t242 = t281 * t299 + (t297 - t314) * t279;
t270 = pkin(5) * t282 + pkin(4);
t346 = r_i_i_C(1) * t282 + t270;
t280 = sin(qJ(2));
t330 = pkin(2) * qJD(2);
t316 = t280 * t330;
t272 = sin(t277);
t335 = pkin(3) * t275;
t319 = t272 * t335;
t278 = -qJ(6) - pkin(10);
t331 = r_i_i_C(3) - t278;
t342 = (-pkin(5) * t322 + t331 * t271) * t269 + (pkin(5) * t279 + pkin(7) + pkin(8) + pkin(9)) * qJD(1) - (t270 * t271 - qJD(6)) * t268 - t316 - t319;
t273 = cos(t277);
t291 = (-r_i_i_C(3) * t268 - t269 * t346) * t271;
t288 = -t273 * t335 + t291;
t336 = pkin(3) * t272;
t332 = r_i_i_C(3) * t269;
t329 = t269 * t271;
t328 = t269 * t278;
t327 = t271 * t268;
t325 = t271 * t284;
t324 = qJD(1) * t281;
t323 = qJD(1) * t284;
t313 = t268 * t325;
t312 = t268 * t323;
t311 = t268 * t324;
t300 = t346 * t271;
t298 = t304 * t279;
t296 = -t268 * t333 - t332;
t295 = -r_i_i_C(2) * t282 - t279 * t341;
t294 = -t268 * t346 - t328;
t293 = t278 * t314 + t347 * t281 + t312 * t333 + t323 * t332;
t292 = t278 * t313 + t347 * t284 + t346 * t311 + t324 * t328;
t290 = t303 * t281 + t313;
t283 = cos(qJ(2));
t289 = -t283 * t330 + t288;
t287 = pkin(5) * t321 + (-pkin(2) * t283 - pkin(3) * t273 - t331 * t268 - t269 * t270 - pkin(1)) * qJD(1);
t240 = t290 * t279 - t284 * t299;
t286 = r_i_i_C(3) * t329 + (t295 * qJD(5) - t271 * t278) * t269 + (-t300 + t348) * t268;
t285 = t286 - t319;
t264 = -pkin(2) * t280 - t336;
t243 = -t282 * t297 + (t282 * t327 + t298) * t281;
t241 = t290 * t282 + t284 * t298;
t1 = [t243 * r_i_i_C(1) + t242 * r_i_i_C(2) - t342 * t281 + t287 * t284 (-t264 + t296) * t324 + t289 * t284 + t292 (t296 + t336) * t324 + t288 * t284 + t292 (-r_i_i_C(3) * t325 - t324 * t333) * t268 + (-r_i_i_C(3) * t324 - t284 * t300) * t269 + t292, t241 * r_i_i_C(2) + t341 * t240, t269 * t325 - t311; -t241 * r_i_i_C(1) + t240 * r_i_i_C(2) + t287 * t281 + t342 * t284, t289 * t281 + (t264 + t294) * t323 + t293, t288 * t281 + (t294 - t336) * t323 + t293, t281 * t291 + t294 * t323 + t293, t243 * r_i_i_C(2) - t341 * t242, t269 * t326 + t312; 0, t285 - t316, t285, t286, t295 * t329 + (-t282 * t341 + t333) * t268 * qJD(5), t327;];
JaD_transl  = t1;
