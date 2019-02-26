% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:47:04
% EndTime: 2019-02-26 22:47:04
% DurationCPUTime: 0.39s
% Computational Cost: add. (1017->84), mult. (662->113), div. (0->0), fcn. (480->12), ass. (0->73)
t284 = qJ(2) + qJ(3);
t281 = qJ(4) + t284;
t277 = qJ(5) + t281;
t272 = sin(t277);
t288 = cos(qJ(6));
t347 = r_i_i_C(1) * t288 + pkin(5);
t306 = t347 * t272;
t285 = sin(qJ(6));
t328 = qJD(6) * t288;
t273 = cos(t277);
t283 = qJD(2) + qJD(3);
t278 = qJD(4) + t283;
t274 = qJD(5) + t278;
t336 = t273 * t274;
t349 = t272 * t328 + t285 * t336;
t344 = pkin(11) + r_i_i_C(3);
t320 = t344 * t273;
t286 = sin(qJ(2));
t337 = pkin(2) * qJD(2);
t322 = t286 * t337;
t275 = sin(t281);
t341 = pkin(4) * t278;
t326 = t275 * t341;
t279 = sin(t284);
t343 = pkin(3) * t283;
t327 = t279 * t343;
t340 = pkin(5) * t272;
t348 = (t320 - t340) * t274 - t322 - t326 - t327;
t276 = cos(t281);
t305 = t347 * t273;
t321 = t344 * t272;
t294 = (-t305 - t321) * t274 - t276 * t341;
t329 = qJD(6) * t285;
t313 = t272 * t329;
t345 = r_i_i_C(1) * t313 + t349 * r_i_i_C(2);
t280 = cos(t284);
t295 = -t280 * t343 + t294;
t342 = pkin(4) * t275;
t338 = r_i_i_C(2) * t285;
t287 = sin(qJ(1));
t335 = t274 * t287;
t334 = t274 * t288;
t290 = cos(qJ(1));
t333 = t274 * t290;
t332 = t288 * t290;
t331 = qJD(1) * t287;
t330 = qJD(1) * t290;
t324 = t272 * t338;
t323 = qJD(1) * t338;
t319 = t344 * t287;
t318 = t272 * t334;
t308 = qJD(6) * t273 - qJD(1);
t307 = qJD(1) * t273 - qJD(6);
t304 = t347 * t290;
t268 = -pkin(3) * t279 - t342;
t303 = t345 * t290 + t331 * t306;
t302 = t308 * t285;
t301 = t290 * t272 * t323 + t345 * t287 + t330 * t320;
t300 = -t320 - t324;
t289 = cos(qJ(2));
t299 = -t289 * pkin(2) - pkin(3) * t280 - pkin(4) * t276 - pkin(5) * t273 - pkin(1) - t321;
t298 = t272 * t333 + t307 * t287;
t296 = -t289 * t337 + t295;
t293 = -t273 * r_i_i_C(2) * t328 + (-t273 * t329 - t318) * r_i_i_C(1) + t344 * t336 + (-t340 + t324) * t274;
t292 = t293 - t326;
t291 = t292 - t327;
t282 = -pkin(10) - pkin(9) - pkin(8) - pkin(7);
t258 = -t286 * pkin(2) + t268;
t251 = -t307 * t332 + (t302 + t318) * t287;
t250 = t308 * t288 * t287 + (-t272 * t335 + t307 * t290) * t285;
t249 = t298 * t288 + t290 * t302;
t248 = t298 * t285 - t308 * t332;
t1 = [t251 * r_i_i_C(1) + t250 * r_i_i_C(2) - t348 * t287 + (t282 * t287 + t299 * t290) * qJD(1) (-t258 + t300) * t331 + t296 * t290 + t303 (-t268 + t300) * t331 + t295 * t290 + t303 (t300 + t342) * t331 + t294 * t290 + t303 (-t287 * t323 - t344 * t333) * t272 + (-qJD(1) * t319 - t274 * t304) * t273 + t303, t248 * r_i_i_C(1) + t249 * r_i_i_C(2); -t249 * r_i_i_C(1) + t248 * r_i_i_C(2) + t348 * t290 + (-t282 * t290 + t299 * t287) * qJD(1) (t258 - t306) * t330 + t296 * t287 + t301 (t268 - t306) * t330 + t295 * t287 + t301 (-t306 - t342) * t330 + t294 * t287 + t301, -t305 * t335 + (-qJD(1) * t304 - t274 * t319) * t272 + t301, -t250 * r_i_i_C(1) + t251 * r_i_i_C(2); 0, t291 - t322, t291, t292, t293 (-t273 * t334 + t313) * r_i_i_C(2) - t349 * r_i_i_C(1);];
JaD_transl  = t1;
