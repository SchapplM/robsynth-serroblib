% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:41
% EndTime: 2019-02-26 21:18:41
% DurationCPUTime: 0.37s
% Computational Cost: add. (528->84), mult. (572->115), div. (0->0), fcn. (439->10), ass. (0->72)
t291 = qJD(5) + qJD(6);
t293 = qJ(5) + qJ(6);
t289 = cos(t293);
t349 = r_i_i_C(2) * t289;
t287 = sin(t293);
t352 = r_i_i_C(1) * t287;
t362 = (t349 + t352) * t291;
t294 = qJ(3) + qJ(4);
t290 = cos(t294);
t288 = sin(t294);
t292 = qJD(3) + qJD(4);
t300 = cos(qJ(1));
t341 = t292 * t300;
t325 = t288 * t341;
t297 = sin(qJ(1));
t336 = qJD(1) * t297;
t361 = t290 * t336 + t325;
t301 = -pkin(10) - pkin(9);
t295 = sin(qJ(5));
t353 = pkin(5) * t295;
t329 = qJD(5) * t353;
t360 = t292 * t301 + t329;
t298 = cos(qJ(5));
t286 = t298 * pkin(5) + pkin(4);
t296 = sin(qJ(3));
t334 = qJD(5) * t298;
t347 = r_i_i_C(3) - t301;
t358 = -pkin(5) * t334 + (-pkin(3) * t296 - t286 * t288 + t347 * t290 - qJ(2)) * qJD(1);
t357 = t298 * (qJD(5) * t288 + qJD(1));
t337 = qJD(1) * t288;
t317 = t291 + t337;
t322 = t290 * t341;
t355 = t317 * t297 - t322;
t343 = t290 * t292;
t323 = t297 * t343;
t354 = t317 * t300 + t323;
t351 = r_i_i_C(1) * t289;
t350 = r_i_i_C(2) * t287;
t348 = r_i_i_C(3) * t288;
t346 = pkin(3) * qJD(3);
t345 = t288 * t291;
t344 = t288 * t292;
t318 = -qJD(1) - t345;
t311 = t318 * t300;
t260 = t355 * t287 + t289 * t311;
t261 = t287 * t311 - t355 * t289;
t339 = -t260 * r_i_i_C(1) + t261 * r_i_i_C(2);
t312 = t318 * t297;
t262 = -t354 * t287 + t289 * t312;
t263 = t287 * t312 + t354 * t289;
t338 = t262 * r_i_i_C(1) - t263 * r_i_i_C(2);
t335 = qJD(1) * t300;
t333 = t290 * t351;
t332 = t290 * t350;
t331 = t288 * t349;
t330 = t296 * t346;
t299 = cos(qJ(3));
t327 = t299 * t346;
t326 = t287 * t344;
t321 = t288 * t336;
t319 = -t286 - t351;
t316 = qJD(1) * t333;
t314 = -qJD(5) - t337;
t313 = t297 * r_i_i_C(2) * t326 + r_i_i_C(3) * t323 + t300 * t316 + (t286 * t290 + t348) * t335;
t310 = -t288 * t301 - t332;
t308 = r_i_i_C(1) * t326 + t292 * t331 + (t332 - t333) * t291;
t307 = r_i_i_C(3) * t321 + t361 * t286 + t297 * t316 + t301 * t322 + t325 * t351 + (t329 + t362) * t290 * t300;
t306 = qJD(1) * (pkin(3) * t299 + t310);
t305 = t291 * t331 + t345 * t352 + (t319 * t290 + t332 - t348) * t292 + t360 * t288;
t304 = t319 * t344 + (-t360 - t362) * t290;
t303 = t327 + t286 * t343 + qJD(2) + (t347 * t292 - t329) * t288 + (-pkin(1) - pkin(8) - pkin(7) - t353) * qJD(1);
t1 = [t261 * r_i_i_C(1) + t260 * r_i_i_C(2) + t358 * t297 + t303 * t300, t335, t300 * t306 + (t304 - t330) * t297 + t313, t304 * t297 + t310 * t335 + t313 (-t297 * t357 + (t314 * t300 - t323) * t295) * pkin(5) + t338, t338; t263 * r_i_i_C(1) + t262 * r_i_i_C(2) + t303 * t297 - t358 * t300, t336 (t330 + (-r_i_i_C(3) * t290 - t288 * t350) * t292) * t300 + t297 * t306 + t307, -r_i_i_C(3) * t322 - t301 * t321 - t361 * t350 + t307 (t300 * t357 + (t314 * t297 + t322) * t295) * pkin(5) + t339, t339; 0, 0, t305 - t327, t305 (-t290 * t334 + t295 * t344) * pkin(5) + t308, t308;];
JaD_transl  = t1;
