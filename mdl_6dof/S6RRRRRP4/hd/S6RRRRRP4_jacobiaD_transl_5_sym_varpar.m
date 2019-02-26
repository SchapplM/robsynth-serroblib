% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP4
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
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:21
% EndTime: 2019-02-26 22:41:21
% DurationCPUTime: 0.42s
% Computational Cost: add. (524->77), mult. (560->110), div. (0->0), fcn. (431->10), ass. (0->69)
t291 = qJ(4) + qJ(5);
t285 = sin(t291);
t292 = qJ(2) + qJ(3);
t286 = sin(t292);
t288 = cos(t292);
t290 = qJD(2) + qJD(3);
t343 = t288 * t290;
t287 = cos(t291);
t289 = qJD(4) + qJD(5);
t344 = t287 * t289;
t364 = t285 * t343 + t286 * t344;
t299 = -pkin(10) - pkin(9);
t293 = sin(qJ(4));
t347 = pkin(4) * qJD(4);
t333 = t293 * t347;
t363 = t290 * t299 + t333;
t296 = cos(qJ(4));
t283 = t296 * pkin(4) + pkin(3);
t361 = r_i_i_C(1) * t287 + t283;
t362 = t286 * t361 + t288 * t299;
t294 = sin(qJ(2));
t348 = pkin(2) * qJD(2);
t330 = t294 * t348;
t345 = t286 * t290;
t349 = r_i_i_C(3) - t299;
t353 = pkin(4) * t293;
t360 = (t349 * t290 - t333) * t288 + (pkin(8) + pkin(7) + t353) * qJD(1) - t283 * t345 - t330;
t346 = t285 * t286;
t329 = t289 * t346;
t359 = r_i_i_C(1) * t329 + t364 * r_i_i_C(2) + t363 * t286;
t298 = cos(qJ(1));
t317 = t288 * t289 - qJD(1);
t358 = t298 * t317;
t337 = qJD(1) * t288;
t316 = -t289 + t337;
t295 = sin(qJ(1));
t326 = t295 * t345;
t356 = t316 * t298 - t326;
t354 = pkin(2) * t294;
t351 = r_i_i_C(2) * t287;
t350 = r_i_i_C(3) * t288;
t341 = t290 * t298;
t325 = t286 * t341;
t304 = t316 * t295 + t325;
t261 = t304 * t285 - t287 * t358;
t262 = t285 * t358 + t304 * t287;
t339 = t261 * r_i_i_C(1) + t262 * r_i_i_C(2);
t311 = t317 * t295;
t263 = t356 * t285 + t287 * t311;
t264 = t285 * t311 - t356 * t287;
t338 = -t263 * r_i_i_C(1) + t264 * r_i_i_C(2);
t336 = qJD(1) * t295;
t335 = qJD(1) * t298;
t334 = r_i_i_C(2) * t346;
t332 = t296 * t347;
t331 = r_i_i_C(2) * qJD(1) * t285;
t314 = -qJD(4) + t337;
t313 = -r_i_i_C(1) * t285 - t351;
t312 = t361 * t290;
t310 = (-qJD(4) * t288 + qJD(1)) * t296;
t309 = t298 * t286 * t331 + t359 * t295 + t335 * t350;
t306 = t359 * t298 + t362 * t336;
t305 = (-r_i_i_C(3) * t286 - t288 * t361) * t290;
t297 = cos(qJ(2));
t303 = t332 + (-t297 * pkin(2) - t283 * t288 - t349 * t286 - pkin(1)) * qJD(1);
t302 = -t297 * t348 + t305;
t301 = t290 * t334 + r_i_i_C(3) * t343 - t286 * t312 + (t313 * t289 - t363) * t288;
t274 = r_i_i_C(2) * t329;
t1 = [t264 * r_i_i_C(1) + t263 * r_i_i_C(2) - t360 * t295 + t303 * t298 (-t334 - t350 + t354) * t336 + t302 * t298 + t306 (-r_i_i_C(3) * t341 - t295 * t331) * t286 + (-r_i_i_C(3) * t336 - t298 * t312) * t288 + t306 (t298 * t310 + (t295 * t314 + t325) * t293) * pkin(4) + t339, t339, 0; -t262 * r_i_i_C(1) + t261 * r_i_i_C(2) + t303 * t295 + t360 * t298, t302 * t295 + (-t362 - t354) * t335 + t309, t295 * t305 - t335 * t362 + t309 (t295 * t310 + (-t298 * t314 + t326) * t293) * pkin(4) + t338, t338, 0; 0, t301 - t330, t301, t274 + (-r_i_i_C(1) * t344 - t332) * t286 + (t313 - t353) * t343, -t364 * r_i_i_C(1) - t343 * t351 + t274, 0;];
JaD_transl  = t1;
