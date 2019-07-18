% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:19:39
% EndTime: 2019-07-18 17:19:39
% DurationCPUTime: 0.39s
% Computational Cost: add. (504->72), mult. (556->111), div. (0->0), fcn. (429->10), ass. (0->70)
t292 = qJ(2) + qJ(3);
t286 = sin(t292);
t296 = cos(qJ(4));
t284 = t296 * pkin(3) + pkin(2);
t291 = qJ(4) + qJ(5);
t287 = cos(t291);
t355 = r_i_i_C(1) * t287 + t284;
t359 = t286 * t355;
t285 = sin(t291);
t288 = cos(t292);
t290 = qJD(2) + qJD(3);
t339 = t288 * t290;
t289 = qJD(4) + qJD(5);
t340 = t287 * t289;
t358 = t285 * t339 + t286 * t340;
t350 = pkin(5) + r_i_i_C(3);
t357 = t350 * t288;
t327 = t350 * t290;
t294 = sin(qJ(2));
t345 = pkin(1) * qJD(2);
t329 = t294 * t345;
t293 = sin(qJ(4));
t344 = pkin(3) * qJD(4);
t332 = t293 * t344;
t342 = t286 * t290;
t348 = pkin(3) * t293;
t356 = qJD(1) * t348 - t284 * t342 + (t327 - t332) * t288 - t329;
t343 = t285 * t286;
t325 = t289 * t343;
t353 = r_i_i_C(1) * t325 + t358 * r_i_i_C(2) + t286 * t332;
t298 = cos(qJ(1));
t314 = t288 * t289 - qJD(1);
t352 = t298 * t314;
t336 = qJD(1) * t288;
t313 = -t289 + t336;
t295 = sin(qJ(1));
t322 = t295 * t342;
t351 = t313 * t298 - t322;
t349 = pkin(1) * t294;
t346 = r_i_i_C(2) * t287;
t341 = t286 * t298;
t321 = t290 * t341;
t303 = t313 * t295 + t321;
t263 = t303 * t285 - t287 * t352;
t264 = t285 * t352 + t303 * t287;
t338 = t263 * r_i_i_C(1) + t264 * r_i_i_C(2);
t307 = t314 * t295;
t265 = t351 * t285 + t287 * t307;
t266 = t285 * t307 - t351 * t287;
t337 = -t265 * r_i_i_C(1) + t266 * r_i_i_C(2);
t335 = qJD(1) * t295;
t334 = qJD(1) * t298;
t333 = r_i_i_C(2) * t343;
t331 = t296 * t344;
t330 = r_i_i_C(2) * qJD(1) * t285;
t328 = t350 * t286;
t326 = t350 * t295;
t311 = -qJD(4) + t336;
t310 = -r_i_i_C(1) * t285 - t346;
t309 = t355 * t290;
t308 = t355 * t298;
t306 = (-qJD(4) * t288 + qJD(1)) * t296;
t305 = t353 * t298 + t335 * t359;
t304 = t353 * t295 + t330 * t341 + t334 * t357;
t297 = cos(qJ(2));
t301 = t331 + (-pkin(1) * t297 - t284 * t288 - t328) * qJD(1);
t300 = -t297 * t345 + (-t288 * t355 - t328) * t290;
t299 = t290 * t333 - t286 * t309 + (t310 * t289 - t332) * t288 + t350 * t339;
t276 = r_i_i_C(2) * t325;
t1 = [t266 * r_i_i_C(1) + t265 * r_i_i_C(2) - t356 * t295 + t301 * t298, (-t333 + t349 - t357) * t335 + t300 * t298 + t305, (-t295 * t330 - t298 * t327) * t286 + (-qJD(1) * t326 - t290 * t308) * t288 + t305, (t298 * t306 + (t311 * t295 + t321) * t293) * pkin(3) + t338, t338; -t264 * r_i_i_C(1) + t263 * r_i_i_C(2) + t301 * t295 + t356 * t298, (-t349 - t359) * t334 + t300 * t295 + t304, -t295 * t288 * t309 + (-qJD(1) * t308 - t290 * t326) * t286 + t304, (t295 * t306 + (-t311 * t298 + t322) * t293) * pkin(3) + t337, t337; 0, t299 - t329, t299, t276 + (-r_i_i_C(1) * t340 - t331) * t286 + (t310 - t348) * t339, -t358 * r_i_i_C(1) - t339 * t346 + t276;];
JaD_transl  = t1;
