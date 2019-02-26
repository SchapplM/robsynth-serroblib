% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR5
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
% Datum: 2019-02-26 21:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:17:11
% EndTime: 2019-02-26 21:17:11
% DurationCPUTime: 0.41s
% Computational Cost: add. (673->81), mult. (566->111), div. (0->0), fcn. (437->11), ass. (0->70)
t295 = pkin(11) + qJ(3);
t291 = qJ(4) + t295;
t286 = sin(t291);
t298 = qJ(5) + qJ(6);
t292 = sin(t298);
t293 = cos(t298);
t296 = qJD(5) + qJD(6);
t345 = t293 * t296;
t287 = cos(t291);
t297 = qJD(3) + qJD(4);
t347 = t287 * t297;
t367 = t286 * t345 + t292 * t347;
t303 = -pkin(10) - pkin(9);
t299 = sin(qJ(5));
t350 = pkin(5) * qJD(5);
t335 = t299 * t350;
t366 = t297 * t303 + t335;
t301 = cos(qJ(5));
t288 = t301 * pkin(5) + pkin(4);
t364 = r_i_i_C(1) * t293 + t288;
t365 = t286 * t364 + t287 * t303;
t289 = sin(t295);
t351 = pkin(3) * qJD(3);
t336 = t289 * t351;
t348 = t286 * t297;
t352 = r_i_i_C(3) - t303;
t356 = pkin(5) * t299;
t363 = (t352 * t297 - t335) * t287 + (pkin(8) + pkin(7) + qJ(2) + t356) * qJD(1) - t288 * t348 - t336;
t349 = t286 * t292;
t332 = t296 * t349;
t362 = r_i_i_C(1) * t332 + t367 * r_i_i_C(2) + t366 * t286;
t302 = cos(qJ(1));
t320 = t287 * t296 - qJD(1);
t361 = t302 * t320;
t340 = qJD(1) * t287;
t319 = -t296 + t340;
t300 = sin(qJ(1));
t330 = t300 * t348;
t359 = t319 * t302 - t330;
t357 = pkin(3) * t289;
t354 = r_i_i_C(2) * t293;
t353 = r_i_i_C(3) * t287;
t344 = t297 * t302;
t329 = t286 * t344;
t307 = t319 * t300 + t329;
t263 = t307 * t292 - t293 * t361;
t264 = t292 * t361 + t307 * t293;
t342 = t263 * r_i_i_C(1) + t264 * r_i_i_C(2);
t314 = t320 * t300;
t265 = t359 * t292 + t293 * t314;
t266 = t292 * t314 - t359 * t293;
t341 = -t265 * r_i_i_C(1) + t266 * r_i_i_C(2);
t339 = qJD(1) * t300;
t338 = qJD(1) * t302;
t337 = r_i_i_C(2) * t349;
t334 = t301 * t350;
t333 = r_i_i_C(2) * qJD(1) * t292;
t317 = -qJD(5) + t340;
t316 = -r_i_i_C(1) * t292 - t354;
t315 = t364 * t297;
t313 = (-qJD(5) * t287 + qJD(1)) * t301;
t312 = t302 * t286 * t333 + t362 * t300 + t338 * t353;
t309 = t362 * t302 + t365 * t339;
t308 = (-r_i_i_C(3) * t286 - t287 * t364) * t297;
t290 = cos(t295);
t306 = -t290 * t351 + t308;
t305 = t334 + qJD(2) + (-t352 * t286 - t287 * t288 - pkin(3) * t290 - cos(pkin(11)) * pkin(2) - pkin(1)) * qJD(1);
t304 = t297 * t337 + r_i_i_C(3) * t347 - t286 * t315 + (t316 * t296 - t366) * t287;
t276 = r_i_i_C(2) * t332;
t1 = [t266 * r_i_i_C(1) + t265 * r_i_i_C(2) - t363 * t300 + t305 * t302, t338 (-t337 - t353 + t357) * t339 + t306 * t302 + t309 (-r_i_i_C(3) * t344 - t300 * t333) * t286 + (-r_i_i_C(3) * t339 - t302 * t315) * t287 + t309 (t302 * t313 + (t317 * t300 + t329) * t299) * pkin(5) + t342, t342; -t264 * r_i_i_C(1) + t263 * r_i_i_C(2) + t305 * t300 + t363 * t302, t339, t306 * t300 + (-t365 - t357) * t338 + t312, t300 * t308 - t338 * t365 + t312 (t300 * t313 + (-t317 * t302 + t330) * t299) * pkin(5) + t341, t341; 0, 0, t304 - t336, t304, t276 + (-r_i_i_C(1) * t345 - t334) * t286 + (t316 - t356) * t347, -t367 * r_i_i_C(1) - t347 * t354 + t276;];
JaD_transl  = t1;
