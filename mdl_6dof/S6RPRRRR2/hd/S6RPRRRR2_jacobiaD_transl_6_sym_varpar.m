% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR2
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

function JaD_transl = S6RPRRRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:15:26
% EndTime: 2019-02-26 21:15:26
% DurationCPUTime: 0.35s
% Computational Cost: add. (676->79), mult. (564->111), div. (0->0), fcn. (433->12), ass. (0->68)
t299 = qJ(3) + qJ(4);
t292 = sin(t299);
t302 = cos(qJ(5));
t287 = t302 * pkin(5) + pkin(4);
t298 = qJ(5) + qJ(6);
t293 = cos(t298);
t362 = r_i_i_C(1) * t293 + t287;
t318 = t362 * t292;
t291 = sin(t298);
t294 = cos(t299);
t296 = qJD(3) + qJD(4);
t345 = t294 * t296;
t295 = qJD(5) + qJD(6);
t346 = t293 * t295;
t364 = t291 * t345 + t292 * t346;
t304 = -pkin(10) - pkin(9);
t300 = sin(qJ(5));
t349 = pkin(5) * qJD(5);
t337 = t300 * t349;
t363 = t296 * t304 + t337;
t301 = sin(qJ(3));
t350 = pkin(3) * qJD(3);
t335 = t301 * t350;
t347 = t292 * t296;
t351 = r_i_i_C(3) - t304;
t361 = -t287 * t347 + (t351 * t296 - t337) * t294 - t335;
t348 = t291 * t292;
t334 = t295 * t348;
t360 = r_i_i_C(1) * t334 + t364 * r_i_i_C(2) + t363 * t292;
t339 = qJD(1) * t294;
t322 = -t295 + t339;
t359 = t293 * t322;
t358 = t300 * (-qJD(5) + t339);
t323 = t294 * t295 - qJD(1);
t333 = t291 * t347;
t357 = t323 * t293 - t333;
t355 = pkin(3) * t301;
t354 = pkin(5) * t300;
t352 = r_i_i_C(2) * t293;
t297 = qJ(1) + pkin(11);
t289 = sin(t297);
t290 = cos(t297);
t317 = t322 * t291;
t265 = t289 * t317 - t357 * t290;
t309 = t323 * t291 + t293 * t347;
t266 = t289 * t359 + t309 * t290;
t343 = t265 * r_i_i_C(1) + t266 * r_i_i_C(2);
t267 = t357 * t289 + t290 * t317;
t268 = t309 * t289 - t290 * t359;
t342 = -t267 * r_i_i_C(1) + t268 * r_i_i_C(2);
t341 = qJD(1) * t289;
t340 = qJD(1) * t290;
t338 = r_i_i_C(2) * t348;
t336 = t302 * t349;
t328 = pkin(8) + pkin(7) + t354;
t319 = -r_i_i_C(1) * t291 - t352;
t316 = -r_i_i_C(3) * t294 - t338;
t315 = t290 * r_i_i_C(3) * t339 + t360 * t289 + t338 * t340;
t303 = cos(qJ(3));
t313 = -t303 * pkin(3) - t287 * t294 - t351 * t292 - pkin(2);
t312 = -t294 * t304 - t318;
t311 = t289 * t304 * t339 + t360 * t290 + t341 * t318;
t310 = (-r_i_i_C(3) * t292 - t294 * t362) * t296;
t308 = t300 * t347 + (-qJD(5) * t294 + qJD(1)) * t302;
t307 = -t303 * t350 + t310;
t306 = r_i_i_C(2) * t333 + r_i_i_C(3) * t345 - t296 * t318 + (t319 * t295 - t363) * t294;
t282 = r_i_i_C(2) * t334;
t1 = [t290 * t336 + t268 * r_i_i_C(1) + t267 * r_i_i_C(2) - t361 * t289 + (-cos(qJ(1)) * pkin(1) - t328 * t289 + t313 * t290) * qJD(1), 0 (t316 + t355) * t341 + t307 * t290 + t311, t290 * t310 + t316 * t341 + t311 (t289 * t358 + t308 * t290) * pkin(5) + t343, t343; t289 * t336 - t266 * r_i_i_C(1) + t265 * r_i_i_C(2) + t361 * t290 + (-sin(qJ(1)) * pkin(1) + t328 * t290 + t313 * t289) * qJD(1), 0, t307 * t289 + (t312 - t355) * t340 + t315, t289 * t310 + t312 * t340 + t315 (t308 * t289 - t290 * t358) * pkin(5) + t342, t342; 0, 0, t306 - t335, t306, t282 + (-r_i_i_C(1) * t346 - t336) * t292 + (t319 - t354) * t345, -t364 * r_i_i_C(1) - t345 * t352 + t282;];
JaD_transl  = t1;
