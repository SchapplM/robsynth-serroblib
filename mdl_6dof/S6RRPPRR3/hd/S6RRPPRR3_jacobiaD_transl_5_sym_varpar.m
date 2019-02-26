% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:35
% EndTime: 2019-02-26 21:29:35
% DurationCPUTime: 0.32s
% Computational Cost: add. (392->88), mult. (992->144), div. (0->0), fcn. (1041->12), ass. (0->62)
t332 = cos(pkin(6));
t329 = sin(pkin(11));
t331 = cos(pkin(11));
t334 = sin(qJ(2));
t336 = cos(qJ(2));
t342 = t336 * t329 + t334 * t331;
t311 = t342 * t332;
t352 = qJD(2) * t336;
t353 = qJD(2) * t334;
t365 = t329 * t353 - t331 * t352;
t305 = t365 * t332;
t315 = t334 * t329 - t336 * t331;
t337 = cos(qJ(1));
t335 = sin(qJ(1));
t354 = qJD(1) * t335;
t313 = -t329 * t352 - t331 * t353;
t357 = t335 * t313;
t296 = t357 - t311 * t354 + (-qJD(1) * t315 - t305) * t337;
t330 = sin(pkin(6));
t360 = t330 * t337;
t345 = qJD(5) * t360;
t366 = t296 - t345;
t364 = pkin(2) * t332;
t363 = pkin(4) * sin(pkin(12));
t362 = r_i_i_C(3) + pkin(9) + qJ(4);
t361 = t330 * t335;
t359 = t334 * t335;
t358 = t334 * t337;
t356 = t335 * t336;
t355 = t336 * t337;
t298 = t337 * t311 - t335 * t315;
t351 = qJD(5) * t298;
t350 = pkin(2) * t353;
t349 = t330 * t354;
t348 = qJD(1) * t360;
t327 = pkin(12) + qJ(5);
t325 = sin(t327);
t326 = cos(t327);
t344 = t325 * r_i_i_C(1) + t326 * r_i_i_C(2);
t310 = t315 * t332;
t297 = -t337 * t310 - t335 * t342;
t343 = t335 * t311 + t337 * t315;
t323 = cos(pkin(12)) * pkin(4) + pkin(3);
t341 = t326 * r_i_i_C(1) - t325 * r_i_i_C(2) + t323;
t340 = qJD(5) * t344;
t339 = t349 - t351;
t309 = t342 * t330;
t338 = t315 * qJD(2);
t293 = -t298 * qJD(1) + t335 * t305 + t337 * t313;
t324 = t336 * pkin(2) + pkin(1);
t317 = t326 * t345;
t314 = -t330 * qJD(3) + t352 * t364;
t312 = t334 * t364 + (-pkin(8) - qJ(3)) * t330;
t306 = qJD(2) * t311;
t304 = qJD(2) * t309;
t303 = t365 * t330;
t299 = t335 * t310 - t337 * t342;
t295 = t310 * t354 + (-qJD(1) * t342 - t306) * t337 + t335 * t338;
t292 = t297 * qJD(1) - t335 * t306 - t337 * t338;
t290 = t325 * t348 + t293 * t326 + (t325 * t343 + t326 * t361) * qJD(5);
t289 = t326 * t348 - t293 * t325 + (-t325 * t361 + t326 * t343) * qJD(5);
t1 = [(-t296 * t326 + t325 * t351 + t317) * r_i_i_C(1) + (t366 * t325 + t326 * t351) * r_i_i_C(2) - t296 * t323 + t297 * qJD(4) + t335 * t350 - t337 * t314 + t362 * t295 + (-t337 * t324 + (t312 + (-t344 - t363) * t330) * t335) * qJD(1), -t343 * qJD(4) + t362 * t293 - t299 * t340 - t341 * t292 + ((t332 * t359 - t355) * qJD(2) + (-t332 * t355 + t359) * qJD(1)) * pkin(2), t348, t292, t289 * r_i_i_C(1) - t290 * r_i_i_C(2), 0; -t337 * t350 + t290 * r_i_i_C(1) + t289 * r_i_i_C(2) - t299 * qJD(4) + t293 * t323 - t335 * t314 + t362 * t292 + (-t324 * t335 + (t330 * t363 - t312) * t337) * qJD(1), t298 * qJD(4) - t362 * (t343 * qJD(1) + t337 * t305 - t357) - t297 * t340 + t341 * t295 + ((-t332 * t358 - t356) * qJD(2) + (-t332 * t356 - t358) * qJD(1)) * pkin(2), t349, -t295, t317 * r_i_i_C(2) + (t339 * r_i_i_C(1) - t296 * r_i_i_C(2)) * t326 + (-t366 * r_i_i_C(1) - t339 * r_i_i_C(2)) * t325, 0; 0, t309 * qJD(4) - t362 * t303 - t341 * t304 + (t315 * t340 - t350) * t330, 0, t304, t344 * t303 + ((-t309 * t326 - t325 * t332) * r_i_i_C(1) + (t309 * t325 - t326 * t332) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
