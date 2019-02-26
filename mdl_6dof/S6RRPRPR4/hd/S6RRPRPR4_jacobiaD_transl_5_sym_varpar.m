% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:39:34
% EndTime: 2019-02-26 21:39:34
% DurationCPUTime: 0.42s
% Computational Cost: add. (445->103), mult. (1155->168), div. (0->0), fcn. (1210->12), ass. (0->64)
t332 = cos(pkin(6));
t329 = sin(pkin(11));
t331 = cos(pkin(11));
t335 = sin(qJ(2));
t338 = cos(qJ(2));
t345 = t338 * t329 + t335 * t331;
t312 = t345 * t332;
t355 = qJD(2) * t338;
t356 = qJD(2) * t335;
t367 = t329 * t356 - t331 * t355;
t306 = t367 * t332;
t316 = t335 * t329 - t338 * t331;
t339 = cos(qJ(1));
t336 = sin(qJ(1));
t357 = qJD(1) * t336;
t314 = -t329 * t355 - t331 * t356;
t360 = t336 * t314;
t297 = t360 - t312 * t357 + (-qJD(1) * t316 - t306) * t339;
t330 = sin(pkin(6));
t363 = t330 * t339;
t348 = qJD(4) * t363;
t368 = t297 - t348;
t366 = pkin(2) * t332;
t365 = r_i_i_C(3) + qJ(5) + pkin(9);
t364 = t330 * t336;
t362 = t335 * t336;
t361 = t335 * t339;
t359 = t336 * t338;
t358 = t338 * t339;
t299 = t339 * t312 - t336 * t316;
t354 = qJD(4) * t299;
t353 = pkin(2) * t356;
t352 = t330 * t357;
t351 = qJD(1) * t363;
t328 = qJ(4) + pkin(12);
t326 = sin(t328);
t327 = cos(t328);
t347 = -t326 * r_i_i_C(1) - t327 * r_i_i_C(2);
t311 = t316 * t332;
t298 = -t339 * t311 - t336 * t345;
t346 = t336 * t312 + t339 * t316;
t337 = cos(qJ(4));
t324 = t337 * pkin(4) + pkin(3);
t344 = t327 * r_i_i_C(1) - t326 * r_i_i_C(2) + t324;
t343 = t352 - t354;
t310 = t345 * t330;
t342 = t316 * qJD(2);
t334 = sin(qJ(4));
t341 = t334 * pkin(4) - t347;
t340 = qJD(4) * t341;
t294 = -t299 * qJD(1) + t336 * t306 + t339 * t314;
t325 = t338 * pkin(2) + pkin(1);
t318 = t327 * t348;
t315 = -t330 * qJD(3) + t355 * t366;
t313 = t335 * t366 + (-pkin(8) - qJ(3)) * t330;
t307 = qJD(2) * t312;
t305 = qJD(2) * t310;
t304 = t367 * t330;
t300 = t336 * t311 - t339 * t345;
t296 = t311 * t357 + (-qJD(1) * t345 - t307) * t339 + t336 * t342;
t293 = t298 * qJD(1) - t336 * t307 - t339 * t342;
t291 = t326 * t351 + t294 * t327 + (t326 * t346 + t327 * t364) * qJD(4);
t290 = t327 * t351 - t294 * t326 + (-t326 * t364 + t327 * t346) * qJD(4);
t1 = [(-t297 * t327 + t326 * t354 + t318) * r_i_i_C(1) + (t326 * t368 + t327 * t354) * r_i_i_C(2) - t297 * t324 + t298 * qJD(5) + t336 * t353 - t339 * t315 + t365 * t296 + (-t339 * t325 + (t347 * t330 + t313) * t336) * qJD(1) + (-t334 * t352 + (t299 * t334 + t337 * t363) * qJD(4)) * pkin(4), -t346 * qJD(5) + t365 * t294 - t344 * t293 - t300 * t340 + ((t332 * t362 - t358) * qJD(2) + (-t332 * t358 + t362) * qJD(1)) * pkin(2), t351, t290 * r_i_i_C(1) - t291 * r_i_i_C(2) + (t337 * t351 - t294 * t334 + (-t334 * t364 + t337 * t346) * qJD(4)) * pkin(4), t293, 0; -t339 * t353 + t291 * r_i_i_C(1) + t290 * r_i_i_C(2) - t300 * qJD(5) + t294 * t324 - t336 * t315 + t365 * t293 + (-t313 * t339 - t325 * t336) * qJD(1) + (t334 * t351 + (t334 * t346 + t337 * t364) * qJD(4)) * pkin(4), t299 * qJD(5) - t365 * (t346 * qJD(1) + t339 * t306 - t360) + t344 * t296 - t298 * t340 + ((-t332 * t361 - t359) * qJD(2) + (-t332 * t359 - t361) * qJD(1)) * pkin(2), t352, t318 * r_i_i_C(2) + (t343 * r_i_i_C(1) - t297 * r_i_i_C(2)) * t327 + (-r_i_i_C(1) * t368 - t343 * r_i_i_C(2)) * t326 + (t337 * t352 - t297 * t334 + (-t299 * t337 + t334 * t363) * qJD(4)) * pkin(4), -t296, 0; 0, t310 * qJD(5) - t365 * t304 - t344 * t305 + (t316 * t340 - t353) * t330, 0, t341 * t304 + ((-t310 * t327 - t326 * t332) * r_i_i_C(1) + (t310 * t326 - t327 * t332) * r_i_i_C(2) + (-t310 * t337 - t332 * t334) * pkin(4)) * qJD(4), t305, 0;];
JaD_transl  = t1;
