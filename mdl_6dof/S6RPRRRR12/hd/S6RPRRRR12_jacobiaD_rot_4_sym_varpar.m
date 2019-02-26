% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR12_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiaD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:09
% EndTime: 2019-02-26 21:21:12
% DurationCPUTime: 1.99s
% Computational Cost: add. (7906->122), mult. (24067->264), div. (442->12), fcn. (30860->17), ass. (0->136)
t333 = cos(pkin(6));
t327 = sin(pkin(14));
t339 = cos(qJ(1));
t396 = t339 * t327;
t331 = cos(pkin(14));
t336 = sin(qJ(1));
t397 = t336 * t331;
t320 = t333 * t396 + t397;
t335 = sin(qJ(3));
t338 = cos(qJ(3));
t329 = sin(pkin(7));
t332 = cos(pkin(7));
t330 = sin(pkin(6));
t400 = t330 * t339;
t325 = t336 * t327;
t395 = t339 * t331;
t431 = -t333 * t395 + t325;
t353 = t329 * t400 + t431 * t332;
t444 = -t320 * t338 + t335 * t353;
t383 = t333 * t325;
t394 = qJD(1) * t339;
t319 = -qJD(1) * t383 + t331 * t394;
t443 = -qJD(3) * t444 + t319 * t335;
t301 = t320 * t335 + t338 * t353;
t314 = -t329 * t431 + t332 * t400;
t328 = sin(pkin(8));
t425 = cos(pkin(8));
t441 = -t301 * t328 + t314 * t425;
t363 = t333 * t397 + t396;
t398 = t336 * t330;
t350 = t329 * t398 - t332 * t363;
t321 = -t383 + t395;
t403 = t321 * t335;
t303 = t338 * t350 - t403;
t337 = cos(qJ(4));
t351 = -t329 * t363 - t332 * t398;
t348 = t351 * t328;
t377 = t337 * t425;
t304 = t321 * t338 + t335 * t350;
t334 = sin(qJ(4));
t409 = t304 * t334;
t274 = -t303 * t377 + t337 * t348 + t409;
t268 = t274 ^ 2;
t346 = t303 * t425 - t348;
t275 = t304 * t337 + t334 * t346;
t270 = 0.1e1 / t275 ^ 2;
t265 = t268 * t270 + 0.1e1;
t263 = 0.1e1 / t265;
t317 = t320 * qJD(1);
t316 = t431 * qJD(1);
t378 = t330 * t394;
t358 = t316 * t332 + t329 * t378;
t283 = -qJD(3) * t403 + t358 * t335 + (qJD(3) * t350 - t317) * t338;
t282 = -t304 * qJD(3) + t317 * t335 + t358 * t338;
t359 = -t316 * t329 + t332 * t378;
t347 = t282 * t425 + t328 * t359;
t262 = t283 * t337 + t347 * t334 + (t337 * t346 - t409) * qJD(4);
t269 = 0.1e1 / t275;
t419 = t262 * t269 * t270;
t261 = qJD(4) * t275 + t283 * t334 - t337 * t347;
t415 = t270 * t274;
t422 = (t261 * t415 - t268 * t419) / t265 ^ 2;
t439 = t263 * t419 + t270 * t422;
t426 = 0.2e1 * t274;
t318 = t363 * qJD(1);
t379 = qJD(1) * t398;
t357 = t318 * t329 + t332 * t379;
t438 = t357 * t425;
t372 = t329 * t379;
t356 = t318 * t332 - t372;
t436 = (t320 * qJD(3) + t356) * t335 + (t353 * qJD(3) - t319) * t338;
t399 = t332 * t338;
t267 = -t438 + (-t318 * t399 + t338 * t372 - t443) * t328;
t289 = t441 ^ 2;
t401 = t329 * t333;
t429 = t330 * (-t327 * t335 + t331 * t399) + t338 * t401;
t299 = -t429 * t328 + (-t330 * t331 * t329 + t333 * t332) * t425;
t297 = 0.1e1 / t299 ^ 2;
t280 = t289 * t297 + 0.1e1;
t296 = 0.1e1 / t299;
t298 = t296 * t297;
t313 = -t335 * t401 + (-t331 * t332 * t335 - t327 * t338) * t330;
t309 = t313 * qJD(3);
t407 = t309 * t328;
t421 = (t267 * t297 * t441 + t289 * t298 * t407) / t280 ^ 2;
t435 = -0.2e1 * t421;
t411 = t441 * t313;
t433 = t328 * (t296 * t444 + t297 * t411);
t391 = -0.2e1 * t422;
t417 = t263 * t270;
t430 = -t262 * t417 + t269 * t391;
t361 = t301 * t425 + t314 * t328;
t428 = -t334 * t444 + t337 * t361;
t389 = t261 * t417;
t418 = t263 * t269;
t427 = qJD(4) * t418 + t426 * t439 - t389;
t281 = atan2(t441, t299);
t276 = sin(t281);
t277 = cos(t281);
t260 = t276 * t441 + t277 * t299;
t257 = 0.1e1 / t260;
t258 = 0.1e1 / t260 ^ 2;
t293 = t303 * t328 + t351 * t425;
t290 = t293 ^ 2;
t255 = t290 * t258 + 0.1e1;
t266 = t282 * t328 - t359 * t425;
t416 = t266 * t258;
t278 = 0.1e1 / t280;
t385 = t297 * t407;
t360 = t267 * t296 + t385 * t441;
t251 = t360 * t278;
t368 = -t276 * t299 + t277 * t441;
t247 = t251 * t368 + t276 * t267 - t277 * t407;
t423 = t247 * t257 * t258;
t424 = (-t290 * t423 + t293 * t416) / t255 ^ 2;
t420 = t258 * t293;
t414 = t276 * t293;
t413 = t277 * t293;
t412 = t441 * t296;
t408 = t309 * t328 ^ 2;
t393 = -0.2e1 * t424;
t392 = -0.2e1 * t423;
t387 = t263 * t415;
t376 = t425 * t283;
t374 = t296 * t435;
t285 = t356 * t338 + t443;
t362 = t285 * t425 - t328 * t357;
t355 = -t303 * t334 - t304 * t377;
t354 = t276 + (t277 * t412 - t276) * t278;
t308 = t429 * qJD(3);
t253 = 0.1e1 / t255;
t252 = t278 * t433;
t250 = t354 * t293;
t248 = (t276 * t444 - t277 * t313) * t328 + t368 * t252;
t246 = t433 * t435 + (0.2e1 * t298 * t408 * t411 + t436 * t296 * t328 + (t444 * t408 + (t267 * t313 - t308 * t441) * t328) * t297) * t278;
t1 = [t293 * t374 + (t266 * t296 + t293 * t385) * t278, 0, t246, 0, 0, 0; t441 * t257 * t393 + ((-t285 * t328 - t438) * t257 + (-t247 * t441 + t250 * t266) * t258) * t253 + ((t250 * t392 + t354 * t416) * t253 + (t250 * t393 + ((-t251 * t278 * t412 + 0.2e1 * t421) * t414 + (t441 * t374 + t251 + (-t251 + t360) * t278) * t413) * t253) * t258) * t293, 0, 0.2e1 * (-t257 * t304 * t328 - t248 * t420) * t424 + ((t368 * t246 + (-t251 * t260 + t267 * t277) * t252) * t420 + (t293 * t392 + t416) * t248 + (t283 * t257 + (-t304 * t247 + (t251 * t444 + t308) * t413 + (t251 * t313 + t252 * t309 + t436) * t414) * t258) * t328) * t253, 0, 0, 0; (t334 * t436 - t337 * t362) * t418 - (t428 * qJD(4) + t362 * t334 + t337 * t436) * t387 - t430 * t428 + t427 * (t334 * t361 + t337 * t444) 0 (t282 * t334 + t337 * t376) * t418 - (qJD(4) * t355 + t282 * t337 - t334 * t376) * t387 - t430 * t355 + t427 * (t303 * t337 - t409 * t425) t391 + (-t274 * t439 + t389) * t426, 0, 0;];
JaD_rot  = t1;
