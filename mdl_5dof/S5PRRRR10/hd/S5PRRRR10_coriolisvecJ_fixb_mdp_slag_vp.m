% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:17
% EndTime: 2019-12-05 17:26:29
% DurationCPUTime: 5.41s
% Computational Cost: add. (2761->387), mult. (7683->581), div. (0->0), fcn. (6239->12), ass. (0->177)
t347 = sin(pkin(6));
t353 = sin(qJ(3));
t451 = t347 * t353;
t342 = pkin(8) * t451;
t349 = cos(pkin(6));
t357 = cos(qJ(3));
t358 = cos(qJ(2));
t443 = t357 * t358;
t354 = sin(qJ(2));
t447 = t353 * t354;
t373 = -t349 * t447 + t443;
t348 = sin(pkin(5));
t435 = qJD(1) * t348;
t448 = t349 * t357;
t439 = t373 * t435 - (pkin(2) * t448 - t342) * qJD(3);
t389 = pkin(3) * t353 - pkin(9) * t357;
t372 = t389 * qJD(3);
t415 = t354 * t435;
t472 = (t372 - t415) * t347;
t352 = sin(qJ(4));
t356 = cos(qJ(4));
t421 = t349 * qJD(2);
t401 = qJD(3) + t421;
t433 = qJD(2) * t347;
t414 = t353 * t433;
t471 = -t352 * t414 + t356 * t401;
t431 = qJD(2) * t357;
t413 = t347 * t431;
t340 = -qJD(4) + t413;
t303 = qJD(5) - t471;
t450 = t347 * t357;
t419 = pkin(8) * t450;
t313 = t419 + (pkin(2) * t353 + pkin(9)) * t349;
t390 = -pkin(3) * t357 - pkin(9) * t353;
t314 = (-pkin(2) + t390) * t347;
t438 = t356 * t313 + t352 * t314;
t470 = t438 * qJD(4) - t439 * t352 - t472 * t356;
t425 = qJD(4) * t356;
t427 = qJD(4) * t352;
t469 = -t313 * t427 + t314 * t425 + t472 * t352 - t439 * t356;
t445 = t354 * t357;
t446 = t353 * t358;
t375 = t349 * t445 + t446;
t449 = t349 * t353;
t437 = -t375 * t435 + (pkin(2) * t449 + t419) * qJD(3);
t326 = pkin(8) * t433 + t415;
t465 = qJD(2) * pkin(2);
t333 = t358 * t435 + t465;
t468 = t357 * t326 + t333 * t449;
t350 = cos(pkin(5));
t434 = qJD(1) * t350;
t416 = t347 * t434;
t277 = t353 * t416 + t468;
t268 = pkin(9) * t401 + t277;
t341 = t349 * t434;
t281 = t341 + (qJD(2) * t390 - t333) * t347;
t250 = t268 * t356 + t281 * t352;
t320 = t353 * t326;
t365 = t373 * qJD(2);
t422 = t347 * qJD(3);
t409 = t357 * t422;
t394 = t350 * t409;
t260 = (t333 * t448 - t320) * qJD(3) + (t348 * t365 + t394) * qJD(1);
t290 = (t372 + t415) * t433;
t362 = -qJD(4) * t250 - t260 * t352 + t356 * t290;
t420 = qJD(2) * qJD(3);
t408 = t347 * t420;
t393 = t353 * t408;
t240 = -pkin(4) * t393 - t362;
t308 = t352 * t401 + t356 * t414;
t467 = t303 * (pkin(4) * t308 + t303 * pkin(10)) + t240;
t455 = t333 * t349;
t276 = t357 * (t416 + t455) - t320;
t344 = t347 ^ 2;
t466 = (-t353 * t357 * MDP(5) + (t353 ^ 2 - t357 ^ 2) * MDP(6)) * t344;
t359 = qJD(2) ^ 2;
t392 = t357 * t408;
t286 = t471 * qJD(4) + t356 * t392;
t351 = sin(qJ(5));
t355 = cos(qJ(5));
t423 = qJD(5) * t355;
t417 = t355 * t286 - t340 * t423 + t351 * t393;
t424 = qJD(5) * t351;
t254 = -t308 * t424 + t417;
t463 = t254 * t351;
t456 = t308 * t351;
t283 = t355 * t340 + t456;
t462 = t283 * t303;
t285 = t308 * t355 - t340 * t351;
t461 = t285 * t303;
t378 = t352 * t392;
t287 = qJD(4) * t308 + t378;
t460 = t287 * t351;
t459 = t287 * t355;
t458 = t471 * t340;
t457 = t308 * t340;
t339 = -pkin(4) * t356 - pkin(10) * t352 - pkin(3);
t454 = t339 * t287;
t453 = t340 * t352;
t452 = t340 * t356;
t444 = t356 * t357;
t410 = t353 * t422;
t442 = -pkin(4) * t410 + t470;
t441 = t277 + t340 * (pkin(4) * t352 - pkin(10) * t356);
t315 = t389 * t433;
t440 = t356 * t276 + t352 * t315;
t432 = qJD(2) * t348;
t430 = qJD(3) * t353;
t429 = qJD(3) * t356;
t428 = qJD(4) * t351;
t426 = qJD(4) * t355;
t418 = t351 * t450;
t412 = t354 * t432;
t368 = t356 * t260 - t268 * t427 + t281 * t425 + t352 * t290;
t239 = pkin(10) * t393 + t368;
t366 = t375 * qJD(2);
t395 = t350 * t410;
t261 = t468 * qJD(3) + (t348 * t366 + t395) * qJD(1);
t244 = pkin(4) * t287 - pkin(10) * t286 + t261;
t406 = -t239 * t351 + t355 * t244;
t404 = t286 * t351 - t355 * t393;
t403 = t303 * t355;
t402 = pkin(10) * t414 - qJD(5) * t339 + t440;
t397 = t347 * t412;
t391 = MDP(16) * t414;
t301 = (t351 * t353 + t355 * t444) * t433;
t387 = t355 * t425 - t301;
t312 = t342 + (-pkin(2) * t357 - pkin(3)) * t349;
t323 = -t356 * t349 + t352 * t451;
t324 = t349 * t352 + t356 * t451;
t271 = pkin(4) * t323 - pkin(10) * t324 + t312;
t386 = -pkin(10) * t410 - qJD(5) * t271 - t469;
t273 = -pkin(10) * t450 + t438;
t293 = -qJD(4) * t323 + t356 * t409;
t294 = qJD(4) * t324 + t352 * t409;
t385 = -pkin(4) * t294 + pkin(10) * t293 + qJD(5) * t273 - t437;
t384 = t239 * t355 + t244 * t351;
t248 = -pkin(10) * t340 + t250;
t267 = -pkin(3) * t401 - t276;
t251 = -pkin(4) * t471 - t308 * pkin(10) + t267;
t242 = t248 * t355 + t251 * t351;
t383 = t248 * t351 - t251 * t355;
t249 = -t268 * t352 + t281 * t356;
t374 = t349 * t446 + t445;
t292 = t348 * t374 + t350 * t451;
t322 = -t347 * t348 * t358 + t349 * t350;
t270 = t292 * t356 + t322 * t352;
t376 = t349 * t443 - t447;
t291 = -t348 * t376 - t350 * t450;
t382 = t270 * t355 + t291 * t351;
t381 = -t270 * t351 + t291 * t355;
t380 = -t276 * t352 + t315 * t356;
t269 = t292 * t352 - t322 * t356;
t379 = -t313 * t352 + t314 * t356;
t302 = -t333 * t347 + t341;
t377 = t302 * t347 - t344 * t465;
t295 = t324 * t351 + t355 * t450;
t370 = -t303 * t423 - t460;
t369 = -t303 * t424 + t459;
t364 = qJD(1) * t349 * t412 + qJD(3) * t326;
t247 = pkin(4) * t340 - t249;
t363 = -pkin(10) * t287 + (t247 + t249) * t303;
t360 = -t302 * t433 - qJD(3) * t455 + (-t350 * t422 - t358 * t432) * qJD(1);
t300 = t351 * t356 * t413 - t355 * t414;
t296 = t324 * t355 - t418;
t272 = pkin(4) * t450 - t379;
t266 = t394 + (qJD(3) * t376 + t365) * t348;
t265 = t395 + (qJD(3) * t374 + t366) * t348;
t263 = -qJD(5) * t418 + t293 * t351 + t324 * t423 - t355 * t410;
t262 = -qJD(5) * t295 + t293 * t355 + t351 * t410;
t257 = -pkin(4) * t414 - t380;
t255 = t285 * qJD(5) + t404;
t246 = -qJD(4) * t269 + t266 * t356 + t352 * t397;
t245 = qJD(4) * t270 + t266 * t352 - t356 * t397;
t238 = -qJD(5) * t242 + t406;
t237 = -qJD(5) * t383 + t384;
t1 = [(-t265 * t401 + t322 * t393) * MDP(10) + (-t266 * t401 + t322 * t392) * MDP(11) + (t245 * t340 - t265 * t471 - t269 * t393 + t287 * t291) * MDP(17) + (t246 * t340 + t265 * t308 - t270 * t393 + t286 * t291) * MDP(18) + ((-qJD(5) * t382 - t246 * t351 + t265 * t355) * t303 + t381 * t287 + t245 * t283 + t269 * t255) * MDP(24) + (-(qJD(5) * t381 + t246 * t355 + t265 * t351) * t303 - t382 * t287 + t245 * t285 + t269 * t254) * MDP(25) + (-MDP(4) * t358 + (-MDP(3) + (-MDP(10) * t357 + MDP(11) * t353) * t344) * t354) * t359 * t348; ((-qJD(2) * t437 - t261) * t349 + (t353 * t377 - t437) * qJD(3)) * MDP(10) + ((qJD(2) * t439 - t260) * t349 + (t357 * t377 + t439) * qJD(3)) * MDP(11) + (t286 * t324 + t293 * t308) * MDP(12) + (-t286 * t323 - t287 * t324 + t293 * t471 - t294 * t308) * MDP(13) + (-t293 * t340 + (-t286 * t357 + (qJD(2) * t324 + t308) * t430) * t347) * MDP(14) + (t294 * t340 + (t287 * t357 + (-qJD(2) * t323 + t471) * t430) * t347) * MDP(15) + (-t340 * t347 - t344 * t431) * MDP(16) * t430 + (t261 * t323 + t267 * t294 + t312 * t287 + t470 * t340 - t437 * t471 + (-t362 * t357 + (qJD(2) * t379 + t249) * t430) * t347) * MDP(17) + (t261 * t324 + t267 * t293 + t312 * t286 + t469 * t340 + t437 * t308 + (t368 * t357 + (-qJD(2) * t438 - t250) * t430) * t347) * MDP(18) + (t254 * t296 + t262 * t285) * MDP(19) + (-t254 * t295 - t255 * t296 - t262 * t283 - t263 * t285) * MDP(20) + (t254 * t323 + t262 * t303 + t285 * t294 + t287 * t296) * MDP(21) + (-t255 * t323 - t263 * t303 - t283 * t294 - t287 * t295) * MDP(22) + (t287 * t323 + t294 * t303) * MDP(23) + ((t271 * t355 - t273 * t351) * t287 + t238 * t323 - t383 * t294 + t272 * t255 + t240 * t295 + t247 * t263 + (t351 * t386 - t355 * t385) * t303 + t442 * t283) * MDP(24) + (-(t271 * t351 + t273 * t355) * t287 - t237 * t323 - t242 * t294 + t272 * t254 + t240 * t296 + t247 * t262 + (t351 * t385 + t355 * t386) * t303 + t442 * t285) * MDP(25) - 0.2e1 * t466 * t420 + (MDP(7) * t409 - MDP(8) * t410) * (qJD(3) + 0.2e1 * t421); (t277 * t401 + t353 * t360 - t357 * t364) * MDP(10) + (t276 * t401 + t353 * t364 + t357 * t360) * MDP(11) + (t286 * t352 - t308 * t452) * MDP(12) + ((t286 - t458) * t356 + (-t287 + t457) * t352) * MDP(13) + (-t340 * t425 + (t340 * t444 + (qJD(3) * t352 - t308) * t353) * t433) * MDP(14) + (t340 * t427 + (-t357 * t453 + (-t471 + t429) * t353) * t433) * MDP(15) + t340 * t391 + (-pkin(3) * t287 - t261 * t356 + t380 * t340 + t277 * t471 + (pkin(9) * t452 + t267 * t352) * qJD(4) + (-t249 * t353 + (-pkin(9) * t430 - t267 * t357) * t352) * t433) * MDP(17) + (-pkin(3) * t286 + t261 * t352 - t440 * t340 - t277 * t308 + (-pkin(9) * t453 + t267 * t356) * qJD(4) + (-t267 * t444 + (-pkin(9) * t429 + t250) * t353) * t433) * MDP(18) + (t254 * t352 * t355 + (-t352 * t424 + t387) * t285) * MDP(19) + (t283 * t301 + t285 * t300 + (-t283 * t355 - t285 * t351) * t425 + (-t463 - t255 * t355 + (t283 * t351 - t285 * t355) * qJD(5)) * t352) * MDP(20) + (-t254 * t356 + t387 * t303 + (-t285 * t340 + t369) * t352) * MDP(21) + (t255 * t356 + (-t351 * t425 + t300) * t303 + (t283 * t340 + t370) * t352) * MDP(22) + (-t287 * t356 - t303 * t453) * MDP(23) + (t355 * t454 - t247 * t300 - t257 * t283 + (t351 * t402 - t355 * t441) * t303 + (t247 * t428 - t238 + (qJD(4) * t283 + t370) * pkin(9)) * t356 + (t247 * t423 + t240 * t351 + t340 * t383 + (t303 * t428 + t255) * pkin(9)) * t352) * MDP(24) + (-t351 * t454 - t247 * t301 - t257 * t285 + (t351 * t441 + t355 * t402) * t303 + (t247 * t426 + t237 + (qJD(4) * t285 - t369) * pkin(9)) * t356 + (-t247 * t424 + t240 * t355 + t340 * t242 + (t303 * t426 + t254) * pkin(9)) * t352) * MDP(25) + ((-MDP(7) * t357 + MDP(8) * t353) * t347 * t349 + t466) * t359; -t471 ^ 2 * MDP(13) + (t286 + t458) * MDP(14) + (-t378 - t457) * MDP(15) + qJD(3) * t391 + (-t250 * t340 + t362) * MDP(17) + (-t249 * t340 - t267 * t471 - t368) * MDP(18) + (t285 * t403 + t463) * MDP(19) + ((t254 - t462) * t355 + (-t255 - t461) * t351) * MDP(20) + (t303 * t403 + t460) * MDP(21) + (-t303 ^ 2 * t351 + t459) * MDP(22) + (-pkin(4) * t255 - t250 * t283 + t363 * t351 - t467 * t355) * MDP(24) + (-pkin(4) * t254 - t250 * t285 + t467 * t351 + t363 * t355) * MDP(25) + (-MDP(12) * t471 + t308 * MDP(13) - MDP(15) * qJD(4) - t267 * MDP(17) - t285 * MDP(21) + t283 * MDP(22) - t303 * MDP(23) + MDP(24) * t383 + t242 * MDP(25)) * t308; t285 * t283 * MDP(19) + (-t283 ^ 2 + t285 ^ 2) * MDP(20) + (t417 + t462) * MDP(21) + (-t404 + t461) * MDP(22) + t287 * MDP(23) + (t242 * t303 - t247 * t285 + t406) * MDP(24) + (t247 * t283 - t303 * t383 - t384) * MDP(25) + (-MDP(21) * t456 - MDP(22) * t285 - MDP(24) * t242 + MDP(25) * t383) * qJD(5);];
tauc = t1;
