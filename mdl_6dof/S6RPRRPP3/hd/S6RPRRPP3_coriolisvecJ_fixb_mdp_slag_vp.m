% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:37:07
% EndTime: 2019-03-09 04:37:14
% DurationCPUTime: 3.94s
% Computational Cost: add. (3049->445), mult. (6767->553), div. (0->0), fcn. (3936->6), ass. (0->175)
t389 = sin(qJ(4));
t390 = sin(qJ(3));
t391 = cos(qJ(4));
t455 = qJD(4) * t391;
t433 = t390 * t455;
t392 = cos(qJ(3));
t457 = qJD(3) * t392;
t512 = t389 * t457 + t433;
t375 = sin(pkin(9)) * pkin(1) + pkin(7);
t357 = t375 * qJD(1);
t381 = t390 * qJD(2);
t325 = t392 * t357 + t381;
t315 = qJD(3) * pkin(8) + t325;
t376 = -cos(pkin(9)) * pkin(1) - pkin(2);
t338 = -pkin(3) * t392 - pkin(8) * t390 + t376;
t318 = t338 * qJD(1);
t283 = t315 * t389 - t391 * t318;
t459 = qJD(3) * t389;
t462 = qJD(1) * t390;
t347 = t391 * t462 + t459;
t404 = pkin(5) * t347 + t283;
t452 = qJD(5) + t404;
t447 = MDP(18) - MDP(21);
t511 = -MDP(20) + MDP(25);
t384 = t390 ^ 2;
t510 = MDP(6) * (-t392 ^ 2 + t384);
t509 = t392 * MDP(5);
t461 = qJD(1) * t392;
t372 = -qJD(4) + t461;
t362 = qJD(5) * t372;
t450 = qJD(1) * qJD(3);
t377 = t390 * t450;
t371 = qJ(5) * t377;
t449 = qJD(3) * qJD(4);
t430 = t389 * t449;
t308 = qJD(1) * t512 + t430;
t505 = qJD(2) * t392 - t390 * t357;
t316 = t505 * qJD(3);
t416 = pkin(3) * t390 - pkin(8) * t392;
t354 = t416 * qJD(3);
t337 = qJD(1) * t354;
t456 = qJD(4) * t389;
t420 = t315 * t456 - t391 * t316 - t318 * t455 - t389 * t337;
t401 = -pkin(5) * t308 - t420;
t262 = -t362 + t371 + t401;
t498 = pkin(4) + qJ(6);
t269 = t372 * t498 + t452;
t508 = -t269 * t372 + t262;
t419 = t315 * t455 + t389 * t316 + t318 * t456 - t391 * t337;
t265 = -pkin(4) * t377 + t419;
t284 = t391 * t315 + t389 * t318;
t278 = qJ(5) * t372 - t284;
t507 = -t278 * t372 + t265;
t502 = qJD(5) * t389 + t325 + (t389 * t461 - t456) * pkin(4);
t448 = -MDP(17) + MDP(20);
t501 = MDP(19) + MDP(23);
t370 = t372 ^ 2;
t500 = pkin(5) + pkin(8);
t454 = t391 * qJD(3);
t345 = t389 * t462 - t454;
t499 = pkin(5) * t345;
t497 = qJ(5) * t308;
t496 = qJ(5) * t345;
t495 = qJ(5) * t391;
t431 = t391 * t450;
t434 = t390 * t456;
t307 = qJD(1) * t434 - t391 * t449 - t392 * t431;
t317 = qJD(3) * t381 + t357 * t457;
t398 = qJ(5) * t307 - qJD(5) * t347 + t317;
t268 = pkin(4) * t308 + t398;
t494 = t268 * t389;
t493 = t268 * t391;
t314 = -qJD(3) * pkin(3) - t505;
t490 = t314 * t389;
t489 = t314 * t391;
t488 = t317 * t389;
t487 = t317 * t391;
t486 = t345 * t347;
t485 = t345 * t372;
t484 = t347 * t372;
t483 = t372 * t391;
t482 = t375 * t389;
t481 = t375 * t391;
t480 = t389 * t390;
t479 = t389 * t392;
t478 = t390 * t391;
t393 = qJD(3) ^ 2;
t477 = t390 * t393;
t476 = t391 * t392;
t475 = t392 * t393;
t361 = t500 * t391;
t351 = t416 * qJD(1);
t424 = t351 * t391 - t389 * t505;
t445 = pkin(5) * t476;
t474 = -(-t390 * t498 + t445) * qJD(1) + t424 + qJD(4) * t361;
t446 = pkin(5) * t479;
t469 = t389 * t351 + t391 * t505;
t473 = (qJ(5) * t390 - t446) * qJD(1) + t469 + t500 * t456;
t413 = qJ(6) * t389 - t495;
t399 = t413 * t392;
t472 = qJD(1) * t399 - qJD(4) * t413 + qJD(6) * t391 + t502;
t471 = qJ(5) * t455 - t461 * t495 + t502;
t436 = t392 * t454;
t470 = -t308 * t478 - t345 * t436;
t468 = t338 * t455 + t389 * t354;
t350 = t375 * t476;
t467 = t389 * t338 + t350;
t466 = pkin(4) * t480 + t390 * t375;
t363 = t384 * t431;
t465 = -t372 * t436 + t363;
t463 = MDP(19) * t391;
t358 = qJD(1) * t376;
t458 = qJD(3) * t390;
t453 = -qJD(5) - t283;
t276 = t284 - t499;
t451 = -qJD(6) - t276;
t444 = pkin(8) * t372 * t389;
t443 = pkin(8) * t483;
t442 = pkin(8) * t458;
t441 = pkin(8) * t454;
t349 = t375 * t479;
t440 = -t307 * t480 + t512 * t347;
t437 = t345 * t458;
t435 = t375 * t456;
t432 = t347 * t458;
t429 = -qJ(5) * t389 - pkin(3);
t428 = -pkin(4) - t482;
t427 = MDP(25) - t448;
t426 = -MDP(24) + t447;
t425 = -qJ(5) + t481;
t423 = t338 * t391 - t349;
t418 = t512 * pkin(4) + qJ(5) * t434 + t375 * t457;
t293 = qJ(5) * t392 - t467;
t415 = -qJD(4) * t350 - t338 * t456 + t354 * t391;
t414 = -t362 - t420;
t264 = -t371 - t414;
t412 = -t264 * t391 + t265 * t389;
t271 = qJD(6) - t278 - t499;
t411 = t269 * t391 - t271 * t389;
t410 = t269 * t389 + t271 * t391;
t277 = pkin(4) * t372 - t453;
t409 = t277 * t391 + t278 * t389;
t408 = t277 * t389 - t278 * t391;
t406 = 0.2e1 * qJD(3) * t358;
t263 = qJD(6) * t345 + t308 * t498 + t398;
t396 = -qJ(5) * t347 + t314;
t274 = t345 * t498 + t396;
t403 = t263 * t389 + t274 * t455;
t402 = -t263 * t391 + t274 * t456;
t400 = pkin(5) * t307 - t419;
t286 = -t307 - t485;
t395 = -t377 * t498 - t400;
t383 = t392 * pkin(4);
t369 = 0.2e1 * t371;
t360 = t500 * t389;
t356 = -pkin(4) * t391 + t429;
t336 = -t391 * t498 + t429;
t309 = -qJ(5) * t478 + t466;
t298 = pkin(4) * t347 + t496;
t295 = t390 * t413 + t466;
t294 = t383 - t423;
t291 = -pkin(5) * t480 - t293;
t290 = -pkin(4) * t462 - t424;
t289 = -qJ(5) * t462 - t469;
t287 = t347 * t498 + t496;
t285 = qJ(6) * t392 + t349 + t383 + (pkin(5) * t390 - t338) * t391;
t282 = pkin(4) * t345 + t396;
t280 = (-qJ(5) * t457 - qJD(5) * t390) * t391 + t418;
t273 = t428 * t458 - t415;
t272 = (qJD(5) + t435) * t392 + t425 * t458 - t468;
t270 = qJD(3) * t399 + (qJD(6) * t389 + (qJ(6) * qJD(4) - qJD(5)) * t391) * t390 + t418;
t267 = -qJD(5) * t392 + (-pkin(5) * t478 - t349) * qJD(4) + (-t390 * t425 - t446) * qJD(3) + t468;
t266 = -pkin(5) * t434 + qJD(6) * t392 + (t445 + (-qJ(6) + t428) * t390) * qJD(3) - t415;
t261 = qJD(6) * t372 + t395;
t1 = [0.2e1 * t377 * t509 - 0.2e1 * t450 * t510 + MDP(7) * t475 - MDP(8) * t477 + (-t375 * t475 + t390 * t406) * MDP(10) + (t375 * t477 + t392 * t406) * MDP(11) + (-t307 * t478 + (-t434 + t436) * t347) * MDP(12) + (t345 * t434 - t440 + t470) * MDP(13) + (t307 * t392 + t372 * t434 + t432 + t465) * MDP(14) + (t372 * t433 + t308 * t392 - t437 + (-qJD(1) * t384 + t372 * t392) * t459) * MDP(15) + (-t372 - t461) * MDP(16) * t458 + (-t415 * t372 + ((t345 * t375 + t490) * qJD(3) + t419) * t392 + (t314 * t455 + t375 * t308 + t488 + (qJD(1) * t423 - t372 * t482 - t283) * qJD(3)) * t390) * MDP(17) + (t468 * t372 + (-t372 * t435 + (t347 * t375 + t489) * qJD(3) - t420) * t392 + (-t314 * t456 - t375 * t307 + t487 + (-qJD(1) * t467 - t372 * t481 - t284) * qJD(3)) * t390) * MDP(18) + (t272 * t345 + t273 * t347 + t293 * t308 - t294 * t307 + t409 * t457 + (-qJD(4) * t408 + t264 * t389 + t265 * t391) * t390) * MDP(19) + (-t273 * t372 - t280 * t345 - t308 * t309 + (-t282 * t459 - t265) * t392 + (-t282 * t455 - t494 + (qJD(1) * t294 + t277) * qJD(3)) * t390) * MDP(20) + (t272 * t372 - t280 * t347 + t307 * t309 + (-t282 * t454 + t264) * t392 + (t282 * t456 - t493 + (-qJD(1) * t293 - t278) * qJD(3)) * t390) * MDP(21) + (t264 * t293 + t265 * t294 + t268 * t309 + t272 * t278 + t273 * t277 + t280 * t282) * MDP(22) + (t266 * t347 - t267 * t345 - t285 * t307 - t291 * t308 + t411 * t457 + (-qJD(4) * t410 + t261 * t391 - t262 * t389) * t390) * MDP(23) + (-t267 * t372 - t270 * t347 + t295 * t307 + (-t274 * t454 - t262) * t392 + ((qJD(1) * t291 + t271) * qJD(3) + t402) * t390) * MDP(24) + (t266 * t372 + t270 * t345 + t295 * t308 + (t274 * t459 + t261) * t392 + ((-qJD(1) * t285 - t269) * qJD(3) + t403) * t390) * MDP(25) + (t261 * t285 + t262 * t291 + t263 * t295 + t266 * t269 + t267 * t271 + t270 * t274) * MDP(26); t440 * MDP(19) + (t440 + t470) * MDP(23) + (-t432 + t465) * MDP(24) + (-t393 * MDP(11) - t268 * MDP(22) - t263 * MDP(26) - t427 * t308 + t426 * t307 + (-t345 * t463 + t408 * MDP(22) + t410 * MDP(26) + (t389 * t427 + t391 * t447) * t372) * qJD(3)) * t392 + (-t393 * MDP(10) - t308 * t463 + (qJD(3) * t282 + t412) * MDP(22) + (qJD(3) * t274 + t261 * t389 + t262 * t391) * MDP(26) + (t409 * MDP(22) + t411 * MDP(26) + t501 * t389 * t345 + (-t389 * t426 + t391 * t427) * t372) * qJD(4)) * t390 + t447 * (t432 - t363) + (MDP(17) + t511) * (-t384 * t389 * t450 + t437); (qJD(3) * t325 - t358 * t462 - t317) * MDP(10) - t358 * t461 * MDP(11) + (-t307 * t389 - t347 * t483) * MDP(12) + ((-t307 + t485) * t391 + (-t308 + t484) * t389) * MDP(13) + (-t372 * t455 + (t372 * t476 + (-t347 + t459) * t390) * qJD(1)) * MDP(14) + (t372 * t456 + (-t372 * t479 + (t345 + t454) * t390) * qJD(1)) * MDP(15) + t372 * MDP(16) * t462 + (-pkin(3) * t308 - t487 + t424 * t372 - t325 * t345 + (t443 + t490) * qJD(4) + (t283 * t390 + (-t314 * t392 - t442) * t389) * qJD(1)) * MDP(17) + (pkin(3) * t307 + t488 - t469 * t372 - t325 * t347 + (-t444 + t489) * qJD(4) + (-t314 * t476 + (t284 - t441) * t390) * qJD(1)) * MDP(18) + (-t289 * t345 - t290 * t347 + (-t264 - t372 * t277 + (qJD(4) * t347 - t308) * pkin(8)) * t391 + ((qJD(4) * t345 - t307) * pkin(8) + t507) * t389) * MDP(19) + (t493 + t290 * t372 - t308 * t356 + t471 * t345 + (-t282 * t389 - t443) * qJD(4) + (-t277 * t390 + (t282 * t392 + t442) * t389) * qJD(1)) * MDP(20) + (-t494 - t289 * t372 + t307 * t356 + t471 * t347 + (-t282 * t391 + t444) * qJD(4) + (t282 * t476 + (t278 + t441) * t390) * qJD(1)) * MDP(21) + (t268 * t356 - t277 * t290 - t278 * t289 - t471 * t282 + (qJD(4) * t409 + t412) * pkin(8)) * MDP(22) + (-t307 * t360 - t308 * t361 + t474 * t347 + t473 * t345 + t508 * t391 + (t271 * t372 + t261) * t389) * MDP(23) + (t307 * t336 + t473 * t372 + t472 * t347 + (t274 * t476 + (qJD(3) * t361 - t271) * t390) * qJD(1) - t403) * MDP(24) + (t308 * t336 + t474 * t372 - t472 * t345 + (-t274 * t479 + (-qJD(3) * t360 + t269) * t390) * qJD(1) + t402) * MDP(25) + (t261 * t360 + t262 * t361 + t263 * t336 + t269 * t474 - t271 * t473 - t274 * t472) * MDP(26) + (-t390 * t509 + t510) * qJD(1) ^ 2; t286 * MDP(14) - MDP(15) * t430 + t420 * MDP(18) + (pkin(4) * t307 - t497) * MDP(19) + (t369 + t414) * MDP(21) + (-pkin(4) * t265 - qJ(5) * t264 - t277 * t284 + t278 * t453 - t282 * t298) * MDP(22) + (t307 * t498 - t497) * MDP(23) + (-0.2e1 * t362 + t369 + t401) * MDP(24) + t400 * MDP(25) + (qJ(5) * t262 - t261 * t498 + t269 * t451 + t271 * t452 - t274 * t287) * MDP(26) + (t283 * MDP(18) + t453 * MDP(21) - t404 * MDP(24) + (-0.2e1 * qJD(6) - t276) * MDP(25) + t448 * t284) * t372 + (-t372 * MDP(15) - t314 * MDP(17) + (-t278 - t284) * MDP(19) + t282 * MDP(20) + t298 * MDP(21) + (t271 + t451) * MDP(23) + t287 * MDP(24) - t274 * MDP(25) + MDP(13) * t347) * t347 + (-MDP(15) * t433 + (-MDP(15) * t479 + (-0.2e1 * pkin(4) * MDP(20) + 0.2e1 * t498 * MDP(25) + MDP(16)) * t390) * qJD(3)) * qJD(1) + (t347 * MDP(12) + t314 * MDP(18) + (t277 + t453) * MDP(19) + t298 * MDP(20) - t282 * MDP(21) + (t269 - t452) * MDP(23) - t274 * MDP(24) - t287 * MDP(25) - MDP(13) * t345) * t345 + t448 * t419; (t282 * t347 + t507) * MDP(22) + (t274 * t347 + (qJD(6) + t271) * t372 + t395) * MDP(26) + t511 * (-t377 + t486) + t501 * t286 + (MDP(21) + MDP(24)) * (-t347 ^ 2 - t370); (-t308 - t484) * MDP(23) + (t377 + t486) * MDP(24) + (-t345 ^ 2 - t370) * MDP(25) + (-t274 * t345 + t508) * MDP(26);];
tauc  = t1;
