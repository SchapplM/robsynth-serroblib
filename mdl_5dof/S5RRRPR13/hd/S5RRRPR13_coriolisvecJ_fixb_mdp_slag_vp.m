% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR13_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR13_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:44
% EndTime: 2019-12-31 21:46:56
% DurationCPUTime: 5.42s
% Computational Cost: add. (3378->453), mult. (9105->616), div. (0->0), fcn. (6719->8), ass. (0->186)
t407 = sin(qJ(3));
t410 = cos(qJ(3));
t517 = cos(pkin(5));
t463 = t517 * qJD(1);
t430 = t463 + qJD(2);
t408 = sin(qJ(2));
t405 = sin(pkin(5));
t491 = qJD(1) * t405;
t470 = t408 * t491;
t352 = t407 * t430 + t410 * t470;
t345 = qJD(5) + t352;
t411 = cos(qJ(2));
t490 = qJD(1) * t411;
t469 = t405 * t490;
t390 = -qJD(3) + t469;
t449 = pkin(1) * t463;
t395 = t408 * t449;
t369 = pkin(7) * t469 + t395;
t332 = pkin(8) * t430 + t369;
t365 = (-pkin(2) * t411 - pkin(8) * t408 - pkin(1)) * t405;
t344 = qJD(1) * t365;
t299 = t332 * t407 - t410 * t344;
t482 = -qJD(4) - t299;
t402 = t405 ^ 2;
t480 = qJD(1) * qJD(2);
t527 = -0.2e1 * t402 * t480;
t526 = (t408 ^ 2 - t411 ^ 2) * MDP(5);
t465 = t405 * t480;
t448 = t408 * t465;
t439 = pkin(3) * t448;
t426 = (pkin(2) * t408 - pkin(8) * t411) * t405;
t368 = qJD(2) * t426;
t360 = qJD(1) * t368;
t471 = pkin(1) * t517;
t501 = t405 * t408;
t524 = -pkin(7) * t501 + t411 * t471;
t370 = t524 * qJD(2);
t361 = qJD(1) * t370;
t486 = qJD(3) * t410;
t487 = qJD(3) * t407;
t454 = t332 * t486 + t344 * t487 - t410 * t360 + t407 * t361;
t277 = -t439 + t454;
t300 = t410 * t332 + t407 * t344;
t297 = qJ(4) * t390 - t300;
t525 = -t297 * t390 + t277;
t447 = t411 * t465;
t488 = qJD(3) * t352;
t318 = t407 * t447 + t488;
t451 = t407 * t469;
t523 = qJD(4) * t407 + t369 + (t451 - t487) * pkin(3);
t481 = pkin(4) * t352 - t482;
t421 = t410 * t430;
t350 = t407 * t470 - t421;
t518 = pkin(4) * t350;
t285 = -t297 - t518;
t475 = t407 * t501;
t450 = qJD(3) * t475;
t317 = qJD(1) * t450 - qJD(3) * t421 - t410 * t447;
t520 = pkin(3) + pkin(9);
t522 = t520 * t317 + (t285 - t300 + t518) * t345;
t453 = t408 * t471;
t500 = t405 * t411;
t364 = pkin(7) * t500 + pkin(8) * t517 + t453;
t422 = t364 * t487 - t365 * t486 - t407 * t368 - t410 * t370;
t483 = t408 * qJD(2);
t278 = -t405 * (qJ(4) * t483 - qJD(4) * t411) + t422;
t521 = t352 ^ 2;
t413 = qJD(1) ^ 2;
t519 = pkin(4) + pkin(8);
t516 = qJ(4) * t350;
t515 = qJ(4) * t410;
t406 = sin(qJ(5));
t409 = cos(qJ(5));
t484 = qJD(5) * t409;
t473 = t406 * t318 + t350 * t484 + t409 * t448;
t485 = qJD(5) * t406;
t283 = t390 * t485 + t473;
t514 = t283 * t409;
t366 = -pkin(7) * t470 + t411 * t449;
t331 = -pkin(2) * t430 - t366;
t414 = -t352 * qJ(4) + t331;
t291 = t350 * pkin(3) + t414;
t513 = t291 * t352;
t504 = t390 * t406;
t314 = -t409 * t350 - t504;
t511 = t314 * t345;
t510 = t314 * t390;
t316 = t350 * t406 - t390 * t409;
t509 = t316 * t345;
t508 = t316 * t390;
t507 = t350 * t352;
t506 = t350 * t390;
t505 = t352 * t390;
t503 = t390 * t410;
t502 = t402 * t413;
t499 = t406 * t317;
t498 = t407 * t411;
t311 = t409 * t317;
t497 = t410 * t411;
t496 = qJ(4) * t486 - t469 * t515 + t523;
t495 = t410 * t364 + t407 * t365;
t367 = qJD(1) * t426;
t494 = t410 * t366 + t407 * t367;
t493 = -t519 * t487 - (-pkin(4) * t498 + qJ(4) * t408) * t491 - t494;
t362 = pkin(7) * t447 + qJD(2) * t395;
t467 = qJD(2) * t500;
t371 = pkin(7) * t467 + qJD(2) * t453;
t489 = qJD(2) * t410;
t479 = pkin(8) * t390 * t407;
t478 = pkin(8) * t503;
t477 = pkin(8) * t483;
t476 = pkin(8) * t489;
t392 = t519 * t410;
t474 = t409 * t500;
t468 = t405 * t483;
t464 = -qJ(4) * t407 - pkin(2);
t437 = t520 * t468;
t270 = -pkin(4) * t317 - qJD(1) * t437 + t454;
t418 = qJ(4) * t317 - qJD(4) * t352 + t362;
t271 = t318 * t520 + t418;
t462 = t409 * t270 - t271 * t406;
t461 = -t407 * t364 + t365 * t410;
t460 = -t407 * t366 + t367 * t410;
t459 = t406 * t345;
t458 = t409 * t345;
t456 = t402 * t408 * t411 * MDP(4);
t455 = t332 * t487 - t344 * t486 - t407 * t360 - t410 * t361;
t445 = pkin(1) * t527;
t303 = pkin(3) * t500 - t461;
t338 = t406 * t470 - t409 * t451;
t444 = t409 * t487 + t338;
t339 = (t406 * t498 + t408 * t409) * t491;
t443 = t406 * t487 - t339;
t375 = -t410 * t520 + t464;
t442 = qJD(5) * t375 + (pkin(4) * t497 - t408 * t520) * t491 - t460 - qJD(3) * t392;
t391 = t519 * t407;
t441 = -qJD(5) * t391 + t523 + t390 * (pkin(9) * t407 - t515);
t438 = -qJD(4) * t390 - t455;
t436 = t270 * t406 + t271 * t409;
t279 = t390 * t520 + t481;
t282 = t350 * t520 + t414;
t266 = t279 * t409 - t282 * t406;
t267 = t279 * t406 + t282 * t409;
t374 = t407 * t517 + t410 * t501;
t288 = pkin(4) * t374 + pkin(9) * t500 + t303;
t373 = -t410 * t517 + t475;
t363 = -pkin(2) * t517 - t524;
t415 = -t374 * qJ(4) + t363;
t292 = t373 * t520 + t415;
t434 = t288 * t409 - t292 * t406;
t433 = t288 * t406 + t292 * t409;
t432 = qJ(4) * t448;
t302 = qJ(4) * t500 - t495;
t428 = -t364 * t486 - t365 * t487 + t368 * t410 - t407 * t370;
t325 = t373 * t409 + t406 * t500;
t423 = -t300 * t390 - t454;
t274 = -t432 - t438;
t269 = -pkin(4) * t318 - t274;
t420 = t269 + (t345 * t520 + t516) * t345;
t419 = -t409 * t318 + t406 * t448;
t324 = -t450 + (qJD(3) * t517 + t467) * t410;
t417 = -qJ(4) * t324 - qJD(4) * t374 + t371;
t416 = -t317 - t506;
t385 = -pkin(3) * t410 + t464;
t326 = t373 * t406 - t474;
t323 = qJD(3) * t374 + t407 * t467;
t313 = t317 * t407;
t308 = pkin(3) * t352 + t516;
t307 = -pkin(3) * t470 - t460;
t305 = -qJ(4) * t470 - t494;
t304 = t317 * t374;
t301 = t373 * pkin(3) + t415;
t295 = pkin(3) * t390 - t482;
t293 = -pkin(4) * t373 - t302;
t290 = qJD(5) * t325 + t323 * t406 + t409 * t468;
t289 = -t323 * t409 - qJD(5) * t474 + (qJD(5) * t373 + t468) * t406;
t284 = t316 * qJD(5) + t419;
t281 = pkin(3) * t323 + t417;
t280 = -pkin(3) * t468 - t428;
t276 = pkin(3) * t318 + t418;
t275 = t323 * t520 + t417;
t273 = -pkin(4) * t323 - t278;
t272 = pkin(4) * t324 - t428 - t437;
t265 = -qJD(5) * t267 + t462;
t264 = qJD(5) * t266 + t436;
t1 = [t526 * t527 + (-t362 * t517 - t371 * t430 + t408 * t445) * MDP(9) + (-t361 * t517 - t370 * t430 + t411 * t445) * MDP(10) + (t324 * t352 - t304) * MDP(11) + (t317 * t373 - t318 * t374 - t323 * t352 - t324 * t350) * MDP(12) + (-t324 * t390 + (t317 * t411 + (qJD(1) * t374 + t352) * t483) * t405) * MDP(13) + (t323 * t390 + (t318 * t411 + (-qJD(1) * t373 - t350) * t483) * t405) * MDP(14) + (-t390 * t405 - t402 * t490) * MDP(15) * t483 + (-t428 * t390 + t371 * t350 + t363 * t318 + t362 * t373 + t331 * t323 + (t454 * t411 + (qJD(1) * t461 - t299) * t483) * t405) * MDP(16) + (-t422 * t390 + t371 * t352 - t363 * t317 + t362 * t374 + t331 * t324 + (-t455 * t411 + (-qJD(1) * t495 - t300) * t483) * t405) * MDP(17) + (t274 * t373 + t277 * t374 + t278 * t350 + t280 * t352 + t295 * t324 + t297 * t323 + t302 * t318 - t303 * t317) * MDP(18) + (-t276 * t373 - t280 * t390 - t281 * t350 - t291 * t323 - t301 * t318 + (-t277 * t411 + (qJD(1) * t303 + t295) * t483) * t405) * MDP(19) + (-t276 * t374 + t278 * t390 - t281 * t352 - t291 * t324 + t301 * t317 + (t274 * t411 + (-qJD(1) * t302 - t297) * t483) * t405) * MDP(20) + (t274 * t302 + t276 * t301 + t277 * t303 + t278 * t297 + t280 * t295 + t281 * t291) * MDP(21) + (t283 * t326 + t290 * t316) * MDP(22) + (t283 * t325 - t284 * t326 - t289 * t316 - t290 * t314) * MDP(23) + (t283 * t374 + t290 * t345 + t316 * t324 - t317 * t326) * MDP(24) + (-t284 * t374 - t289 * t345 - t314 * t324 - t317 * t325) * MDP(25) + (t324 * t345 - t304) * MDP(26) + ((-qJD(5) * t433 + t272 * t409 - t275 * t406) * t345 - t434 * t317 + t265 * t374 + t266 * t324 + t273 * t314 + t293 * t284 - t269 * t325 + t285 * t289) * MDP(27) + (-(qJD(5) * t434 + t272 * t406 + t275 * t409) * t345 + t433 * t317 - t264 * t374 - t267 * t324 + t273 * t316 + t293 * t283 + t269 * t326 + t285 * t290) * MDP(28) + 0.2e1 * t456 * t480 + (MDP(6) * t467 - MDP(7) * t468) * (0.2e1 * t463 + qJD(2)); t502 * t526 + (pkin(1) * t408 * t502 + t369 * t430 - t362) * MDP(9) + (pkin(7) * t448 + t366 * t430 + (-qJD(2) * t463 + t502) * t411 * pkin(1)) * MDP(10) + (-t352 * t503 - t313) * MDP(11) + ((-t317 + t506) * t410 + (-t318 + t505) * t407) * MDP(12) + (-t390 * t486 + (t390 * t497 + (qJD(2) * t407 - t352) * t408) * t491) * MDP(13) + (t390 * t487 + (-t390 * t498 + (t350 + t489) * t408) * t491) * MDP(14) + (-pkin(2) * t318 - t362 * t410 + t460 * t390 - t369 * t350 + (t331 * t407 + t478) * qJD(3) + (t299 * t408 + (-t331 * t411 - t477) * t407) * t491) * MDP(16) + (pkin(2) * t317 + t362 * t407 - t494 * t390 - t369 * t352 + (t331 * t410 - t479) * qJD(3) + (-t331 * t497 + (t300 - t476) * t408) * t491) * MDP(17) + (-t305 * t350 - t307 * t352 + (-t274 - t390 * t295 + (-t318 + t488) * pkin(8)) * t410 + ((qJD(3) * t350 - t317) * pkin(8) + t525) * t407) * MDP(18) + (t276 * t410 + t307 * t390 - t318 * t385 + t496 * t350 + (-t291 * t407 - t478) * qJD(3) + (-t295 * t408 + (t291 * t411 + t477) * t407) * t491) * MDP(19) + (-t276 * t407 - t305 * t390 + t317 * t385 + t496 * t352 + (-t291 * t410 + t479) * qJD(3) + (t291 * t497 + (t297 + t476) * t408) * t491) * MDP(20) + (t276 * t385 - t295 * t307 - t297 * t305 - t496 * t291 + (-t274 * t410 + t277 * t407 + (t295 * t410 + t297 * t407) * qJD(3)) * pkin(8)) * MDP(21) + (-t283 * t406 * t410 + (-t410 * t484 + t443) * t316) * MDP(22) + (t314 * t339 + t316 * t338 + (-t314 * t406 + t316 * t409) * t487 + (-t514 + t284 * t406 + (t314 * t409 + t316 * t406) * qJD(5)) * t410) * MDP(23) + (t283 * t407 + t443 * t345 + (-t345 * t484 + t499 - t508) * t410) * MDP(24) + (-t284 * t407 + t444 * t345 + (t345 * t485 + t311 + t510) * t410) * MDP(25) + (-t345 * t503 - t313) * MDP(26) + (-(-t375 * t406 + t391 * t409) * t317 + t265 * t407 + t392 * t284 + (t406 * t441 - t409 * t442) * t345 + t493 * t314 - t444 * t285 + (-t266 * t390 + t269 * t409 - t285 * t485) * t410) * MDP(27) + ((t375 * t409 + t391 * t406) * t317 - t264 * t407 + t392 * t283 + (t406 * t442 + t409 * t441) * t345 + t493 * t316 + t443 * t285 + (t267 * t390 - t269 * t406 - t285 * t484) * t410) * MDP(28) + t390 * MDP(15) * t470 + (-t456 + (-t411 * MDP(6) + t408 * MDP(7)) * t405 * t517) * t413; MDP(11) * t507 + (-t350 ^ 2 + t521) * MDP(12) + t416 * MDP(13) + (-t505 - t318) * MDP(14) + MDP(15) * t448 + (-t331 * t352 + t423) * MDP(16) + (t299 * t390 + t331 * t350 + t455) * MDP(17) + (pkin(3) * t317 - qJ(4) * t318 + (-t297 - t300) * t352 + (t295 + t482) * t350) * MDP(18) + (t308 * t350 - t423 - 0.2e1 * t439 + t513) * MDP(19) + (-t291 * t350 + t308 * t352 + t390 * t482 + 0.2e1 * t432 + t438) * MDP(20) + (-pkin(3) * t277 - qJ(4) * t274 - t291 * t308 - t295 * t300 + t297 * t482) * MDP(21) + (-t316 * t459 + t514) * MDP(22) + ((-t284 - t509) * t409 + (-t283 + t511) * t406) * MDP(23) + (t316 * t350 - t345 * t459 - t311) * MDP(24) + (-t314 * t350 - t345 * t458 + t499) * MDP(25) + t345 * t350 * MDP(26) + (qJ(4) * t284 + t266 * t350 + t481 * t314 + t420 * t406 + t409 * t522) * MDP(27) + (qJ(4) * t283 - t267 * t350 + t481 * t316 - t406 * t522 + t420 * t409) * MDP(28); t416 * MDP(18) + (t448 - t507) * MDP(19) + (-t390 ^ 2 - t521) * MDP(20) + (t513 + t525) * MDP(21) + (-t311 + t510) * MDP(27) + (t499 + t508) * MDP(28) + (-MDP(27) * t459 - MDP(28) * t458) * t345; t316 * t314 * MDP(22) + (-t314 ^ 2 + t316 ^ 2) * MDP(23) + (t473 + t511) * MDP(24) + (-t419 + t509) * MDP(25) - t317 * MDP(26) + (t267 * t345 - t285 * t316 + t462) * MDP(27) + (t266 * t345 + t285 * t314 - t436) * MDP(28) + (MDP(24) * t504 - MDP(25) * t316 - MDP(27) * t267 - MDP(28) * t266) * qJD(5);];
tauc = t1;
