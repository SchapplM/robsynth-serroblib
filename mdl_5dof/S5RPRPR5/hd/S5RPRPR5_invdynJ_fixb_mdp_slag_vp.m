% Calculate vector of inverse dynamics joint torques for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:17
% EndTime: 2020-01-03 11:43:24
% DurationCPUTime: 3.89s
% Computational Cost: add. (2387->353), mult. (5936->482), div. (0->0), fcn. (4375->12), ass. (0->194)
t411 = sin(pkin(8));
t406 = t411 ^ 2;
t521 = 0.2e1 * t406;
t410 = sin(pkin(9));
t412 = cos(pkin(9));
t416 = sin(qJ(3));
t419 = cos(qJ(3));
t368 = t410 * t419 + t412 * t416;
t433 = qJD(1) * t368;
t342 = t411 * t433;
t418 = cos(qJ(5));
t334 = t418 * t342;
t478 = qJD(1) * t419;
t458 = t411 * t478;
t480 = qJD(1) * t411;
t459 = t416 * t480;
t345 = -t410 * t459 + t412 * t458;
t415 = sin(qJ(5));
t507 = t345 * t415;
t302 = t334 + t507;
t413 = cos(pkin(8));
t479 = qJD(1) * t413;
t388 = -qJD(3) + t479;
t382 = -qJD(5) + t388;
t508 = t302 * t382;
t510 = qJDD(1) * pkin(1);
t397 = qJDD(2) - t510;
t417 = sin(qJ(1));
t420 = cos(qJ(1));
t524 = -g(2) * t420 - g(3) * t417;
t529 = t524 - t397;
t371 = -pkin(2) * t413 - pkin(6) * t411 - pkin(1);
t355 = qJDD(1) * t371 + qJDD(2);
t351 = t419 * t355;
t466 = qJDD(1) * t413;
t386 = -qJDD(3) + t466;
t511 = qJ(2) * t419;
t387 = t413 * t511;
t474 = qJD(4) * t411;
t477 = qJD(2) * t413;
t429 = -t416 * t477 - t419 * t474;
t502 = t411 * t419;
t461 = qJ(4) * t502;
t512 = qJ(2) * t416;
t463 = t413 * t512;
t430 = -t461 - t463;
t503 = t411 * t416;
t462 = qJ(4) * t503;
t356 = qJD(1) * t371 + qJD(2);
t476 = qJD(3) * t356;
t280 = -t416 * t476 - pkin(3) * t386 + t351 + t430 * qJDD(1) + ((-t387 + t462) * qJD(3) + t429) * qJD(1);
t423 = qJD(3) * t430 - t416 * t474;
t470 = qJD(1) * qJD(2);
t457 = t419 * t470;
t475 = qJD(3) * t419;
t449 = qJDD(1) * t387 + t416 * t355 + t356 * t475 + t413 * t457;
t465 = qJDD(1) * t416;
t455 = t411 * t465;
t286 = -qJ(4) * t455 + qJD(1) * t423 + t449;
t269 = t412 * t280 - t286 * t410;
t432 = t368 * qJD(3);
t464 = qJDD(1) * t419;
t315 = (qJD(1) * t432 + t410 * t465 - t412 * t464) * t411;
t267 = -pkin(4) * t386 + pkin(7) * t315 + t269;
t270 = t410 * t280 + t412 * t286;
t439 = t410 * t416 - t412 * t419;
t431 = t439 * qJD(3);
t314 = (qJD(1) * t431 - qJDD(1) * t368) * t411;
t268 = pkin(7) * t314 + t270;
t352 = t419 * t356;
t323 = qJD(1) * t430 + t352;
t313 = -pkin(3) * t388 + t323;
t324 = -qJ(4) * t459 + qJD(1) * t387 + t356 * t416;
t501 = t412 * t324;
t291 = t410 * t313 + t501;
t519 = pkin(7) * t342;
t277 = t291 - t519;
t473 = qJD(5) * t415;
t276 = t277 * t473;
t357 = pkin(3) * t459 + qJ(2) * t480 + qJD(4);
t322 = pkin(4) * t342 + t357;
t400 = qJ(3) + pkin(9) + qJ(5);
t394 = cos(t400);
t393 = sin(t400);
t494 = t420 * t393;
t500 = t413 * t417;
t338 = t394 * t500 - t494;
t506 = t394 * t420;
t340 = t393 * t417 + t413 * t506;
t517 = g(1) * t411;
t528 = g(2) * t338 - g(3) * t340 - t415 * t267 - t418 * t268 + t322 * t302 + t394 * t517 + t276;
t377 = -qJDD(5) + t386;
t369 = t377 * MDP(21);
t441 = -t342 * t415 + t418 * t345;
t527 = t302 * MDP(17) * t441 + (-t302 ^ 2 + t441 ^ 2) * MDP(18) - t369;
t468 = qJDD(1) * qJ(2);
t436 = t468 + t470;
t525 = t413 * t436;
t509 = t441 * t382;
t495 = t419 * t420;
t499 = t416 * t417;
t360 = -t413 * t499 - t495;
t497 = t417 * t419;
t498 = t416 * t420;
t362 = t413 * t498 - t497;
t523 = -g(2) * t360 - g(3) * t362;
t337 = -t393 * t500 - t506;
t339 = -t394 * t417 + t413 * t494;
t453 = t418 * t267 - t415 * t268;
t522 = -g(2) * t337 - g(3) * t339 - t322 * t441 + t393 * t517 + t453;
t452 = t418 * t314 + t315 * t415;
t272 = qJD(5) * t441 - t452;
t407 = t413 ^ 2;
t520 = pkin(3) * t410;
t518 = pkin(7) * t345;
t515 = g(2) * t417;
t513 = g(3) * t420;
t421 = qJD(1) ^ 2;
t505 = t406 * t421;
t318 = t410 * t324;
t504 = t411 * (-qJ(4) - pkin(6));
t290 = t412 * t313 - t318;
t275 = -pkin(4) * t388 + t290 - t518;
t496 = t418 * t275;
t491 = t371 * t475 + t419 * t477;
t311 = t423 + t491;
t312 = (-t387 + (qJ(4) * t411 - t371) * t416) * qJD(3) + t429;
t284 = t412 * t311 + t410 * t312;
t294 = t412 * t323 - t318;
t366 = t419 * t371;
t328 = -t461 + t366 + (-pkin(3) - t512) * t413;
t490 = t416 * t371 + t387;
t332 = -t462 + t490;
t297 = t410 * t328 + t412 * t332;
t493 = -t413 * t433 + t432;
t492 = t439 * t479 - t431;
t384 = t411 * pkin(3) * t475;
t489 = t411 * qJD(2) + t384;
t488 = pkin(3) * t503 + t411 * qJ(2);
t487 = t420 * pkin(1) + t417 * qJ(2);
t485 = t406 + t407;
t409 = t419 ^ 2;
t484 = t416 ^ 2 - t409;
t483 = MDP(10) * t411;
t482 = MDP(11) * t411;
t472 = t386 * MDP(12);
t471 = qJD(3) + t388;
t469 = qJD(1) * qJD(3);
t467 = qJDD(1) * t411;
t460 = -qJD(5) * t334 + t415 * t314 - t418 * t315;
t456 = t416 * t469;
t454 = t485 * t421;
t283 = -t311 * t410 + t412 * t312;
t293 = -t323 * t410 - t501;
t296 = t412 * t328 - t332 * t410;
t451 = qJD(1) * t471;
t450 = t386 + t466;
t448 = 0.2e1 * t485;
t447 = qJD(3) * t463;
t445 = -t513 + t515;
t444 = qJD(5) * t368 + t493;
t443 = -qJD(5) * t439 + t492;
t442 = -t415 * t275 - t418 * t277;
t353 = t368 * t411;
t354 = t439 * t411;
t440 = -t418 * t353 + t354 * t415;
t317 = -t353 * t415 - t354 * t418;
t438 = pkin(3) * t455 + qJ(2) * t467 + qJD(1) * t384 + t411 * t470 + qJDD(4);
t437 = qJD(3) * (t388 + t479);
t395 = pkin(3) * t412 + pkin(4);
t435 = t395 * t415 + t418 * t520;
t434 = t395 * t418 - t415 * t520;
t271 = -t345 * t473 + t460;
t427 = -t388 ^ 2 - t505;
t425 = t448 * t470 + t513;
t402 = t417 * pkin(1);
t396 = pkin(3) * t419 + pkin(2);
t363 = t413 * t495 + t499;
t361 = t413 * t497 - t498;
t348 = t411 * t432;
t344 = t411 * t431;
t331 = pkin(3) * t458 + pkin(4) * t345;
t329 = pkin(4) * t353 + t488;
t325 = -pkin(4) * t344 + t489;
t295 = -pkin(4) * t314 + t438;
t292 = -pkin(7) * t353 + t297;
t289 = -pkin(4) * t413 + pkin(7) * t354 + t296;
t288 = qJD(5) * t317 - t418 * t344 - t348 * t415;
t287 = qJD(5) * t440 + t344 * t415 - t348 * t418;
t282 = t294 - t518;
t281 = t293 + t519;
t274 = pkin(7) * t344 + t284;
t273 = pkin(7) * t348 + t283;
t1 = [qJDD(1) * MDP(1) + t524 * MDP(2) + t445 * MDP(3) + (t448 * t468 + t425 - t515) * MDP(6) + (-t397 * pkin(1) - g(2) * t487 - g(3) * t402 + (t485 * t468 + t425) * qJ(2)) * MDP(7) + (qJDD(1) * t409 - 0.2e1 * t419 * t456) * t406 * MDP(8) + (-t416 * t464 + t469 * t484) * MDP(9) * t521 + (t416 * t437 - t419 * t450) * t483 + (t416 * t450 + t419 * t437) * t482 + t413 * t472 + (-g(2) * t363 - g(3) * t361 - t351 * t413 - t366 * t386 + (t388 * t413 + (t521 + t407) * qJD(1)) * qJ(2) * t475 + (qJD(3) * t371 * t388 + t436 * t521 + (qJ(2) * t386 + qJD(2) * t388 + t476 + t525) * t413) * t416) * MDP(13) + ((-t447 + t491) * t388 + t490 * t386 + (-qJD(1) * t447 + t449) * t413 + g(2) * t362 - g(3) * t360 + (t457 + (-t456 + t464) * qJ(2)) * t521) * MDP(14) + (t269 * t354 - t270 * t353 - t283 * t345 - t284 * t342 + t290 * t348 + t291 * t344 + t296 * t315 + t297 * t314 + t411 * t524) * MDP(15) + (t270 * t297 + t291 * t284 + t269 * t296 + t290 * t283 + t438 * t488 + t357 * t489 - g(2) * (pkin(3) * t499 + t487) - g(3) * (t396 * t500 - t417 * t504 + t402) + (-g(2) * (t396 * t413 - t504) - g(3) * (-pkin(3) * t416 - qJ(2))) * t420) * MDP(16) + (t271 * t317 + t287 * t441) * MDP(17) + (t271 * t440 - t272 * t317 - t287 * t302 - t288 * t441) * MDP(18) + (-t271 * t413 - t287 * t382 - t317 * t377) * MDP(19) + (t272 * t413 + t288 * t382 - t377 * t440) * MDP(20) + t413 * t369 + (-(t273 * t418 - t274 * t415) * t382 - (t289 * t418 - t292 * t415) * t377 - t453 * t413 + t325 * t302 + t329 * t272 - t295 * t440 + t322 * t288 - g(2) * t340 - g(3) * t338 + (-(-t289 * t415 - t292 * t418) * t382 - t442 * t413) * qJD(5)) * MDP(22) + (g(2) * t339 - g(3) * t337 + t329 * t271 - t276 * t413 + t322 * t287 + t295 * t317 + t325 * t441 + ((-qJD(5) * t292 + t273) * t382 + t289 * t377 + t267 * t413) * t415 + ((qJD(5) * t289 + t274) * t382 + t292 * t377 + (qJD(5) * t275 + t268) * t413) * t418) * MDP(23) + (MDP(4) * t413 - MDP(5) * t411) * (t510 + t529); -MDP(4) * t466 + MDP(5) * t467 - MDP(6) * t454 + (-qJ(2) * t454 - t529) * MDP(7) + (-t386 * t419 + t416 * t427) * MDP(13) + (t386 * t416 + t419 * t427) * MDP(14) + (t314 * t368 - t315 * t439 - t342 * t492 + t345 * t493) * MDP(15) + (-t269 * t439 + t270 * t368 - t290 * t493 + t291 * t492 - t357 * t480 - t524) * MDP(16) + (-(-t368 * t415 - t418 * t439) * t377 - t302 * t480 + (t415 * t443 + t418 * t444) * t382) * MDP(22) + ((t368 * t418 - t415 * t439) * t377 - t441 * t480 + (-t415 * t444 + t418 * t443) * t382) * MDP(23); t419 * t416 * MDP(8) * t505 - t484 * MDP(9) * t505 + (-t416 * t451 + t464) * t483 + (-t419 * t451 - t465) * t482 - t472 + (t351 + (-t413 * t451 - t505) * t511 + (-t356 * t471 + t517 - t525) * t416 + t523) * MDP(13) + (g(1) * t502 + g(2) * t361 - g(3) * t363 - t352 * t388 + (t471 * t479 + t505) * t512 - t449) * MDP(14) + ((t291 + t293) * t345 - (t290 - t294) * t342 + (t314 * t410 + t315 * t412) * pkin(3)) * MDP(15) + (-t290 * t293 - t291 * t294 + (t270 * t410 + t269 * t412 + (g(1) * t416 - t357 * t478) * t411 + t523) * pkin(3)) * MDP(16) + (t271 - t508) * MDP(19) + (-t272 - t509) * MDP(20) + (-t434 * t377 + (t281 * t418 - t282 * t415) * t382 - t331 * t302 + (t382 * t435 + t442) * qJD(5) + t522) * MDP(22) + (t435 * t377 - (t281 * t415 + t282 * t418) * t382 - t331 * t441 + (t382 * t434 - t496) * qJD(5) + t528) * MDP(23) + t527; (-t342 ^ 2 - t345 ^ 2) * MDP(15) + (g(1) * t413 + t290 * t345 + t291 * t342 - t411 * t445 + t438) * MDP(16) + (t272 - t509) * MDP(22) + (t271 + t508) * MDP(23); (t460 - t508) * MDP(19) + (t452 - t509) * MDP(20) + (t382 * t442 + t522) * MDP(22) + (-(-t277 * t415 + t496) * t382 + t528) * MDP(23) + (-MDP(19) * t507 - t441 * MDP(20) + t442 * MDP(22) - MDP(23) * t496) * qJD(5) + t527;];
tau = t1;
