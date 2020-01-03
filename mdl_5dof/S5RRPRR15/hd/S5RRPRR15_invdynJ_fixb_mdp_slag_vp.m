% Calculate vector of inverse dynamics joint torques for
% S5RRPRR15
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR15_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR15_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR15_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:19
% EndTime: 2019-12-31 20:43:26
% DurationCPUTime: 5.16s
% Computational Cost: add. (2352->439), mult. (4912->577), div. (0->0), fcn. (3167->10), ass. (0->212)
t419 = sin(qJ(2));
t484 = qJD(1) * qJD(2);
t472 = t419 * t484;
t423 = cos(qJ(2));
t482 = qJDD(1) * t423;
t551 = t472 - t482;
t538 = pkin(3) + pkin(6);
t422 = cos(qJ(4));
t418 = sin(qJ(4));
t497 = qJD(2) * t418;
t500 = qJD(1) * t423;
t364 = t422 * t500 + t497;
t421 = cos(qJ(5));
t473 = t418 * t500;
t495 = qJD(2) * t422;
t366 = -t473 + t495;
t417 = sin(qJ(5));
t522 = t366 * t417;
t311 = t421 * t364 + t522;
t501 = qJD(1) * t419;
t395 = qJD(4) + t501;
t388 = qJD(5) + t395;
t550 = t311 * t388;
t453 = t364 * t417 - t421 * t366;
t549 = t388 * t453;
t367 = t417 * t422 + t418 * t421;
t443 = t367 * t419;
t541 = qJD(4) + qJD(5);
t508 = -qJD(1) * t443 - t541 * t367;
t530 = qJ(3) * t419;
t457 = pkin(2) * t423 + t530;
t448 = pkin(1) + t457;
t349 = t448 * qJD(1);
t548 = t448 * qJDD(1);
t425 = -pkin(2) - pkin(7);
t359 = t425 * t423 - pkin(1) - t530;
t325 = t359 * qJD(1);
t400 = pkin(6) * t501;
t542 = qJD(3) + t400;
t485 = pkin(3) * t501 + t542;
t331 = t425 * qJD(2) + t485;
t295 = t325 * t422 + t331 * t418;
t288 = -pkin(8) * t364 + t295;
t487 = qJD(5) * t417;
t286 = t288 * t487;
t401 = pkin(6) * t500;
t374 = pkin(3) * t500 + t401;
t413 = qJD(2) * qJ(3);
t348 = t413 + t374;
t317 = pkin(4) * t364 + t348;
t416 = qJ(4) + qJ(5);
t406 = sin(t416);
t407 = cos(t416);
t420 = sin(qJ(1));
t424 = cos(qJ(1));
t515 = t419 * t424;
t333 = t406 * t515 + t407 * t420;
t517 = t419 * t420;
t335 = -t406 * t517 + t407 * t424;
t410 = g(3) * t423;
t547 = g(1) * t333 - g(2) * t335 + t311 * t317 - t406 * t410 + t286;
t332 = -t406 * t420 + t407 * t515;
t334 = t406 * t424 + t407 * t517;
t305 = -t364 * qJD(4) + t422 * qJDD(2) + t551 * t418;
t471 = t423 * t484;
t483 = qJDD(1) * t419;
t441 = t471 + t483;
t363 = qJDD(4) + t441;
t394 = pkin(2) * t472;
t529 = qJ(3) * t423;
t455 = pkin(7) * t419 - t529;
t493 = qJD(3) * t419;
t432 = t455 * qJD(2) - t493;
t293 = t432 * qJD(1) + t359 * qJDD(1) + t394;
t397 = pkin(6) * t483;
t481 = qJDD(3) + t397;
t470 = pkin(6) * t471 + t481;
t308 = t441 * pkin(3) + t425 * qJDD(2) + t470;
t466 = -t293 * t418 + t422 * t308;
t430 = -t295 * qJD(4) + t466;
t276 = pkin(4) * t363 - pkin(8) * t305 + t430;
t306 = -qJD(4) * t473 + qJDD(2) * t418 + (qJD(2) * qJD(4) - t551) * t422;
t490 = qJD(4) * t422;
t480 = -t422 * t293 - t418 * t308 - t331 * t490;
t491 = qJD(4) * t418;
t277 = -pkin(8) * t306 - t325 * t491 - t480;
t467 = t421 * t276 - t417 * t277;
t546 = -g(1) * t332 - g(2) * t334 + t317 * t453 + t407 * t410 + t467;
t358 = qJDD(5) + t363;
t545 = t358 * MDP(26) + (-t311 ^ 2 + t453 ^ 2) * MDP(23) - t311 * MDP(22) * t453;
t544 = -t349 * t501 + t410;
t382 = t538 * t419;
t369 = t418 * t382;
t505 = t422 * t359 + t369;
t344 = t422 * t363;
t543 = -t395 * t491 + t344;
t461 = g(1) * t420 - g(2) * t424;
t462 = g(1) * t424 + g(2) * t420;
t540 = t395 * t348 + t425 * t363;
t465 = t305 * t417 + t421 * t306;
t281 = -t453 * qJD(5) + t465;
t539 = t541 * t423;
t533 = g(3) * t419;
t532 = pkin(8) - t425;
t531 = pkin(6) * qJDD(2);
t528 = qJDD(2) * pkin(2);
t294 = -t325 * t418 + t422 * t331;
t287 = -pkin(8) * t366 + t294;
t285 = pkin(4) * t395 + t287;
t527 = t285 * t421;
t526 = t288 * t421;
t525 = t305 * t422;
t524 = t364 * t395;
t523 = t366 * t395;
t521 = t417 * t418;
t520 = t418 * t363;
t519 = t418 * t419;
t518 = t418 * t423;
t516 = t419 * t422;
t514 = t420 * t422;
t513 = t421 * t422;
t512 = t422 * t423;
t511 = t422 * t424;
t427 = qJD(1) ^ 2;
t510 = t423 * t427;
t476 = t422 * t501;
t507 = -t417 * t491 - t418 * t487 + t421 * t476 - t501 * t521 + t541 * t513;
t404 = pkin(2) * t501;
t337 = t455 * qJD(1) + t404;
t506 = t422 * t337 + t418 * t374;
t477 = -pkin(4) * t422 - pkin(3);
t504 = pkin(4) * t490 - t477 * t501 + t542;
t383 = t538 * t423;
t414 = t419 ^ 2;
t415 = t423 ^ 2;
t503 = t414 - t415;
t499 = qJD(2) * t364;
t498 = qJD(2) * t366;
t496 = qJD(2) * t419;
t494 = qJD(2) * t423;
t492 = qJD(4) * t325;
t489 = qJD(4) * t423;
t488 = qJD(4) * t425;
t486 = qJD(5) * t421;
t479 = t421 * t305 - t417 * t306 - t364 * t486;
t398 = pkin(6) * t482;
t411 = qJDD(2) * qJ(3);
t412 = qJD(2) * qJD(3);
t478 = t398 + t411 + t412;
t474 = t418 * t489;
t378 = t532 * t422;
t469 = pkin(8) * t423 - t359;
t468 = -qJD(2) * pkin(2) + qJD(3);
t403 = pkin(2) * t496;
t322 = t403 + t432;
t375 = t538 * t494;
t464 = -t322 * t418 + t422 * t375;
t463 = qJD(5) * t285 + t277;
t373 = t538 * t496;
t451 = -t513 + t521;
t460 = -t451 * t358 + t388 * t508;
t356 = t422 * t374;
t377 = t532 * t418;
t449 = pkin(4) * t423 - pkin(8) * t519;
t459 = t449 * qJD(1) - qJD(5) * t377 - t337 * t418 - t532 * t491 + t356;
t458 = pkin(8) * t476 + t378 * t541 + t506;
t456 = pkin(2) * t419 - t529;
t279 = t285 * t417 + t526;
t376 = t400 + t468;
t381 = -t401 - t413;
t452 = t376 * t423 + t381 * t419;
t450 = t395 * t418;
t336 = t470 - t528;
t446 = -0.2e1 * pkin(1) * t484 - t531;
t445 = -t395 * t490 - t520;
t444 = -qJ(3) * t494 - t493;
t442 = -t358 * t367 - t507 * t388;
t440 = t422 * t322 - t359 * t491 + t418 * t375 + t382 * t490;
t280 = -t366 * t487 + t479;
t439 = pkin(1) * t427 + t462;
t426 = qJD(2) ^ 2;
t438 = pkin(6) * t426 - t461;
t437 = 0.2e1 * qJD(2) * t349 + t531;
t435 = -t462 * t423 - t533;
t433 = 0.2e1 * qJDD(1) * pkin(1) - t438;
t309 = pkin(3) * t482 - qJD(1) * t373 + t478;
t431 = t309 + t435;
t307 = t444 * qJD(1) + t394 - t548;
t342 = t403 + t444;
t429 = qJD(1) * t342 + t307 + t438 - t548;
t328 = pkin(6) * t472 - t478;
t428 = t452 * qJD(2) - t328 * t423 + t336 * t419 - t462;
t396 = pkin(4) * t418 + qJ(3);
t371 = -qJ(3) * t500 + t404;
t370 = t422 * t382;
t353 = -t418 * t517 + t511;
t352 = t418 * t424 + t419 * t514;
t351 = t418 * t515 + t514;
t350 = -t418 * t420 + t419 * t511;
t347 = pkin(4) * t512 + t383;
t339 = t367 * t423;
t338 = t451 * t423;
t318 = -pkin(4) * t474 + (-pkin(6) + t477) * t496;
t302 = -pkin(8) * t512 + t505;
t297 = pkin(4) * t419 + t469 * t418 + t370;
t290 = t367 * t539 - t451 * t496;
t289 = qJD(2) * t443 + t451 * t539;
t284 = pkin(4) * t306 + t309;
t283 = (t419 * t495 + t474) * pkin(8) + t440;
t282 = t449 * qJD(2) + (t469 * t422 - t369) * qJD(4) + t464;
t278 = -t288 * t417 + t527;
t1 = [(t280 * t338 + t281 * t339 - t289 * t311 - t290 * t453) * MDP(23) + (-t280 * t339 - t289 * t453) * MDP(22) + (t280 * t419 + t289 * t388 - t339 * t358 - t453 * t494) * MDP(24) + (-t279 * t494 + g(1) * t334 - g(2) * t332 + t347 * t280 - t284 * t339 + t286 * t419 + t317 * t289 - t318 * t453 + (-(-qJD(5) * t302 + t282) * t388 - t297 * t358 - t276 * t419) * t417 + (-(qJD(5) * t297 + t283) * t388 - t302 * t358 - t463 * t419) * t421) * MDP(28) + (t446 * t419 + t433 * t423) * MDP(9) + (-t433 * t419 + t446 * t423) * MDP(10) + (t428 * pkin(6) - t349 * t342 + (-t307 + t461) * t448) * MDP(14) + (-t429 * t419 + t437 * t423) * MDP(13) + (t437 * t419 + t429 * t423) * MDP(12) + ((t414 + t415) * qJDD(1) * pkin(6) + t428) * MDP(11) + (-t440 * t395 - t505 * t363 - t373 * t366 + t383 * t305 + g(1) * t352 - g(2) * t350 + ((qJD(2) * t348 + t492) * t418 + t480) * t419 + (-t295 * qJD(2) - t309 * t418 - t348 * t490) * t423) * MDP(21) + qJDD(1) * MDP(1) + (t464 * t395 + (-t359 * t418 + t370) * t363 + t466 * t419 - t373 * t364 + t383 * t306 + t309 * t512 - g(1) * t353 - g(2) * t351 + (t294 * t423 - t348 * t516) * qJD(2) + (-t295 * t419 - t348 * t518 - t505 * t395) * qJD(4)) * MDP(20) + ((t395 * t495 - t306) * t419 + (-t499 - t543) * t423) * MDP(18) + ((-t364 * t418 + t366 * t422) * t496 + (-t525 + t306 * t418 + (t364 * t422 + t366 * t418) * qJD(4)) * t423) * MDP(16) + (-t305 * t518 + (t418 * t496 - t422 * t489) * t366) * MDP(15) + 0.2e1 * (t419 * t482 - t503 * t484) * MDP(5) + (t363 * t419 + t395 * t494) * MDP(19) + (t358 * t419 + t388 * t494) * MDP(26) + ((t395 * t497 + t305) * t419 + (t445 + t498) * t423) * MDP(17) + (-t281 * t419 + t290 * t388 - t311 * t494 + t338 * t358) * MDP(25) + ((t282 * t421 - t283 * t417) * t388 + (t297 * t421 - t302 * t417) * t358 + t467 * t419 + t278 * t494 + t318 * t311 + t347 * t281 - t284 * t338 - t317 * t290 - g(1) * t335 - g(2) * t333 + ((-t297 * t417 - t302 * t421) * t388 - t279 * t419) * qJD(5)) * MDP(27) + (qJDD(1) * t414 + 0.2e1 * t419 * t471) * MDP(4) + t461 * MDP(2) + t462 * MDP(3) + (qJDD(2) * t419 + t423 * t426) * MDP(6) + (qJDD(2) * t423 - t419 * t426) * MDP(7); (t439 * t419 - t397 - t410) * MDP(9) + t503 * MDP(5) * t427 - t419 * MDP(4) * t510 + MDP(7) * t482 + MDP(6) * t483 + qJDD(2) * MDP(8) + (-t452 * qJD(1) * pkin(6) - t336 * pkin(2) - g(3) * t457 - t328 * qJ(3) - t381 * qJD(3) + t349 * t371 + t462 * t456) * MDP(14) + (t439 * t423 - t398 + t533) * MDP(10) + (-t366 * t450 + t525) * MDP(15) + (-t462 * t419 + t481 - 0.2e1 * t528 + t544) * MDP(12) + ((-t306 - t523) * t422 + (-t305 + t524) * t418) * MDP(16) + ((-t366 * t423 - t395 * t519) * qJD(1) + t543) * MDP(17) + ((t364 * t423 - t395 * t516) * qJD(1) + t445) * MDP(18) + ((t377 * t417 - t378 * t421) * t358 + t396 * t281 + t284 * t367 + (t458 * t417 - t459 * t421) * t388 + t507 * t317 + t504 * t311 + t435 * t406) * MDP(27) + (-t280 * t451 - t453 * t508) * MDP(22) + (-(-t377 * t421 - t378 * t417) * t358 + t396 * t280 - t284 * t451 + (t459 * t417 + t458 * t421) * t388 + t508 * t317 - t504 * t453 + t435 * t407) * MDP(28) + (-t280 * t367 + t281 * t451 - t508 * t311 + t453 * t507) * MDP(23) + (qJ(3) * t305 + t506 * t395 + t485 * t366 - t540 * t418 + (-t395 * t488 + t431) * t422) * MDP(21) + (qJ(3) * t306 - t356 * t395 + t485 * t364 + t540 * t422 + ((t337 - t488) * t395 + t431) * t418) * MDP(20) + (-t456 * qJDD(1) + ((-t381 - t413) * t419 + (-t376 + t468) * t423) * qJD(1)) * MDP(11) + t460 * MDP(24) + t442 * MDP(25) + (t398 + 0.2e1 * t411 + 0.2e1 * t412 + (qJD(1) * t371 - g(3)) * t419 + (-qJD(1) * t349 - t462) * t423) * MDP(13) + (-t371 * MDP(12) - t395 * MDP(19) - t294 * MDP(20) + t295 * MDP(21) + MDP(24) * t453 + t311 * MDP(25) - t388 * MDP(26) - t278 * MDP(27) + t279 * MDP(28)) * t500; qJDD(2) * MDP(12) + (-t414 * t427 - t426) * MDP(13) + (qJD(2) * t381 + t336 + t544) * MDP(14) + (t344 - t499) * MDP(20) + (-t498 - t520) * MDP(21) + (-qJD(2) * t311 + t460) * MDP(27) + (qJD(2) * t453 + t442) * MDP(28) + (qJDD(1) * MDP(11) + MDP(12) * t510 - t462 * MDP(14)) * t419 + (-t395 * MDP(21) * t422 - MDP(20) * t450) * t395; t366 * t364 * MDP(15) + (-t364 ^ 2 + t366 ^ 2) * MDP(16) + (t305 + t524) * MDP(17) + (-t306 + t523) * MDP(18) + t363 * MDP(19) + (-g(1) * t350 - g(2) * t352 + g(3) * t512 + t295 * t395 - t348 * t366 + t430) * MDP(20) + (g(1) * t351 - g(2) * t353 + t294 * t395 + t348 * t364 + (t492 - t410) * t418 + t480) * MDP(21) + (t280 + t550) * MDP(24) + (-t281 - t549) * MDP(25) + (-(-t287 * t417 - t526) * t388 - t279 * qJD(5) + (-t311 * t366 + t421 * t358 - t388 * t487) * pkin(4) + t546) * MDP(27) + ((-t288 * t388 - t276) * t417 + (t287 * t388 - t463) * t421 + (-t417 * t358 + t366 * t453 - t388 * t486) * pkin(4) + t547) * MDP(28) + t545; (t479 + t550) * MDP(24) + (-t465 - t549) * MDP(25) + (t279 * t388 + t546) * MDP(27) + (-t417 * t276 - t421 * t277 + t278 * t388 + t547) * MDP(28) + (-MDP(24) * t522 + t453 * MDP(25) - t279 * MDP(27) - MDP(28) * t527) * qJD(5) + t545;];
tau = t1;
