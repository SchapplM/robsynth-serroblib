% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:24:01
% EndTime: 2019-03-09 02:24:06
% DurationCPUTime: 4.77s
% Computational Cost: add. (2719->447), mult. (5121->589), div. (0->0), fcn. (3419->14), ass. (0->201)
t417 = cos(pkin(10));
t394 = -pkin(1) * t417 - pkin(2);
t390 = -pkin(7) + t394;
t360 = qJDD(1) * t390 + qJDD(3);
t424 = cos(qJ(4));
t363 = qJD(1) * t390 + qJD(3);
t420 = sin(qJ(4));
t485 = qJD(2) * qJD(4);
t496 = qJD(4) * t420;
t473 = -t420 * qJDD(2) - t363 * t496 - t424 * t485;
t298 = -qJDD(4) * pkin(4) - t360 * t424 - t473;
t500 = qJD(1) * t420;
t391 = qJD(5) + t500;
t410 = qJ(1) + pkin(10);
t400 = sin(t410);
t401 = cos(t410);
t453 = g(1) * t400 - g(2) * t401;
t535 = g(3) * t420;
t432 = t424 * t453 - t535;
t552 = qJD(5) * pkin(8) * t391 + t298 + t432;
t416 = sin(pkin(10));
t392 = pkin(1) * t416 + qJ(3);
t504 = qJDD(1) * t392;
t419 = sin(qJ(5));
t423 = cos(qJ(5));
t487 = t423 * qJD(4);
t499 = qJD(1) * t424;
t368 = t419 * t499 - t487;
t422 = cos(qJ(6));
t370 = qJD(4) * t419 + t423 * t499;
t418 = sin(qJ(6));
t521 = t370 * t418;
t315 = t422 * t368 + t521;
t388 = qJD(6) + t391;
t551 = t315 * t388;
t447 = t368 * t418 - t422 * t370;
t550 = t388 * t447;
t332 = -t420 * qJD(2) + t363 * t424;
t549 = qJD(4) * t332;
t333 = t424 * qJD(2) + t420 * t363;
t326 = qJD(4) * pkin(8) + t333;
t454 = pkin(4) * t420 - pkin(8) * t424;
t359 = t392 + t454;
t334 = t359 * qJD(1);
t296 = t326 * t423 + t334 * t419;
t290 = -pkin(9) * t368 + t296;
t490 = qJD(6) * t418;
t287 = t290 * t490;
t325 = -qJD(4) * pkin(4) - t332;
t303 = pkin(5) * t368 + t325;
t415 = qJ(5) + qJ(6);
t408 = sin(t415);
t409 = cos(t415);
t517 = t409 * t420;
t329 = t400 * t517 + t401 * t408;
t331 = -t400 * t408 + t401 * t517;
t534 = g(3) * t424;
t548 = g(1) * t329 - g(2) * t331 + t303 * t315 + t409 * t534 + t287;
t518 = t408 * t420;
t328 = -t400 * t518 + t401 * t409;
t330 = t400 * t409 + t401 * t518;
t297 = qJDD(4) * pkin(8) + qJDD(2) * t424 + t360 * t420 + t549;
t455 = pkin(4) * t424 + pkin(8) * t420;
t366 = qJD(4) * t455 + qJD(3);
t308 = qJD(1) * t366 + qJDD(1) * t454 + t504;
t302 = t423 * t308;
t469 = t420 * t487;
t491 = qJD(5) * t424;
t437 = -t419 * t491 - t469;
t482 = qJDD(1) * t424;
t311 = qJD(1) * t437 + qJD(5) * t487 + t419 * qJDD(4) + t423 * t482;
t486 = qJD(1) * qJD(4);
t465 = t424 * t486;
t483 = qJDD(1) * t420;
t365 = qJDD(5) + t465 + t483;
t276 = pkin(5) * t365 - pkin(9) * t311 - qJD(5) * t296 - t297 * t419 + t302;
t470 = t419 * t496;
t312 = -qJD(1) * t470 + qJD(5) * t370 - t423 * qJDD(4) + t419 * t482;
t492 = qJD(5) * t423;
t476 = -t423 * t297 - t419 * t308 - t334 * t492;
t494 = qJD(5) * t419;
t439 = -t326 * t494 - t476;
t277 = -pkin(9) * t312 + t439;
t463 = t422 * t276 - t418 * t277;
t547 = -g(1) * t328 - g(2) * t330 + t303 * t447 + t408 * t534 + t463;
t358 = qJDD(6) + t365;
t546 = t358 * MDP(26) + (-t315 ^ 2 + t447 ^ 2) * MDP(23) - t315 * MDP(22) * t447;
t372 = t418 * t423 + t419 * t422;
t345 = t372 * t424;
t378 = qJD(1) * t392;
t544 = -qJD(1) * t378 - t453;
t543 = -qJD(6) * t423 - t492;
t542 = qJDD(1) * t394;
t541 = qJD(5) + qJD(6);
t540 = 0.2e1 * qJD(4) * t378 + qJDD(4) * t390;
t462 = t311 * t418 + t422 * t312;
t282 = -qJD(6) * t447 + t462;
t539 = pkin(8) + pkin(9);
t537 = pkin(9) * t424;
t533 = t282 * t420;
t295 = -t326 * t419 + t423 * t334;
t289 = -pkin(9) * t370 + t295;
t286 = pkin(5) * t391 + t289;
t532 = t286 * t422;
t531 = t290 * t422;
t530 = t311 * t419;
t529 = t312 * t420;
t371 = t418 * t419 - t422 * t423;
t528 = t358 * t371;
t527 = t358 * t372;
t525 = t365 * t419;
t524 = t365 * t423;
t523 = t368 * t391;
t522 = t370 * t391;
t520 = t370 * t423;
t519 = t390 * t419;
t516 = t419 * t420;
t515 = t420 * t423;
t514 = t423 * t424;
t513 = t424 * t311;
t512 = qJDD(2) - g(3);
t489 = qJD(6) * t422;
t475 = t422 * t311 - t418 * t312 - t368 * t489;
t281 = -t370 * t490 + t475;
t495 = qJD(4) * t424;
t511 = t281 * t420 - t447 * t495;
t466 = t419 * t490;
t292 = -t424 * t466 + (t514 * t541 - t470) * t422 + t437 * t418;
t510 = -t292 * t388 - t345 * t358;
t509 = t311 * t420 + t370 * t495;
t441 = t371 * t420;
t508 = -qJD(1) * t441 - t371 * t541;
t438 = qJD(1) * t372;
t507 = t372 * t541 + t420 * t438;
t373 = t455 * qJD(1);
t506 = t423 * t332 + t419 * t373;
t367 = t390 * t515;
t505 = t419 * t359 + t367;
t413 = t424 ^ 2;
t503 = t420 ^ 2 - t413;
t426 = qJD(4) ^ 2;
t427 = qJD(1) ^ 2;
t502 = -t426 - t427;
t498 = qJD(4) * t315;
t497 = qJD(4) * t368;
t493 = qJD(5) * t420;
t484 = qJD(3) * qJD(1);
t480 = qJDD(4) * t420;
t479 = qJDD(4) * t424;
t478 = pkin(9) * t515;
t468 = t424 * t487;
t474 = t359 * t492 + t419 * t366 + t390 * t468;
t472 = qJD(5) * t539;
t471 = t419 * t500;
t464 = pkin(5) - t519;
t461 = t390 * t391 + t326;
t459 = -qJD(5) * t334 - t297;
t458 = qJD(6) * t286 + t277;
t457 = qJD(1) + t493;
t456 = -t333 + (t471 + t494) * pkin(5);
t421 = sin(qJ(1));
t425 = cos(qJ(1));
t452 = g(1) * t421 - g(2) * t425;
t357 = t423 * t373;
t384 = t539 * t423;
t451 = qJD(6) * t384 - t332 * t419 + t357 + (pkin(5) * t424 + t478) * qJD(1) + t423 * t472;
t383 = t539 * t419;
t450 = pkin(9) * t471 + qJD(6) * t383 + t419 * t472 + t506;
t279 = t286 * t418 + t531;
t291 = qJD(4) * t441 - t345 * t541;
t346 = t371 * t424;
t448 = -t291 * t388 + t346 * t358;
t445 = t391 * t492 + t525;
t444 = t391 * t494 - t524;
t442 = t388 * t371;
t440 = qJDD(3) + t542;
t436 = -t423 * t491 + t470;
t435 = -g(1) * t401 - g(2) * t400 + t504;
t434 = -pkin(8) * t365 + t325 * t391;
t433 = t360 + t544;
t364 = t484 + t504;
t429 = -t390 * t426 + t364 + t435 + t484;
t399 = -pkin(5) * t423 - pkin(4);
t377 = -t420 * t426 + t479;
t376 = -t424 * t426 - t480;
t362 = t391 * t470;
t353 = t423 * t366;
t347 = (pkin(5) * t419 - t390) * t424;
t344 = t423 * t359;
t340 = -t400 * t419 + t401 * t515;
t339 = t400 * t423 + t401 * t516;
t338 = t400 * t515 + t401 * t419;
t337 = -t400 * t516 + t401 * t423;
t322 = -pkin(5) * t436 + t390 * t496;
t309 = -t419 * t537 + t505;
t300 = -pkin(9) * t514 + t420 * t464 + t344;
t285 = pkin(9) * t436 - t493 * t519 + t474;
t284 = t353 + (-t367 + (-t359 + t537) * t419) * qJD(5) + (t424 * t464 + t478) * qJD(4);
t283 = pkin(5) * t312 + t298;
t278 = -t290 * t418 + t532;
t1 = [(qJDD(3) + 0.2e1 * t542 - t453) * MDP(5) + (t429 * t420 + t424 * t540) * MDP(13) + (-t420 * t540 + t429 * t424) * MDP(14) + (-t281 * t345 + t282 * t346 - t291 * t315 + t292 * t447) * MDP(23) + (-t281 * t346 - t291 * t447) * MDP(22) + (-t279 * t495 + g(1) * t330 - g(2) * t328 + t347 * t281 - t283 * t346 + t287 * t420 + t303 * t291 - t322 * t447 + (-(-qJD(6) * t309 + t284) * t388 - t300 * t358 - t276 * t420) * t418 + (-(qJD(6) * t300 + t285) * t388 - t309 * t358 - t458 * t420) * t422) * MDP(28) + (t365 * t514 + t391 * t437 + t509) * MDP(17) + (g(1) * t425 + g(2) * t421) * MDP(3) + (-t448 + t511) * MDP(24) + (t370 * t437 + t423 * t513) * MDP(15) + 0.2e1 * (-t420 * t482 + t486 * t503) * MDP(9) + (t435 + 0.2e1 * t484 + t504) * MDP(6) + (-t474 * t391 - t505 * t365 + g(1) * t339 - g(2) * t337 + (t461 * t494 + (-t325 * t423 + t370 * t390) * qJD(4) + t476) * t420 + (-t296 * qJD(4) + t298 * t423 - t390 * t311 - t325 * t494) * t424) * MDP(21) + (-g(1) * t340 - g(2) * t338 + t344 * t365 + t353 * t391 + (t390 * t497 - t461 * t492 + t302) * t420 + (t295 * qJD(4) - t390 * t312 + t325 * t492) * t424 + ((-qJD(5) * t359 - t390 * t495) * t391 + t298 * t424 + (-qJD(4) * t325 - t365 * t390 + t459) * t420) * t419) * MDP(20) + ((t284 * t422 - t285 * t418) * t388 + (t300 * t422 - t309 * t418) * t358 + t463 * t420 + t278 * t495 + t322 * t315 + t347 * t282 + t283 * t345 + t303 * t292 - g(1) * t331 - g(2) * t329 + ((-t300 * t418 - t309 * t422) * t388 - t279 * t420) * qJD(6)) * MDP(27) + (t365 * t420 + t391 * t495) * MDP(19) + (t358 * t420 + t388 * t495) * MDP(26) + t376 * MDP(11) + t377 * MDP(10) + (-t315 * t495 + t510 - t533) * MDP(25) + (-t529 + t362 + (-t445 - t497) * t424) * MDP(18) + ((t368 * t423 + t370 * t419) * t496 + (-t530 - t312 * t423 + (t368 * t419 - t520) * qJD(5)) * t424) * MDP(16) + qJDD(1) * MDP(1) + (qJDD(1) * t413 - 0.2e1 * t420 * t465) * MDP(8) + (t452 + (t416 ^ 2 + t417 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t364 * t392 + t378 * qJD(3) + t440 * t394 - g(1) * (-pkin(1) * t421 - pkin(2) * t400 + qJ(3) * t401) - g(2) * (pkin(1) * t425 + pkin(2) * t401 + qJ(3) * t400)) * MDP(7) + t452 * MDP(2); t376 * MDP(13) - t377 * MDP(14) + (t362 + t529) * MDP(20) + (t391 * t469 + t509) * MDP(21) + (t510 + t533) * MDP(27) + (t448 + t511) * MDP(28) + (MDP(4) + MDP(7)) * t512 + ((-t445 + t497) * MDP(20) + t444 * MDP(21) + MDP(27) * t498) * t424; qJDD(1) * MDP(5) - t427 * MDP(6) + (t440 + t544) * MDP(7) + (t420 * t502 + t479) * MDP(13) + (t424 * t502 - t480) * MDP(14) + (-t424 * t312 + (t497 - t525) * t420 + (-t419 * t495 - t423 * t457) * t391) * MDP(20) + (-t513 + (qJD(4) * t370 - t524) * t420 + (t419 * t457 - t468) * t391) * MDP(21) + (qJD(1) * t442 + (-qJD(4) * t372 * t388 - t282) * t424 + ((t418 * t494 + t422 * t543 + t466) * t388 - t527 + t498) * t420) * MDP(27) + (t388 * t438 + (qJD(4) * t442 - t281) * t424 + (-(t418 * t543 - t419 * t489 - t422 * t494) * t388 + t528 - qJD(4) * t447) * t420) * MDP(28); MDP(10) * t482 - MDP(11) * t483 + qJDD(4) * MDP(12) + (qJD(4) * t333 + t424 * t433 + t473 + t535) * MDP(13) + (t549 + (-qJD(4) * t363 - t512) * t424 + (-t433 + t485) * t420) * MDP(14) + (t391 * t520 + t530) * MDP(15) + ((t311 - t523) * t423 + (-t312 - t522) * t419) * MDP(16) + ((-t370 * t424 + t391 * t515) * qJD(1) + t445) * MDP(17) + ((t368 * t424 - t391 * t516) * qJD(1) - t444) * MDP(18) + (-pkin(4) * t312 - t333 * t368 - t357 * t391 + (t332 * t391 + t434) * t419 - t552 * t423) * MDP(20) + (-pkin(4) * t311 - t333 * t370 + t506 * t391 + t419 * t552 + t434 * t423) * MDP(21) + (t281 * t372 - t447 * t508) * MDP(22) + (-t281 * t371 - t282 * t372 - t315 * t508 + t447 * t507) * MDP(23) + (t388 * t508 + t527) * MDP(24) + (-t388 * t507 - t528) * MDP(25) + ((-t383 * t422 - t384 * t418) * t358 + t399 * t282 + t283 * t371 + (t418 * t450 - t422 * t451) * t388 + t456 * t315 + t507 * t303 - t432 * t409) * MDP(27) + (-(-t383 * t418 + t384 * t422) * t358 + t399 * t281 + t283 * t372 + (t418 * t451 + t422 * t450) * t388 - t456 * t447 + t508 * t303 + t432 * t408) * MDP(28) + (-t391 * MDP(19) - MDP(20) * t295 + t296 * MDP(21) + MDP(24) * t447 + t315 * MDP(25) - t388 * MDP(26) - t278 * MDP(27) + t279 * MDP(28)) * t499 + (MDP(8) * t420 * t424 - MDP(9) * t503) * t427; t370 * t368 * MDP(15) + (-t368 ^ 2 + t370 ^ 2) * MDP(16) + (t311 + t523) * MDP(17) + (-t312 + t522) * MDP(18) + t365 * MDP(19) + (-t326 * t492 - g(1) * t337 - g(2) * t339 + t296 * t391 - t325 * t370 + t302 + (t459 + t534) * t419) * MDP(20) + (g(1) * t338 - g(2) * t340 + g(3) * t514 + t295 * t391 + t325 * t368 - t439) * MDP(21) + (t281 + t551) * MDP(24) + (-t282 - t550) * MDP(25) + (-(-t289 * t418 - t531) * t388 - t279 * qJD(6) + (-t315 * t370 + t422 * t358 - t388 * t490) * pkin(5) + t547) * MDP(27) + ((-t290 * t388 - t276) * t418 + (t289 * t388 - t458) * t422 + (-t418 * t358 + t370 * t447 - t388 * t489) * pkin(5) + t548) * MDP(28) + t546; (t475 + t551) * MDP(24) + (-t462 - t550) * MDP(25) + (t279 * t388 + t547) * MDP(27) + (-t418 * t276 - t422 * t277 + t278 * t388 + t548) * MDP(28) + (-MDP(24) * t521 + MDP(25) * t447 - MDP(27) * t279 - MDP(28) * t532) * qJD(6) + t546;];
tau  = t1;
