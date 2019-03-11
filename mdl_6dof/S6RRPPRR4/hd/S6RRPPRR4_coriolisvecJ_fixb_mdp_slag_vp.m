% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:13
% EndTime: 2019-03-09 09:05:23
% DurationCPUTime: 5.44s
% Computational Cost: add. (5548->500), mult. (16794->654), div. (0->0), fcn. (13419->10), ass. (0->221)
t458 = sin(pkin(11));
t460 = cos(pkin(11));
t467 = cos(qJ(2));
t459 = sin(pkin(6));
t542 = qJD(1) * t459;
t524 = t467 * t542;
t504 = t460 * t524;
t464 = sin(qJ(2));
t520 = t464 * t542;
t425 = t458 * t520 - t504;
t461 = cos(pkin(6));
t541 = qJD(1) * t461;
t448 = qJD(2) + t541;
t463 = sin(qJ(5));
t466 = cos(qJ(5));
t389 = -t466 * t425 + t448 * t463;
t388 = qJD(6) + t389;
t540 = qJD(2) * t459;
t523 = t464 * t540;
t584 = pkin(2) * t523;
t583 = MDP(4) * t464;
t582 = MDP(5) * (t464 ^ 2 - t467 ^ 2);
t488 = t458 * t467 + t460 * t464;
t432 = t488 * t459;
t420 = qJD(1) * t432 + qJD(5);
t465 = cos(qJ(6));
t391 = t425 * t463 + t448 * t466;
t462 = sin(qJ(6));
t558 = t391 * t462;
t360 = -t465 * t420 + t558;
t429 = t488 * t542;
t531 = qJD(5) + t429;
t581 = t360 * t531;
t509 = t466 * t531;
t566 = pkin(8) + qJ(3);
t519 = t566 * t464;
t503 = t459 * t519;
t570 = pkin(1) * t467;
t407 = (pkin(2) + t570) * t461 - t503;
t571 = pkin(1) * t464;
t530 = t461 * t571;
t552 = t459 * t467;
t422 = t552 * t566 + t530;
t370 = t407 * t460 - t458 * t422;
t572 = pkin(3) + pkin(9);
t339 = pkin(4) * t432 - t461 * t572 - t370;
t553 = t459 * t464;
t431 = t458 * t553 - t460 * t552;
t500 = (-pkin(2) * t467 - pkin(1)) * t459;
t476 = -qJ(4) * t432 + t500;
t353 = t431 * t572 + t476;
t580 = t463 * t339 + t466 * t353;
t556 = t429 * t463;
t377 = -t425 * t462 + t465 * t556;
t538 = qJD(5) * t465;
t579 = t463 * t538 + t377;
t577 = MDP(14) * t429;
t409 = t422 * qJD(1);
t396 = t458 * t409;
t529 = t461 * t570;
t446 = qJD(1) * t529;
t408 = -qJD(1) * t503 + t446;
t367 = t408 * t460 - t396;
t532 = qJD(4) - t367;
t501 = qJD(2) * t520;
t419 = qJD(2) * t504 - t458 * t501;
t392 = pkin(2) * t448 + t408;
t356 = t392 * t460 - t396;
t498 = qJD(4) - t356;
t567 = pkin(4) * t429;
t323 = -t448 * t572 + t498 + t567;
t487 = qJD(1) * t500;
t433 = qJD(3) + t487;
t469 = -qJ(4) * t429 + t433;
t338 = t425 * t572 + t469;
t303 = t323 * t463 + t338 * t466;
t440 = qJD(2) * t446;
t470 = (-qJD(2) * t519 + qJD(3) * t467) * t459;
t382 = qJD(1) * t470 + t440;
t394 = -qJD(2) * t422 - qJD(3) * t553;
t383 = t394 * qJD(1);
t335 = t382 * t458 - t460 * t383;
t321 = pkin(4) * t419 + t335;
t427 = qJD(2) * t432;
t418 = qJD(1) * t427;
t439 = pkin(2) * t501;
t516 = -qJ(4) * t419 + t439;
t482 = -qJD(4) * t429 + t516;
t322 = t418 * t572 + t482;
t513 = -t466 * t321 + t322 * t463;
t574 = -t303 * qJD(5) - t513;
t292 = -pkin(5) * t419 - t574;
t576 = t388 * (pkin(5) * t391 + pkin(10) * t388) + t292;
t551 = t460 * t409;
t357 = t458 * t392 + t551;
t347 = -t448 * qJ(4) - t357;
t568 = pkin(4) * t425;
t327 = -t347 - t568;
t452 = -pkin(2) * t460 - pkin(3);
t450 = -pkin(9) + t452;
t575 = -t327 * t531 - t450 * t419;
t421 = t429 ^ 2;
t569 = pkin(3) * t418;
t537 = qJD(5) * t466;
t539 = qJD(5) * t463;
t358 = t463 * t418 + t425 * t537 - t448 * t539;
t535 = qJD(6) * t465;
t525 = t465 * t358 + t462 * t419 + t420 * t535;
t536 = qJD(6) * t462;
t310 = -t391 * t536 + t525;
t565 = t310 * t462;
t362 = t391 * t465 + t420 * t462;
t512 = t358 * t462 - t465 * t419;
t311 = qJD(6) * t362 + t512;
t564 = t311 * t463;
t563 = t360 * t388;
t562 = t362 * t388;
t561 = t388 * t462;
t507 = t388 * t465;
t560 = t389 * t425;
t559 = t391 * t425;
t557 = t425 * t448;
t455 = t459 ^ 2;
t468 = qJD(1) ^ 2;
t554 = t455 * t468;
t405 = t466 * t418;
t359 = qJD(5) * t391 - t405;
t550 = t462 * t359;
t549 = t463 * t419;
t548 = t465 * t359;
t547 = t466 * t310;
t499 = pkin(5) * t466 + pkin(10) * t463;
t546 = (-pkin(4) - t499) * t429 - qJD(5) * t499 - t532;
t366 = t408 * t458 + t551;
t340 = t366 - t568;
t515 = pkin(2) * t520 + qJ(4) * t425;
t345 = t429 * t572 + t515;
t545 = t463 * t340 + t466 * t345;
t336 = t460 * t382 + t458 * t383;
t447 = qJD(2) * t529;
t393 = t447 + t470;
t351 = t460 * t393 + t458 * t394;
t371 = t458 * t407 + t460 * t422;
t533 = t567 + t532;
t528 = t455 * t571;
t527 = t467 * t554;
t526 = t310 * t463 + (t429 * t466 + t537) * t362;
t330 = -t448 * qJD(4) - t336;
t344 = -t461 * qJD(4) - t351;
t364 = -t461 * qJ(4) - t371;
t522 = t467 * t540;
t518 = qJD(1) * qJD(2) * t455;
t451 = pkin(2) * t458 + qJ(4);
t478 = t463 * t321 + t466 * t322 + t323 * t537 - t338 * t539;
t291 = pkin(10) * t419 + t478;
t316 = -pkin(4) * t418 - t330;
t298 = pkin(5) * t359 - pkin(10) * t358 + t316;
t514 = -t291 * t462 + t465 * t298;
t376 = t465 * t425 + t462 * t556;
t511 = t376 * t388 + t539 * t561;
t350 = t393 * t458 - t460 * t394;
t510 = t420 * t531;
t508 = t531 * t391;
t434 = pkin(5) * t463 - pkin(10) * t466 + t451;
t505 = -pkin(10) * t425 - qJD(6) * t434 + t545;
t502 = t467 * t518;
t497 = t291 * t465 + t298 * t462;
t300 = pkin(10) * t420 + t303;
t308 = pkin(5) * t389 - pkin(10) * t391 + t327;
t294 = t300 * t465 + t308 * t462;
t496 = t300 * t462 - t308 * t465;
t307 = pkin(10) * t432 + t580;
t343 = -pkin(4) * t431 - t364;
t400 = t431 * t463 + t461 * t466;
t489 = t431 * t466 - t461 * t463;
t312 = -pkin(5) * t489 - pkin(10) * t400 + t343;
t495 = t307 * t465 + t312 * t462;
t494 = -t307 * t462 + t312 * t465;
t302 = t323 * t466 - t338 * t463;
t428 = -t458 * t523 + t460 * t522;
t481 = -qJ(4) * t428 - qJD(4) * t432 + t584;
t328 = t427 * t572 + t481;
t329 = pkin(4) * t428 + t350;
t493 = -t328 * t463 + t329 * t466;
t492 = t335 * t432 + t350 * t429;
t491 = t339 * t466 - t353 * t463;
t374 = t400 * t465 + t432 * t462;
t373 = t400 * t462 - t432 * t465;
t324 = -pkin(4) * t427 - t344;
t363 = pkin(3) * t425 + t469;
t486 = t363 * t429 + t335;
t485 = t466 * t419 - t463 * t510;
t484 = -t388 * t535 - t550;
t483 = t388 * t536 - t548;
t480 = -pkin(8) * t552 - t530;
t479 = -pkin(8) * t501 + t440;
t477 = t466 * t328 + t463 * t329 + t339 * t537 - t353 * t539;
t475 = -t420 * t509 - t549;
t474 = t480 * t448;
t473 = -t466 * t536 - t579;
t472 = qJD(5) * t360 + t484;
t299 = -pkin(5) * t420 - t302;
t471 = -pkin(10) * t359 + (t299 + t302) * t388;
t375 = pkin(3) * t431 + t476;
t372 = pkin(3) * t429 + t515;
t369 = qJD(5) * t400 - t427 * t466;
t368 = qJD(5) * t489 + t427 * t463;
t365 = -pkin(3) * t461 - t370;
t349 = pkin(3) * t427 + t481;
t346 = -pkin(3) * t448 + t498;
t337 = t482 + t569;
t314 = -qJD(6) * t373 + t368 * t465 + t428 * t462;
t313 = qJD(6) * t374 + t368 * t462 - t428 * t465;
t306 = -pkin(5) * t432 - t491;
t304 = pkin(5) * t425 - t340 * t466 + t345 * t463;
t301 = pkin(5) * t369 - pkin(10) * t368 + t324;
t296 = -pkin(5) * t428 + qJD(5) * t580 - t493;
t295 = pkin(10) * t428 + t477;
t290 = -qJD(6) * t294 + t514;
t289 = -qJD(6) * t496 + t497;
t1 = [(-0.2e1 * pkin(1) * t502 - (-pkin(8) * t523 + t447) * t448 - t479 * t461) * MDP(10) + (t474 + (t461 * t480 - 0.2e1 * t528) * qJD(1)) * qJD(2) * MDP(9) + (t335 * t461 - t337 * t431 - t349 * t425 + t350 * t448 - t363 * t427 - t375 * t418) * MDP(14) + (-t330 * t461 - t337 * t432 - t344 * t448 - t349 * t429 - t363 * t428 - t375 * t419) * MDP(15) + (t358 * t432 + t368 * t420 + t391 * t428 + t400 * t419) * MDP(19) + (t419 * t432 + t420 * t428) * MDP(21) + (t358 * t400 + t368 * t391) * MDP(17) + (-t310 * t373 - t311 * t374 - t313 * t362 - t314 * t360) * MDP(25) + (t310 * t374 + t314 * t362) * MDP(24) + (t330 * t364 + t335 * t365 + t337 * t375 + t344 * t347 + t346 * t350 + t349 * t363) * MDP(16) + (-t335 * t370 + t336 * t371 - t356 * t350 + t357 * t351 + (t433 + t487) * t584) * MDP(12) - 0.2e1 * t518 * t582 + 0.2e1 * t502 * t583 + (-t303 * t428 + t316 * t400 + t324 * t391 + t327 * t368 + t343 * t358 - t419 * t580 - t420 * t477 - t432 * t478) * MDP(23) + (-t336 * t431 - t351 * t425 - t356 * t428 - t357 * t427 - t370 * t419 - t371 * t418 + t492) * MDP(11) + (t330 * t431 + t344 * t425 + t346 * t428 + t347 * t427 + t364 * t418 + t365 * t419 + t492) * MDP(13) + (MDP(6) * t522 - MDP(7) * t523) * (t448 + t541) + (-t359 * t432 - t369 * t420 - t389 * t428 + t419 * t489) * MDP(20) + (-t359 * t489 + t369 * t388) * MDP(28) + (t311 * t489 - t313 * t388 - t359 * t373 - t360 * t369) * MDP(27) + (-t310 * t489 + t314 * t388 + t359 * t374 + t362 * t369) * MDP(26) + (t358 * t489 - t359 * t400 - t368 * t389 - t369 * t391) * MDP(18) + (-(qJD(6) * t494 + t295 * t465 + t301 * t462) * t388 - t495 * t359 + t289 * t489 - t294 * t369 + t296 * t362 + t306 * t310 + t292 * t374 + t299 * t314) * MDP(30) + ((-qJD(6) * t495 - t295 * t462 + t301 * t465) * t388 + t494 * t359 - t290 * t489 - t496 * t369 + t296 * t360 + t306 * t311 + t292 * t373 + t299 * t313) * MDP(29) + (t493 * t420 + t491 * t419 - t513 * t432 + t302 * t428 + t324 * t389 + t343 * t359 - t316 * t489 + t327 * t369 + (-t303 * t432 - t420 * t580) * qJD(5)) * MDP(22); (t362 * t473 + t465 * t547) * MDP(24) + (t388 * t473 + t466 * t548 + t526) * MDP(26) + (-t330 * t451 + t335 * t452 - t346 * t366 - t347 * t532 - t363 * t372) * MDP(16) + (t468 * t528 + (qJD(2) * t480 - t474) * qJD(1)) * MDP(9) + (pkin(1) * t527 + (-pkin(8) * t520 + t446) * t448 - t479) * MDP(10) + (t356 * t366 - t357 * t367 + (-t335 * t460 + t336 * t458 - t433 * t520) * pkin(2)) * MDP(12) + (-t366 * t448 + t372 * t425 + t486) * MDP(14) + (MDP(6) * t524 - MDP(7) * t520) * (qJD(2) - t448) + (-t418 * t451 + t419 * t452 + (-t347 - t366) * t429 + (t346 - t532) * t425) * MDP(13) + (-t363 * t425 + t372 * t429 + t448 * t532 - t330) * MDP(15) + t420 * t425 * MDP(21) - t527 * t583 + (-t564 + (t484 - t581) * t466 + t511) * MDP(27) + (t359 * t463 + t388 * t509) * MDP(28) + (t302 * t425 + t451 * t359 + t533 * t389 + (t316 + (-qJD(5) * t450 + t345) * t420) * t463 + (-t340 * t420 - t575) * t466) * MDP(22) + (-t303 * t425 + t316 * t466 + t451 * t358 + (-t450 * t537 + t545) * t420 + t533 * t391 + t575 * t463) * MDP(23) + ((t357 - t366) * t429 + (-t356 + t367) * t425 + (-t418 * t458 - t419 * t460) * pkin(2)) * MDP(11) + (-t434 * t550 - t299 * t377 - t304 * t362 + (t462 * t546 + t465 * t505) * t388 + (-t299 * t538 - t289 + (qJD(5) * t362 + t483) * t450) * t463 + (-t299 * t536 + t292 * t465 - t294 * t429 - t450 * t310 + (-t450 * t507 - t294) * qJD(5)) * t466) * MDP(30) + ((-t359 - t508) * t466 + (t389 * t531 - t358) * t463) * MDP(18) + (t358 * t466 - t463 * t508) * MDP(17) + (t434 * t548 - t299 * t376 - t304 * t360 + (t462 * t505 - t465 * t546) * t388 + (-t299 * t462 * qJD(5) + t450 * t472 + t290) * t463 + (t299 * t535 + t292 * t462 - t496 * t429 - t450 * t311 + (-t450 * t561 - t496) * qJD(5)) * t466) * MDP(29) + (t485 + t559) * MDP(19) + (t475 - t560) * MDP(20) + (t360 * t377 + t362 * t376 + (t360 * t465 + t362 * t462) * t539 + (-t565 - t311 * t465 + (t360 * t462 - t362 * t465) * qJD(6)) * t466) * MDP(25) + t554 * t582; (t356 * t429 + t357 * t425 + t439) * MDP(12) - (qJD(2) + t448) * t577 + (-t419 + t557) * MDP(15) + (t569 - t347 * t425 + (-qJD(4) - t346) * t429 + t516) * MDP(16) + (-t549 + t560) * MDP(22) + (t559 + (t539 + t556) * t420) * MDP(23) + (t511 + t564) * MDP(29) + (t388 * t579 + t526) * MDP(30) + (-t419 * MDP(23) + (t360 * t429 + t472) * MDP(29) + t483 * MDP(30) - MDP(22) * t510) * t466 + (MDP(11) + MDP(13)) * (-t425 ^ 2 - t421); (t419 + t557) * MDP(13) - t425 * t577 + (-t448 ^ 2 - t421) * MDP(15) + (t347 * t448 + t486) * MDP(16) + (-t389 * t448 + t485) * MDP(22) + (-t391 * t448 + t475) * MDP(23) + (-t466 * t311 + (-t465 * t448 - t462 * t509) * t388 + (t484 + t581) * t463) * MDP(29) + (-t547 + (t462 * t448 - t465 * t509) * t388 + (t362 * t531 + t483) * t463) * MDP(30); -t389 ^ 2 * MDP(18) + (t389 * t420 + t358) * MDP(19) + t405 * MDP(20) + t419 * MDP(21) + (t303 * t420 + t574) * MDP(22) + (t302 * t420 + t327 * t389 - t478) * MDP(23) + (t362 * t507 + t565) * MDP(24) + ((t310 - t563) * t465 + (-t311 - t562) * t462) * MDP(25) + (t388 * t507 + t550) * MDP(26) + (-t388 ^ 2 * t462 + t548) * MDP(27) + (-pkin(5) * t311 - t303 * t360 + t471 * t462 - t465 * t576) * MDP(29) + (-pkin(5) * t310 - t303 * t362 + t462 * t576 + t471 * t465) * MDP(30) + (t389 * MDP(17) + (-qJD(5) + t420) * MDP(20) - t327 * MDP(22) - t362 * MDP(26) + t360 * MDP(27) - t388 * MDP(28) + t496 * MDP(29) + t294 * MDP(30) + t391 * MDP(18)) * t391; t362 * t360 * MDP(24) + (-t360 ^ 2 + t362 ^ 2) * MDP(25) + (t525 + t563) * MDP(26) + (-t512 + t562) * MDP(27) + t359 * MDP(28) + (t294 * t388 - t299 * t362 + t514) * MDP(29) + (t299 * t360 - t388 * t496 - t497) * MDP(30) + (-MDP(26) * t558 - MDP(27) * t362 - MDP(29) * t294 + MDP(30) * t496) * qJD(6);];
tauc  = t1;
