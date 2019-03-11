% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:17
% EndTime: 2019-03-08 20:25:24
% DurationCPUTime: 5.43s
% Computational Cost: add. (3423->407), mult. (7689->566), div. (0->0), fcn. (6444->16), ass. (0->200)
t437 = qJD(4) + qJD(5);
t448 = sin(qJ(5));
t452 = cos(qJ(4));
t559 = cos(qJ(5));
t507 = qJD(2) * t559;
t449 = sin(qJ(4));
t525 = qJD(2) * t449;
t567 = -t448 * t525 + t452 * t507;
t568 = t437 * t567;
t436 = qJDD(4) + qJDD(5);
t446 = cos(pkin(6));
t421 = qJD(1) * t446 + qJD(3);
t453 = cos(qJ(2));
t443 = sin(pkin(6));
t527 = qJD(1) * t443;
t509 = t453 * t527;
t398 = qJD(2) * pkin(2) + t509;
t444 = cos(pkin(12));
t450 = sin(qJ(2));
t510 = t450 * t527;
t404 = t444 * t510;
t441 = sin(pkin(12));
t368 = t441 * t398 + t404;
t500 = t368 + (pkin(8) + pkin(9)) * qJD(2);
t343 = t421 * t449 + t452 * t500;
t414 = t446 * qJDD(1) + qJDD(3);
t402 = t452 * t414;
t517 = t443 * qJDD(1);
t413 = t453 * t517;
t526 = qJD(2) * t443;
t505 = qJD(1) * t526;
t379 = qJDD(2) * pkin(2) - t450 * t505 + t413;
t565 = t450 * t517 + t453 * t505;
t345 = t441 * t379 + t444 * t565;
t338 = qJDD(2) * pkin(8) + t345;
t499 = pkin(9) * qJDD(2) + t338;
t299 = qJDD(4) * pkin(4) - qJD(4) * t343 - t499 * t449 + t402;
t342 = t452 * t421 - t500 * t449;
t300 = qJD(4) * t342 + t449 * t414 + t499 * t452;
t334 = qJD(4) * pkin(4) + t342;
t511 = t559 * t343;
t306 = t448 * t334 + t511;
t561 = qJD(5) * t306 - t559 * t299 + t448 * t300;
t287 = -t436 * pkin(5) + t561;
t478 = t441 * t453 + t444 * t450;
t381 = t478 * t443;
t440 = qJ(4) + qJ(5);
t434 = sin(t440);
t435 = cos(t440);
t382 = t478 * t446;
t534 = t453 * t444;
t392 = t441 * t450 - t534;
t442 = sin(pkin(11));
t445 = cos(pkin(11));
t480 = -t382 * t442 - t392 * t445;
t481 = t382 * t445 - t392 * t442;
t542 = t443 * t445;
t543 = t442 * t443;
t469 = -g(3) * (-t381 * t434 + t435 * t446) - g(2) * (-t434 * t481 - t435 * t542) - g(1) * (-t434 * t480 + t435 * t543);
t466 = -t287 + t469;
t476 = (g(1) * t442 - g(2) * t445) * t443;
t566 = -t414 + t476;
t536 = t448 * t452;
t395 = t449 * t559 + t536;
t362 = t437 * t395;
t501 = qJDD(2) * t559;
t519 = qJDD(2) * t449;
t488 = t448 * t519 - t452 * t501;
t340 = qJD(2) * t362 + t488;
t335 = qJDD(6) + t340;
t512 = -pkin(4) * t452 - pkin(3);
t558 = pkin(2) * t444;
t406 = t512 - t558;
t473 = -t448 * t449 + t452 * t559;
t353 = -pkin(5) * t473 - pkin(10) * t395 + t406;
t470 = t392 * t446;
t347 = -t442 * t478 - t445 * t470;
t350 = t442 * t470 - t445 * t478;
t541 = t443 * t450;
t380 = t441 * t541 - t443 * t534;
t468 = g(1) * t350 + g(2) * t347 - g(3) * t380;
t464 = t468 * t435;
t564 = t353 * t335 - t464;
t363 = -t381 * t449 + t446 * t452;
t364 = t381 * t452 + t446 * t449;
t319 = t448 * t363 + t364 * t559;
t451 = cos(qJ(6));
t374 = t380 * t451;
t447 = sin(qJ(6));
t563 = -t319 * t447 + t374;
t375 = t441 * t509 + t404;
t524 = qJD(4) * t449;
t491 = pkin(4) * t524 - t375;
t506 = qJD(5) * t559;
t523 = qJD(5) * t448;
t462 = t448 * t299 + t300 * t559 + t334 * t506 - t343 * t523;
t286 = t436 * pkin(10) + t462;
t537 = t448 * t343;
t305 = t334 * t559 - t537;
t303 = -t437 * pkin(5) - t305;
t403 = t441 * t510;
t367 = t398 * t444 - t403;
t357 = qJD(2) * t512 - t367;
t388 = -qJD(2) * t536 - t449 * t507;
t317 = -pkin(5) * t567 + pkin(10) * t388 + t357;
t427 = pkin(2) * t441 + pkin(8);
t554 = pkin(9) + t427;
t390 = t554 * t449;
t391 = t554 * t452;
t356 = -t448 * t390 + t391 * t559;
t361 = t437 * t473;
t385 = qJD(6) - t567;
t378 = t444 * t509 - t403;
t502 = qJD(4) * t554;
t383 = t449 * t502;
t384 = t452 * t502;
t474 = -t390 * t559 - t448 * t391;
t531 = -t474 * qJD(5) + t473 * t378 + t383 * t559 + t448 * t384;
t560 = -g(1) * t480 - g(2) * t481;
t562 = (qJD(6) * t317 + t286) * t473 + t287 * t395 + t303 * t361 + (-qJD(6) * t353 + t531) * t385 - g(3) * t381 - t356 * t335 + t560;
t553 = t303 * t567;
t552 = t303 * t395;
t518 = qJDD(2) * t452;
t339 = t448 * t518 + t449 * t501 + t568;
t521 = qJD(6) * t451;
t513 = t451 * t339 + t447 * t436 + t437 * t521;
t522 = qJD(6) * t447;
t310 = t388 * t522 + t513;
t551 = t310 * t447;
t545 = t388 * t447;
t369 = -t451 * t437 - t545;
t548 = t369 * t385;
t479 = t388 * t451 - t437 * t447;
t547 = t479 * t385;
t546 = t380 * t447;
t540 = t446 * t450;
t539 = t446 * t453;
t538 = t447 * t335;
t535 = t451 * t335;
t533 = qJDD(1) - g(3);
t532 = -t310 * t473 - t362 * t479;
t530 = t356 * qJD(5) - t395 * t378 - t448 * t383 + t384 * t559;
t529 = pkin(5) * t362 - pkin(10) * t361 + t491;
t438 = t449 ^ 2;
t528 = -t452 ^ 2 + t438;
t520 = qJD(2) * qJD(4);
t515 = t395 * t538;
t514 = t395 * t535;
t504 = t449 * t520;
t304 = t437 * pkin(10) + t306;
t485 = t304 * t447 - t317 * t451;
t498 = t303 * t522 - t388 * t485;
t496 = t339 * t447 - t451 * t436;
t495 = t385 * t451;
t358 = -pkin(5) * t388 - pkin(10) * t567;
t431 = pkin(4) * t448 + pkin(10);
t493 = pkin(4) * t525 + qJD(6) * t431 + t358;
t307 = t448 * t342 + t511;
t490 = pkin(4) * t523 - t307;
t344 = t379 * t444 - t441 * t565;
t486 = -t431 * t335 - t553;
t292 = t304 * t451 + t317 * t447;
t311 = -qJD(6) * t479 + t496;
t484 = t311 * t473 - t362 * t369;
t483 = t319 * t451 + t546;
t482 = t361 * t437 + t395 * t436;
t308 = t342 * t559 - t537;
t477 = -pkin(4) * t506 + t308;
t475 = t363 * t559 - t448 * t364;
t472 = -t361 * t447 - t395 * t521;
t471 = -t361 * t451 + t395 * t522;
t467 = -t292 * t388 + t303 * t521 - t447 * t466;
t365 = -qJD(2) * pkin(3) - t367;
t465 = -qJD(2) * t365 - t338 - t560;
t428 = -pkin(3) - t558;
t463 = -qJDD(4) * t427 + (qJD(2) * t428 + t365 + t378) * qJD(4);
t320 = pkin(4) * t504 + qJDD(2) * t512 - t344;
t461 = ((t310 - t548) * t451 + (-t311 + t547) * t447) * MDP(21) + (-t479 * t495 + t551) * MDP(20) + (-t385 ^ 2 * t447 - t369 * t388 + t535) * MDP(23) + (t385 * t495 - t388 * t479 + t538) * MDP(22) + (t339 - t568) * MDP(15) + (-t488 + (-qJD(2) * t395 - t388) * t437) * MDP(16) + (t388 ^ 2 - t567 ^ 2) * MDP(14) + t436 * MDP(17) + (MDP(13) * t567 + MDP(24) * t385) * t388;
t454 = qJD(4) ^ 2;
t460 = -qJD(2) * t375 + t427 * t454 - t344 + t468 + (-pkin(3) + t428) * qJDD(2);
t459 = -g(1) * (-t442 * t539 - t445 * t450) - g(2) * (-t442 * t450 + t445 * t539) - g(3) * t443 * t453;
t457 = t357 * t388 + t469 - t561;
t326 = -t434 * t542 + t435 * t481;
t328 = t434 * t543 + t435 * t480;
t360 = t381 * t435 + t434 * t446;
t456 = g(1) * t328 + g(2) * t326 + g(3) * t360 - t357 * t567 - t462;
t455 = qJD(2) ^ 2;
t432 = -pkin(4) * t559 - pkin(5);
t408 = qJDD(4) * t452 - t449 * t454;
t407 = qJDD(4) * t449 + t452 * t454;
t377 = t392 * t526;
t376 = qJD(2) * t381;
t336 = -t362 * t437 + t436 * t473;
t324 = qJD(4) * t363 - t377 * t452;
t323 = -qJD(4) * t364 + t377 * t449;
t295 = pkin(5) * t340 - pkin(10) * t339 + t320;
t294 = t451 * t295;
t290 = t319 * qJD(5) - t323 * t559 + t448 * t324;
t289 = t475 * qJD(5) + t448 * t323 + t324 * t559;
t1 = [t533 * MDP(1) + (-t344 * t380 + t345 * t381 - t367 * t376 - t368 * t377 + t414 * t446 - g(3)) * MDP(5) + (qJD(4) * t323 + qJDD(4) * t363 - t380 * t518) * MDP(11) + (-qJD(4) * t324 - qJDD(4) * t364 + t380 * t519) * MDP(12) + (-t290 * t437 + t340 * t380 - t376 * t567 + t436 * t475) * MDP(18) + (-t289 * t437 - t319 * t436 + t339 * t380 - t376 * t388) * MDP(19) + ((-qJD(6) * t483 - t289 * t447 + t376 * t451) * t385 + t563 * t335 + t290 * t369 - t475 * t311) * MDP(25) + (-(qJD(6) * t563 + t289 * t451 + t376 * t447) * t385 - t483 * t335 - t290 * t479 - t475 * t310) * MDP(26) + ((-t376 * t452 + t380 * t524) * MDP(11) + (qJD(4) * t380 * t452 + t376 * t449) * MDP(12)) * qJD(2) + ((qJDD(2) * t453 - t450 * t455) * MDP(3) + (-qJDD(2) * t450 - t453 * t455) * MDP(4)) * t443; qJDD(2) * MDP(2) + (t413 + t459) * MDP(3) + (-g(1) * (t442 * t540 - t445 * t453) - g(2) * (-t442 * t453 - t445 * t540) - t533 * t541) * MDP(4) + (t367 * t375 - t368 * t378 + (t344 * t444 + t345 * t441 + t459) * pkin(2)) * MDP(5) + (qJDD(2) * t438 + 0.2e1 * t452 * t504) * MDP(6) + 0.2e1 * (t449 * t518 - t520 * t528) * MDP(7) + t407 * MDP(8) + t408 * MDP(9) + (t449 * t463 - t452 * t460) * MDP(11) + (t449 * t460 + t452 * t463) * MDP(12) + (t339 * t395 - t361 * t388) * MDP(13) + (t339 * t473 - t340 * t395 + t361 * t567 + t362 * t388) * MDP(14) + t482 * MDP(15) + t336 * MDP(16) + (-t320 * t473 + t340 * t406 + t357 * t362 + t436 * t474 - t437 * t530 - t491 * t567 - t464) * MDP(18) + (t320 * t395 + t339 * t406 - t356 * t436 + t357 * t361 - t388 * t491 + t434 * t468 + t437 * t531) * MDP(19) + (t310 * t395 * t451 + t471 * t479) * MDP(20) + ((-t369 * t451 + t447 * t479) * t361 + (-t551 - t311 * t451 + (t369 * t447 + t451 * t479) * qJD(6)) * t395) * MDP(21) + (-t385 * t471 + t514 + t532) * MDP(22) + (t385 * t472 + t484 - t515) * MDP(23) + (-t335 * t473 + t362 * t385) * MDP(24) + (-t485 * t362 - t294 * t473 - t474 * t311 + t530 * t369 + (t529 * t385 + (t304 * t473 - t356 * t385 + t552) * qJD(6) + t564) * t451 + t562 * t447) * MDP(25) + (-t292 * t362 - t474 * t310 - t530 * t479 + ((-qJD(6) * t304 + t295) * t473 - qJD(6) * t552 + (qJD(6) * t356 - t529) * t385 - t564) * t447 + t562 * t451) * MDP(26); (-g(3) * t446 - t566) * MDP(5) + t408 * MDP(11) - t407 * MDP(12) + t336 * MDP(18) - t482 * MDP(19) + (-t484 - t515) * MDP(25) + (-t514 + t532) * MDP(26) + (MDP(25) * t472 + MDP(26) * t471) * t385; (t432 * t311 + t490 * t369 + (t385 * t477 + t486) * t447 + (-t385 * t493 + t466) * t451 + t498) * MDP(25) + MDP(9) * t518 + MDP(8) * t519 + (g(3) * t364 + t449 * t566 + t465 * t452) * MDP(12) + qJDD(4) * MDP(10) + t461 + (-g(3) * t363 + t449 * t465 - t452 * t476 + t402) * MDP(11) + (t432 * t310 + t486 * t451 - t490 * t479 + (t447 * t493 + t451 * t477) * t385 + t467) * MDP(26) + (t308 * t437 + (t388 * t525 - t436 * t448 - t437 * t506) * pkin(4) + t456) * MDP(19) + (t307 * t437 + (t436 * t559 - t437 * t523 + t525 * t567) * pkin(4) + t457) * MDP(18) + (-MDP(6) * t449 * t452 + MDP(7) * t528) * t455; t461 + (t305 * t437 + t456) * MDP(19) + (t306 * t437 + t457) * MDP(18) + (-pkin(5) * t310 + (t305 * t451 + t358 * t447) * t385 + t306 * t479 - t451 * t553 + (t385 * t522 - t535) * pkin(10) + t467) * MDP(26) + (-pkin(5) * t311 - t306 * t369 + (-pkin(10) * t335 + t305 * t385 - t553) * t447 + ((-pkin(10) * qJD(6) - t358) * t385 + t466) * t451 + t498) * MDP(25); -t479 * t369 * MDP(20) + (-t369 ^ 2 + t479 ^ 2) * MDP(21) + (t513 + t548) * MDP(22) + (-t496 - t547) * MDP(23) + t335 * MDP(24) + (-t447 * t286 + t294 + t292 * t385 + t303 * t479 - g(1) * (-t328 * t447 - t350 * t451) - g(2) * (-t326 * t447 - t347 * t451) - g(3) * (-t360 * t447 + t374)) * MDP(25) + (-t451 * t286 - t447 * t295 - t485 * t385 + t303 * t369 - g(1) * (-t328 * t451 + t350 * t447) - g(2) * (-t326 * t451 + t347 * t447) - g(3) * (-t360 * t451 - t546)) * MDP(26) + (MDP(22) * t545 + MDP(23) * t479 - MDP(25) * t292 + MDP(26) * t485) * qJD(6);];
tau  = t1;
