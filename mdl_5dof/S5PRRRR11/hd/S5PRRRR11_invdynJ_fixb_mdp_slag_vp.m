% Calculate vector of inverse dynamics joint torques for
% S5PRRRR11
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR11_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR11_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:46:00
% EndTime: 2024-09-27 21:46:04
% DurationCPUTime: 1.70s
% Computational Cost: add. (2835->376), mult. (7428->509), div. (0->0), fcn. (5857->22), ass. (0->204)
t461 = cos(qJ(3));
t526 = qJD(2) * qJD(3);
t512 = t461 * t526;
t458 = sin(qJ(3));
t524 = qJDD(2) * t458;
t576 = t512 + t524;
t454 = sin(pkin(5));
t457 = sin(qJ(4));
t460 = cos(qJ(4));
t478 = t457 * t461 + t458 * t460;
t568 = qJD(3) + qJD(4);
t465 = t568 * t478;
t523 = qJDD(2) * t461;
t510 = t454 * t523;
t323 = (qJD(2) * t465 + t457 * t524) * t454 - t460 * t510;
t533 = qJD(2) * t461;
t514 = t454 * t533;
t535 = qJD(2) * t454;
t515 = t458 * t535;
t575 = -t457 * t515 + t460 * t514;
t459 = cos(qJ(5));
t379 = -t457 * t514 - t460 * t515;
t456 = sin(qJ(5));
t560 = t379 * t456;
t331 = t459 * t575 + t560;
t455 = cos(pkin(5));
t534 = qJD(2) * t455;
t432 = qJD(3) + t534;
t422 = qJD(4) + t432;
t411 = qJD(5) + t422;
t573 = t331 * t411;
t479 = t459 * t379 - t456 * t575;
t572 = t411 * t479;
t386 = (-t457 * t458 + t460 * t461) * t454;
t571 = t576 * t455;
t567 = pkin(7) + pkin(8);
t518 = t567 * t458;
t494 = t454 * t518;
t556 = t454 * t461;
t427 = qJD(1) * t556;
t554 = t455 * t461;
t435 = pkin(2) * t554;
t540 = qJD(2) * t435 + t427;
t350 = -qJD(2) * t494 + t540;
t339 = pkin(3) * t432 + t350;
t555 = t455 * t458;
t434 = pkin(2) * t555;
t536 = qJD(1) * t458;
t351 = t454 * t536 + (t556 * t567 + t434) * qJD(2);
t345 = t460 * t351;
t482 = -t339 * t457 - t345;
t564 = pkin(9) * t575;
t305 = -t482 + t564;
t529 = qJD(5) * t456;
t303 = t305 * t529;
t396 = (-pkin(3) * t461 - pkin(2)) * t454;
t442 = t455 * qJD(1);
t380 = qJD(2) * t396 + t442;
t338 = -pkin(4) * t575 + t380;
t453 = qJ(3) + qJ(4);
t506 = pkin(5) - t453;
t488 = -qJ(5) + t506;
t475 = sin(t488);
t445 = pkin(5) + t453;
t437 = qJ(5) + t445;
t520 = sin(t437) / 0.2e1;
t389 = t520 - t475 / 0.2e1;
t448 = qJ(5) + t453;
t439 = cos(t448);
t450 = pkin(10) + qJ(2);
t443 = sin(t450);
t444 = cos(t450);
t347 = -t389 * t444 - t439 * t443;
t349 = t389 * t443 - t439 * t444;
t476 = cos(t488);
t521 = cos(t437) / 0.2e1;
t469 = -t338 * t331 - g(1) * t349 - g(2) * t347 - g(3) * (t521 - t476 / 0.2e1) + t303;
t438 = sin(t448);
t471 = t476 / 0.2e1 + t521;
t346 = t438 * t443 - t444 * t471;
t348 = -t444 * t438 - t443 * t471;
t511 = t454 * t524;
t322 = t457 * t510 + t460 * t511 + t575 * t568;
t428 = qJDD(2) * t455 + qJDD(3);
t414 = qJDD(4) + t428;
t513 = t458 * t526;
t490 = t455 * t513;
t517 = t567 * t461;
t508 = t455 * t523;
t525 = qJDD(1) * t454;
t541 = pkin(2) * t508 + t461 * t525;
t313 = -pkin(2) * t490 + pkin(3) * t428 + (-qJDD(2) * t518 + (-qJD(2) * t517 - t536) * qJD(3)) * t454 + t541;
t485 = t571 * pkin(2) + pkin(7) * t510 + qJD(3) * t427 + t458 * t525;
t317 = (-pkin(7) * t513 + (-t513 + t523) * pkin(8)) * t454 + t485;
t504 = t460 * t313 - t457 * t317;
t468 = t482 * qJD(4) + t504;
t289 = pkin(4) * t414 - pkin(9) * t322 + t468;
t530 = qJD(4) * t460;
t531 = qJD(4) * t457;
t493 = -t457 * t313 - t460 * t317 - t339 * t530 + t351 * t531;
t290 = -pkin(9) * t323 - t493;
t505 = t459 * t289 - t456 * t290;
t467 = t338 * t479 - g(1) * t348 + g(2) * t346 - g(3) * (t520 + t475 / 0.2e1) + t505;
t404 = qJDD(5) + t414;
t395 = t404 * MDP(23);
t570 = t395 + t331 * MDP(19) * t479 + (-t331 ^ 2 + t479 ^ 2) * MDP(20);
t366 = pkin(3) * t455 + t435 - t494;
t433 = pkin(7) * t556;
t539 = t433 + t434;
t375 = pkin(8) * t556 + t539;
t543 = t457 * t366 + t460 * t375;
t373 = t379 * pkin(9);
t343 = t457 * t351;
t502 = t460 * t339 - t343;
t304 = t373 + t502;
t503 = t322 * t456 + t459 * t323;
t296 = -qJD(5) * t479 + t503;
t566 = pkin(3) * t457;
t565 = pkin(7) * t458;
t563 = pkin(2) * qJDD(2);
t562 = t296 * t455;
t561 = t323 * t455;
t399 = -pkin(2) * t535 + t442;
t559 = t399 * t461;
t558 = t428 * MDP(9);
t553 = t456 * t404;
t551 = t457 * t459;
t302 = pkin(4) * t422 + t304;
t550 = t459 * t302;
t549 = t459 * t305;
t548 = t459 * t404;
t387 = t478 * t454;
t337 = t386 * t456 + t387 * t459;
t340 = t568 * t386;
t341 = t465 * t454;
t301 = qJD(5) * t337 + t340 * t456 + t459 * t341;
t336 = -t459 * t386 + t387 * t456;
t546 = -t301 * t411 - t336 * t404;
t545 = -t341 * t422 + t386 * t414;
t544 = t460 * t350 - t343;
t542 = t571 * t454;
t451 = t458 ^ 2;
t538 = -t461 ^ 2 + t451;
t532 = qJD(3) * t458;
t528 = qJD(5) * t459;
t402 = t414 * MDP(16);
t527 = qJD(3) - t432;
t431 = cos(t445) / 0.2e1;
t436 = cos(t506);
t522 = t436 / 0.2e1 + t431;
t519 = sin(t445) / 0.2e1;
t516 = t459 * t322 - t456 * t323 + t528 * t575;
t501 = -t350 * t457 - t345;
t500 = t460 * t366 - t375 * t457;
t424 = qJD(3) * t435;
t367 = -qJD(3) * t494 + t424;
t368 = (-t454 * t517 - t434) * qJD(3);
t499 = -t367 * t457 + t460 * t368;
t498 = t432 + t534;
t497 = qJD(5) * t302 + t290;
t496 = t527 * qJD(2);
t495 = pkin(3) * t454 * t532;
t487 = sin(t506);
t486 = t458 * t496;
t300 = -qJD(5) * t336 + t340 * t459 - t341 * t456;
t484 = -t300 * t411 - t337 * t404;
t483 = -t456 * t302 - t549;
t481 = -t340 * t422 - t387 * t414;
t392 = t519 - t487 / 0.2e1;
t447 = cos(t453);
t357 = -t392 * t444 - t443 * t447;
t359 = t392 * t443 - t444 * t447;
t474 = qJD(3) * t432 * t461 + t428 * t458;
t473 = t366 * t530 + t460 * t367 + t457 * t368 - t375 * t531;
t295 = t379 * t529 + t516;
t441 = t455 * qJDD(1);
t355 = qJD(2) * t495 + qJDD(2) * t396 + t441;
t466 = -g(1) * t359 - g(2) * t357 - g(3) * (t431 - t436 / 0.2e1) - t380 * t575 + t493;
t464 = t379 * t575 * MDP(12) + (t295 - t573) * MDP(21) + (-t296 - t572) * MDP(22) + (-t422 * t575 + t322) * MDP(14) + (-t379 * t422 - t323) * MDP(15) + (t379 ^ 2 - t575 ^ 2) * MDP(13) + t402 + t570;
t446 = sin(t453);
t356 = t443 * t446 - t444 * t522;
t358 = -t443 * t522 - t444 * t446;
t463 = -g(1) * t358 + g(2) * t356 - g(3) * (t519 + t487 / 0.2e1) + t380 * t379 + t468;
t449 = t454 ^ 2;
t440 = pkin(3) * t460 + pkin(4);
t397 = -t454 * t563 + t441;
t394 = t428 * t556;
t384 = -t443 * t555 + t444 * t461;
t383 = -t443 * t554 - t444 * t458;
t382 = -t443 * t461 - t444 * t555;
t381 = t443 * t458 - t444 * t554;
t353 = pkin(3) * t515 - pkin(4) * t379;
t352 = -pkin(4) * t386 + t396;
t326 = pkin(4) * t341 + t495;
t318 = t322 * t455;
t316 = pkin(9) * t386 + t543;
t314 = pkin(4) * t455 - pkin(9) * t387 + t500;
t310 = t373 + t544;
t309 = t501 - t564;
t308 = pkin(4) * t323 + t355;
t298 = -pkin(9) * t340 - qJD(4) * t543 + t499;
t297 = -pkin(9) * t341 + t473;
t294 = t295 * t455;
t1 = [(qJDD(1) - g(3)) * MDP(1) + t394 * MDP(10) + t542 * MDP(11) + (t545 + t561) * MDP(17) + (t318 + t481) * MDP(18) + (t546 + t562) * MDP(24) + (t294 + t484) * MDP(25) + ((-t432 * t532 + t490 - t508) * MDP(10) - t474 * MDP(11)) * t454; t455 * t558 + t455 * t402 + t455 * t395 + ((-t454 * t565 + t435) * t428 + (-pkin(7) * t511 + t541) * t455 - g(1) * t382 - g(2) * t384 + (-t397 * t454 + t449 * t563) * t461 + (-t498 * t433 + ((t399 - t442) * t454 + (-t455 * t432 + (-t455 ^ 2 - t449) * qJD(2)) * pkin(2)) * t458) * qJD(3)) * MDP(10) + (-t424 * t432 - t539 * t428 - t485 * t455 - g(1) * t381 - g(2) * t383 + (t397 * t458 + (t498 * t565 + t559) * qJD(3)) * t454 - t576 * t449 * pkin(2)) * MDP(11) + (t322 * t387 - t340 * t379) * MDP(12) + (t322 * t386 - t323 * t387 + t340 * t575 + t341 * t379) * MDP(13) + (t318 - t481) * MDP(14) + (t545 - t561) * MDP(15) + (t499 * t422 + t500 * t414 + t504 * t455 - t575 * t495 + t396 * t323 - t355 * t386 + t380 * t341 - g(1) * t357 + g(2) * t359 + (-t422 * t543 + t455 * t482) * qJD(4)) * MDP(17) + (-g(1) * t356 - g(2) * t358 + t396 * t322 + t380 * t340 + t355 * t387 - t379 * t495 - t414 * t543 - t422 * t473 + t455 * t493) * MDP(18) + (t295 * t337 - t300 * t479) * MDP(19) + (-t295 * t336 - t296 * t337 + t300 * t331 + t301 * t479) * MDP(20) + (t294 - t484) * MDP(21) + (t546 - t562) * MDP(22) + ((-t297 * t456 + t298 * t459) * t411 + (t314 * t459 - t316 * t456) * t404 + t505 * t455 - t326 * t331 + t352 * t296 + t308 * t336 + t338 * t301 - g(1) * t347 + g(2) * t349 + ((-t314 * t456 - t316 * t459) * t411 + t483 * t455) * qJD(5)) * MDP(24) + (-g(1) * t346 - g(2) * t348 + t352 * t295 + t338 * t300 + t303 * t455 + t308 * t337 - t326 * t479 + (-(-qJD(5) * t316 + t298) * t411 - t314 * t404 - t289 * t455) * t456 + (-(qJD(5) * t314 + t297) * t411 - t316 * t404 - t497 * t455) * t459) * MDP(25) + qJDD(2) * MDP(2) + (qJDD(2) * t451 + 0.2e1 * t458 * t512) * t449 * MDP(5) + 0.2e1 * (t458 * t523 - t526 * t538) * t449 * MDP(6) + (t454 * t474 + t542) * MDP(7) + (t394 + (-t498 * t532 + t508) * t454) * MDP(8) + (g(1) * t443 - g(2) * t444) * MDP(3) + (g(1) * t444 + g(2) * t443) * MDP(4); (-t501 * t422 + (t414 * t460 - t422 * t531 + t515 * t575) * pkin(3) + t463) * MDP(17) + t464 + (t544 * t422 + (t379 * t515 - t457 * t414 - t422 * t530) * pkin(3) + t466) * MDP(18) + t558 + (t353 * t479 + (-t440 * t404 - t289 + (t309 - (-qJD(4) - qJD(5)) * t566) * t411) * t456 + (-t404 * t566 + (-pkin(3) * t530 - qJD(5) * t440 + t310) * t411 - t497) * t459 + t469) * MDP(25) + (t440 * t548 - (t309 * t459 - t310 * t456) * t411 + t353 * t331 + (-t457 * t553 + (-t456 * t460 - t551) * t411 * qJD(4)) * pkin(3) + ((-pkin(3) * t551 - t440 * t456) * t411 + t483) * qJD(5) + t467) * MDP(24) + (-t455 * pkin(2) * t486 - g(1) * t383 + g(2) * t381 + t541) * MDP(10) + (g(1) * t384 - g(2) * t382 + t540 * t432 - t485) * MDP(11) + (-t458 * t461 * MDP(5) + t538 * MDP(6)) * t449 * qJD(2) ^ 2 + ((-t486 + t523) * MDP(8) + (t527 * t533 + t524) * MDP(7) + ((-pkin(7) * t496 - g(3)) * t461 + (-qJDD(2) * pkin(7) - qJD(1) * t527 - t399 * qJD(2)) * t458) * MDP(10) + (g(3) * t458 + (t527 * t565 - t559) * qJD(2)) * MDP(11)) * t454; (-t422 * t482 + t463) * MDP(17) + t464 + (t422 * t502 + t466) * MDP(18) + ((-t305 * t411 - t289) * t456 + (t304 * t411 - t497) * t459 + (-t379 * t479 - t411 * t528 - t553) * pkin(4) + t469) * MDP(25) + (-(-t304 * t456 - t549) * t411 + t483 * qJD(5) + (-t331 * t379 - t411 * t529 + t548) * pkin(4) + t467) * MDP(24); (t516 - t573) * MDP(21) + (-t503 - t572) * MDP(22) + (-t411 * t483 + t467) * MDP(24) + (-t459 * t290 - t456 * t289 + (-t305 * t456 + t550) * t411 + t469) * MDP(25) + (MDP(21) * t560 + MDP(22) * t479 + t483 * MDP(24) - MDP(25) * t550) * qJD(5) + t570;];
tau = t1;
