% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPPRRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:34:00
% EndTime: 2019-03-09 02:34:08
% DurationCPUTime: 6.00s
% Computational Cost: add. (4802->423), mult. (10031->530), div. (0->0), fcn. (7823->14), ass. (0->192)
t459 = sin(pkin(10));
t460 = cos(pkin(10));
t468 = cos(qJ(4));
t534 = qJD(1) * t468;
t464 = sin(qJ(4));
t535 = qJD(1) * t464;
t398 = -t459 * t534 - t460 * t535;
t515 = t459 * t535;
t517 = t460 * t534;
t399 = -t515 + t517;
t463 = sin(qJ(5));
t467 = cos(qJ(5));
t363 = -t467 * t398 + t399 * t463;
t466 = cos(qJ(6));
t528 = qJD(6) * t466;
t600 = t363 * t466 + t528;
t407 = t459 * t468 + t460 * t464;
t400 = t407 * qJD(4);
t532 = qJD(4) * t468;
t516 = t460 * t532;
t533 = qJD(4) * t464;
t401 = -t459 * t533 + t516;
t548 = t460 * t468;
t408 = -t459 * t464 + t548;
t485 = t467 * t407 + t408 * t463;
t344 = qJD(5) * t485 + t467 * t400 + t401 * t463;
t451 = qJDD(4) + qJDD(5);
t455 = qJD(4) + qJD(5);
t582 = -t407 * t463 + t467 * t408;
t599 = -t344 * t455 + t582 * t451;
t486 = t398 * t463 + t467 * t399;
t422 = qJDD(1) * t548;
t523 = qJDD(1) * t459;
t491 = -t464 * t523 + t422;
t370 = -qJD(1) * t400 + t491;
t419 = qJD(4) * t515;
t571 = -t407 * qJDD(1) + t419;
t371 = qJD(1) * t516 - t571;
t530 = qJD(5) * t467;
t531 = qJD(5) * t463;
t333 = t467 * t370 - t463 * t371 + t398 * t530 - t399 * t531;
t462 = sin(qJ(6));
t518 = t466 * t333 + t462 * t451 + t455 * t528;
t529 = qJD(6) * t462;
t315 = -t486 * t529 + t518;
t314 = t315 * t466;
t357 = t455 * t462 + t466 * t486;
t431 = t466 * t451;
t316 = qJD(6) * t357 + t333 * t462 - t431;
t355 = -t466 * t455 + t462 * t486;
t598 = -t462 * t316 - t600 * t355 + t314;
t313 = t315 * t462;
t334 = qJD(5) * t486 + t370 * t463 + t467 * t371;
t332 = qJDD(6) + t334;
t329 = t462 * t332;
t555 = t363 * t455;
t557 = t486 * t455;
t559 = t357 * t486;
t590 = qJD(6) + t363;
t597 = t451 * MDP(22) + (-t334 + t557) * MDP(21) - t363 ^ 2 * MDP(19) + (MDP(18) * t363 + MDP(19) * t486 - MDP(29) * t590) * t486 + (t333 + t555) * MDP(20) + (t600 * t357 + t313) * MDP(25) + (t600 * t590 + t329 - t559) * MDP(27);
t461 = -pkin(1) - qJ(3);
t420 = qJD(1) * t461 + qJD(2);
t513 = -pkin(7) * qJD(1) + t420;
t392 = t513 * t459;
t393 = t513 * t460;
t583 = -t392 * t464 + t468 * t393;
t350 = -pkin(8) * t399 + t583;
t349 = qJD(4) * pkin(4) + t350;
t487 = -t392 * t468 - t393 * t464;
t351 = pkin(8) * t398 - t487;
t562 = t351 * t463;
t325 = t349 * t467 - t562;
t321 = -pkin(5) * t455 - t325;
t596 = t321 * t363;
t570 = -qJD(1) * qJD(3) + qJDD(1) * t461;
t411 = qJDD(2) + t570;
t509 = -pkin(7) * qJDD(1) + t411;
t378 = t509 * t459;
t379 = t509 * t460;
t505 = -t464 * t378 + t468 * t379;
t323 = qJDD(4) * pkin(4) - pkin(8) * t370 + qJD(4) * t487 + t505;
t488 = t468 * t378 + t464 * t379;
t324 = -pkin(8) * t371 + qJD(4) * t583 + t488;
t561 = t351 * t467;
t326 = t349 * t463 + t561;
t572 = qJD(5) * t326 - t467 * t323 + t324 * t463;
t305 = -pkin(5) * t451 + t572;
t454 = pkin(10) + qJ(4);
t444 = qJ(5) + t454;
t435 = cos(t444);
t465 = sin(qJ(1));
t450 = g(1) * t465;
t522 = t435 * t450;
t593 = t305 + t522;
t434 = sin(t444);
t469 = cos(qJ(1));
t449 = g(2) * t469;
t592 = g(3) * t434 + t435 * t449;
t580 = t450 - t449;
t441 = qJD(1) * qJ(2) + qJD(3);
t445 = t459 * pkin(3);
t414 = qJD(1) * t445 + t441;
t376 = -pkin(4) * t398 + t414;
t429 = g(3) * t435;
t573 = (qJD(5) * t349 + t324) * t467 + t323 * t463 - t351 * t531;
t589 = t363 * t376 + t580 * t434 + t429 - t573;
t586 = pkin(5) * t486;
t560 = t355 * t486;
t585 = t590 * (t590 * pkin(9) + t586);
t537 = t459 ^ 2 + t460 ^ 2;
t584 = t420 * t537;
t566 = -pkin(7) + t461;
t412 = t566 * t459;
t413 = t566 * t460;
t540 = t468 * t412 + t464 * t413;
t456 = qJDD(1) * qJ(2);
t457 = qJD(1) * qJD(2);
t579 = t456 + t457;
t418 = qJDD(3) + t579;
t493 = g(1) * t469 + g(2) * t465;
t581 = t418 - t493;
t578 = t459 * MDP(7) + t460 * MDP(8);
t322 = pkin(9) * t455 + t326;
t335 = pkin(5) * t363 - pkin(9) * t486 + t376;
t308 = -t322 * t462 + t335 * t466;
t577 = -t308 * t486 + t321 * t529 + t592 * t466;
t309 = t322 * t466 + t335 * t462;
t576 = t309 * t486 + t321 * t528 + t593 * t462;
t575 = -t486 * t376 - t522 - t572 + t592;
t473 = qJD(5) * t582 - t400 * t463 + t467 * t401;
t569 = -t451 * t485 - t455 * t473;
t475 = -qJD(3) * t407 - t412 * t533 + t413 * t532;
t347 = -pkin(8) * t401 + t475;
t474 = -t408 * qJD(3) - qJD(4) * t540;
t348 = pkin(8) * t400 + t474;
t503 = -t412 * t464 + t468 * t413;
t358 = -pkin(8) * t408 + t503;
t359 = -pkin(8) * t407 + t540;
t489 = t358 * t467 - t359 * t463;
t310 = qJD(5) * t489 + t347 * t467 + t348 * t463;
t338 = t358 * t463 + t359 * t467;
t433 = qJ(2) + t445;
t381 = pkin(4) * t407 + t433;
t339 = pkin(5) * t485 - pkin(9) * t582 + t381;
t304 = pkin(9) * t451 + t573;
t500 = qJD(6) * t335 + t304;
t568 = t305 * t582 - t321 * t344 - t338 * t332 - (qJD(6) * t339 + t310) * t590 - t485 * t500;
t567 = 0.2e1 * t457;
t565 = pkin(1) * qJDD(1);
t564 = t321 * t582;
t563 = t339 * t332;
t558 = t357 * t462;
t547 = t462 * t465;
t546 = t462 * t469;
t545 = t465 * t466;
t330 = t466 * t332;
t544 = t466 * t469;
t541 = -t400 * qJD(4) + t408 * qJDD(4);
t539 = t469 * pkin(1) + t465 * qJ(2);
t514 = qJDD(2) - t580;
t383 = pkin(4) * t401 + qJD(2);
t511 = t537 * MDP(9);
t510 = t537 * t411;
t502 = t590 * t462;
t403 = pkin(3) * t523 + t418;
t354 = pkin(4) * t371 + t403;
t307 = pkin(5) * t334 - pkin(9) * t333 + t354;
t499 = qJD(6) * t322 - t307;
t439 = pkin(4) * t463 + pkin(9);
t497 = pkin(4) * t399 + pkin(9) * t363 + qJD(6) * t439 + t586;
t327 = t350 * t463 + t561;
t495 = pkin(4) * t531 - t327;
t328 = t350 * t467 - t562;
t494 = -pkin(4) * t530 + t328;
t490 = -t332 * t439 + t596;
t484 = -qJD(4) * t401 - qJDD(4) * t407;
t482 = t330 - (t363 * t462 + t529) * t590;
t481 = -t344 * t466 - t529 * t582;
t479 = -pkin(9) * t332 + t325 * t590 + t596;
t470 = qJD(1) ^ 2;
t447 = t469 * qJ(2);
t443 = cos(t454);
t442 = sin(t454);
t440 = -pkin(4) * t467 - pkin(5);
t391 = t434 * t544 - t547;
t390 = t434 * t546 + t545;
t389 = t434 * t545 + t546;
t388 = -t434 * t547 + t544;
t317 = pkin(5) * t473 + pkin(9) * t344 + t383;
t311 = qJD(5) * t338 + t347 * t463 - t348 * t467;
t306 = t466 * t307;
t1 = [(-t310 * t455 + t333 * t381 - t338 * t451 - t344 * t376 + t354 * t582 + t383 * t486 - t435 * t493) * MDP(24) + (-t333 * t485 - t334 * t582 + t344 * t363 - t473 * t486) * MDP(19) + (-(-t355 * t466 - t558) * t344 + (-t313 - t316 * t466 + (t355 * t462 - t357 * t466) * qJD(6)) * t582) * MDP(26) + (t333 * t582 - t344 * t486) * MDP(18) + (t314 * t582 + t357 * t481) * MDP(25) + t599 * MDP(20) + (-t582 * t329 - t316 * t485 - t473 * t355 + (t344 * t462 - t528 * t582) * t590) * MDP(28) + (t315 * t485 + t330 * t582 + t357 * t473 + t481 * t590) * MDP(27) + (-g(1) * t391 - g(2) * t389 + t306 * t485 + t308 * t473 + t311 * t355 - t489 * t316 + (t317 * t590 + t563 + (-t322 * t485 - t338 * t590 + t564) * qJD(6)) * t466 + t568 * t462) * MDP(30) + (g(1) * t390 - g(2) * t388 - t309 * t473 + t311 * t357 - t489 * t315 + (-(-qJD(6) * t338 + t317) * t590 - t563 + t499 * t485 - qJD(6) * t564) * t462 + t568 * t466) * MDP(31) + (t332 * t485 + t473 * t590) * MDP(29) + (t418 * qJ(2) + t441 * qJD(2) - g(1) * (t461 * t465 + t447) - g(2) * (qJ(3) * t469 + t539) + t461 * t510 - qJD(3) * t584) * MDP(10) + t578 * (t579 + t581) + t569 * MDP(21) + t484 * MDP(14) + (t514 - 0.2e1 * t565) * MDP(4) + (-(qJDD(2) - t565) * pkin(1) - g(1) * (-pkin(1) * t465 + t447) - g(2) * t539 + (t456 + t567) * qJ(2)) * MDP(6) + (0.2e1 * t456 + t567 - t493) * MDP(5) + t541 * MDP(13) + (-t311 * t455 + t334 * t381 + t354 * t485 + t363 * t383 + t376 * t473 - t434 * t493 + t451 * t489) * MDP(23) + (t580 + t537 * (-t411 - t570)) * MDP(9) + t580 * MDP(2) + qJDD(1) * MDP(1) + (-qJD(2) * t398 + qJD(4) * t474 + qJDD(4) * t503 + t433 * t371 + t414 * t401 + t403 * t407 - t442 * t493) * MDP(16) + (qJD(2) * t399 - qJD(4) * t475 - qJDD(4) * t540 + t433 * t370 - t414 * t400 + t403 * t408 - t443 * t493) * MDP(17) + (-t370 * t407 - t371 * t408 - t398 * t400 - t399 * t401) * MDP(12) + (t370 * t408 - t399 * t400) * MDP(11) + t493 * MDP(3); t514 * MDP(6) + (-qJD(1) * t441 + t510 - t580) * MDP(10) + (qJD(1) * t398 + t541) * MDP(16) + (-qJD(1) * t399 + t484) * MDP(17) + (-qJD(1) * t363 + t599) * MDP(23) + (-qJD(1) * t486 + t569) * MDP(24) + (-t316 * t582 - t329 * t485 + t344 * t355) * MDP(30) + (-t315 * t582 - t330 * t485 + t344 * t357) * MDP(31) + (-MDP(6) * qJ(2) - MDP(5) - t578) * t470 + (-pkin(1) * MDP(6) + MDP(4) - t511) * qJDD(1) + ((-qJD(1) * t466 - t462 * t473 - t485 * t528) * MDP(30) + (qJD(1) * t462 - t466 * t473 + t485 * t529) * MDP(31)) * t590; (qJD(1) * t584 + t581) * MDP(10) - t419 * MDP(16) + t422 * MDP(17) + (t334 + t557) * MDP(23) + (t333 - t555) * MDP(24) + (t482 - t560) * MDP(30) + (-t466 * t590 ^ 2 - t329 - t559) * MDP(31) - t470 * t511 + ((MDP(16) * t464 + MDP(8)) * t460 + (MDP(16) * t468 - MDP(17) * t464 + MDP(7)) * t459) * qJDD(1) + ((t399 + t517) * MDP(16) + 0.2e1 * t398 * MDP(17)) * qJD(4); (t327 * t455 + (-t363 * t399 + t451 * t467 - t455 * t531) * pkin(4) + t575) * MDP(23) + ((t399 - t517) * qJD(4) + t571) * MDP(14) + ((-qJD(1) * t407 - t398) * qJD(4) + t491) * MDP(13) - t399 * t398 * MDP(11) + (t440 * t316 - t593 * t466 + t490 * t462 + t495 * t355 + (t462 * t494 - t466 * t497) * t590 + t577) * MDP(30) + (t440 * t315 + t490 * t466 - t592 * t462 + t495 * t357 + (t462 * t497 + t466 * t494) * t590 + t576) * MDP(31) + (t328 * t455 + (-t399 * t486 - t451 * t463 - t455 * t530) * pkin(4) + t589) * MDP(24) + (g(3) * t442 - t414 * t399 - t443 * t580 + t505) * MDP(16) + (g(3) * t443 - t414 * t398 + t442 * t580 - t488) * MDP(17) + (-t398 ^ 2 + t399 ^ 2) * MDP(12) + qJDD(4) * MDP(15) + (t482 + t560) * MDP(28) + (-t558 * t590 + t598) * MDP(26) + t597; (t326 * t455 + t575) * MDP(23) + (t325 * t455 + t589) * MDP(24) + (-t357 * t502 + t598) * MDP(26) + (-t502 * t590 + t330 + t560) * MDP(28) + (-pkin(5) * t316 - t326 * t355 + t479 * t462 + (-t593 - t585) * t466 + t577) * MDP(30) + (-pkin(5) * t315 - t326 * t357 + t479 * t466 + (-t592 + t585) * t462 + t576) * MDP(31) + t597; t357 * t355 * MDP(25) + (-t355 ^ 2 + t357 ^ 2) * MDP(26) + (t355 * t590 + t518) * MDP(27) + (t357 * t590 + t431) * MDP(28) + t332 * MDP(29) + (-g(1) * t388 - g(2) * t390 + t309 * t590 - t321 * t357 + t306) * MDP(30) + (g(1) * t389 - g(2) * t391 + t308 * t590 + t321 * t355) * MDP(31) + ((-t304 + t429) * MDP(31) + (-MDP(28) * t486 - MDP(30) * t322 - MDP(31) * t335) * qJD(6)) * t466 + (-qJD(6) * t486 * MDP(27) + (-qJD(6) * t455 - t333) * MDP(28) + (-t500 + t429) * MDP(30) + t499 * MDP(31)) * t462;];
tau  = t1;
