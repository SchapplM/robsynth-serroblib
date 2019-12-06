% Calculate vector of inverse dynamics joint torques for
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:45
% EndTime: 2019-12-05 19:00:50
% DurationCPUTime: 3.19s
% Computational Cost: add. (3626->357), mult. (5481->453), div. (0->0), fcn. (3946->16), ass. (0->195)
t477 = qJD(1) + qJD(2);
t488 = cos(qJ(4));
t489 = cos(qJ(3));
t568 = t488 * t489;
t545 = t477 * t568;
t483 = sin(qJ(4));
t484 = sin(qJ(3));
t575 = t483 * t484;
t546 = t477 * t575;
t392 = -t545 + t546;
t487 = cos(qJ(5));
t413 = t483 * t489 + t484 * t488;
t394 = t413 * t477;
t482 = sin(qJ(5));
t585 = t394 * t482;
t353 = -t487 * t392 - t585;
t550 = -qJD(4) - qJD(5);
t466 = qJD(3) - t550;
t602 = t353 * t466;
t514 = t392 * t482 - t487 * t394;
t601 = t466 * t514;
t481 = qJ(1) + qJ(2);
t468 = sin(t481);
t470 = cos(t481);
t600 = g(2) * t470 + g(3) * t468;
t537 = -g(2) * t468 + g(3) * t470;
t412 = -t568 + t575;
t485 = sin(qJ(2));
t594 = pkin(1) * t485;
t549 = qJD(1) * t594;
t426 = pkin(7) * t477 + t549;
t538 = pkin(8) * t477 + t426;
t383 = t538 * t489;
t376 = t488 * t383;
t382 = t538 * t484;
t377 = qJD(3) * pkin(3) - t382;
t515 = -t483 * t377 - t376;
t589 = pkin(9) * t392;
t330 = -t515 - t589;
t460 = -pkin(3) * t489 - pkin(2);
t490 = cos(qJ(2));
t561 = qJD(1) * t490;
t548 = pkin(1) * t561;
t395 = t460 * t477 - t548;
t356 = pkin(4) * t392 + t395;
t480 = qJ(3) + qJ(4);
t471 = qJ(5) + t480;
t455 = sin(t471);
t554 = qJD(5) * t482;
t456 = cos(t471);
t581 = t468 * t456;
t584 = t456 * t470;
t506 = g(1) * t455 - g(2) * t581 + g(3) * t584 + t330 * t554 - t356 * t353;
t475 = qJDD(1) + qJDD(2);
t476 = qJD(3) + qJD(4);
t557 = qJD(3) * t489;
t539 = t477 * t557;
t336 = qJD(4) * t545 + t413 * t475 - t476 * t546 + t488 * t539;
t474 = qJDD(3) + qJDD(4);
t551 = qJDD(1) * t485;
t559 = qJD(2) * t490;
t397 = pkin(7) * t475 + (qJD(1) * t559 + t551) * pkin(1);
t578 = t475 * t484;
t340 = -t426 * t557 + qJDD(3) * pkin(3) - t397 * t484 + (-t539 - t578) * pkin(8);
t558 = qJD(3) * t484;
t540 = t477 * t558;
t577 = t475 * t489;
t341 = -t426 * t558 + t397 * t489 + (-t540 + t577) * pkin(8);
t501 = t515 * qJD(4) + t488 * t340 - t483 * t341;
t304 = pkin(4) * t474 - pkin(9) * t336 + t501;
t367 = t476 * t413;
t521 = t412 * t475;
t337 = t367 * t477 + t521;
t556 = qJD(4) * t483;
t596 = (qJD(4) * t377 + t341) * t488 + t483 * t340 - t383 * t556;
t305 = -pkin(9) * t337 + t596;
t497 = -g(1) * t456 + t487 * t304 - t482 * t305 + t356 * t514 + t537 * t455;
t465 = qJDD(5) + t474;
t599 = t465 * MDP(25) + t353 * t514 * MDP(21) + (-t353 ^ 2 + t514 ^ 2) * MDP(22);
t595 = -pkin(7) - pkin(8);
t542 = qJD(3) * t595;
t419 = t484 * t542;
t420 = t489 * t542;
t441 = t595 * t484;
t472 = t489 * pkin(8);
t442 = pkin(7) * t489 + t472;
t564 = t483 * t441 + t488 * t442;
t598 = -t564 * qJD(4) + t413 * t548 - t483 * t419 + t488 * t420;
t555 = qJD(4) * t488;
t597 = -t412 * t548 - t488 * t419 - t483 * t420 - t441 * t555 + t442 * t556;
t389 = t394 * pkin(9);
t374 = t483 * t383;
t533 = t488 * t377 - t374;
t329 = -t389 + t533;
t457 = pkin(7) + t594;
t586 = -pkin(8) - t457;
t409 = t586 * t484;
t410 = t457 * t489 + t472;
t565 = t483 * t409 + t488 * t410;
t535 = t336 * t482 + t487 * t337;
t310 = -t514 * qJD(5) + t535;
t593 = pkin(1) * t490;
t592 = pkin(2) * t475;
t591 = pkin(2) * t477;
t366 = t476 * t412;
t590 = pkin(9) * t366;
t588 = pkin(9) * t413;
t583 = t465 * t483;
t582 = t465 * t487;
t469 = cos(t480);
t580 = t468 * t469;
t579 = t469 * t470;
t576 = t477 * t484;
t574 = t483 * t487;
t571 = t484 * t489;
t327 = pkin(4) * t476 + t329;
t570 = t487 * t327;
t569 = t487 * t330;
t567 = -t488 * t382 - t374;
t560 = qJD(2) * t485;
t464 = pkin(1) * t560;
t563 = -qJD(1) * t464 + qJDD(1) * t593;
t396 = -t563 - t592;
t427 = -t548 - t591;
t566 = t396 * t484 + t427 * t557;
t478 = t484 ^ 2;
t562 = -t489 ^ 2 + t478;
t553 = qJD(5) * t487;
t547 = pkin(1) * t559;
t463 = pkin(3) * t558;
t544 = t487 * t336 - t482 * t337 - t392 * t553;
t543 = t427 * t558 + t600 * t489;
t541 = t477 * t560;
t359 = pkin(4) * t367 + t463;
t536 = qJD(3) * t586;
t532 = t382 * t483 - t376;
t531 = t488 * t409 - t410 * t483;
t530 = t488 * t441 - t442 * t483;
t528 = -qJD(5) * t327 - t305;
t527 = t477 * t549;
t363 = -t412 * t482 + t413 * t487;
t321 = t363 * qJD(5) - t366 * t482 + t487 * t367;
t360 = pkin(3) * t540 + t460 * t475 - t563;
t322 = pkin(4) * t337 + t360;
t362 = t487 * t412 + t413 * t482;
t526 = g(2) * t584 + g(3) * t581 + t356 * t321 + t322 * t362;
t525 = g(2) * t579 + g(3) * t580 + t360 * t412 + t395 * t367;
t524 = t563 + t600;
t523 = t359 - t549;
t357 = t530 - t588;
t364 = t367 * pkin(9);
t520 = -qJD(5) * t357 + t364 + t597;
t408 = t412 * pkin(9);
t358 = -t408 + t564;
t519 = qJD(5) * t358 - t590 - t598;
t518 = -t482 * t327 - t569;
t346 = t531 - t588;
t347 = -t408 + t565;
t517 = t346 * t487 - t347 * t482;
t516 = t346 * t482 + t347 * t487;
t391 = pkin(4) * t412 + t460;
t378 = t484 * t536 + t489 * t547;
t379 = -t484 * t547 + t489 * t536;
t512 = t488 * t378 + t483 * t379 + t409 * t555 - t410 * t556;
t309 = -t394 * t554 + t544;
t510 = -t463 + t549;
t509 = -t427 * t477 - t397 + t537;
t320 = -t362 * qJD(5) - t366 * t487 - t367 * t482;
t508 = t356 * t320 + t322 * t363 - t455 * t600;
t467 = sin(t480);
t507 = t360 * t413 - t395 * t366 - t467 * t600;
t492 = qJD(3) ^ 2;
t505 = -pkin(7) * t492 + t527 + t592;
t459 = -pkin(2) - t593;
t504 = -pkin(1) * t541 - t457 * t492 - t459 * t475;
t502 = -pkin(7) * qJDD(3) + (t548 - t591) * qJD(3);
t500 = -t565 * qJD(4) - t378 * t483 + t488 * t379;
t498 = -qJDD(3) * t457 + (t459 * t477 - t547) * qJD(3);
t496 = (-t309 * t362 - t310 * t363 + t320 * t353 + t321 * t514) * MDP(22) + (t309 * t363 - t320 * t514) * MDP(21) + (-t336 * t412 - t337 * t413 + t366 * t392 - t367 * t394) * MDP(15) + (t320 * t466 + t363 * t465) * MDP(23) + (-t321 * t466 - t362 * t465) * MDP(24) + (t336 * t413 - t366 * t394) * MDP(14) + (-t366 * t476 + t413 * t474) * MDP(16) + (-t367 * t476 - t412 * t474) * MDP(17) + 0.2e1 * (-t562 * t477 * qJD(3) + t475 * t571) * MDP(8) + (t475 * t478 + 0.2e1 * t484 * t539) * MDP(7) + (qJDD(3) * t489 - t484 * t492) * MDP(10) + (qJDD(3) * t484 + t489 * t492) * MDP(9) + t475 * MDP(4);
t495 = t394 * t392 * MDP(14) + (t309 - t602) * MDP(23) + (-t310 - t601) * MDP(24) + (t392 * t476 + t336) * MDP(16) - t521 * MDP(17) + (-t392 ^ 2 + t394 ^ 2) * MDP(15) + t474 * MDP(18) + t599;
t494 = g(1) * t467 - g(2) * t580 + g(3) * t579 + t395 * t392 - t596;
t493 = -g(1) * t469 - t395 * t394 + t537 * t467 + t501;
t491 = cos(qJ(1));
t486 = sin(qJ(1));
t458 = pkin(3) * t488 + pkin(4);
t434 = t460 - t593;
t421 = t464 + t463;
t381 = t391 - t593;
t368 = pkin(3) * t576 + pkin(4) * t394;
t355 = t359 + t464;
t332 = -t389 + t567;
t331 = t532 + t589;
t319 = t500 + t590;
t318 = -t364 + t512;
t1 = [(t498 * t489 + (-t504 - t600) * t484 + t566) * MDP(13) + (t498 * t484 + (-t396 + t504) * t489 + t543) * MDP(12) + (((-qJDD(1) - t475) * t485 + (-qJD(1) - t477) * t559) * pkin(1) + t537) * MDP(6) + ((t475 * t490 - t541) * pkin(1) + t524) * MDP(5) + (-t355 * t353 + t381 * t310 + (-t516 * qJD(5) - t318 * t482 + t319 * t487) * t466 + t517 * t465 + t526) * MDP(26) + (-t355 * t514 + t381 * t309 - (t517 * qJD(5) + t318 * t487 + t319 * t482) * t466 - t516 * t465 + t508) * MDP(27) + (t434 * t337 + t421 * t392 + t531 * t474 + t500 * t476 + t525) * MDP(19) + t496 + (t434 * t336 + t421 * t394 - t565 * t474 - t512 * t476 + t507) * MDP(20) + qJDD(1) * MDP(1) + (g(2) * t491 + g(3) * t486) * MDP(2) + (-g(2) * t486 + g(3) * t491) * MDP(3); (t502 * t489 + (-t505 - t600) * t484 + t566) * MDP(13) + (t502 * t484 + (-t396 + t505) * t489 + t543) * MDP(12) + ((-t551 + (-qJD(2) + t477) * t561) * pkin(1) + t537) * MDP(6) + t496 + (t391 * t310 + (t357 * t487 - t358 * t482) * t465 + (t520 * t482 - t519 * t487) * t466 - t523 * t353 + t526) * MDP(26) + (t391 * t309 - (t357 * t482 + t358 * t487) * t465 + (t519 * t482 + t520 * t487) * t466 - t523 * t514 + t508) * MDP(27) + (t460 * t336 - t510 * t394 - t564 * t474 + t597 * t476 + t507) * MDP(20) + (t460 * t337 - t510 * t392 + t530 * t474 + t598 * t476 + t525) * MDP(19) + (t524 + t527) * MDP(5); (-g(1) * t489 + t509 * t484) * MDP(12) + (g(1) * t484 + t509 * t489) * MDP(13) + t495 + (-t532 * t476 + (-t392 * t576 + t474 * t488 - t476 * t556) * pkin(3) + t493) * MDP(19) + (t567 * t476 + (-t394 * t576 - t474 * t483 - t476 * t555) * pkin(3) + t494) * MDP(20) + qJDD(3) * MDP(11) + (t368 * t514 + (-t458 * t465 - t304 + (-pkin(3) * t483 * t550 + t331) * t466) * t482 + (-pkin(3) * t583 + (-pkin(3) * t555 - qJD(5) * t458 + t332) * t466 + t528) * t487 + t506) * MDP(27) + MDP(10) * t577 + MDP(9) * t578 + (t458 * t582 + t368 * t353 - (t331 * t487 - t332 * t482) * t466 + (-t482 * t583 + (-t482 * t488 - t574) * t466 * qJD(4)) * pkin(3) + ((-pkin(3) * t574 - t482 * t458) * t466 + t518) * qJD(5) + t497) * MDP(26) + (-MDP(7) * t571 + t562 * MDP(8)) * t477 ^ 2; (-(-t329 * t482 - t569) * t466 + t518 * qJD(5) + (t353 * t394 - t466 * t554 + t582) * pkin(4) + t497) * MDP(26) + t495 + (-t515 * t476 + t493) * MDP(19) + (t533 * t476 + t494) * MDP(20) + ((-t330 * t466 - t304) * t482 + (t329 * t466 + t528) * t487 + (t394 * t514 - t465 * t482 - t466 * t553) * pkin(4) + t506) * MDP(27); (t544 - t602) * MDP(23) + (-t535 - t601) * MDP(24) + (-t518 * t466 + t497) * MDP(26) + (-t487 * t305 - t482 * t304 + (-t330 * t482 + t570) * t466 + t506) * MDP(27) + (-MDP(23) * t585 + t514 * MDP(24) + t518 * MDP(26) - MDP(27) * t570) * qJD(5) + t599;];
tau = t1;
