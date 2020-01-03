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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:15:28
% EndTime: 2020-01-03 12:15:33
% DurationCPUTime: 3.51s
% Computational Cost: add. (3626->357), mult. (5481->451), div. (0->0), fcn. (3946->16), ass. (0->197)
t483 = sin(qJ(2));
t556 = qJD(2) * t483;
t462 = pkin(1) * t556;
t488 = cos(qJ(2));
t593 = pkin(1) * t488;
t560 = -qJD(1) * t462 + qJDD(1) * t593;
t473 = qJDD(1) + qJDD(2);
t592 = pkin(2) * t473;
t396 = -t560 - t592;
t479 = qJ(1) + qJ(2);
t466 = sin(t479);
t468 = cos(t479);
t520 = -g(2) * t468 - g(3) * t466;
t603 = t396 - t520;
t475 = qJD(1) + qJD(2);
t486 = cos(qJ(4));
t487 = cos(qJ(3));
t564 = t486 * t487;
t542 = t475 * t564;
t481 = sin(qJ(4));
t482 = sin(qJ(3));
t571 = t481 * t482;
t543 = t475 * t571;
t392 = -t542 + t543;
t485 = cos(qJ(5));
t413 = t481 * t487 + t482 * t486;
t394 = t413 * t475;
t480 = sin(qJ(5));
t580 = t394 * t480;
t353 = -t392 * t485 - t580;
t547 = -qJD(4) - qJD(5);
t464 = qJD(3) - t547;
t602 = t353 * t464;
t512 = t392 * t480 - t394 * t485;
t601 = t464 * t512;
t535 = g(2) * t466 - g(3) * t468;
t412 = -t564 + t571;
t594 = pkin(1) * t483;
t546 = qJD(1) * t594;
t426 = pkin(7) * t475 + t546;
t536 = pkin(8) * t475 + t426;
t383 = t536 * t487;
t376 = t486 * t383;
t382 = t536 * t482;
t377 = qJD(3) * pkin(3) - t382;
t513 = -t377 * t481 - t376;
t587 = pkin(9) * t392;
t330 = -t513 - t587;
t458 = -pkin(3) * t487 - pkin(2);
t557 = qJD(1) * t488;
t545 = pkin(1) * t557;
t395 = t458 * t475 - t545;
t356 = pkin(4) * t392 + t395;
t478 = qJ(3) + qJ(4);
t469 = qJ(5) + t478;
t453 = sin(t469);
t454 = cos(t469);
t550 = qJD(5) * t480;
t504 = g(1) * t453 + t330 * t550 - t356 * t353 + t454 * t535;
t474 = qJD(3) + qJD(4);
t553 = qJD(3) * t487;
t537 = t475 * t553;
t336 = qJD(4) * t542 + t413 * t473 - t474 * t543 + t486 * t537;
t472 = qJDD(3) + qJDD(4);
t548 = qJDD(1) * t483;
t555 = qJD(2) * t488;
t397 = pkin(7) * t473 + (qJD(1) * t555 + t548) * pkin(1);
t575 = t473 * t482;
t340 = -t426 * t553 + qJDD(3) * pkin(3) - t397 * t482 + (-t537 - t575) * pkin(8);
t554 = qJD(3) * t482;
t538 = t475 * t554;
t574 = t473 * t487;
t341 = -t426 * t554 + t397 * t487 + (-t538 + t574) * pkin(8);
t499 = qJD(4) * t513 + t340 * t486 - t481 * t341;
t304 = pkin(4) * t472 - pkin(9) * t336 + t499;
t367 = t474 * t413;
t519 = t412 * t473;
t337 = t367 * t475 + t519;
t552 = qJD(4) * t481;
t596 = (qJD(4) * t377 + t341) * t486 + t481 * t340 - t383 * t552;
t305 = -pkin(9) * t337 + t596;
t578 = t453 * t468;
t579 = t453 * t466;
t495 = -g(1) * t454 + g(2) * t579 - g(3) * t578 + t304 * t485 - t480 * t305 + t356 * t512;
t463 = qJDD(5) + t472;
t600 = t463 * MDP(25) + t353 * MDP(21) * t512 + (-t353 ^ 2 + t512 ^ 2) * MDP(22);
t595 = -pkin(7) - pkin(8);
t540 = qJD(3) * t595;
t419 = t482 * t540;
t420 = t487 * t540;
t441 = t595 * t482;
t470 = t487 * pkin(8);
t442 = pkin(7) * t487 + t470;
t561 = t441 * t481 + t442 * t486;
t599 = -qJD(4) * t561 + t413 * t545 - t481 * t419 + t420 * t486;
t551 = qJD(4) * t486;
t598 = -t412 * t545 - t419 * t486 - t420 * t481 - t441 * t551 + t442 * t552;
t389 = t394 * pkin(9);
t374 = t481 * t383;
t531 = t377 * t486 - t374;
t329 = -t389 + t531;
t455 = pkin(7) + t594;
t582 = -pkin(8) - t455;
t409 = t582 * t482;
t410 = t455 * t487 + t470;
t562 = t409 * t481 + t410 * t486;
t490 = qJD(3) ^ 2;
t597 = pkin(7) * t490 - t592;
t533 = t336 * t480 + t337 * t485;
t310 = -qJD(5) * t512 + t533;
t591 = pkin(2) * t475;
t590 = pkin(3) * t481;
t366 = t474 * t412;
t588 = pkin(9) * t366;
t586 = pkin(9) * t413;
t327 = pkin(4) * t474 + t329;
t581 = t327 * t485;
t465 = sin(t478);
t577 = t465 * t466;
t576 = t465 * t468;
t573 = t475 * t482;
t572 = t480 * t463;
t570 = t481 * t485;
t567 = t482 * t487;
t566 = t485 * t330;
t565 = t485 * t463;
t563 = -t382 * t486 - t374;
t476 = t482 ^ 2;
t559 = -t487 ^ 2 + t476;
t549 = qJD(5) * t485;
t544 = pkin(1) * t555;
t461 = pkin(3) * t554;
t541 = t336 * t485 - t337 * t480 - t392 * t549;
t539 = t475 * t556;
t359 = pkin(4) * t367 + t461;
t534 = qJD(3) * t582;
t530 = t382 * t481 - t376;
t529 = t409 * t486 - t410 * t481;
t528 = t441 * t486 - t442 * t481;
t526 = -qJD(5) * t327 - t305;
t525 = t475 * t546;
t362 = t412 * t485 + t413 * t480;
t320 = -qJD(5) * t362 - t366 * t485 - t367 * t480;
t360 = pkin(3) * t538 + t458 * t473 - t560;
t322 = pkin(4) * t337 + t360;
t363 = -t412 * t480 + t413 * t485;
t524 = g(2) * t578 + g(3) * t579 + t320 * t356 + t322 * t363;
t523 = g(2) * t576 + g(3) * t577 + t360 * t413 - t366 * t395;
t427 = -t545 - t591;
t522 = t427 * t553 + t482 * t603;
t521 = t359 - t546;
t357 = t528 - t586;
t364 = t367 * pkin(9);
t518 = -qJD(5) * t357 + t364 + t598;
t408 = t412 * pkin(9);
t358 = -t408 + t561;
t517 = qJD(5) * t358 - t588 - t599;
t516 = -t327 * t480 - t566;
t346 = t529 - t586;
t347 = -t408 + t562;
t515 = t346 * t485 - t347 * t480;
t514 = t346 * t480 + t347 * t485;
t391 = pkin(4) * t412 + t458;
t378 = t482 * t534 + t487 * t544;
t379 = -t482 * t544 + t487 * t534;
t510 = t378 * t486 + t379 * t481 + t409 * t551 - t410 * t552;
t309 = -t394 * t550 + t541;
t508 = -t461 + t546;
t507 = -t427 * t475 - t397 + t535;
t321 = qJD(5) * t363 - t366 * t480 + t367 * t485;
t506 = t321 * t356 + t322 * t362 + t454 * t520;
t467 = cos(t478);
t505 = t360 * t412 + t367 * t395 + t467 * t520;
t503 = t520 + t525;
t457 = -pkin(2) - t593;
t502 = pkin(1) * t539 + t455 * t490 + t457 * t473;
t500 = -pkin(7) * qJDD(3) + (t545 - t591) * qJD(3);
t498 = -qJD(4) * t562 - t378 * t481 + t379 * t486;
t496 = -qJDD(3) * t455 + (t457 * t475 - t544) * qJD(3);
t494 = (-t309 * t362 - t310 * t363 + t320 * t353 + t321 * t512) * MDP(22) + (t309 * t363 - t320 * t512) * MDP(21) + (-t336 * t412 - t337 * t413 + t366 * t392 - t367 * t394) * MDP(15) + (t320 * t464 + t363 * t463) * MDP(23) + (-t321 * t464 - t362 * t463) * MDP(24) + (t336 * t413 - t366 * t394) * MDP(14) + (-t366 * t474 + t413 * t472) * MDP(16) + (-t367 * t474 - t412 * t472) * MDP(17) + 0.2e1 * (-qJD(3) * t475 * t559 + t473 * t567) * MDP(8) + (t473 * t476 + 0.2e1 * t482 * t537) * MDP(7) + (qJDD(3) * t487 - t482 * t490) * MDP(10) + (qJDD(3) * t482 + t487 * t490) * MDP(9) + t473 * MDP(4);
t493 = t394 * t392 * MDP(14) + (t309 - t602) * MDP(23) + (-t310 - t601) * MDP(24) + (t392 * t474 + t336) * MDP(16) - t519 * MDP(17) + (-t392 ^ 2 + t394 ^ 2) * MDP(15) + t472 * MDP(18) + t600;
t492 = g(1) * t465 + t395 * t392 + t467 * t535 - t596;
t491 = -g(1) * t467 + g(2) * t577 - g(3) * t576 - t395 * t394 + t499;
t489 = cos(qJ(1));
t484 = sin(qJ(1));
t456 = pkin(3) * t486 + pkin(4);
t434 = t458 - t593;
t421 = t462 + t461;
t402 = t427 * t554;
t381 = t391 - t593;
t368 = pkin(3) * t573 + pkin(4) * t394;
t355 = t359 + t462;
t332 = -t389 + t563;
t331 = t530 + t587;
t319 = t498 + t588;
t318 = -t364 + t510;
t1 = [(-g(2) * t489 - g(3) * t484) * MDP(2) + (g(2) * t484 - g(3) * t489) * MDP(3) + (t434 * t336 + t421 * t394 - t472 * t562 - t474 * t510 + t523) * MDP(20) + ((t473 * t488 - t539) * pkin(1) + t520 + t560) * MDP(5) + (t402 + t496 * t482 + (-t502 - t603) * t487) * MDP(12) + qJDD(1) * MDP(1) + (t482 * t502 + t487 * t496 + t522) * MDP(13) + (t434 * t337 + t421 * t392 + t472 * t529 + t474 * t498 + t505) * MDP(19) + (((-qJDD(1) - t473) * t483 + (-qJD(1) - t475) * t555) * pkin(1) + t535) * MDP(6) + (-t355 * t353 + t381 * t310 + (-qJD(5) * t514 - t318 * t480 + t319 * t485) * t464 + t515 * t463 + t506) * MDP(26) + (-t355 * t512 + t381 * t309 - (qJD(5) * t515 + t318 * t485 + t319 * t480) * t464 - t514 * t463 + t524) * MDP(27) + t494; ((-t548 + (-qJD(2) + t475) * t557) * pkin(1) + t535) * MDP(6) + (t458 * t337 - t508 * t392 + t528 * t472 + t474 * t599 + t505) * MDP(19) + (t391 * t309 - (t357 * t480 + t358 * t485) * t463 + (t480 * t517 + t485 * t518) * t464 - t521 * t512 + t524) * MDP(27) + (t500 * t487 + (-t525 + t597) * t482 + t522) * MDP(13) + (t391 * t310 + (t357 * t485 - t358 * t480) * t463 + (t480 * t518 - t485 * t517) * t464 - t521 * t353 + t506) * MDP(26) + (t503 + t560) * MDP(5) + (t402 + t500 * t482 + (-t396 + t503 - t597) * t487) * MDP(12) + t494 + (t458 * t336 - t508 * t394 - t561 * t472 + t474 * t598 + t523) * MDP(20); (t456 * t565 + t368 * t353 - (t331 * t485 - t332 * t480) * t464 + (-t481 * t572 + (-t480 * t486 - t570) * t464 * qJD(4)) * pkin(3) + ((-pkin(3) * t570 - t456 * t480) * t464 + t516) * qJD(5) + t495) * MDP(26) + (-t530 * t474 + (-t392 * t573 + t472 * t486 - t474 * t552) * pkin(3) + t491) * MDP(19) + (t563 * t474 + (-t394 * t573 - t472 * t481 - t474 * t551) * pkin(3) + t492) * MDP(20) + qJDD(3) * MDP(11) + (t368 * t512 + (-t456 * t463 - t304 + (-t547 * t590 + t331) * t464) * t480 + (-t463 * t590 + (-pkin(3) * t551 - qJD(5) * t456 + t332) * t464 + t526) * t485 + t504) * MDP(27) + MDP(10) * t574 + MDP(9) * t575 + (g(1) * t482 + t487 * t507) * MDP(13) + t493 + (-g(1) * t487 + t482 * t507) * MDP(12) + (-MDP(7) * t567 + MDP(8) * t559) * t475 ^ 2; (-t474 * t513 + t491) * MDP(19) + (t474 * t531 + t492) * MDP(20) + (-(-t329 * t480 - t566) * t464 + t516 * qJD(5) + (t353 * t394 - t464 * t550 + t565) * pkin(4) + t495) * MDP(26) + t493 + ((-t330 * t464 - t304) * t480 + (t329 * t464 + t526) * t485 + (t394 * t512 - t464 * t549 - t572) * pkin(4) + t504) * MDP(27); (t541 - t602) * MDP(23) + (-t533 - t601) * MDP(24) + (-t464 * t516 + t495) * MDP(26) + (-t485 * t305 - t480 * t304 + (-t330 * t480 + t581) * t464 + t504) * MDP(27) + (-MDP(23) * t580 + MDP(24) * t512 + MDP(26) * t516 - MDP(27) * t581) * qJD(5) + t600;];
tau = t1;
