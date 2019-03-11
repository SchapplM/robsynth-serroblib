% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:51:52
% EndTime: 2019-03-09 02:52:03
% DurationCPUTime: 9.28s
% Computational Cost: add. (4546->506), mult. (10542->617), div. (0->0), fcn. (8053->14), ass. (0->219)
t508 = sin(pkin(9));
t510 = cos(pkin(9));
t514 = sin(qJ(3));
t629 = cos(qJ(3));
t463 = t508 * t629 + t514 * t510;
t637 = t463 * qJD(1);
t444 = qJD(6) + t637;
t587 = qJD(3) * t514;
t571 = t508 * t587;
t573 = t629 * t510;
t561 = qJD(1) * t573;
t568 = qJDD(1) * t629;
t578 = qJDD(1) * t514;
t576 = qJD(3) * t561 + t508 * t568 + t510 * t578;
t412 = qJD(1) * t571 - t576;
t410 = -qJDD(6) + t412;
t507 = sin(pkin(10));
t509 = cos(pkin(10));
t513 = sin(qJ(6));
t516 = cos(qJ(6));
t460 = t507 * t516 + t509 * t513;
t593 = t444 * t460;
t638 = -t507 * t513 + t509 * t516;
t646 = -t638 * t410 - t444 * t593;
t598 = t508 * t514;
t461 = -t573 + t598;
t404 = t638 * t461;
t585 = qJD(6) * t516;
t586 = qJD(6) * t513;
t592 = -t507 * t586 + t509 * t585 + t638 * t637;
t557 = t460 * t410 - t444 * t592;
t572 = qJD(1) * t598;
t451 = -t561 + t572;
t422 = qJD(3) * t507 - t509 * t451;
t424 = qJD(3) * t509 + t451 * t507;
t547 = t422 * t513 - t424 * t516;
t645 = t444 * t547;
t615 = qJDD(1) * pkin(1);
t515 = sin(qJ(1));
t517 = cos(qJ(1));
t641 = g(1) * t515 - g(2) * t517;
t541 = -qJDD(2) + t615 + t641;
t503 = pkin(9) + qJ(3);
t496 = sin(t503);
t498 = cos(t503);
t624 = g(1) * t517;
t559 = g(2) * t515 + t624;
t531 = -g(3) * t496 - t498 * t559;
t579 = qJD(1) * qJD(2);
t619 = pkin(7) + qJ(2);
t632 = qJDD(1) * t619 + t579;
t435 = t632 * t508;
t436 = t632 * t510;
t470 = t619 * t508;
t464 = qJD(1) * t470;
t471 = t619 * t510;
t465 = qJD(1) * t471;
t570 = qJD(3) * t629;
t564 = t514 * t435 - t629 * t436 + t464 * t570 + t465 * t587;
t521 = t531 - t564;
t594 = -t629 * t464 - t514 * t465;
t644 = qJD(3) * t594 - t521;
t396 = (qJD(2) * t508 + qJD(3) * t471) * t514 - qJD(2) * t573 + t470 * t570;
t421 = -t514 * t470 + t471 * t629;
t539 = t641 * t496;
t643 = qJD(3) * t396 - qJDD(3) * t421 - t539;
t636 = qJD(4) - t594;
t590 = t498 * pkin(3) + t496 * qJ(4);
t504 = qJDD(3) * qJ(4);
t505 = qJD(3) * qJD(4);
t640 = -t504 - t505;
t639 = t510 * MDP(4) - t508 * MDP(5);
t635 = qJ(2) * qJDD(1);
t458 = t463 * qJD(3);
t556 = t508 * t578 - t510 * t568;
t413 = qJD(1) * t458 + t556;
t634 = -pkin(4) * t413 + qJDD(5);
t630 = t637 ^ 2;
t633 = t412 * t507 - t509 * t630;
t631 = t451 ^ 2;
t628 = pkin(3) * t413;
t626 = pkin(8) * t509;
t491 = g(3) * t498;
t621 = pkin(3) + qJ(5);
t620 = pkin(4) + t619;
t618 = -pkin(8) - t621;
t617 = qJ(4) * t451;
t616 = qJ(5) * t498;
t614 = qJDD(3) * pkin(3);
t608 = t424 * t513;
t375 = t516 * t422 + t608;
t613 = t375 * t444;
t612 = t375 * t451;
t611 = t547 * t451;
t609 = t412 * t621;
t607 = t451 * t637;
t606 = t496 * t515;
t605 = t496 * t517;
t502 = pkin(10) + qJ(6);
t497 = cos(t502);
t604 = t497 * t515;
t603 = t497 * t517;
t602 = t498 * t515;
t601 = t498 * t517;
t493 = pkin(2) * t510 + pkin(1);
t468 = -qJDD(1) * t493 + qJDD(2);
t529 = qJ(4) * t412 + t468;
t522 = -qJD(4) * t637 + t529;
t336 = qJD(5) * t451 + t413 * t621 + t522;
t535 = t435 * t629 + t514 * t436 - t464 * t587 + t465 * t570;
t533 = qJDD(4) + t535;
t348 = -t412 * pkin(4) - qJD(3) * qJD(5) - qJDD(3) * t621 + t533;
t329 = t509 * t336 + t507 * t348;
t457 = -t510 * t570 + t571;
t545 = qJ(4) * t457 - qJD(4) * t463;
t361 = qJD(5) * t461 + t458 * t621 + t545;
t397 = qJD(2) * t463 + qJD(3) * t421;
t370 = -t457 * pkin(4) + t397;
t338 = t509 * t361 + t507 * t370;
t469 = -qJD(1) * t493 + qJD(2);
t530 = -qJ(4) * t637 + t469;
t368 = t451 * t621 + t530;
t581 = pkin(4) * t637 + t636;
t374 = -qJD(3) * t621 + t581;
t347 = t509 * t368 + t507 * t374;
t383 = t621 * t637 + t617;
t417 = -t514 * t464 + t629 * t465;
t387 = -pkin(4) * t451 + t417;
t352 = t509 * t383 + t507 * t387;
t542 = -qJ(4) * t463 - t493;
t385 = t461 * t621 + t542;
t420 = t470 * t629 + t514 * t471;
t402 = t463 * pkin(4) + t420;
t355 = t509 * t385 + t507 * t402;
t589 = t508 ^ 2 + t510 ^ 2;
t584 = t417 * qJD(3);
t574 = -pkin(5) * t509 - pkin(4);
t582 = -t574 * t637 + t636;
t506 = qJD(3) * qJ(4);
t379 = qJD(5) + t387 + t506;
t580 = qJD(5) - t379;
t400 = qJDD(3) * t507 - t509 * t413;
t401 = qJDD(3) * t509 + t413 * t507;
t577 = -t513 * t400 + t516 * t401 - t422 * t585;
t575 = -g(1) * t605 - g(2) * t606 + t491;
t567 = t589 * qJD(1) ^ 2;
t328 = -t336 * t507 + t509 * t348;
t324 = -pkin(5) * t412 - pkin(8) * t401 + t328;
t327 = -pkin(8) * t400 + t329;
t566 = t516 * t324 - t327 * t513;
t346 = -t368 * t507 + t509 * t374;
t565 = t516 * t400 + t513 * t401;
t563 = 0.2e1 * t589;
t562 = g(2) * (pkin(3) * t601 + qJ(4) * t605 + t517 * t493);
t560 = g(1) * t602 - g(2) * t601;
t555 = t324 * t513 + t327 * t516;
t554 = t328 * t509 + t329 * t507;
t553 = -t328 * t507 + t329 * t509;
t332 = pkin(5) * t637 - pkin(8) * t424 + t346;
t335 = -pkin(8) * t422 + t347;
t325 = t332 * t516 - t335 * t513;
t326 = t332 * t513 + t335 * t516;
t395 = t509 * t402;
t342 = pkin(5) * t463 + t395 + (-pkin(8) * t461 - t385) * t507;
t350 = t461 * t626 + t355;
t552 = t342 * t516 - t350 * t513;
t551 = t342 * t513 + t350 * t516;
t550 = -t346 * t509 - t347 * t507;
t549 = -t346 * t507 + t347 * t509;
t548 = -t412 * t509 - t507 * t630;
t540 = -t493 - t590;
t356 = t564 + t640;
t382 = t509 * t387;
t466 = t618 * t507;
t538 = qJD(5) * t509 + qJD(6) * t466 - pkin(5) * t451 + t382 + (-pkin(8) * t637 - t383) * t507;
t467 = t618 * t509;
t537 = qJD(5) * t507 - qJD(6) * t467 + t626 * t637 + t352;
t405 = t460 * t461;
t339 = -t424 * t586 + t577;
t349 = -t356 + t634;
t527 = -t535 - t575;
t526 = t349 * t461 + t379 * t458 + t559;
t525 = -qJD(3) * t397 - qJDD(3) * t420 + t560;
t524 = t349 + t531;
t340 = -qJD(6) * t547 + t565;
t523 = t563 * t579 - t559;
t393 = pkin(3) * t451 + t530;
t520 = t393 * t637 + qJDD(4) - t527;
t495 = sin(t502);
t489 = pkin(5) * t507 + qJ(4);
t474 = qJ(4) * t601;
t472 = qJ(4) * t602;
t439 = qJD(3) * t451;
t434 = -t495 * t606 + t603;
t433 = t495 * t517 + t496 * t604;
t432 = t495 * t605 + t604;
t431 = -t495 * t515 + t496 * t603;
t411 = pkin(3) * t461 + t542;
t409 = pkin(3) * t637 + t617;
t408 = -t506 - t417;
t407 = -qJD(3) * pkin(3) + t636;
t403 = -t461 * pkin(4) + t421;
t384 = pkin(3) * t458 + t545;
t373 = t461 * t574 + t421;
t369 = -pkin(4) * t458 - t396;
t367 = t509 * t370;
t363 = pkin(5) * t422 + t379;
t362 = t458 * t574 - t396;
t360 = t533 - t614;
t359 = qJD(6) * t405 - t458 * t638;
t358 = qJD(6) * t404 + t458 * t460;
t354 = -t385 * t507 + t395;
t353 = t522 + t628;
t351 = -t383 * t507 + t382;
t337 = -t361 * t507 + t367;
t333 = pkin(5) * t400 + t349;
t331 = t458 * t626 + t338;
t330 = -pkin(5) * t457 + t367 + (-pkin(8) * t458 - t361) * t507;
t1 = [t639 * (t541 + t615) + (pkin(1) * t541 + (t589 * t635 + t523) * qJ(2)) * MDP(7) + (t563 * t635 + t523) * MDP(6) + (-(t330 * t513 + t331 * t516) * t444 + t551 * t410 - t555 * t463 + t326 * t457 - t362 * t547 + t373 * t339 + t333 * t405 + t363 * t358 + g(1) * t433 - g(2) * t431 + (-t325 * t463 - t444 * t552) * qJD(6)) * MDP(29) + (t339 * t463 + t358 * t444 - t405 * t410 + t457 * t547) * MDP(25) + (t339 * t405 - t358 * t547) * MDP(23) + (t356 * t461 + t360 * t463 + t396 * t451 + t397 * t637 - t407 * t457 + t408 * t458 - t412 * t420 - t413 * t421 - t559) * MDP(15) + (-t329 * t463 - t338 * t637 + t347 * t457 + t355 * t412 + t369 * t424 + t403 * t401 + t507 * t526 + t509 * t539) * MDP(20) + (t328 * t463 + t337 * t637 - t346 * t457 - t354 * t412 + t369 * t422 + t403 * t400 + t507 * t539 - t509 * t526) * MDP(19) + (t412 * t461 - t413 * t463 + t451 * t457 - t458 * t637) * MDP(9) + (-t412 * t463 - t457 * t637) * MDP(8) + (t353 * t411 + t393 * t384 - t356 * t421 + t408 * t396 + t360 * t420 + t407 * t397 - t619 * t624 - t562 + (-g(1) * t540 - g(2) * t619) * t515) * MDP(18) + t641 * MDP(2) + (t339 * t404 - t340 * t405 - t358 * t375 + t359 * t547) * MDP(24) + ((t330 * t516 - t331 * t513) * t444 - t552 * t410 + t566 * t463 - t325 * t457 + t362 * t375 + t373 * t340 - t333 * t404 + t363 * t359 - g(1) * t434 - g(2) * t432 + (-t326 * t463 - t444 * t551) * qJD(6)) * MDP(28) + (-t340 * t463 - t359 * t444 + t375 * t457 - t404 * t410) * MDP(26) + (t329 * t355 + t347 * t338 + t328 * t354 + t346 * t337 + t349 * t403 + t379 * t369 - t562 + (-g(1) * t620 - g(2) * t616) * t517 + (-g(1) * (t540 - t616) - g(2) * t620) * t515) * MDP(22) + t559 * MDP(3) + qJDD(1) * MDP(1) + (-t353 * t461 - t384 * t451 - t393 * t458 - t411 * t413 - t525) * MDP(16) + (-t413 * t493 + t458 * t469 + t461 * t468 + t525) * MDP(13) + (-t353 * t463 - t384 * t637 + t393 * t457 + t411 * t412 - t643) * MDP(17) + (t412 * t493 - t457 * t469 + t463 * t468 + t643) * MDP(14) + (-t337 * t424 - t338 * t422 - t354 * t401 - t355 * t400 + t458 * t549 + t461 * t553 + t560) * MDP(21) + (-t410 * t463 - t444 * t457) * MDP(27) + (-qJD(3) * t458 - qJDD(3) * t461) * MDP(11) + (-qJD(3) * t457 + qJDD(3) * t463) * MDP(10); -MDP(6) * t567 + (-qJ(2) * t567 - t541) * MDP(7) + ((-t451 - t572) * qJD(3) + t576) * MDP(14) + (-t630 - t631) * MDP(15) + (t412 + t439) * MDP(17) + (t628 - t408 * t451 + (-qJD(4) - t407) * t637 + t529 - t641) * MDP(18) + (t422 * t451 + t633) * MDP(19) + (t424 * t451 - t548) * MDP(20) + (-t400 * t509 + t401 * t507 + (t422 * t507 + t424 * t509) * t637) * MDP(21) + (t379 * t451 + t550 * t637 + t553 - t641) * MDP(22) + (t557 + t612) * MDP(28) + (-t611 - t646) * MDP(29) - t639 * qJDD(1) + (MDP(13) - MDP(16)) * (0.2e1 * t637 * qJD(3) + t556); MDP(8) * t607 + (t630 - t631) * MDP(9) + ((t451 - t572) * qJD(3) + t576) * MDP(10) - t556 * MDP(11) + qJDD(3) * MDP(12) + (-t469 * t637 + t527 + t584) * MDP(13) + (t451 * t469 + t644) * MDP(14) + (pkin(3) * t412 - qJ(4) * t413 + (-t408 - t417) * t637 + (t407 - t636) * t451) * MDP(15) + (t409 * t451 + t520 - t584 - 0.2e1 * t614) * MDP(16) + (-t393 * t451 + t409 * t637 + 0.2e1 * t504 + 0.2e1 * t505 - t644) * MDP(17) + (-t356 * qJ(4) - t360 * pkin(3) - t393 * t409 - t407 * t417 - g(1) * (-pkin(3) * t605 + t474) - g(2) * (-pkin(3) * t606 + t472) - g(3) * t590 - t636 * t408) * MDP(18) + (t509 * t609 + qJ(4) * t400 + t346 * t451 + t581 * t422 + (-t509 * t580 - t351) * t637 + t524 * t507) * MDP(19) + (-t507 * t609 + qJ(4) * t401 - t347 * t451 + t581 * t424 + (t507 * t580 + t352) * t637 + t524 * t509) * MDP(20) + (t351 * t424 + t352 * t422 + (qJD(5) * t424 - t347 * t637 + t401 * t621 - t328) * t509 + (qJD(5) * t422 + t346 * t637 + t400 * t621 - t329) * t507 - t575) * MDP(21) + (t349 * qJ(4) - t347 * t352 - t346 * t351 - g(1) * t474 - g(2) * t472 - g(3) * (t590 + t616) + t581 * t379 + t550 * qJD(5) + (t496 * t559 - t554) * t621) * MDP(22) + (t339 * t638 + t547 * t593) * MDP(23) + (-t339 * t460 - t340 * t638 + t375 * t593 + t547 * t592) * MDP(24) + (-t611 + t646) * MDP(25) + (t557 - t612) * MDP(26) + t444 * t451 * MDP(27) + (-(-t466 * t513 + t467 * t516) * t410 + t489 * t340 + t333 * t460 + t325 * t451 + (t513 * t537 - t516 * t538) * t444 + t582 * t375 + t592 * t363 + t531 * t495) * MDP(28) + ((t466 * t516 + t467 * t513) * t410 + t489 * t339 + t333 * t638 - t326 * t451 + (t513 * t538 + t516 * t537) * t444 - t582 * t547 - t593 * t363 + t531 * t497) * MDP(29); (-t412 + t439) * MDP(15) + (qJDD(3) - t607) * MDP(16) + (-qJD(3) ^ 2 - t630) * MDP(17) + (t408 * qJD(3) + t520 - t614) * MDP(18) + (-qJD(3) * t422 + t548) * MDP(19) + (-qJD(3) * t424 + t633) * MDP(20) + (-t400 * t507 - t401 * t509 + (-t422 * t509 + t424 * t507) * t637) * MDP(21) + (-qJD(3) * t379 + t549 * t637 + t554 + t575) * MDP(22) + (-qJD(3) * t375 + t646) * MDP(28) + (qJD(3) * t547 + t557) * MDP(29); (t424 * t637 + t400) * MDP(19) + (-t422 * t637 + t401) * MDP(20) + (-t422 ^ 2 - t424 ^ 2) * MDP(21) + (t346 * t424 + t347 * t422 + t521 + t634 - t640) * MDP(22) + (t340 - t645) * MDP(28) + (t339 - t613) * MDP(29); -t547 * t375 * MDP(23) + (-t375 ^ 2 + t547 ^ 2) * MDP(24) + (t577 + t613) * MDP(25) + (-t565 - t645) * MDP(26) - t410 * MDP(27) + (-g(1) * t431 - g(2) * t433 + t326 * t444 + t363 * t547 + t491 * t497 + t566) * MDP(28) + (g(1) * t432 - g(2) * t434 + t325 * t444 + t363 * t375 - t491 * t495 - t555) * MDP(29) + (-MDP(25) * t608 + MDP(26) * t547 - MDP(28) * t326 - MDP(29) * t325) * qJD(6);];
tau  = t1;
