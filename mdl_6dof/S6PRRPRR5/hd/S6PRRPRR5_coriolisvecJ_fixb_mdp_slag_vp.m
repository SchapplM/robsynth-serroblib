% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:21:06
% EndTime: 2019-03-08 22:21:18
% DurationCPUTime: 8.31s
% Computational Cost: add. (4443->489), mult. (11407->685), div. (0->0), fcn. (9030->12), ass. (0->208)
t521 = cos(qJ(3));
t581 = qJD(2) * t521;
t503 = -qJD(5) + t581;
t497 = -qJD(6) + t503;
t511 = sin(pkin(12));
t513 = cos(pkin(12));
t571 = t513 * qJD(3);
t517 = sin(qJ(3));
t582 = qJD(2) * t517;
t473 = t511 * t582 - t571;
t580 = qJD(3) * t511;
t475 = t513 * t582 + t580;
t516 = sin(qJ(5));
t520 = cos(qJ(5));
t413 = t520 * t473 + t475 * t516;
t519 = cos(qJ(6));
t412 = t473 * t516 - t475 * t520;
t515 = sin(qJ(6));
t610 = t412 * t515;
t632 = -t519 * t413 + t610;
t630 = t497 * t632;
t543 = pkin(3) * t517 - qJ(4) * t521;
t457 = qJD(3) * t543 - qJD(4) * t517;
t579 = qJD(3) * t517;
t569 = pkin(8) * t579;
t495 = t511 * t569;
t518 = sin(qJ(2));
t512 = sin(pkin(6));
t585 = qJD(1) * t512;
t522 = cos(qJ(2));
t600 = t521 * t522;
t592 = t513 * t457 - (-t511 * t600 + t513 * t518) * t585 + t495;
t636 = t511 * t457 - (t511 * t518 + t513 * t600) * t585;
t482 = t543 * qJD(2);
t567 = t518 * t585;
t486 = qJD(2) * pkin(8) + t567;
t514 = cos(pkin(6));
t584 = qJD(1) * t514;
t619 = -t517 * t486 + t521 * t584;
t396 = t511 * t482 + t513 * t619;
t564 = t511 * t581;
t382 = -pkin(9) * t564 + t396;
t635 = qJD(4) * t513 - t382;
t603 = t513 * t521;
t536 = pkin(4) * t517 - pkin(9) * t603;
t530 = t536 * qJD(3);
t634 = t530 + t592;
t604 = t513 * t517;
t607 = t511 * t521;
t633 = (-pkin(8) * t604 - pkin(9) * t607) * qJD(3) + t636;
t558 = -qJD(3) * pkin(3) + qJD(4);
t434 = -t619 + t558;
t397 = pkin(4) * t473 + t434;
t353 = pkin(5) * t413 + t397;
t631 = t353 * t632;
t629 = t412 * t503;
t628 = t413 * t503;
t539 = t412 * t519 + t413 * t515;
t627 = t497 * t539;
t601 = t520 * t513;
t609 = t511 * t516;
t480 = -t601 + t609;
t533 = t480 * t521;
t590 = qJD(2) * t533 - t480 * qJD(5);
t481 = t511 * t520 + t513 * t516;
t534 = t521 * t481;
t589 = -qJD(2) * t534 + t481 * qJD(5);
t570 = qJD(2) * qJD(3);
t560 = t521 * t570;
t547 = t511 * t560;
t574 = qJD(5) * t520;
t372 = -t473 * t574 + t560 * t601 + (-qJD(5) * t475 - t547) * t516;
t498 = t517 * t584;
t443 = t521 * t486 + t498;
t438 = qJD(3) * qJ(4) + t443;
t488 = -pkin(3) * t521 - qJ(4) * t517 - pkin(2);
t566 = t522 * t585;
t444 = qJD(2) * t488 - t566;
t376 = -t438 * t511 + t513 * t444;
t354 = -pkin(4) * t581 - pkin(9) * t475 + t376;
t377 = t513 * t438 + t511 * t444;
t362 = -pkin(9) * t473 + t377;
t336 = t354 * t516 + t362 * t520;
t583 = qJD(2) * t512;
t562 = t522 * t583;
t548 = qJD(1) * t562;
t403 = t521 * t548 + (qJD(4) + t619) * qJD(3);
t425 = (t457 + t567) * qJD(2);
t351 = -t403 * t511 + t513 * t425;
t344 = qJD(2) * t530 + t351;
t352 = t513 * t403 + t511 * t425;
t348 = -pkin(9) * t547 + t352;
t555 = t520 * t344 - t348 * t516;
t525 = -qJD(5) * t336 + t555;
t561 = t517 * t570;
t321 = pkin(5) * t561 - pkin(10) * t372 + t525;
t529 = qJD(3) * t534;
t617 = qJD(5) * t412;
t373 = qJD(2) * t529 - t617;
t576 = qJD(5) * t516;
t532 = t516 * t344 + t520 * t348 + t354 * t574 - t362 * t576;
t322 = -pkin(10) * t373 + t532;
t557 = t519 * t321 - t515 * t322;
t625 = t353 * t539 + t557;
t624 = MDP(27) * t561 + (t539 ^ 2 - t632 ^ 2) * MDP(24) + t632 * MDP(23) * t539;
t623 = MDP(5) * t517;
t622 = MDP(6) * (t517 ^ 2 - t521 ^ 2);
t621 = t516 * t633 - t520 * t634;
t472 = t513 * t488;
t421 = -pkin(9) * t604 + t472 + (-pkin(8) * t511 - pkin(4)) * t521;
t446 = pkin(8) * t603 + t511 * t488;
t608 = t511 * t517;
t433 = -pkin(9) * t608 + t446;
t620 = t421 * t574 - t433 * t576 + t516 * t634 + t520 * t633;
t593 = t516 * t421 + t520 * t433;
t616 = pkin(9) + qJ(4);
t490 = t616 * t511;
t491 = t616 * t513;
t588 = -t516 * t490 + t520 * t491;
t395 = t513 * t482 - t511 * t619;
t375 = qJD(2) * t536 + t395;
t537 = qJD(4) * t511 + qJD(5) * t491;
t618 = -t490 * t574 + t635 * t520 + (-t375 - t537) * t516;
t554 = t515 * t372 + t519 * t373;
t327 = -qJD(6) * t539 + t554;
t615 = qJD(2) * pkin(2);
t335 = t520 * t354 - t362 * t516;
t331 = pkin(10) * t412 + t335;
t329 = -pkin(5) * t503 + t331;
t614 = t329 * t519;
t332 = -pkin(10) * t413 + t336;
t613 = t332 * t519;
t578 = qJD(3) * t521;
t407 = qJD(3) * t498 + t486 * t578 + t517 * t548;
t612 = t407 * t511;
t611 = t407 * t513;
t606 = t512 * t518;
t605 = t512 * t522;
t523 = qJD(3) ^ 2;
t602 = t517 * t523;
t599 = t521 * t523;
t575 = qJD(5) * t517;
t400 = -qJD(3) * t533 - t481 * t575;
t598 = -pkin(5) * t579 + pkin(10) * t400 + qJD(5) * t593 + t621;
t401 = t574 * t604 - t575 * t609 + t529;
t597 = -pkin(10) * t401 + t620;
t417 = t519 * t480 + t481 * t515;
t596 = -qJD(6) * t417 - t515 * t589 + t519 * t590;
t418 = -t480 * t515 + t481 * t519;
t595 = qJD(6) * t418 + t515 * t590 + t519 * t589;
t549 = t513 * t569;
t591 = -t549 + t636;
t467 = (pkin(4) * t511 + pkin(8)) * t578;
t483 = pkin(4) * t608 + t517 * pkin(8);
t573 = qJD(6) * t515;
t572 = qJD(6) * t519;
t568 = t519 * t372 - t515 * t373 - t413 * t572;
t424 = pkin(4) * t564 + t443;
t505 = -pkin(4) * t513 - pkin(3);
t563 = t518 * t583;
t559 = pkin(5) * t589 - t424;
t330 = t332 * t573;
t556 = t515 * t321 - t330;
t552 = t520 * t421 - t433 * t516;
t551 = -t520 * t490 - t491 * t516;
t550 = qJD(6) * t329 + t322;
t385 = pkin(4) * t547 + t407;
t487 = -t566 - t615;
t546 = -t487 - t566;
t371 = t520 * t375;
t399 = -pkin(10) * t480 + t588;
t545 = pkin(5) * t582 + pkin(10) * t590 + t481 * qJD(4) + qJD(5) * t588 + qJD(6) * t399 - t382 * t516 + t371;
t398 = -pkin(10) * t481 + t551;
t544 = -pkin(10) * t589 + qJD(6) * t398 + t618;
t324 = t329 * t515 + t613;
t455 = t480 * t517;
t345 = -pkin(5) * t521 + pkin(10) * t455 + t552;
t454 = t481 * t517;
t347 = -pkin(10) * t454 + t593;
t542 = t345 * t515 + t347 * t519;
t464 = t514 * t517 + t521 * t606;
t429 = -t464 * t511 - t513 * t605;
t430 = t464 * t513 - t511 * t605;
t363 = t429 * t520 - t430 * t516;
t364 = t429 * t516 + t430 * t520;
t541 = t363 * t519 - t364 * t515;
t540 = t363 * t515 + t364 * t519;
t390 = t519 * t454 - t455 * t515;
t391 = -t454 * t515 - t455 * t519;
t463 = -t514 * t521 + t517 * t606;
t326 = t412 * t573 + t568;
t528 = qJD(3) * (-t546 - t615);
t526 = -qJ(4) * t579 + (-t434 + t558) * t521;
t524 = qJD(2) ^ 2;
t450 = pkin(5) * t480 + t505;
t445 = -pkin(8) * t607 + t472;
t431 = -qJD(3) * t463 + t521 * t562;
t423 = pkin(5) * t454 + t483;
t394 = t431 * t513 + t511 * t563;
t393 = -t431 * t511 + t513 * t563;
t374 = pkin(5) * t401 + t467;
t341 = pkin(5) * t373 + t385;
t340 = qJD(6) * t391 + t400 * t515 + t519 * t401;
t339 = -qJD(6) * t390 + t400 * t519 - t401 * t515;
t334 = -qJD(5) * t364 + t393 * t520 - t394 * t516;
t333 = qJD(5) * t363 + t393 * t516 + t394 * t520;
t323 = -t332 * t515 + t614;
t1 = [-t431 * qJD(3) * MDP(11) + (-t393 * t475 - t394 * t473) * MDP(14) + (t351 * t429 + t352 * t430 + t376 * t393 + t377 * t394 + t407 * t463) * MDP(15) + (-t334 * t503 + t373 * t463) * MDP(21) + (t333 * t503 + t372 * t463) * MDP(22) + (-(-qJD(6) * t540 - t333 * t515 + t334 * t519) * t497 + t463 * t327) * MDP(28) + ((qJD(6) * t541 + t333 * t519 + t334 * t515) * t497 + t463 * t326) * MDP(29) + (-MDP(4) * t522 + (-MDP(10) * t521 + MDP(11) * t517 - MDP(3)) * t518) * t524 * t512 + (-MDP(10) * qJD(3) + MDP(12) * t473 + MDP(13) * t475 + MDP(15) * t434 + MDP(21) * t413 - MDP(22) * t412 - MDP(28) * t632 - MDP(29) * t539) * (qJD(3) * t464 + t517 * t562) + ((-MDP(12) * t393 + MDP(13) * t394) * t521 + ((-MDP(11) * t605 + (-t429 * t513 - t430 * t511) * MDP(14) + (t511 * MDP(12) + t513 * MDP(13)) * t463) * t521 + (-MDP(10) * t605 + t429 * MDP(12) - t430 * MDP(13) + t363 * MDP(21) - t364 * MDP(22) + MDP(28) * t541 - MDP(29) * t540) * t517) * qJD(3)) * qJD(2); -0.2e1 * t570 * t622 + (-pkin(8) * t599 + t517 * t528) * MDP(10) + (pkin(8) * t602 + t521 * t528) * MDP(11) + ((-t473 * t566 + t612 + (qJD(2) * t445 + t376) * qJD(3)) * t517 + (-t351 + (pkin(8) * t473 + t434 * t511) * qJD(3) + (t495 - t592) * qJD(2)) * t521) * MDP(12) + ((-t475 * t566 + t611 + (-qJD(2) * t446 - t377) * qJD(3)) * t517 + (t352 + (pkin(8) * t475 + t434 * t513) * qJD(3) + (t549 + t591) * qJD(2)) * t521) * MDP(13) + ((-t351 * t513 - t352 * t511) * t517 - t592 * t475 - t591 * t473 + (-t376 * t513 - t377 * t511 + (-t445 * t513 - t446 * t511) * qJD(2)) * t578) * MDP(14) + (-t434 * t517 * t566 + t351 * t445 + t352 * t446 + t591 * t377 + t592 * t376 + (t407 * t517 + t434 * t578) * pkin(8)) * MDP(15) + (-t372 * t455 - t400 * t412) * MDP(16) + (-t372 * t454 + t373 * t455 - t400 * t413 + t401 * t412) * MDP(17) + (-t372 * t521 - t400 * t503) * MDP(18) + (t373 * t521 + t401 * t503) * MDP(19) + (-t555 * t521 + t467 * t413 + t483 * t373 + t385 * t454 + t397 * t401 + t621 * t503 + (t336 * t521 + t503 * t593) * qJD(5) + (-t413 * t566 + (qJD(2) * t552 + t335) * qJD(3)) * t517) * MDP(21) + (t532 * t521 - t467 * t412 + t483 * t372 - t385 * t455 + t397 * t400 + t620 * t503 + (t412 * t566 + (-qJD(2) * t593 - t336) * qJD(3)) * t517) * MDP(22) + (t326 * t391 - t339 * t539) * MDP(23) + (-t326 * t390 - t327 * t391 + t339 * t632 + t340 * t539) * MDP(24) + (-t326 * t521 - t339 * t497) * MDP(25) + (t327 * t521 + t340 * t497) * MDP(26) + (-t557 * t521 - t374 * t632 + t423 * t327 + t341 * t390 + t353 * t340 + (t515 * t597 + t519 * t598) * t497 + (t324 * t521 + t497 * t542) * qJD(6) + (t632 * t566 + ((t345 * t519 - t347 * t515) * qJD(2) + t323) * qJD(3)) * t517) * MDP(28) + ((t550 * t519 + t556) * t521 - t374 * t539 + t423 * t326 + t341 * t391 + t353 * t339 + ((qJD(6) * t345 + t597) * t519 + (-qJD(6) * t347 - t598) * t515) * t497 + (t539 * t566 + (-qJD(2) * t542 - t324) * qJD(3)) * t517) * MDP(29) + MDP(7) * t599 - MDP(8) * t602 + 0.2e1 * t560 * t623 + ((-qJD(2) * t455 - t412) * MDP(18) + (-qJD(2) * t454 - t413) * MDP(19) + (-t503 - t581) * MDP(20) + (qJD(2) * t391 - t539) * MDP(25) + (-qJD(2) * t390 + t632) * MDP(26) + (-t497 - t581) * MDP(27)) * t579; (qJD(3) * t443 - t407) * MDP(10) + t546 * t581 * MDP(11) + (-t611 - t443 * t473 + (-t376 * t517 + t395 * t521 + t511 * t526) * qJD(2)) * MDP(12) + (t612 - t443 * t475 + (t377 * t517 - t396 * t521 + t513 * t526) * qJD(2)) * MDP(13) + (t395 * t475 + t396 * t473 + (-qJD(4) * t473 + t376 * t581 + t352) * t513 + (qJD(4) * t475 + t377 * t581 - t351) * t511) * MDP(14) + (-pkin(3) * t407 - t376 * t395 - t377 * t396 - t434 * t443 + (-t376 * t511 + t377 * t513) * qJD(4) + (-t351 * t511 + t352 * t513) * qJ(4)) * MDP(15) + (t372 * t481 - t412 * t590) * MDP(16) + (-t372 * t480 - t373 * t481 + t412 * t589 - t413 * t590) * MDP(17) + (t505 * t373 + t385 * t480 + t589 * t397 - t424 * t413) * MDP(21) + (t505 * t372 + t385 * t481 + t590 * t397 + t412 * t424) * MDP(22) + (t326 * t418 - t539 * t596) * MDP(23) + (-t326 * t417 - t327 * t418 + t539 * t595 + t596 * t632) * MDP(24) + (t450 * t327 + t341 * t417 + t595 * t353 - t559 * t632) * MDP(28) + (t450 * t326 + t341 * t418 + t596 * t353 - t539 * t559) * MDP(29) + (-t590 * MDP(18) + t589 * MDP(19) + (t371 + t537 * t520 + (-qJD(5) * t490 + t635) * t516) * MDP(21) + t618 * MDP(22)) * t503 + (-t596 * MDP(25) + t595 * MDP(26) + (t515 * t544 + t519 * t545) * MDP(28) + (-t515 * t545 + t519 * t544) * MDP(29)) * t497 + (-t487 * MDP(10) + (qJD(3) * t481 + t412) * MDP(18) + (-qJD(3) * t480 + t413) * MDP(19) + t503 * MDP(20) + (qJD(3) * t551 - t335) * MDP(21) + (-qJD(3) * t588 + t336) * MDP(22) + (qJD(3) * t418 + t539) * MDP(25) + (-qJD(3) * t417 - t632) * MDP(26) + t497 * MDP(27) + ((t398 * t519 - t399 * t515) * qJD(3) - t323) * MDP(28) + (-(t398 * t515 + t399 * t519) * qJD(3) + t324) * MDP(29)) * t582 + (-t521 * t623 + t622) * t524; (-t473 ^ 2 - t475 ^ 2) * MDP(14) + (t376 * t475 + t377 * t473 + t407) * MDP(15) + (t373 + t629) * MDP(21) + (t372 + t628) * MDP(22) + (t327 + t627) * MDP(28) + (t326 - t630) * MDP(29) + ((-t475 + t580) * MDP(12) + (t473 + t571) * MDP(13)) * t581; -t412 * t413 * MDP(16) + (t412 ^ 2 - t413 ^ 2) * MDP(17) + (t372 - t628) * MDP(18) + (-t481 * t560 + t617 + t629) * MDP(19) + MDP(20) * t561 + (-t336 * t503 + t397 * t412 + t525) * MDP(21) + (-t335 * t503 + t397 * t413 - t532) * MDP(22) + (t326 + t630) * MDP(25) + (-t327 + t627) * MDP(26) + ((-t331 * t515 - t613) * t497 - t324 * qJD(6) + (-t412 * t632 + t497 * t573 + t519 * t561) * pkin(5) + t625) * MDP(28) + (-t631 + t330 + (t332 * t497 - t321) * t515 + (-t331 * t497 - t550) * t519 + (-t412 * t539 + t497 * t572 - t515 * t561) * pkin(5)) * MDP(29) + t624; (t568 + t630) * MDP(25) + (-t554 + t627) * MDP(26) + (-t324 * t497 + t625) * MDP(28) + (-t519 * t322 - t323 * t497 - t556 - t631) * MDP(29) + (MDP(25) * t610 + MDP(26) * t539 - MDP(28) * t324 - MDP(29) * t614) * qJD(6) + t624;];
tauc  = t1;
