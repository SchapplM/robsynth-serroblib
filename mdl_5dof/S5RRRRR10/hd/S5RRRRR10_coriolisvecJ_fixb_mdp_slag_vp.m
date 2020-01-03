% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:36:15
% EndTime: 2019-12-31 22:36:34
% DurationCPUTime: 10.22s
% Computational Cost: add. (6820->484), mult. (17802->683), div. (0->0), fcn. (14149->10), ass. (0->204)
t505 = sin(pkin(5));
t514 = cos(qJ(2));
t591 = qJD(1) * t514;
t569 = t505 * t591;
t487 = -qJD(3) + t569;
t649 = qJD(4) - t487;
t511 = cos(qJ(5));
t582 = qJD(5) * t511;
t506 = cos(pkin(5));
t592 = qJD(1) * t506;
t495 = qJD(2) + t592;
t513 = cos(qJ(3));
t509 = sin(qJ(3));
t510 = sin(qJ(2));
t593 = qJD(1) * t505;
t570 = t510 * t593;
t546 = t509 * t570;
t445 = t495 * t513 - t546;
t446 = t495 * t509 + t513 * t570;
t508 = sin(qJ(4));
t512 = cos(qJ(4));
t403 = -t445 * t512 + t446 * t508;
t644 = t403 * t511;
t650 = t582 + t644;
t577 = pkin(1) * t592;
t460 = -pkin(7) * t570 + t514 * t577;
t526 = (pkin(2) * t510 - pkin(8) * t514) * t505;
t461 = qJD(1) * t526;
t554 = -t460 * t509 + t461 * t513;
t623 = pkin(8) + pkin(9);
t572 = qJD(3) * t623;
t606 = t513 * t514;
t648 = (pkin(3) * t510 - pkin(9) * t606) * t593 + t554 + t513 * t572;
t547 = t509 * t569;
t596 = t460 * t513 + t461 * t509;
t647 = pkin(9) * t547 - t509 * t572 - t596;
t528 = t445 * t508 + t446 * t512;
t578 = qJD(1) * qJD(2);
t563 = t505 * t578;
t544 = t514 * t563;
t586 = qJD(3) * t513;
t417 = -qJD(3) * t546 + t495 * t586 + t513 * t544;
t588 = qJD(2) * t514;
t566 = t509 * t588;
t587 = qJD(3) * t509;
t418 = (t510 * t586 + t566) * t593 + t495 * t587;
t584 = qJD(4) * t512;
t585 = qJD(4) * t508;
t362 = t417 * t512 - t418 * t508 + t445 * t584 - t446 * t585;
t507 = sin(qJ(5));
t545 = t510 * t563;
t573 = t362 * t511 + t507 * t545 + t582 * t649;
t583 = qJD(5) * t507;
t336 = -t528 * t583 + t573;
t334 = t336 * t507;
t386 = t507 * t649 + t511 * t528;
t559 = t362 * t507 - t511 * t545;
t337 = qJD(5) * t386 + t559;
t363 = qJD(4) * t528 + t417 * t508 + t418 * t512;
t358 = t511 * t363;
t616 = t528 * t507;
t384 = -t511 * t649 + t616;
t579 = -qJD(5) - t403;
t356 = t507 * t363;
t632 = t579 * t582 - t356;
t643 = t579 * t507;
t646 = MDP(22) * t545 - t363 * MDP(21) - t403 ^ 2 * MDP(19) + (t403 * t649 + t362) * MDP(20) + (-t579 * t644 - t632) * MDP(27) + (MDP(18) * t403 + MDP(19) * t528 + MDP(21) * t649 - MDP(27) * t386 + MDP(29) * t579) * t528 + (t386 * t650 + t334) * MDP(25) + (t384 * t528 - t579 * t643 + t358) * MDP(28) + (t336 * t511 - t507 * t337 - t384 * t650 + t386 * t643) * MDP(26);
t492 = t510 * t577;
t463 = pkin(7) * t569 + t492;
t429 = pkin(8) * t495 + t463;
t458 = (-pkin(2) * t514 - pkin(8) * t510 - pkin(1)) * t505;
t441 = qJD(1) * t458;
t394 = -t429 * t509 + t441 * t513;
t376 = -pkin(9) * t446 + t394;
t373 = -pkin(3) * t487 + t376;
t395 = t429 * t513 + t441 * t509;
t377 = pkin(9) * t445 + t395;
t620 = t377 * t508;
t342 = t373 * t512 - t620;
t340 = -pkin(4) * t649 - t342;
t645 = t340 * t403;
t470 = t508 * t509 - t512 * t513;
t642 = t649 * t470;
t471 = t508 * t513 + t509 * t512;
t598 = t649 * t471;
t371 = pkin(4) * t528 + pkin(10) * t403;
t428 = -pkin(2) * t495 - t460;
t407 = -pkin(3) * t445 + t428;
t462 = qJD(2) * t526;
t453 = qJD(1) * t462;
t610 = t505 * t510;
t496 = pkin(7) * t610;
t621 = pkin(1) * t514;
t464 = (t506 * t621 - t496) * qJD(2);
t454 = qJD(1) * t464;
t517 = -qJD(3) * t395 + t453 * t513 - t454 * t509;
t350 = pkin(3) * t545 - pkin(9) * t417 + t517;
t525 = -t429 * t587 + t441 * t586 + t453 * t509 + t454 * t513;
t355 = -pkin(9) * t418 + t525;
t548 = -t350 * t508 - t355 * t512 - t373 * t584 + t377 * t585;
t641 = t403 * t407 + t548;
t502 = t505 ^ 2;
t638 = -0.2e1 * t502 * t578;
t636 = (t510 ^ 2 - t514 ^ 2) * MDP(5);
t467 = t506 * t509 + t513 * t610;
t609 = t505 * t514;
t622 = pkin(1) * t510;
t457 = pkin(7) * t609 + (pkin(8) + t622) * t506;
t555 = -t457 * t509 + t458 * t513;
t382 = -pkin(3) * t609 - pkin(9) * t467 + t555;
t466 = -t506 * t513 + t509 * t610;
t597 = t457 * t513 + t458 * t509;
t389 = -pkin(9) * t466 + t597;
t635 = t382 * t508 + t389 * t512;
t488 = t623 * t509;
t489 = t623 * t513;
t527 = -t488 * t512 - t489 * t508;
t634 = qJD(4) * t527 - t508 * t648 + t512 * t647;
t434 = -t488 * t508 + t489 * t512;
t633 = qJD(4) * t434 + t508 * t647 + t512 * t648;
t543 = -t463 + (-t547 + t587) * pkin(3);
t619 = t377 * t512;
t343 = t373 * t508 + t619;
t518 = -qJD(4) * t343 + t350 * t512 - t355 * t508;
t325 = -pkin(4) * t545 - t518;
t341 = pkin(10) * t649 + t343;
t351 = pkin(4) * t403 - pkin(10) * t528 + t407;
t535 = t341 * t507 - t351 * t511;
t630 = -t325 * t511 + t340 * t583 + t528 * t535;
t323 = t325 * t507;
t329 = t341 * t511 + t351 * t507;
t629 = t329 * t528 + t340 * t582 + t323;
t626 = -t407 * t528 + t518;
t567 = t505 * t588;
t424 = -qJD(3) * t466 + t513 * t567;
t519 = -qJD(3) * t597 + t462 * t513 - t464 * t509;
t590 = qJD(2) * t510;
t568 = t505 * t590;
t364 = pkin(3) * t568 - pkin(9) * t424 + t519;
t423 = qJD(3) * t467 + t505 * t566;
t524 = -t457 * t587 + t458 * t586 + t462 * t509 + t464 * t513;
t367 = -pkin(9) * t423 + t524;
t625 = -qJD(4) * t635 + t364 * t512 - t367 * t508;
t615 = t445 * t487;
t614 = t446 * t487;
t613 = t471 * t511;
t612 = t487 * t513;
t515 = qJD(1) ^ 2;
t611 = t502 * t515;
t608 = t509 * t487;
t600 = pkin(4) * t570 + t633;
t455 = pkin(7) * t544 + qJD(2) * t492;
t465 = pkin(1) * t506 * t590 + pkin(7) * t567;
t589 = qJD(2) * t513;
t580 = qJD(2) - t495;
t575 = t507 * t609;
t501 = -pkin(3) * t513 - pkin(2);
t571 = t502 * t591;
t324 = pkin(10) * t545 - t548;
t393 = pkin(3) * t418 + t455;
t331 = pkin(4) * t363 - pkin(10) * t362 + t393;
t562 = -t324 * t507 + t331 * t511;
t558 = t507 * t642 - t511 * t570;
t557 = t507 * t570 + t511 * t642;
t499 = pkin(3) * t508 + pkin(10);
t550 = pkin(3) * t446 + qJD(5) * t499 + t371;
t549 = t502 * t510 * t514 * MDP(4);
t408 = pkin(3) * t423 + t465;
t346 = t376 * t508 + t619;
t542 = pkin(3) * t585 - t346;
t347 = t376 * t512 - t620;
t541 = -pkin(3) * t584 + t347;
t540 = pkin(1) * t638;
t419 = pkin(4) * t470 - pkin(10) * t471 + t501;
t539 = pkin(10) * t570 - qJD(5) * t419 - t634;
t538 = -pkin(4) * t598 - pkin(10) * t642 + qJD(5) * t434 - t543;
t537 = t324 * t511 + t331 * t507;
t536 = -t363 * t499 + t645;
t353 = -pkin(10) * t609 + t635;
t414 = t466 * t512 + t467 * t508;
t415 = -t466 * t508 + t467 * t512;
t456 = t496 + (-pkin(2) - t621) * t506;
t416 = pkin(3) * t466 + t456;
t368 = pkin(4) * t414 - pkin(10) * t415 + t416;
t534 = t353 * t511 + t368 * t507;
t533 = -t353 * t507 + t368 * t511;
t531 = t382 * t512 - t389 * t508;
t399 = t415 * t507 + t511 * t609;
t523 = t364 * t508 + t367 * t512 + t382 * t584 - t389 * t585;
t522 = t471 * t582 - t558;
t521 = -t471 * t583 - t557;
t500 = -pkin(3) * t512 - pkin(4);
t400 = t415 * t511 - t575;
t370 = qJD(4) * t415 + t423 * t512 + t424 * t508;
t369 = -qJD(4) * t414 - t423 * t508 + t424 * t512;
t352 = pkin(4) * t609 - t531;
t345 = -qJD(5) * t575 + t369 * t507 + t415 * t582 - t511 * t568;
t344 = -qJD(5) * t399 + t369 * t511 + t507 * t568;
t332 = pkin(4) * t370 - pkin(10) * t369 + t408;
t327 = -pkin(4) * t568 - t625;
t326 = pkin(10) * t568 + t523;
t322 = -qJD(5) * t329 + t562;
t321 = -qJD(5) * t535 + t537;
t1 = [(-t370 * t649 + (t363 * t514 + (-qJD(1) * t414 - t403) * t590) * t505) * MDP(21) + (t625 * t649 + t408 * t403 + t416 * t363 + t393 * t414 + t407 * t370 + (-t518 * t514 + (qJD(1) * t531 + t342) * t590) * t505) * MDP(23) + (t369 * t649 + (-t362 * t514 + (qJD(1) * t415 + t528) * t590) * t505) * MDP(20) + (-t523 * t649 + t408 * t528 + t416 * t362 + t393 * t415 + t407 * t369 + (-t548 * t514 + (-qJD(1) * t635 - t343) * t590) * t505) * MDP(24) + (t524 * t487 + t465 * t446 + t456 * t417 + t455 * t467 + t428 * t424 + (t525 * t514 + (-qJD(1) * t597 - t395) * t590) * t505) * MDP(17) + (-t519 * t487 - t465 * t445 + t456 * t418 + t455 * t466 + t428 * t423 + (-t517 * t514 + (qJD(1) * t555 + t394) * t590) * t505) * MDP(16) + (-t424 * t487 + (-t417 * t514 + (qJD(1) * t467 + t446) * t590) * t505) * MDP(13) + (t423 * t487 + (t418 * t514 + (-qJD(1) * t466 + t445) * t590) * t505) * MDP(14) + 0.2e1 * t549 * t578 + (-t454 * t506 - t464 * t495 + t514 * t540) * MDP(10) + (-t455 * t506 - t465 * t495 + t510 * t540) * MDP(9) + (-t362 * t414 - t363 * t415 - t369 * t403 - t370 * t528) * MDP(19) + (t362 * t415 + t369 * t528) * MDP(18) + (MDP(6) * t567 - MDP(7) * t568) * (t495 + t592) + ((qJD(5) * t533 + t326 * t511 + t332 * t507) * t579 - t534 * t363 - t321 * t414 - t329 * t370 + t327 * t386 + t352 * t336 + t325 * t400 + t340 * t344) * MDP(31) + (-t337 * t414 + t345 * t579 - t363 * t399 - t370 * t384) * MDP(28) + (t363 * t414 - t370 * t579) * MDP(29) + (t336 * t414 - t344 * t579 + t363 * t400 + t370 * t386) * MDP(27) + (-(-qJD(5) * t534 - t326 * t507 + t332 * t511) * t579 + t533 * t363 + t322 * t414 - t535 * t370 + t327 * t384 + t352 * t337 + t325 * t399 + t340 * t345) * MDP(30) + (t505 * t649 - t571) * MDP(22) * t590 + t636 * t638 + (-t487 * t505 - t571) * MDP(15) * t590 + (t417 * t467 + t424 * t446) * MDP(11) + (-t417 * t466 - t418 * t467 - t423 * t446 + t424 * t445) * MDP(12) + (t336 * t400 + t344 * t386) * MDP(25) + (-t336 * t399 - t337 * t400 - t344 * t384 - t345 * t386) * MDP(26); t580 * MDP(6) * t569 + (pkin(7) * t545 + t460 * t495 + (-t506 * t578 + t611) * t621) * MDP(10) + (t463 * t495 + t611 * t622 - t455) * MDP(9) + (t558 * t386 + t557 * t384 + (-t334 - t337 * t511 + (t384 * t507 - t386 * t511) * qJD(5)) * t471) * MDP(26) + ((t419 * t511 - t434 * t507) * t363 + t322 * t470 - t527 * t337 + t471 * t323 - (t507 * t539 - t511 * t538) * t579 + t600 * t384 - t598 * t535 + t522 * t340) * MDP(30) + ((t417 - t615) * t513 + (-t418 + t614) * t509) * MDP(12) + (t417 * t509 - t446 * t612) * MDP(11) + (-pkin(2) * t418 - t455 * t513 + t554 * t487 + t463 * t445 + (pkin(8) * t612 + t428 * t509) * qJD(3) + (-t394 * t510 + (-pkin(8) * t590 - t428 * t514) * t509) * t593) * MDP(16) + (t336 * t613 + t386 * t521) * MDP(25) + (-(t419 * t507 + t434 * t511) * t363 - t321 * t470 - t527 * t336 + t325 * t613 - (t507 * t538 + t511 * t539) * t579 + t600 * t386 - t598 * t329 + t521 * t340) * MDP(31) + (t336 * t470 + t358 * t471 + t386 * t598 - t521 * t579) * MDP(27) + (t487 * t587 + (-t514 * t608 + (-t445 + t589) * t510) * t593) * MDP(14) + (-pkin(2) * t417 + t455 * t509 - t596 * t487 - t463 * t446 + (-pkin(8) * t608 + t428 * t513) * qJD(3) + (-t428 * t606 + (-pkin(8) * t589 + t395) * t510) * t593) * MDP(17) + (-t337 * t470 - t356 * t471 - t384 * t598 + t522 * t579) * MDP(28) + (-t487 * t586 + (t487 * t606 + (qJD(2) * t509 - t446) * t510) * t593) * MDP(13) + (t363 * t470 - t579 * t598) * MDP(29) + (t501 * t363 + t393 * t470 + t543 * t403 + t598 * t407) * MDP(23) + (t362 * t471 - t528 * t642) * MDP(18) + (-t362 * t470 - t363 * t471 + t403 * t642 - t528 * t598) * MDP(19) + (t501 * t362 + t393 * t471 - t407 * t642 + t528 * t543) * MDP(24) - t515 * t549 + t611 * t636 - (MDP(20) * t642 + MDP(21) * t598 + MDP(23) * t633 + MDP(24) * t634) * t649 + (-t580 * MDP(7) + (-qJD(2) * t470 + t403) * MDP(21) + (qJD(2) * t527 - t342) * MDP(23) + (qJD(2) * t471 - t528) * MDP(20) + (-qJD(2) * t434 + t343) * MDP(24) + t487 * MDP(15) - t649 * MDP(22)) * t570; (-t395 * t487 - t428 * t446 + t517) * MDP(16) + MDP(15) * t545 + (t417 + t615) * MDP(13) + (t346 * t649 + (-t403 * t446 + t512 * t545 - t585 * t649) * pkin(3) + t626) * MDP(23) + (-t418 - t614) * MDP(14) + (-t394 * t487 - t428 * t445 - t525) * MDP(17) + (t500 * t337 + t536 * t507 + t542 * t384 - (t507 * t541 - t511 * t550) * t579 + t630) * MDP(30) + (t500 * t336 + t536 * t511 + t542 * t386 - (t507 * t550 + t511 * t541) * t579 + t629) * MDP(31) + (t347 * t649 + (-t446 * t528 - t508 * t545 - t584 * t649) * pkin(3) + t641) * MDP(24) + (-t445 ^ 2 + t446 ^ 2) * MDP(12) - t446 * t445 * MDP(11) + t646; (t343 * t649 + t626) * MDP(23) + (t342 * t649 + t641) * MDP(24) + (-pkin(4) * t337 + (-t342 * t507 + t371 * t511) * t579 - t343 * t384 + t507 * t645 + t632 * pkin(10) + t630) * MDP(30) + (-pkin(4) * t336 - (t342 * t511 + t371 * t507) * t579 - t343 * t386 + t340 * t644 + (-t579 * t583 - t358) * pkin(10) + t629) * MDP(31) + t646; t386 * t384 * MDP(25) + (-t384 ^ 2 + t386 ^ 2) * MDP(26) + (-t384 * t579 + t573) * MDP(27) + (-t386 * t579 - t559) * MDP(28) + t363 * MDP(29) + (-t329 * t579 - t340 * t386 + t562) * MDP(30) + (t340 * t384 + t535 * t579 - t537) * MDP(31) + (-MDP(27) * t616 - MDP(28) * t386 - MDP(30) * t329 + MDP(31) * t535) * qJD(5);];
tauc = t1;
