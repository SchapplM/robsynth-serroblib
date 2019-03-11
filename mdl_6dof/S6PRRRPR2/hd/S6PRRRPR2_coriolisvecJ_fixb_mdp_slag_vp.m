% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:09:08
% EndTime: 2019-03-08 23:09:19
% DurationCPUTime: 6.16s
% Computational Cost: add. (5383->424), mult. (13477->608), div. (0->0), fcn. (10452->12), ass. (0->196)
t493 = sin(qJ(3));
t494 = sin(qJ(2));
t488 = sin(pkin(6));
t562 = qJD(1) * t488;
t546 = t494 * t562;
t590 = qJD(3) * pkin(3);
t509 = t493 * t590 - t546;
t492 = sin(qJ(4));
t496 = cos(qJ(3));
t593 = cos(qJ(4));
t540 = qJD(2) * t593;
t559 = qJD(2) * t493;
t602 = -t492 * t559 + t496 * t540;
t442 = qJD(6) - t602;
t484 = qJD(3) + qJD(4);
t516 = -t492 * t493 + t593 * t496;
t423 = t484 * t516;
t576 = t492 * t496;
t456 = t593 * t493 + t576;
t424 = t484 * t456;
t604 = pkin(4) * t424 - qJ(5) * t423 - qJD(5) * t456 + t509;
t594 = -pkin(9) - pkin(8);
t548 = qJD(3) * t594;
t458 = t493 * t548;
t459 = t496 * t548;
t497 = cos(qJ(2));
t545 = t497 * t562;
t465 = t594 * t493;
t466 = t594 * t496;
t596 = t593 * t465 + t492 * t466;
t603 = qJD(4) * t596 + t593 * t458 + t492 * t459 - t516 * t545;
t450 = -qJD(2) * t576 - t493 * t540;
t487 = sin(pkin(12));
t489 = cos(pkin(12));
t433 = t450 * t487 + t489 * t484;
t495 = cos(qJ(6));
t601 = t433 * t495;
t491 = sin(qJ(6));
t454 = t487 * t495 + t489 * t491;
t600 = t442 * t454;
t574 = t495 * t489;
t579 = t487 * t491;
t453 = -t574 + t579;
t565 = t442 * t453;
t431 = t450 * t489 - t484 * t487;
t599 = t431 * t495 - t433 * t491;
t598 = MDP(5) * t496;
t597 = (t493 ^ 2 - t496 ^ 2) * MDP(6);
t536 = -qJD(2) * t594 + t546;
t490 = cos(pkin(6));
t561 = qJD(1) * t490;
t428 = -t536 * t493 + t496 * t561;
t560 = qJD(2) * t488;
t538 = qJD(1) * t560;
t531 = t497 * t538;
t391 = qJD(3) * t428 + t496 * t531;
t429 = t493 * t561 + t536 * t496;
t392 = -qJD(3) * t429 - t493 * t531;
t420 = t428 + t590;
t539 = qJD(4) * t593;
t558 = qJD(4) * t492;
t502 = t593 * t391 + t492 * t392 + t420 * t539 - t429 * t558;
t324 = t484 * qJD(5) + t502;
t411 = t602 * t484;
t412 = t424 * qJD(2);
t552 = qJD(2) * qJD(3);
t537 = t493 * t552;
t443 = pkin(3) * t537 + t494 * t538;
t352 = pkin(4) * t412 - qJ(5) * t411 + qJD(5) * t450 + t443;
t313 = -t324 * t487 + t489 * t352;
t314 = t489 * t324 + t487 * t352;
t525 = -t313 * t487 + t314 * t489;
t570 = -t603 * t487 + t489 * t604;
t569 = t487 * t604 + t603 * t489;
t418 = t492 * t429;
t374 = t593 * t428 - t418;
t413 = -pkin(4) * t450 - qJ(5) * t602;
t401 = pkin(3) * t559 + t413;
t346 = -t374 * t487 + t489 * t401;
t470 = pkin(3) * t539 + qJD(5);
t534 = t470 * t487 + t346;
t347 = t489 * t374 + t487 * t401;
t533 = t470 * t489 - t347;
t419 = t593 * t429;
t373 = t492 * t428 + t419;
t528 = pkin(3) * t558 - t373;
t437 = t492 * t465 - t593 * t466;
t567 = t437 * qJD(4) - t456 * t545 + t492 * t458 - t593 * t459;
t595 = t493 * MDP(10) + t496 * MDP(11);
t592 = t489 * pkin(5);
t483 = t489 * pkin(10);
t591 = qJD(2) * pkin(2);
t588 = t411 * t487;
t587 = t411 * t489;
t586 = t423 * t487;
t585 = t602 * t487;
t584 = t602 * t489;
t583 = t456 * t487;
t582 = t456 * t489;
t578 = t488 * t494;
t577 = t488 * t497;
t498 = qJD(3) ^ 2;
t575 = t493 * t498;
t573 = t496 * t498;
t572 = pkin(5) * t424 - t423 * t483 + t570;
t571 = pkin(10) * t586 - t569;
t568 = pkin(5) * t586 + t567;
t372 = t492 * t420 + t419;
t366 = qJ(5) * t484 + t372;
t482 = -pkin(3) * t496 - pkin(2);
t440 = t482 * qJD(2) - t545;
t387 = -pkin(4) * t602 + qJ(5) * t450 + t440;
t338 = t489 * t366 + t487 * t387;
t371 = t593 * t420 - t418;
t351 = t489 * t371 + t487 * t413;
t555 = qJD(6) * t495;
t564 = t411 * t574 + t433 * t555;
t415 = -pkin(4) * t516 - qJ(5) * t456 + t482;
t370 = t487 * t415 + t489 * t437;
t327 = pkin(10) * t433 + t338;
t557 = qJD(6) * t327;
t556 = qJD(6) * t456;
t551 = pkin(10) * t585;
t543 = t494 * t560;
t542 = t497 * t560;
t326 = t492 * t391 - t593 * t392 + t420 * t558 + t429 * t539;
t535 = t326 * t487 - t338 * t450;
t337 = -t366 * t487 + t489 * t387;
t350 = -t371 * t487 + t489 * t413;
t369 = t489 * t415 - t437 * t487;
t311 = -pkin(10) * t588 + t314;
t322 = -pkin(5) * t602 + pkin(10) * t431 + t337;
t532 = -qJD(6) * t322 - t311;
t481 = -t593 * pkin(3) - pkin(4);
t530 = -t450 * pkin(5) - pkin(10) * t584;
t439 = pkin(5) * t585;
t529 = -t439 + t528;
t309 = t322 * t495 - t327 * t491;
t310 = t322 * t491 + t327 * t495;
t524 = -t326 * t489 + t337 * t450;
t356 = -pkin(5) * t516 - pkin(10) * t582 + t369;
t361 = -pkin(10) * t583 + t370;
t523 = t356 * t495 - t361 * t491;
t522 = t356 * t491 + t361 * t495;
t446 = t490 * t496 - t493 * t578;
t447 = t490 * t493 + t496 * t578;
t403 = t492 * t446 + t593 * t447;
t385 = -t403 * t487 - t489 * t577;
t386 = t403 * t489 - t487 * t577;
t521 = t385 * t495 - t386 * t491;
t520 = t385 * t491 + t386 * t495;
t519 = qJD(6) * t431 - t588;
t518 = t337 * t584 + t338 * t585 + t525;
t517 = t593 * t446 - t492 * t447;
t478 = pkin(3) * t492 + qJ(5);
t452 = t478 * t489 + t483;
t515 = qJD(6) * t452 + t530 + t534;
t451 = (-pkin(10) - t478) * t487;
t514 = -qJD(6) * t451 - t533 - t551;
t464 = qJ(5) * t489 + t483;
t512 = qJD(5) * t487 + qJD(6) * t464 + t350 + t530;
t463 = (-pkin(10) - qJ(5)) * t487;
t511 = -qJD(5) * t489 - qJD(6) * t463 + t351 - t551;
t510 = t440 * t450 - t326;
t320 = pkin(5) * t588 + t326;
t365 = -t484 * pkin(4) + qJD(5) - t371;
t355 = -pkin(5) * t433 + t365;
t508 = t309 * t450 + t320 * t453 + t355 * t600;
t507 = -t310 * t450 + t320 * t454 - t355 * t565;
t506 = t326 * t456 + t365 * t423 - t411 * t596;
t505 = -0.2e1 * qJD(3) * t591;
t504 = -pkin(4) * t411 - qJ(5) * t412 - (-qJD(5) + t365) * t602;
t503 = t411 * t481 - t412 * t478 - (t365 - t470) * t602;
t333 = t519 * t491 + t564;
t334 = -qJD(6) * t599 + t454 * t411;
t376 = -t431 * t491 - t601;
t501 = (-t333 * t453 - t334 * t454 + t565 * t376 + t599 * t600) * MDP(24) + (t333 * t454 + t565 * t599) * MDP(23) + (t412 * t454 - t565 * t442 - t450 * t599) * MDP(25) + (-t376 * t450 - t412 * t453 - t442 * t600) * MDP(26) + t411 * MDP(14) + (t450 ^ 2 - t602 ^ 2) * MDP(13) + (MDP(12) * t602 + t442 * MDP(27)) * t450 + (-t602 * MDP(14) + (-qJD(2) * t456 - t450) * MDP(15)) * t484;
t500 = -t440 * t602 - t502;
t499 = qJD(2) ^ 2;
t479 = -pkin(4) - t592;
t462 = t481 - t592;
t427 = -t447 * qJD(3) - t493 * t542;
t426 = t446 * qJD(3) + t496 * t542;
t407 = t453 * t456;
t406 = t454 * t456;
t397 = pkin(5) * t583 - t596;
t360 = t439 + t372;
t354 = t403 * qJD(4) + t492 * t426 - t593 * t427;
t353 = t517 * qJD(4) + t593 * t426 + t492 * t427;
t345 = t454 * t423 + t555 * t582 - t556 * t579;
t344 = -t453 * t423 - t454 * t556;
t341 = t353 * t489 + t487 * t543;
t340 = -t353 * t487 + t489 * t543;
t307 = pkin(5) * t412 - pkin(10) * t587 + t313;
t306 = t495 * t307;
t1 = [(-t340 * t602 - t354 * t433 + t385 * t412 - t517 * t588) * MDP(19) + (t341 * t602 - t354 * t431 - t386 * t412 - t517 * t587) * MDP(20) + (t340 * t431 + t341 * t433 + (-t385 * t489 - t386 * t487) * t411) * MDP(21) + (t313 * t385 + t314 * t386 - t326 * t517 + t337 * t340 + t338 * t341 + t354 * t365) * MDP(22) + ((-t520 * qJD(6) + t340 * t495 - t341 * t491) * t442 + t521 * t412 + t354 * t376 - t517 * t334) * MDP(28) + (-(t521 * qJD(6) + t340 * t491 + t341 * t495) * t442 - t520 * t412 - t354 * t599 - t517 * t333) * MDP(29) + (-MDP(17) * t354 - MDP(18) * t353) * t484 + (MDP(10) * t427 - MDP(11) * t426) * qJD(3) + ((-t412 * MDP(17) - t411 * MDP(18)) * t497 + ((-MDP(17) * t602 - t450 * MDP(18)) * t494 - t595 * t497 * qJD(3)) * qJD(2) + (-t497 * MDP(4) + (-MDP(10) * t496 + MDP(11) * t493 - MDP(3)) * t494) * t499) * t488; 0.2e1 * t537 * t598 - 0.2e1 * t552 * t597 + MDP(7) * t573 - MDP(8) * t575 + (-pkin(8) * t573 + t493 * t505) * MDP(10) + (pkin(8) * t575 + t496 * t505) * MDP(11) + (t411 * t456 - t423 * t450) * MDP(12) + (t411 * t516 - t412 * t456 + t423 * t602 + t424 * t450) * MDP(13) + (t412 * t482 + t424 * t440 - t443 * t516 - t509 * t602) * MDP(17) + (t411 * t482 + t423 * t440 + t443 * t456 - t509 * t450) * MDP(18) + (-t313 * t516 + t337 * t424 + t369 * t412 - t433 * t567 + t506 * t487 - t570 * t602) * MDP(19) + (t314 * t516 - t338 * t424 - t370 * t412 - t431 * t567 + t506 * t489 + t569 * t602) * MDP(20) + (t570 * t431 + t569 * t433 + (-t313 * t456 - t337 * t423 - t369 * t411) * t489 + (-t314 * t456 - t338 * t423 - t370 * t411) * t487) * MDP(21) + (t313 * t369 + t314 * t370 - t326 * t596 + t570 * t337 + t569 * t338 + t567 * t365) * MDP(22) + (-t333 * t407 - t344 * t599) * MDP(23) + (-t333 * t406 + t334 * t407 - t344 * t376 + t345 * t599) * MDP(24) + (-t333 * t516 + t344 * t442 - t407 * t412 - t424 * t599) * MDP(25) + (t334 * t516 - t345 * t442 - t376 * t424 - t406 * t412) * MDP(26) + (-t412 * t516 + t424 * t442) * MDP(27) + (t523 * t412 - (-t311 * t491 + t306) * t516 + t309 * t424 + t397 * t334 + t320 * t406 + t355 * t345 + (t571 * t491 + t572 * t495) * t442 + t568 * t376 + (t310 * t516 - t522 * t442) * qJD(6)) * MDP(28) + (-t522 * t412 + (t307 * t491 + t311 * t495) * t516 - t310 * t424 + t397 * t333 - t320 * t407 + t355 * t344 + (-t572 * t491 + t571 * t495) * t442 - t568 * t599 + (t309 * t516 - t523 * t442) * qJD(6)) * MDP(29) + (t423 * MDP(14) - t424 * MDP(15) - t567 * MDP(17) - MDP(18) * t603) * t484; (-t431 * t534 + t433 * t533 + t518) * MDP(21) + (t346 * t602 - t433 * t528 + t503 * t487 + t524) * MDP(19) + (-t347 * t602 - t431 * t528 + t503 * t489 + t535) * MDP(20) + (t326 * t481 - t534 * t337 + t533 * t338 + t528 * t365 + t525 * t478) * MDP(22) + (t373 * t484 + (-t484 * t558 + t559 * t602) * pkin(3) + t510) * MDP(17) + (t374 * t484 + (t450 * t559 - t484 * t539) * pkin(3) + t500) * MDP(18) + t501 + (-(t451 * t491 + t452 * t495) * t412 + t462 * t333 + (t515 * t491 + t514 * t495) * t442 - t529 * t599 + t507) * MDP(29) + ((t451 * t495 - t452 * t491) * t412 + t462 * t334 + (t514 * t491 - t515 * t495) * t442 + t529 * t376 + t508) * MDP(28) + t595 * qJD(2) * t591 + (-t493 * t598 + t597) * t499; (t372 * t484 + t510) * MDP(17) + (-t350 * t431 - t351 * t433 + (-t431 * t487 + t433 * t489) * qJD(5) + t518) * MDP(21) + (t350 * t602 + t372 * t433 + t504 * t487 + t524) * MDP(19) + (-t351 * t602 + t372 * t431 + t504 * t489 + t535) * MDP(20) + (-pkin(4) * t326 - t337 * t350 - t338 * t351 - t365 * t372 + (-t337 * t487 + t338 * t489) * qJD(5) + t525 * qJ(5)) * MDP(22) + (t371 * t484 + t500) * MDP(18) + t501 + ((t463 * t495 - t464 * t491) * t412 + t479 * t334 - t360 * t376 + (t511 * t491 - t512 * t495) * t442 + t508) * MDP(28) + (-(t463 * t491 + t464 * t495) * t412 + t479 * t333 + t360 * t599 + (t512 * t491 + t511 * t495) * t442 + t507) * MDP(29); (t431 * t602 + t588) * MDP(19) + (-t433 * t602 + t587) * MDP(20) + (-t431 ^ 2 - t433 ^ 2) * MDP(21) + (-t337 * t431 - t338 * t433 + t326) * MDP(22) + (-t442 * t599 + t334) * MDP(28) + (t442 * t601 + (t431 * t442 + t519) * t491 + t564) * MDP(29); -t376 ^ 2 * MDP(24) + (t376 * t442 + t564) * MDP(25) + t412 * MDP(27) + (t310 * t442 + t306) * MDP(28) + (t309 * t442 + t355 * t376) * MDP(29) - (MDP(23) * t376 - MDP(24) * t599 + t442 * MDP(26) - t355 * MDP(28)) * t599 + (MDP(26) * t519 - MDP(28) * t557 + MDP(29) * t532) * t495 + (t519 * MDP(25) + (-qJD(6) * t433 - t587) * MDP(26) + t532 * MDP(28) + (-t307 + t557) * MDP(29)) * t491;];
tauc  = t1;
