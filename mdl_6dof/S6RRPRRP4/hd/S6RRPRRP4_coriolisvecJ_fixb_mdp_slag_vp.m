% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:51
% EndTime: 2019-03-09 11:56:02
% DurationCPUTime: 6.56s
% Computational Cost: add. (9251->507), mult. (22930->636), div. (0->0), fcn. (16992->8), ass. (0->215)
t502 = sin(pkin(10));
t505 = sin(qJ(2));
t507 = cos(qJ(2));
t598 = cos(pkin(10));
t477 = t502 * t507 + t505 * t598;
t467 = t477 * qJD(1);
t504 = sin(qJ(4));
t506 = cos(qJ(4));
t436 = qJD(2) * t506 - t504 * t467;
t437 = qJD(2) * t504 + t467 * t506;
t503 = sin(qJ(5));
t602 = cos(qJ(5));
t376 = -t602 * t436 + t437 * t503;
t522 = t503 * t436 + t437 * t602;
t628 = t376 * t522;
t537 = t598 * t507;
t489 = qJD(1) * t537;
t561 = qJD(1) * t505;
t465 = -t502 * t561 + t489;
t406 = pkin(2) * t561 + pkin(3) * t467 - pkin(8) * t465;
t399 = t506 * t406;
t600 = -qJ(3) - pkin(7);
t486 = t600 * t507;
t482 = qJD(1) * t486;
t470 = t502 * t482;
t485 = t600 * t505;
t481 = qJD(1) * t485;
t425 = t481 * t598 + t470;
t494 = pkin(2) * t502 + pkin(8);
t601 = pkin(9) + t494;
t539 = qJD(4) * t601;
t587 = t465 * t506;
t627 = pkin(4) * t467 - pkin(9) * t587 - t425 * t504 + t506 * t539 + t399;
t568 = t504 * t406 + t506 * t425;
t588 = t465 * t504;
t626 = -pkin(9) * t588 + t504 * t539 + t568;
t556 = qJD(1) * qJD(2);
t541 = t505 * t556;
t487 = t502 * t541;
t513 = qJD(2) * t489 - t487;
t511 = qJD(4) * t436 + t506 * t513;
t542 = t602 * qJD(5);
t559 = qJD(4) * t506;
t560 = qJD(4) * t504;
t553 = qJD(2) * t560 + t467 * t559 + t504 * t513;
t558 = qJD(5) * t503;
t342 = -t436 * t542 + t437 * t558 + t503 * t553 - t602 * t511;
t459 = qJD(4) - t465;
t455 = qJD(5) + t459;
t333 = t376 * t455 - t342;
t343 = t436 * t558 + t437 * t542 + t503 * t511 + t602 * t553;
t466 = t477 * qJD(2);
t456 = qJD(1) * t466;
t603 = t522 ^ 2;
t625 = t456 * MDP(24) + (t455 * t522 - t343) * MDP(23) + MDP(20) * t628 + (-t376 ^ 2 + t603) * MDP(21) + t333 * MDP(22);
t599 = qJD(2) * pkin(2);
t473 = t481 + t599;
t421 = t473 * t598 + t470;
t414 = -qJD(2) * pkin(3) - t421;
t374 = -t436 * pkin(4) + t414;
t340 = t376 * pkin(5) - qJ(6) * t522 + t374;
t624 = t340 * t376;
t623 = t374 * t376;
t581 = t503 * t504;
t520 = t602 * t506 - t581;
t402 = t520 * t465;
t609 = qJD(4) + qJD(5);
t611 = t602 * qJD(4) + t542;
t429 = -t506 * t611 + t581 * t609;
t565 = -t429 - t402;
t550 = t602 * t504;
t480 = t503 * t506 + t550;
t401 = t480 * t465;
t430 = t609 * t480;
t564 = t430 - t401;
t518 = -t502 * t505 + t537;
t469 = t518 * qJD(2);
t622 = t469 * t504 + t477 * t559;
t358 = pkin(5) * t522 + qJ(6) * t376;
t619 = -0.2e1 * t556;
t617 = MDP(5) * (t505 ^ 2 - t507 ^ 2);
t579 = t504 * t456;
t552 = -pkin(2) * t507 - pkin(1);
t420 = -pkin(3) * t518 - pkin(8) * t477 + t552;
t413 = t506 * t420;
t435 = t502 * t485 - t486 * t598;
t583 = t477 * t506;
t362 = -pkin(4) * t518 - pkin(9) * t583 - t435 * t504 + t413;
t428 = t506 * t435;
t566 = t504 * t420 + t428;
t584 = t477 * t504;
t369 = -pkin(9) * t584 + t566;
t616 = t503 * t362 + t602 * t369;
t474 = t601 * t504;
t475 = t601 * t506;
t521 = -t474 * t602 - t503 * t475;
t615 = -qJD(5) * t521 + t627 * t503 + t626 * t602;
t418 = -t503 * t474 + t475 * t602;
t614 = -qJD(5) * t418 + t626 * t503 - t627 * t602;
t538 = t598 * t482;
t424 = t481 * t502 - t538;
t533 = -t424 + (t560 - t588) * pkin(4);
t440 = t506 * t456;
t613 = t459 * t560 - t440;
t612 = -t437 * t560 + t511 * t506;
t547 = t477 * t560;
t585 = t469 * t506;
t610 = -t547 + t585;
t451 = t456 * pkin(5);
t530 = t552 * qJD(1);
t483 = qJD(3) + t530;
t395 = -pkin(3) * t465 - pkin(8) * t467 + t483;
t492 = pkin(2) * t541;
t394 = t456 * pkin(3) - pkin(8) * t513 + t492;
t387 = t506 * t394;
t422 = t502 * t473 - t538;
t415 = qJD(2) * pkin(8) + t422;
t529 = -t415 * t559 + t387;
t540 = qJD(2) * t600;
t463 = qJD(3) * t507 + t505 * t540;
t442 = t463 * qJD(1);
t464 = -qJD(3) * t505 + t507 * t540;
t443 = t464 * qJD(1);
t390 = t442 * t598 + t502 * t443;
t580 = t504 * t390;
t323 = t456 * pkin(4) - pkin(9) * t511 - t395 * t560 + t529 - t580;
t516 = t506 * t390 + t504 * t394 + t395 * t559 - t415 * t560;
t328 = -pkin(9) * t553 + t516;
t365 = t506 * t395 - t415 * t504;
t356 = -pkin(9) * t437 + t365;
t348 = pkin(4) * t459 + t356;
t366 = t504 * t395 + t506 * t415;
t357 = pkin(9) * t436 + t366;
t534 = -t602 * t323 + t503 * t328 + t348 * t558 + t357 * t542;
t314 = -t451 + t534;
t512 = t340 * t522 + t314;
t608 = -t374 * t522 - t534;
t607 = -t459 ^ 2 * t506 - t579;
t606 = -t342 * t520 - t522 * t564;
t555 = t505 * t599;
t407 = pkin(3) * t466 - pkin(8) * t469 + t555;
t400 = t506 * t407;
t405 = t463 * t598 + t502 * t464;
t335 = -pkin(9) * t585 + pkin(4) * t466 - t405 * t504 + t400 + (-t428 + (pkin(9) * t477 - t420) * t504) * qJD(4);
t515 = t506 * t405 + t504 * t407 + t420 * t559 - t435 * t560;
t339 = -pkin(9) * t622 + t515;
t605 = -qJD(5) * t616 + t335 * t602 - t503 * t339;
t551 = t602 * t357;
t325 = t503 * t348 + t551;
t597 = t325 * t455;
t595 = t376 * t467;
t594 = t522 * t467;
t593 = t521 * t456;
t592 = t418 * t456;
t591 = t436 * t467;
t590 = t437 * t459;
t589 = t437 * t467;
t426 = t520 * t456;
t427 = t480 * t456;
t582 = t503 * t357;
t508 = qJD(2) ^ 2;
t578 = t505 * t508;
t577 = t507 * t508;
t509 = qJD(1) ^ 2;
t576 = t507 * t509;
t575 = qJ(6) * t467 + t615;
t574 = -t467 * pkin(5) + t614;
t332 = t356 * t602 - t582;
t573 = -pkin(4) * t542 - qJD(6) + t332;
t572 = -pkin(5) * t564 + qJ(6) * t565 + qJD(6) * t480 - t533;
t569 = -t430 * t455 + t426;
t567 = -t429 * t455 + t427;
t324 = t348 * t602 - t582;
t557 = qJD(6) - t324;
t545 = t414 * t559;
t536 = pkin(1) * t619;
t389 = t502 * t442 - t598 * t443;
t404 = t463 * t502 - t598 * t464;
t434 = -t598 * t485 - t486 * t502;
t535 = t503 * t323 + t602 * t328 + t348 * t542 - t357 * t558;
t496 = -t598 * pkin(2) - pkin(3);
t331 = t503 * t356 + t551;
t532 = pkin(4) * t558 - t331;
t531 = -t480 * t343 - t376 * t565;
t403 = pkin(4) * t584 + t434;
t528 = t389 * t477 - t435 * t456;
t371 = pkin(4) * t622 + t404;
t527 = t459 * t588 - t613;
t444 = t455 * qJD(6);
t447 = t456 * qJ(6);
t313 = t447 + t444 + t535;
t525 = t362 * t602 - t503 * t369;
t517 = t324 * t455 - t535;
t484 = -t506 * pkin(4) + t496;
t514 = t503 * t335 + t602 * t339 + t362 * t542 - t369 * t558;
t363 = pkin(4) * t553 + t389;
t499 = -pkin(4) * t602 - pkin(5);
t495 = pkin(4) * t503 + qJ(6);
t423 = t456 * t518;
t411 = -pkin(5) * t520 - t480 * qJ(6) + t484;
t410 = t520 * t477;
t409 = t480 * t477;
t351 = pkin(5) * t409 - qJ(6) * t410 + t403;
t350 = t469 * t550 - t503 * t547 - t558 * t584 + (t469 * t503 + t477 * t611) * t506;
t349 = t430 * t477 - t469 * t520;
t345 = pkin(4) * t437 + t358;
t337 = pkin(5) * t518 - t525;
t336 = -qJ(6) * t518 + t616;
t320 = t455 * qJ(6) + t325;
t319 = -t455 * pkin(5) + t557;
t318 = pkin(5) * t350 + qJ(6) * t349 - qJD(6) * t410 + t371;
t317 = t343 * pkin(5) + t342 * qJ(6) - qJD(6) * t522 + t363;
t316 = -t466 * pkin(5) - t605;
t315 = qJ(6) * t466 - qJD(6) * t518 + t514;
t1 = [0.2e1 * t507 * MDP(4) * t541 + (-t342 * t410 - t349 * t522) * MDP(20) + (t342 * t409 - t343 * t410 + t349 * t376 - t350 * t522) * MDP(21) + (-t313 * t409 + t314 * t410 - t315 * t376 + t316 * t522 - t319 * t349 - t320 * t350 - t336 * t343 - t337 * t342) * MDP(28) + (t436 * t466 - t459 * t622 - t477 * t579 + t518 * t553) * MDP(16) + (t436 * t610 - t437 * t622 - t511 * t584 - t553 * t583) * MDP(14) + (t390 * t518 + t404 * t467 + t405 * t465 - t421 * t469 - t422 * t466 + t434 * t513 + t528) * MDP(11) + (t342 * t518 - t349 * t455 + t410 * t456 + t466 * t522) * MDP(22) + (-t313 * t518 + t315 * t455 - t317 * t410 - t318 * t522 + t320 * t466 + t336 * t456 + t340 * t349 + t342 * t351) * MDP(29) + ((-t435 * t559 + t400) * t459 + t413 * t456 - t529 * t518 + t365 * t466 - t404 * t436 + t434 * t553 + t477 * t545 + ((-qJD(4) * t420 - t405) * t459 - (-qJD(4) * t395 - t390) * t518 + t414 * t469 + t528) * t504) * MDP(18) + (t343 * t518 - t350 * t455 - t376 * t466 - t409 * t456) * MDP(23) + (t314 * t518 - t316 * t455 + t317 * t409 + t318 * t376 - t319 * t466 - t337 * t456 + t340 * t350 + t343 * t351) * MDP(27) + (t324 * t466 + t403 * t343 + t374 * t350 + t363 * t409 + t371 * t376 + t455 * t605 + t525 * t456 + t518 * t534) * MDP(25) + (-t366 * t466 + t389 * t583 + t404 * t437 + t414 * t610 + t434 * t511 - t456 * t566 - t459 * t515 + t516 * t518) * MDP(19) + (t437 * t466 + t440 * t477 + t459 * t610 - t511 * t518) * MDP(15) + (-t325 * t466 - t403 * t342 - t374 * t349 + t363 * t410 + t371 * t522 - t455 * t514 - t456 * t616 + t518 * t535) * MDP(26) + (t455 * t466 - t423) * MDP(24) + (t459 * t466 - t423) * MDP(17) + MDP(6) * t577 + t617 * t619 + (t313 * t336 + t314 * t337 + t315 * t320 + t316 * t319 + t317 * t351 + t318 * t340) * MDP(30) + (t437 * t585 + t477 * t612) * MDP(13) + (t389 * t434 + t390 * t435 - t421 * t404 + t422 * t405 + (t483 + t530) * t555) * MDP(12) + (-pkin(7) * t577 + t505 * t536) * MDP(9) - MDP(7) * t578 + (pkin(7) * t578 + t507 * t536) * MDP(10); t509 * t617 + ((-t425 + t421) * t465 + (-t502 * t456 - t513 * t598) * pkin(2)) * MDP(11) + (t421 * t424 - t422 * t425 + (-t389 * t598 + t390 * t502 - t483 * t561) * pkin(2)) * MDP(12) + ((-t487 + (t489 + qJD(4)) * qJD(2)) * t504 + t590) * t506 * MDP(13) + (t437 * t588 - t504 * t553 + (t559 - t587) * t436 + t612) * MDP(14) + (-t589 - t607) * MDP(15) + (t527 - t591) * MDP(16) + (-t494 * t579 - t389 * t506 + t424 * t436 + t496 * t553 + (-t494 * t559 - t399 + (t414 + t425) * t504) * t459) * MDP(18) + (t389 * t504 - t414 * t587 - t424 * t437 + t459 * t568 + t494 * t613 + t496 * t511 + t545) * MDP(19) + (-t342 * t480 + t522 * t565) * MDP(20) + (t531 + t606) * MDP(21) + (-t402 * t455 + t567 - t594) * MDP(22) + (t401 * t455 + t569 + t595) * MDP(23) + (t484 * t343 - t363 * t520 + t564 * t374 + t533 * t376 + t455 * t614 + t593) * MDP(25) + (-t484 * t342 + t363 * t480 + t565 * t374 + t455 * t615 + t533 * t522 - t592) * MDP(26) + (-t317 * t520 + t340 * t564 + t343 * t411 - t376 * t572 + t455 * t574 + t593) * MDP(27) + (t313 * t520 + t314 * t480 + t319 * t565 - t320 * t564 + t342 * t521 - t343 * t418 + t376 * t575 - t522 * t574) * MDP(28) + (-t317 * t480 - t340 * t565 + t342 * t411 - t455 * t575 + t522 * t572 + t592) * MDP(29) + (t313 * t418 - t314 * t521 + t317 * t411 - t319 * t574 - t320 * t575 - t340 * t572) * MDP(30) - t505 * MDP(4) * t576 + ((t422 - t424) * MDP(11) - qJD(4) * t504 ^ 2 * MDP(13) - t365 * MDP(18) + t366 * MDP(19) - t324 * MDP(25) + t325 * MDP(26) + t319 * MDP(27) - t320 * MDP(29) - t459 * MDP(17) - t455 * MDP(24)) * t467 + (MDP(9) * t505 * t509 + MDP(10) * t576) * pkin(1); (-t465 ^ 2 - t467 ^ 2) * MDP(11) + (t421 * t467 - t422 * t465 + t492) * MDP(12) + (t527 + t591) * MDP(18) + (-t589 + t607) * MDP(19) + (t569 - t595) * MDP(25) + (-t427 - t594) * MDP(26) + (t426 - t595) * MDP(27) + (t531 - t606) * MDP(28) + (t567 + t594) * MDP(29) + (t313 * t480 - t314 * t520 + t319 * t564 + t320 * t565 - t340 * t467) * MDP(30) + (t401 * MDP(25) - MDP(26) * t565 - MDP(27) * t564 - t402 * MDP(29)) * t455; -t437 * t436 * MDP(13) + (-t436 ^ 2 + t437 ^ 2) * MDP(14) + (-t436 * t459 + t511) * MDP(15) + (-t553 + t590) * MDP(16) + t456 * MDP(17) + (-t414 * t437 + t387 - t580 + (-qJD(4) + t459) * t366) * MDP(18) + (t365 * t459 - t414 * t436 - t516) * MDP(19) + (t331 * t455 + (-t376 * t437 - t455 * t558 + t456 * t602) * pkin(4) + t608) * MDP(25) + (t332 * t455 + t623 + (-t437 * t522 - t455 * t542 - t456 * t503) * pkin(4) - t535) * MDP(26) + (-t345 * t376 - t455 * t532 - t456 * t499 - t512) * MDP(27) + (-t342 * t499 - t343 * t495 + (t320 + t532) * t522 + (t319 + t573) * t376) * MDP(28) + (t345 * t522 - t455 * t573 + t456 * t495 + t313 - t624) * MDP(29) + (t313 * t495 + t314 * t499 + t319 * t532 - t320 * t573 - t340 * t345) * MDP(30) + t625; (t597 + t608) * MDP(25) + (t517 + t623) * MDP(26) + (-t358 * t376 + t451 - t512 + t597) * MDP(27) + (pkin(5) * t342 - qJ(6) * t343 + (t320 - t325) * t522 + (t319 - t557) * t376) * MDP(28) + (t358 * t522 + 0.2e1 * t444 + 0.2e1 * t447 - t517 - t624) * MDP(29) + (-pkin(5) * t314 + qJ(6) * t313 - t319 * t325 + t320 * t557 - t340 * t358) * MDP(30) + t625; (-qJD(2) * t467 + t628) * MDP(27) + t333 * MDP(28) + (-t455 ^ 2 - t603) * MDP(29) + (-t320 * t455 + t512) * MDP(30);];
tauc  = t1;
