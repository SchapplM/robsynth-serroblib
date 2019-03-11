% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:13:34
% EndTime: 2019-03-09 18:13:44
% DurationCPUTime: 6.33s
% Computational Cost: add. (6008->428), mult. (14459->550), div. (0->0), fcn. (10606->8), ass. (0->193)
t501 = sin(qJ(6));
t563 = qJD(6) * t501;
t497 = qJD(2) + qJD(3);
t503 = sin(qJ(3));
t504 = sin(qJ(2));
t583 = t503 * t504;
t530 = t497 * t583;
t507 = cos(qJ(2));
t604 = cos(qJ(3));
t552 = t604 * t507;
t532 = qJD(1) * t552;
t571 = t497 * t532;
t407 = qJD(1) * t530 - t571;
t466 = t503 * t507 + t504 * t604;
t421 = t497 * t466;
t408 = t421 * qJD(1);
t567 = qJD(1) * t504;
t551 = t503 * t567;
t448 = -t532 + t551;
t502 = sin(qJ(5));
t506 = cos(qJ(5));
t564 = qJD(5) * t506;
t565 = qJD(5) * t502;
t568 = qJD(1) * t466;
t344 = -t506 * t407 + t502 * t408 + t448 * t564 - t565 * t568;
t505 = cos(qJ(6));
t559 = qJD(5) - t497;
t562 = qJD(6) * t505;
t577 = t505 * t344 + t559 * t562;
t609 = t448 * t502 + t506 * t568;
t328 = -t563 * t609 + t577;
t383 = t501 * t559 + t505 * t609;
t597 = t344 * t501;
t329 = qJD(6) * t383 + t597;
t345 = qJD(5) * t609 - t407 * t502 - t506 * t408;
t591 = t609 * t501;
t381 = -t505 * t559 + t591;
t581 = t505 * t345;
t585 = t501 * t345;
t610 = -t506 * t448 + t502 * t568;
t619 = qJD(6) + t610;
t592 = t383 * t619;
t593 = t381 * t619;
t599 = t328 * t501;
t626 = t619 * t505;
t634 = ((t329 + t592) * t501 - (t328 - t593) * t505) * MDP(30) - (t383 * t626 + t599) * MDP(29) - (-t383 * t609 + t619 * t626 + t585) * MDP(31) + (t501 * t619 ^ 2 - t381 * t609 - t581) * MDP(32) - (t609 ^ 2 - t610 ^ 2) * MDP(23) - (t559 * t610 + t344) * MDP(24) - (t559 * t609 - t345) * MDP(25) + (-MDP(22) * t610 + t619 * MDP(33)) * t609;
t443 = t568 * pkin(9);
t605 = -pkin(8) - pkin(7);
t478 = t605 * t504;
t472 = qJD(1) * t478;
t601 = qJD(2) * pkin(2);
t458 = t472 + t601;
t479 = t605 * t507;
t474 = qJD(1) * t479;
t584 = t503 * t474;
t411 = t604 * t458 + t584;
t608 = qJD(4) - t411;
t620 = -t443 + t608;
t457 = t604 * t474;
t418 = t503 * t472 - t457;
t566 = qJD(3) * t503;
t531 = pkin(2) * t566 - t418;
t419 = t604 * t472 + t584;
t549 = qJD(3) * t604;
t573 = pkin(2) * t549 + qJD(4) - t419;
t623 = pkin(5) * t609;
t491 = -t507 * pkin(2) - pkin(1);
t477 = t491 * qJD(1);
t388 = t448 * pkin(3) - qJ(4) * t568 + t477;
t370 = -pkin(4) * t448 - t388;
t333 = pkin(5) * t610 - pkin(10) * t609 + t370;
t508 = -pkin(3) - pkin(4);
t368 = t497 * t508 + t620;
t412 = t503 * t458 - t457;
t602 = t448 * pkin(9);
t380 = t412 + t602;
t494 = t497 * qJ(4);
t376 = t380 + t494;
t340 = t368 * t502 + t376 * t506;
t337 = pkin(10) * t559 + t340;
t325 = t333 * t501 + t337 * t505;
t600 = t325 * t609;
t622 = t602 - t531;
t621 = t443 - t573;
t493 = t497 * qJD(4);
t553 = qJD(2) * t605;
t533 = qJD(1) * t553;
t459 = t504 * t533;
t460 = t507 * t533;
t534 = t458 * t549 + t604 * t459 + t503 * t460 + t474 * t566;
t364 = t493 + t534;
t348 = pkin(9) * t408 + t364;
t365 = t458 * t566 + t503 * t459 - t604 * t460 - t474 * t549;
t350 = pkin(9) * t407 + t365;
t319 = t502 * t348 - t350 * t506 + t368 * t565 + t376 * t564;
t324 = t333 * t505 - t337 * t501;
t545 = t319 * t505 + t324 * t609;
t318 = t506 * t348 + t502 * t350 + t368 * t564 - t376 * t565;
t514 = -t370 * t610 + t318;
t516 = t370 * t609 + t319;
t557 = qJD(1) * qJD(2);
t612 = -0.2e1 * t557;
t611 = MDP(5) * (t504 ^ 2 - t507 ^ 2);
t423 = t503 * t478 - t604 * t479;
t473 = t504 * t553;
t475 = t507 * t553;
t371 = t604 * t473 + t503 * t475 + t478 * t549 + t479 * t566;
t358 = pkin(9) * t421 + t371;
t372 = qJD(3) * t423 + t503 * t473 - t604 * t475;
t420 = -qJD(2) * t552 - t507 * t549 + t530;
t359 = t420 * pkin(9) + t372;
t422 = -t478 * t604 - t503 * t479;
t393 = -t466 * pkin(9) + t422;
t465 = -t552 + t583;
t394 = pkin(9) * t465 + t423;
t527 = t393 * t506 - t394 * t502;
t326 = qJD(5) * t527 + t358 * t506 + t359 * t502;
t339 = t368 * t506 - t376 * t502;
t336 = -pkin(5) * t559 - t339;
t410 = t465 * pkin(3) - t466 * qJ(4) + t491;
t385 = -pkin(4) * t465 - t410;
t415 = t465 * t502 + t466 * t506;
t525 = t506 * t465 - t466 * t502;
t346 = -pkin(5) * t525 - pkin(10) * t415 + t385;
t355 = qJD(5) * t525 - t420 * t506 + t421 * t502;
t361 = t393 * t502 + t394 * t506;
t607 = t319 * t415 + t336 * t355 - t361 * t345 - (qJD(6) * t346 + t326) * t619 + (qJD(6) * t333 + t318) * t525;
t606 = t568 ^ 2;
t603 = pkin(4) * t568;
t598 = t336 * t415;
t596 = t346 * t345;
t595 = t371 * t497;
t594 = t372 * t497;
t590 = t411 * t497;
t589 = t412 * t497;
t588 = t448 * t568;
t509 = qJD(2) ^ 2;
t582 = t504 * t509;
t579 = t507 * t509;
t510 = qJD(1) ^ 2;
t578 = t507 * t510;
t528 = -qJ(4) * t502 + t506 * t508;
t576 = -qJD(5) * t528 + t380 * t502 - t620 * t506;
t490 = -pkin(2) * t604 - pkin(3);
t487 = -pkin(4) + t490;
t488 = pkin(2) * t503 + qJ(4);
t524 = t487 * t506 - t488 * t502;
t575 = -qJD(5) * t524 + t622 * t502 + t621 * t506;
t523 = t487 * t502 + t488 * t506;
t574 = qJD(5) * t523 - t621 * t502 + t622 * t506;
t529 = qJ(4) * t506 + t502 * t508;
t572 = qJD(5) * t529 + t380 * t506 + t620 * t502;
t409 = pkin(3) * t568 + t448 * qJ(4);
t556 = pkin(2) * t567;
t555 = t504 * t601;
t548 = t504 * t557;
t547 = t619 * pkin(10) + t623;
t546 = pkin(1) * t612;
t543 = t559 ^ 2;
t541 = t506 * t559;
t540 = t619 * t336;
t378 = -t409 - t603;
t338 = -pkin(10) * t610 + t378 - t623;
t425 = -pkin(10) + t523;
t536 = qJD(6) * t425 + t338 - t556;
t471 = -pkin(10) + t529;
t535 = qJD(6) * t471 + t338;
t392 = t556 + t409;
t522 = t355 * t505 - t415 * t563;
t362 = t421 * pkin(3) + t420 * qJ(4) - t466 * qJD(4) + t555;
t521 = -t388 * t568 - t365;
t520 = -t388 * t448 + t534;
t519 = -t477 * t568 - t365;
t518 = t477 * t448 - t534;
t357 = pkin(2) * t548 + t408 * pkin(3) + t407 * qJ(4) - qJD(4) * t568;
t347 = -pkin(4) * t421 - t362;
t515 = -pkin(10) * t345 + (t336 + t339) * t619;
t334 = -pkin(4) * t408 - t357;
t513 = -t425 * t345 + t575 * t619 - t540;
t512 = -t471 * t345 + t576 * t619 - t540;
t377 = t571 + (t448 - t551) * t497;
t511 = MDP(11) * t588 + t377 * MDP(13) + (-t448 ^ 2 + t606) * MDP(12) + t634;
t470 = pkin(5) - t528;
t424 = pkin(5) - t524;
t396 = t494 + t412;
t391 = -pkin(3) * t497 + t608;
t373 = -t392 - t603;
t356 = qJD(5) * t415 - t420 * t502 - t506 * t421;
t327 = qJD(5) * t361 + t358 * t502 - t359 * t506;
t323 = pkin(5) * t356 - pkin(10) * t355 + t347;
t316 = pkin(5) * t345 - pkin(10) * t344 + t334;
t315 = t505 * t316;
t1 = [(-MDP(13) * t420 - MDP(14) * t421) * t497 + (-t345 * t525 + t356 * t619) * MDP(33) + (-t328 * t525 + t356 * t383 + t415 * t581 + t522 * t619) * MDP(31) + (-t415 * t585 + t329 * t525 - t356 * t381 + (-t355 * t501 - t415 * t562) * t619) * MDP(32) + (-t325 * t356 + t327 * t383 - t527 * t328 + (-(-qJD(6) * t361 + t323) * t619 - t596 + (-qJD(6) * t337 + t316) * t525 - qJD(6) * t598) * t501 + t607 * t505) * MDP(35) + (-t315 * t525 + t324 * t356 + t327 * t381 - t527 * t329 + (t323 * t619 + t596 + (t337 * t525 - t361 * t619 + t598) * qJD(6)) * t505 + t607 * t501) * MDP(34) + (t328 * t415 * t505 + t383 * t522) * MDP(29) + 0.2e1 * t507 * MDP(4) * t548 + (t334 * t415 + t344 * t385 + t347 * t609 + t355 * t370) * MDP(28) + (t344 * t415 + t355 * t609) * MDP(22) + t611 * t612 + (-t334 * t525 + t345 * t385 + t347 * t610 + t356 * t370) * MDP(27) + (t344 * t525 - t345 * t415 - t355 * t610 - t356 * t609) * MDP(23) + (-t407 * t491 - t420 * t477 + 0.2e1 * t555 * t568 - t595) * MDP(17) + (-t357 * t466 - t362 * t568 + t388 * t420 + t407 * t410 + t595) * MDP(20) + (-t407 * t466 - t420 * t568) * MDP(11) + (t407 * t465 - t408 * t466 + t420 * t448 - t421 * t568) * MDP(12) + (-t364 * t465 + t365 * t466 - t371 * t448 + t372 * t568 - t391 * t420 - t396 * t421 - t407 * t422 - t408 * t423) * MDP(19) - (-MDP(24) * t355 + MDP(25) * t356 + MDP(27) * t327 + MDP(28) * t326) * t559 + (t357 * t410 + t362 * t388 + t364 * t423 + t365 * t422 + t371 * t396 + t372 * t391) * MDP(21) + (-pkin(7) * t579 + t504 * t546) * MDP(9) - MDP(7) * t582 + (pkin(7) * t582 + t507 * t546) * MDP(10) + MDP(6) * t579 + (-t594 + t408 * t491 + t421 * t477 + (qJD(1) * t465 + t448) * t555) * MDP(16) + (t357 * t465 + t362 * t448 + t388 * t421 + t408 * t410 - t594) * MDP(18) + ((-t381 * t505 - t383 * t501) * t355 + (-t599 - t329 * t505 + (t381 * t501 - t383 * t505) * qJD(6)) * t415) * MDP(30); (t392 * t568 + t497 * t573 + t493 + t520) * MDP(20) + (-t407 * t490 - t408 * t488 + (t396 + t531) * t568 + (t391 - t573) * t448) * MDP(19) + (t418 * t497 + (-t448 * t567 - t497 * t566) * pkin(2) + t519) * MDP(16) + (-t392 * t448 - t497 * t531 + t521) * MDP(18) + (t419 * t497 + (-t497 * t549 - t567 * t568) * pkin(2) + t518) * MDP(17) + (-t373 * t609 + t559 * t575 + t514) * MDP(28) + (-t373 * t610 - t559 * t574 + t516) * MDP(27) + (t364 * t488 + t365 * t490 - t388 * t392 + t391 * t531 + t396 * t573) * MDP(21) + (t424 * t329 + t381 * t574 + t501 * t513 - t536 * t626 + t545) * MDP(34) - t504 * MDP(4) * t578 + t510 * t611 + (-t600 + t424 * t328 + t574 * t383 + (t536 * t619 - t319) * t501 + t513 * t505) * MDP(35) + t511 + (MDP(9) * t504 * t510 + MDP(10) * t578) * pkin(1); (t519 + t589) * MDP(16) + (-t409 * t448 + t521 + t589) * MDP(18) + (pkin(3) * t407 - qJ(4) * t408 + (t396 - t412) * t568 + (t391 - t608) * t448) * MDP(19) + (t409 * t568 + 0.2e1 * t493 + t520 - t590) * MDP(20) + (-t378 * t609 + t559 * t576 + t514) * MDP(28) + (-t378 * t610 - t559 * t572 + t516) * MDP(27) + (-pkin(3) * t365 + qJ(4) * t364 - t388 * t409 - t391 * t412 + t396 * t608) * MDP(21) + (t470 * t329 + t381 * t572 + t501 * t512 - t535 * t626 + t545) * MDP(34) + (-t600 + t470 * t328 + t572 * t383 + (t535 * t619 - t319) * t501 + t512 * t505) * MDP(35) + t511 + (t518 + t590) * MDP(17); MDP(18) * t588 + t377 * MDP(19) + (-t497 ^ 2 - t606) * MDP(20) + (-t396 * t497 - t521) * MDP(21) + (-t502 * t543 - t568 * t610) * MDP(27) + (-t506 * t543 - t568 * t609) * MDP(28) + (-t506 * t329 + (-t501 * t541 - t505 * t568) * t619 + (t381 * t559 - t562 * t619 - t585) * t502) * MDP(34) + (-t506 * t328 + (t501 * t568 - t505 * t541) * t619 + (t383 * t559 + t563 * t619 - t581) * t502) * MDP(35); (t340 * t559 - t516) * MDP(27) + (t339 * t559 - t514) * MDP(28) + (-pkin(5) * t329 - t340 * t381 + t501 * t515 - t547 * t626 - t545) * MDP(34) + (-pkin(5) * t328 + t600 - t340 * t383 + (t547 * t619 + t319) * t501 + t515 * t505) * MDP(35) - t634; t383 * t381 * MDP(29) + (-t381 ^ 2 + t383 ^ 2) * MDP(30) + (t577 + t593) * MDP(31) + (t592 - t597) * MDP(32) + t345 * MDP(33) + (-t318 * t501 + t325 * t619 - t336 * t383 + t315) * MDP(34) + (-t316 * t501 - t318 * t505 + t324 * t619 + t336 * t381) * MDP(35) + (-MDP(31) * t591 - MDP(32) * t383 - MDP(34) * t325 - MDP(35) * t324) * qJD(6);];
tauc  = t1;
