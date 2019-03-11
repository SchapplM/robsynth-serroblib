% Calculate vector of inverse dynamics joint torques for
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:12:21
% EndTime: 2019-03-09 08:12:32
% DurationCPUTime: 9.88s
% Computational Cost: add. (4802->533), mult. (11008->657), div. (0->0), fcn. (8073->14), ass. (0->231)
t502 = sin(pkin(9));
t504 = cos(pkin(9));
t507 = sin(qJ(2));
t510 = cos(qJ(2));
t458 = t502 * t510 + t504 * t507;
t447 = t458 * qJD(1);
t437 = qJD(6) + t447;
t567 = qJD(1) * qJD(2);
t557 = t510 * t567;
t558 = t507 * t567;
t533 = -t502 * t558 + t504 * t557;
t412 = qJDD(1) * t458 + t533;
t407 = qJDD(6) + t412;
t501 = sin(pkin(10));
t503 = cos(pkin(10));
t506 = sin(qJ(6));
t509 = cos(qJ(6));
t457 = t501 * t509 + t503 * t506;
t578 = t437 * t457;
t620 = -t501 * t506 + t503 * t509;
t626 = t620 * t407 - t437 * t578;
t581 = t504 * t510;
t584 = t502 * t507;
t456 = -t581 + t584;
t399 = t620 * t456;
t572 = qJD(6) * t509;
t573 = qJD(6) * t506;
t577 = t620 * t447 - t501 * t573 + t503 * t572;
t546 = -t457 * t407 - t577 * t437;
t575 = qJD(1) * t507;
t444 = -qJD(1) * t581 + t502 * t575;
t422 = qJD(2) * t501 - t503 * t444;
t424 = qJD(2) * t503 + t444 * t501;
t536 = t422 * t506 - t424 * t509;
t625 = t437 * t536;
t492 = t510 * pkin(2);
t485 = t492 + pkin(1);
t495 = qJ(2) + pkin(9);
t491 = cos(t495);
t508 = sin(qJ(1));
t511 = cos(qJ(1));
t621 = g(1) * t508 - g(2) * t511;
t624 = t491 * t621;
t489 = sin(t495);
t548 = g(1) * t511 + g(2) * t508;
t623 = t548 * t489;
t622 = -t491 * pkin(3) - t489 * qJ(4);
t600 = qJ(3) + pkin(7);
t466 = t600 * t510;
t462 = qJD(1) * t466;
t452 = t502 * t462;
t465 = t600 * t507;
t461 = qJD(1) * t465;
t417 = -t461 * t504 - t452;
t618 = qJD(4) - t417;
t446 = t458 * qJD(2);
t565 = qJDD(1) * t510;
t473 = t504 * t565;
t566 = qJDD(1) * t507;
t411 = qJD(1) * t446 + t502 * t566 - t473;
t617 = -pkin(4) * t411 + qJDD(5);
t482 = g(3) * t491;
t525 = t482 - t623;
t593 = t412 * t501;
t614 = t447 ^ 2;
t616 = -t503 * t614 - t593;
t613 = pkin(2) * t507;
t612 = pkin(3) * t411;
t610 = pkin(4) * t444;
t609 = pkin(4) * t447;
t608 = pkin(8) * t503;
t604 = g(3) * t510;
t603 = pkin(3) + qJ(5);
t602 = pkin(4) + t600;
t484 = -pkin(2) * t504 - pkin(3);
t477 = -qJ(5) + t484;
t601 = -pkin(8) + t477;
t599 = qJ(4) * t412;
t598 = qJ(5) * t491;
t597 = qJDD(2) * pkin(3);
t591 = t424 * t506;
t370 = t509 * t422 + t591;
t596 = t370 * t437;
t595 = t370 * t444;
t594 = t536 * t444;
t592 = t412 * t503;
t494 = pkin(10) + qJ(6);
t488 = sin(t494);
t590 = t488 * t508;
t589 = t489 * t511;
t490 = cos(t494);
t588 = t490 * t508;
t587 = t490 * t511;
t586 = t491 * t511;
t582 = t504 * t462;
t580 = t600 * t511;
t564 = pkin(2) * t558 + qJDD(3);
t524 = -qJDD(1) * t485 + t564;
t516 = -qJD(4) * t447 + t524 - t599;
t326 = qJD(5) * t444 + t411 * t603 + t516;
t555 = qJD(2) * t600;
t439 = -qJD(3) * t507 - t510 * t555;
t406 = qJDD(2) * pkin(2) + qJD(1) * t439 - qJDD(1) * t465;
t438 = qJD(3) * t510 - t507 * t555;
t415 = qJD(1) * t438 + qJDD(1) * t466;
t355 = t406 * t504 - t502 * t415;
t544 = qJDD(4) - t355;
t339 = pkin(4) * t412 - qJD(2) * qJD(5) - qJDD(2) * t603 + t544;
t320 = t503 * t326 + t501 * t339;
t574 = qJD(2) * t507;
t449 = qJD(2) * t581 - t502 * t574;
t487 = pkin(2) * t574;
t529 = -qJ(4) * t449 - qJD(4) * t458 + t487;
t348 = qJD(5) * t456 + t446 * t603 + t529;
t394 = t438 * t502 - t504 * t439;
t368 = pkin(4) * t449 + t394;
t329 = t503 * t348 + t501 * t368;
t464 = -qJD(1) * t485 + qJD(3);
t522 = -qJ(4) * t447 + t464;
t361 = t444 * t603 + t522;
t455 = qJD(2) * pkin(2) - t461;
t409 = t455 * t504 - t452;
t545 = qJD(4) - t409;
t365 = -qJD(2) * t603 + t545 + t609;
t336 = t503 * t361 + t501 * t365;
t554 = pkin(2) * t575 + qJ(4) * t444;
t366 = t447 * t603 + t554;
t416 = -t461 * t502 + t582;
t379 = t416 - t610;
t343 = t503 * t366 + t501 * t379;
t550 = -qJ(4) * t458 - t485;
t378 = t456 * t603 + t550;
t420 = t504 * t465 + t466 * t502;
t396 = pkin(4) * t458 + t420;
t346 = t503 * t378 + t501 * t396;
t356 = t502 * t406 + t504 * t415;
t410 = t502 * t455 + t582;
t499 = t507 ^ 2;
t576 = -t510 ^ 2 + t499;
t559 = -pkin(5) * t503 - pkin(4);
t570 = -t447 * t559 + t618;
t569 = t609 + t618;
t402 = -qJD(2) * qJ(4) - t410;
t373 = qJD(5) - t402 - t610;
t568 = qJD(5) - t373;
t391 = qJDD(2) * t501 - t503 * t411;
t392 = qJDD(2) * t503 + t411 * t501;
t563 = -t506 * t391 + t509 * t392 - t422 * t572;
t562 = qJDD(2) * qJ(4) + t356;
t561 = t492 - t622;
t480 = pkin(2) * t502 + qJ(4);
t319 = -t326 * t501 + t503 * t339;
t315 = pkin(5) * t412 - pkin(8) * t392 + t319;
t318 = -pkin(8) * t391 + t320;
t553 = t509 * t315 - t318 * t506;
t335 = -t361 * t501 + t503 * t365;
t552 = t509 * t391 + t506 * t392;
t497 = qJD(2) * qJD(4);
t352 = -t497 - t562;
t472 = t511 * t485;
t551 = g(2) * (pkin(3) * t586 + qJ(4) * t589 + t472);
t549 = -pkin(3) * t489 - t613;
t543 = t315 * t506 + t318 * t509;
t542 = t319 * t503 + t320 * t501;
t541 = -t319 * t501 + t320 * t503;
t323 = pkin(5) * t447 - pkin(8) * t424 + t335;
t327 = -pkin(8) * t422 + t336;
t316 = t323 * t509 - t327 * t506;
t317 = t323 * t506 + t327 * t509;
t390 = t503 * t396;
t334 = pkin(5) * t458 + t390 + (-pkin(8) * t456 - t378) * t501;
t341 = t456 * t608 + t346;
t540 = t334 * t509 - t341 * t506;
t539 = t334 * t506 + t341 * t509;
t538 = -t335 * t503 - t336 * t501;
t537 = -t335 * t501 + t336 * t503;
t395 = t438 * t504 + t439 * t502;
t421 = -t465 * t502 + t466 * t504;
t534 = -t621 + t564;
t532 = -t485 + t622;
t531 = t621 * t489;
t530 = -0.2e1 * pkin(1) * t567 - pkin(7) * qJDD(2);
t377 = t503 * t379;
t442 = t601 * t501;
t528 = qJD(5) * t503 + qJD(6) * t442 - pkin(5) * t444 + t377 + (-pkin(8) * t447 - t366) * t501;
t443 = t601 * t503;
t527 = qJD(5) * t501 - qJD(6) * t443 + t447 * t608 + t343;
t400 = t457 * t456;
t331 = -t424 * t573 + t563;
t340 = -t352 + t617;
t523 = -g(3) * t489 - t491 * t548;
t512 = qJD(2) ^ 2;
t521 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t512 + t621;
t513 = qJD(1) ^ 2;
t520 = pkin(1) * t513 - pkin(7) * qJDD(1) + t548;
t519 = t340 * t456 + t373 * t446 + t548;
t518 = t340 + t523;
t332 = -qJD(6) * t536 + t552;
t517 = t523 + t562;
t515 = t394 * t447 - t395 * t444 - t411 * t421 + t412 * t420 - t548;
t386 = pkin(3) * t444 + t522;
t514 = t386 * t447 + t525 + t544;
t469 = qJ(4) * t586;
t467 = t508 * t491 * qJ(4);
t463 = pkin(5) * t501 + t480;
t436 = qJD(2) * t444;
t431 = -t489 * t590 + t587;
t430 = t488 * t511 + t489 * t588;
t429 = t488 * t589 + t588;
t428 = t489 * t587 - t590;
t408 = pkin(3) * t456 + t550;
t398 = -qJD(2) * pkin(3) + t545;
t397 = -pkin(4) * t456 + t421;
t393 = pkin(3) * t447 + t554;
t374 = pkin(3) * t446 + t529;
t369 = -pkin(4) * t446 + t395;
t367 = t456 * t559 + t421;
t364 = t503 * t368;
t354 = t446 * t559 + t395;
t353 = t544 - t597;
t351 = pkin(5) * t422 + t373;
t350 = qJD(6) * t400 - t446 * t620;
t349 = qJD(6) * t399 + t446 * t457;
t345 = -t378 * t501 + t390;
t344 = t516 + t612;
t342 = -t366 * t501 + t377;
t328 = -t348 * t501 + t364;
t324 = pkin(5) * t391 + t340;
t322 = t446 * t608 + t329;
t321 = pkin(5) * t449 + t364 + (-pkin(8) * t446 - t348) * t501;
t1 = [(qJDD(2) * t507 + t510 * t512) * MDP(6) + (qJDD(2) * t510 - t507 * t512) * MDP(7) + (t407 * t458 + t437 * t449) * MDP(25) + (qJD(2) * t394 + qJDD(2) * t420 - t344 * t456 - t374 * t444 - t386 * t446 - t408 * t411 - t624) * MDP(14) + (-t328 * t424 - t329 * t422 - t345 * t392 - t346 * t391 + t446 * t537 + t456 * t541 + t624) * MDP(19) + (qJDD(1) * t499 + 0.2e1 * t507 * t557) * MDP(4) + t548 * MDP(3) + (-t507 * t521 + t510 * t530) * MDP(10) + (-t320 * t458 - t329 * t447 - t336 * t449 - t346 * t412 + t369 * t424 + t397 * t392 + t501 * t519 + t503 * t531) * MDP(18) + (t319 * t458 + t328 * t447 + t335 * t449 + t345 * t412 + t369 * t422 + t397 * t391 + t501 * t531 - t503 * t519) * MDP(17) + (qJD(2) * t395 + qJDD(2) * t421 - t344 * t458 - t374 * t447 - t386 * t449 - t408 * t412 + t531) * MDP(15) + (-t355 * t458 - t356 * t456 - t409 * t449 - t410 * t446 + t515) * MDP(11) + (t352 * t456 + t353 * t458 + t398 * t449 + t402 * t446 + t515) * MDP(13) + (t344 * t408 + t386 * t374 - t352 * t421 - t402 * t395 + t353 * t420 + t398 * t394 - g(1) * t580 - t551 + (-g(1) * t532 - g(2) * t600) * t508) * MDP(16) + (t356 * t421 + t410 * t395 - t355 * t420 - t409 * t394 - t524 * t485 + t464 * t487 - g(1) * (-t485 * t508 + t580) - g(2) * (t508 * t600 + t472)) * MDP(12) + t621 * MDP(2) + (-t332 * t458 - t350 * t437 - t370 * t449 + t399 * t407) * MDP(24) + ((t321 * t509 - t322 * t506) * t437 + t540 * t407 + t553 * t458 + t316 * t449 + t354 * t370 + t367 * t332 - t324 * t399 + t351 * t350 - g(1) * t431 - g(2) * t429 + (-t317 * t458 - t437 * t539) * qJD(6)) * MDP(26) + (t331 * t399 - t332 * t400 - t349 * t370 + t350 * t536) * MDP(22) + (t331 * t458 + t349 * t437 + t400 * t407 - t449 * t536) * MDP(23) + (t331 * t400 - t349 * t536) * MDP(21) + (-(t321 * t506 + t322 * t509) * t437 - t539 * t407 - t543 * t458 - t317 * t449 - t354 * t536 + t367 * t331 + t324 * t400 + t351 * t349 + g(1) * t430 - g(2) * t428 + (-t316 * t458 - t437 * t540) * qJD(6)) * MDP(27) + (t507 * t530 + t510 * t521) * MDP(9) + 0.2e1 * (t507 * t565 - t567 * t576) * MDP(5) + (t320 * t346 + t336 * t329 + t319 * t345 + t335 * t328 + t340 * t397 + t373 * t369 - t551 + (-g(1) * t602 - g(2) * t598) * t511 + (-g(1) * (t532 - t598) - g(2) * t602) * t508) * MDP(20) + qJDD(1) * MDP(1); t437 * t444 * MDP(25) + MDP(6) * t566 + MDP(7) * t565 + qJDD(2) * MDP(8) + (t507 * t520 - t604) * MDP(9) + (g(3) * t507 + t510 * t520) * MDP(10) + ((t410 - t416) * t447 + (-t409 + t417) * t444 + (-t411 * t502 - t412 * t504) * pkin(2)) * MDP(11) + (t409 * t416 - t410 * t417 + (-t604 + t355 * t504 + t356 * t502 + (-qJD(1) * t464 + t548) * t507) * pkin(2)) * MDP(12) + (-t411 * t480 + t412 * t484 + (-t402 - t416) * t447 + (t398 - t618) * t444) * MDP(13) + (-qJD(2) * t416 + t393 * t444 + (-pkin(3) + t484) * qJDD(2) + t514) * MDP(14) + (-qJD(2) * t417 + qJDD(2) * t480 - t386 * t444 + t393 * t447 + 0.2e1 * t497 + t517) * MDP(15) + (-t352 * t480 + t353 * t484 - t386 * t393 - t398 * t416 - g(1) * (t511 * t549 + t469) - g(2) * (t508 * t549 + t467) - g(3) * t561 - t618 * t402) * MDP(16) + (t477 * t592 + t335 * t444 + t391 * t480 + t569 * t422 + (-t503 * t568 - t342) * t447 + t518 * t501) * MDP(17) + (-t477 * t593 - t336 * t444 + t392 * t480 + t569 * t424 + (t501 * t568 + t343) * t447 + t518 * t503) * MDP(18) + (t342 * t424 + t343 * t422 + (qJD(5) * t424 - t336 * t447 - t392 * t477 - t319) * t503 + (qJD(5) * t422 + t335 * t447 - t391 * t477 - t320) * t501 - t525) * MDP(19) + (t340 * t480 - t336 * t343 - t335 * t342 - g(1) * (-t511 * t613 + t469) - g(2) * (-t508 * t613 + t467) - g(3) * (t561 + t598) + t542 * t477 + t569 * t373 + t538 * qJD(5) + t603 * t623) * MDP(20) + (t331 * t620 + t536 * t578) * MDP(21) + (-t331 * t457 - t332 * t620 + t370 * t578 + t536 * t577) * MDP(22) + (-t594 + t626) * MDP(23) + (t546 - t595) * MDP(24) + ((-t442 * t506 + t443 * t509) * t407 + t463 * t332 + t324 * t457 + t316 * t444 + (t506 * t527 - t509 * t528) * t437 + t570 * t370 + t577 * t351 + t523 * t488) * MDP(26) + (-(t442 * t509 + t443 * t506) * t407 + t463 * t331 + t324 * t620 - t317 * t444 + (t506 * t528 + t509 * t527) * t437 - t570 * t536 - t578 * t351 + t523 * t490) * MDP(27) + (-MDP(4) * t507 * t510 + MDP(5) * t576) * t513; (t410 * t444 + t534) * MDP(12) + (-t502 * t557 - t504 * t558 + t473) * MDP(14) + (t436 - t533) * MDP(15) + (-t402 * t444 + t534 - t599 + t612) * MDP(16) + (t422 * t444 + t616) * MDP(17) + (t424 * t444 - t592) * MDP(18) + (-t391 * t503 + t392 * t501) * MDP(19) + (t373 * t444 + t541 - t621) * MDP(20) + (t546 + t595) * MDP(26) + (-t594 - t626) * MDP(27) + (t409 * MDP(12) - qJD(2) * MDP(14) + (-qJD(4) - t398) * MDP(16) + (t422 * t501 + t424 * t503) * MDP(19) + t538 * MDP(20) + t501 * MDP(18) * t447) * t447 + (-MDP(14) * t584 - t458 * MDP(15) - (MDP(12) + MDP(16)) * t485) * qJDD(1) + (MDP(11) + MDP(13)) * (-t444 ^ 2 - t614); (t436 + t412) * MDP(13) + (-t444 * t447 + qJDD(2)) * MDP(14) + (-t614 - t512) * MDP(15) + (qJD(2) * t402 + t514 - t597) * MDP(16) + (-qJD(2) * t422 - t501 * t614 + t592) * MDP(17) + (-qJD(2) * t424 + t616) * MDP(18) + (-t391 * t501 - t392 * t503 + (-t422 * t503 + t424 * t501) * t447) * MDP(19) + (-qJD(2) * t373 + t447 * t537 + t525 + t542) * MDP(20) + (-qJD(2) * t370 + t626) * MDP(26) + (qJD(2) * t536 + t546) * MDP(27); (t424 * t447 + t391) * MDP(17) + (-t422 * t447 + t392) * MDP(18) + (-t422 ^ 2 - t424 ^ 2) * MDP(19) + (t335 * t424 + t336 * t422 + t497 + t517 + t617) * MDP(20) + (t332 - t625) * MDP(26) + (t331 - t596) * MDP(27); -t536 * t370 * MDP(21) + (-t370 ^ 2 + t536 ^ 2) * MDP(22) + (t563 + t596) * MDP(23) + (-t552 - t625) * MDP(24) + t407 * MDP(25) + (-g(1) * t428 - g(2) * t430 + t317 * t437 + t351 * t536 + t482 * t490 + t553) * MDP(26) + (g(1) * t429 - g(2) * t431 + t316 * t437 + t351 * t370 - t482 * t488 - t543) * MDP(27) + (-MDP(23) * t591 + MDP(24) * t536 - MDP(26) * t317 - MDP(27) * t316) * qJD(6);];
tau  = t1;
