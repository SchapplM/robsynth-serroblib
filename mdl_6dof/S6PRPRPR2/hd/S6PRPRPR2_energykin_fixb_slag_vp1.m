% Calculate kinetic energy for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:31
% EndTime: 2019-03-08 19:29:34
% DurationCPUTime: 3.31s
% Computational Cost: add. (3542->369), mult. (8591->555), div. (0->0), fcn. (10934->14), ass. (0->167)
t588 = Icges(5,2) + Icges(6,3);
t533 = sin(pkin(10));
t536 = cos(pkin(10));
t540 = sin(qJ(2));
t571 = sin(pkin(11));
t572 = cos(pkin(11));
t576 = cos(qJ(2));
t518 = -t540 * t571 + t576 * t572;
t537 = cos(pkin(6));
t545 = t537 * t518;
t546 = t540 * t572 + t571 * t576;
t495 = -t533 * t546 + t536 * t545;
t511 = t546 * t537;
t496 = t511 * t536 + t518 * t533;
t534 = sin(pkin(6));
t566 = t534 * t536;
t442 = Icges(4,5) * t496 + Icges(4,6) * t495 - Icges(4,3) * t566;
t497 = -t533 * t545 - t536 * t546;
t498 = -t511 * t533 + t518 * t536;
t567 = t533 * t534;
t443 = Icges(4,5) * t498 + Icges(4,6) * t497 + Icges(4,3) * t567;
t509 = t518 * t534;
t510 = t546 * t534;
t478 = Icges(4,5) * t510 + Icges(4,6) * t509 + Icges(4,3) * t537;
t557 = t537 * t576;
t513 = -t533 * t540 + t536 * t557;
t564 = t537 * t540;
t514 = t533 * t576 + t536 * t564;
t482 = Icges(3,5) * t514 + Icges(3,6) * t513 - Icges(3,3) * t566;
t515 = -t533 * t557 - t536 * t540;
t516 = -t533 * t564 + t536 * t576;
t483 = Icges(3,5) * t516 + Icges(3,6) * t515 + Icges(3,3) * t567;
t505 = Icges(3,3) * t537 + (Icges(3,5) * t540 + Icges(3,6) * t576) * t534;
t587 = -(t478 + t505) * t537 + (t482 + t442) * t566 - (t483 + t443) * t567;
t539 = sin(qJ(4));
t565 = t534 * t539;
t575 = cos(qJ(4));
t470 = t496 * t575 - t536 * t565;
t532 = sin(pkin(12));
t535 = cos(pkin(12));
t438 = -t470 * t532 - t495 * t535;
t570 = t495 * t532;
t439 = t470 * t535 - t570;
t558 = t534 * t575;
t469 = t496 * t539 + t536 * t558;
t586 = -Icges(5,4) * t470 + Icges(6,5) * t439 + Icges(5,6) * t495 + Icges(6,6) * t438 + t588 * t469;
t472 = t498 * t575 + t533 * t565;
t440 = -t472 * t532 - t497 * t535;
t569 = t497 * t532;
t441 = t472 * t535 - t569;
t471 = t498 * t539 - t533 * t558;
t585 = -Icges(5,4) * t472 + Icges(6,5) * t441 + Icges(5,6) * t497 + Icges(6,6) * t440 + t588 * t471;
t501 = t510 * t575 + t537 * t539;
t465 = -t501 * t532 - t509 * t535;
t568 = t509 * t532;
t466 = t501 * t535 - t568;
t500 = t510 * t539 - t537 * t575;
t584 = -Icges(5,4) * t501 + Icges(6,5) * t466 + Icges(5,6) * t509 + Icges(6,6) * t465 + t588 * t500;
t580 = qJD(2) ^ 2;
t574 = pkin(2) * t576;
t573 = pkin(5) * t535;
t519 = pkin(2) * t534 * t540 + t537 * qJ(3);
t562 = -pkin(3) * t510 + pkin(8) * t509 - t519;
t561 = qJD(2) * t534;
t523 = t533 * t561;
t473 = -qJD(4) * t497 + t523;
t528 = qJD(2) * t537;
t502 = -qJD(4) * t509 + t528;
t560 = qJD(3) * t536;
t556 = t536 * t561;
t555 = pkin(2) * t564 - qJ(3) * t534;
t553 = (-rSges(4,1) * t510 - rSges(4,2) * t509 - rSges(4,3) * t537 - t519) * t534;
t490 = t533 * t574 + t536 * t555;
t491 = -t533 * t555 + t536 * t574;
t552 = qJD(3) * t537 + t490 * t523 + t491 * t556 + qJD(1);
t474 = -qJD(4) * t495 - t556;
t453 = pkin(3) * t496 - pkin(8) * t495;
t454 = pkin(3) * t498 - pkin(8) * t497;
t548 = t453 * t523 + t454 * t556 + t552;
t429 = pkin(4) * t470 + qJ(5) * t469;
t547 = qJD(5) * t500 + t473 * t429 + t548;
t477 = t491 * t528;
t544 = t454 * t528 + t477 + (qJD(2) * t533 * t562 - t560) * t534;
t522 = qJD(3) * t567;
t543 = t522 + ((-t453 - t490) * t537 + t562 * t566) * qJD(2);
t430 = pkin(4) * t472 + qJ(5) * t471;
t542 = qJD(5) * t469 + t502 * t430 + t544;
t461 = pkin(4) * t501 + qJ(5) * t500;
t541 = qJD(5) * t471 + t474 * t461 + t543;
t531 = pkin(12) + qJ(6);
t530 = cos(t531);
t529 = sin(t531);
t508 = t537 * rSges(3,3) + (rSges(3,1) * t540 + rSges(3,2) * t576) * t534;
t507 = Icges(3,5) * t537 + (Icges(3,1) * t540 + Icges(3,4) * t576) * t534;
t506 = Icges(3,6) * t537 + (Icges(3,4) * t540 + Icges(3,2) * t576) * t534;
t489 = rSges(3,1) * t516 + rSges(3,2) * t515 + rSges(3,3) * t567;
t488 = rSges(3,1) * t514 + rSges(3,2) * t513 - rSges(3,3) * t566;
t487 = Icges(3,1) * t516 + Icges(3,4) * t515 + Icges(3,5) * t567;
t486 = Icges(3,1) * t514 + Icges(3,4) * t513 - Icges(3,5) * t566;
t485 = Icges(3,4) * t516 + Icges(3,2) * t515 + Icges(3,6) * t567;
t484 = Icges(3,4) * t514 + Icges(3,2) * t513 - Icges(3,6) * t566;
t480 = Icges(4,1) * t510 + Icges(4,4) * t509 + Icges(4,5) * t537;
t479 = Icges(4,4) * t510 + Icges(4,2) * t509 + Icges(4,6) * t537;
t464 = qJD(6) * t500 + t502;
t463 = t501 * t530 - t509 * t529;
t462 = -t501 * t529 - t509 * t530;
t460 = (-t488 * t537 - t508 * t566) * qJD(2);
t459 = (t489 * t537 - t508 * t567) * qJD(2);
t458 = rSges(5,1) * t501 - rSges(5,2) * t500 - rSges(5,3) * t509;
t457 = Icges(5,1) * t501 - Icges(5,4) * t500 - Icges(5,5) * t509;
t455 = Icges(5,5) * t501 - Icges(5,6) * t500 - Icges(5,3) * t509;
t449 = rSges(4,1) * t498 + rSges(4,2) * t497 + rSges(4,3) * t567;
t448 = rSges(4,1) * t496 + rSges(4,2) * t495 - rSges(4,3) * t566;
t447 = Icges(4,1) * t498 + Icges(4,4) * t497 + Icges(4,5) * t567;
t446 = Icges(4,1) * t496 + Icges(4,4) * t495 - Icges(4,5) * t566;
t445 = Icges(4,4) * t498 + Icges(4,2) * t497 + Icges(4,6) * t567;
t444 = Icges(4,4) * t496 + Icges(4,2) * t495 - Icges(4,6) * t566;
t437 = t472 * t530 - t497 * t529;
t436 = -t472 * t529 - t497 * t530;
t435 = t470 * t530 - t495 * t529;
t434 = -t470 * t529 - t495 * t530;
t433 = qJD(1) + (t488 * t533 + t489 * t536) * t561;
t432 = qJD(6) * t469 + t474;
t431 = qJD(6) * t471 + t473;
t426 = rSges(6,1) * t466 + rSges(6,2) * t465 + rSges(6,3) * t500;
t425 = Icges(6,1) * t466 + Icges(6,4) * t465 + Icges(6,5) * t500;
t424 = Icges(6,4) * t466 + Icges(6,2) * t465 + Icges(6,6) * t500;
t422 = rSges(5,1) * t472 - rSges(5,2) * t471 - rSges(5,3) * t497;
t421 = rSges(5,1) * t470 - rSges(5,2) * t469 - rSges(5,3) * t495;
t420 = Icges(5,1) * t472 - Icges(5,4) * t471 - Icges(5,5) * t497;
t419 = Icges(5,1) * t470 - Icges(5,4) * t469 - Icges(5,5) * t495;
t416 = Icges(5,5) * t472 - Icges(5,6) * t471 - Icges(5,3) * t497;
t415 = Icges(5,5) * t470 - Icges(5,6) * t469 - Icges(5,3) * t495;
t414 = -pkin(5) * t568 + pkin(9) * t500 + t501 * t573;
t413 = rSges(7,1) * t463 + rSges(7,2) * t462 + rSges(7,3) * t500;
t412 = Icges(7,1) * t463 + Icges(7,4) * t462 + Icges(7,5) * t500;
t411 = Icges(7,4) * t463 + Icges(7,2) * t462 + Icges(7,6) * t500;
t410 = Icges(7,5) * t463 + Icges(7,6) * t462 + Icges(7,3) * t500;
t408 = t522 + ((-t448 - t490) * t537 + t536 * t553) * qJD(2);
t407 = -t534 * t560 + t477 + (t449 * t537 + t533 * t553) * qJD(2);
t406 = rSges(6,1) * t441 + rSges(6,2) * t440 + rSges(6,3) * t471;
t405 = rSges(6,1) * t439 + rSges(6,2) * t438 + rSges(6,3) * t469;
t404 = Icges(6,1) * t441 + Icges(6,4) * t440 + Icges(6,5) * t471;
t403 = Icges(6,1) * t439 + Icges(6,4) * t438 + Icges(6,5) * t469;
t402 = Icges(6,4) * t441 + Icges(6,2) * t440 + Icges(6,6) * t471;
t401 = Icges(6,4) * t439 + Icges(6,2) * t438 + Icges(6,6) * t469;
t398 = rSges(7,1) * t437 + rSges(7,2) * t436 + rSges(7,3) * t471;
t397 = rSges(7,1) * t435 + rSges(7,2) * t434 + rSges(7,3) * t469;
t396 = Icges(7,1) * t437 + Icges(7,4) * t436 + Icges(7,5) * t471;
t395 = Icges(7,1) * t435 + Icges(7,4) * t434 + Icges(7,5) * t469;
t394 = Icges(7,4) * t437 + Icges(7,2) * t436 + Icges(7,6) * t471;
t393 = Icges(7,4) * t435 + Icges(7,2) * t434 + Icges(7,6) * t469;
t392 = Icges(7,5) * t437 + Icges(7,6) * t436 + Icges(7,3) * t471;
t391 = Icges(7,5) * t435 + Icges(7,6) * t434 + Icges(7,3) * t469;
t390 = -pkin(5) * t569 + pkin(9) * t471 + t472 * t573;
t389 = -pkin(5) * t570 + pkin(9) * t469 + t470 * t573;
t388 = (t448 * t533 + t449 * t536) * t561 + t552;
t387 = -t421 * t502 + t458 * t474 + t543;
t386 = t422 * t502 - t458 * t473 + t544;
t385 = t421 * t473 - t422 * t474 + t548;
t384 = t426 * t474 + (-t405 - t429) * t502 + t541;
t383 = t406 * t502 + (-t426 - t461) * t473 + t542;
t382 = t405 * t473 + (-t406 - t430) * t474 + t547;
t381 = -t397 * t464 + t413 * t432 + t414 * t474 + (-t389 - t429) * t502 + t541;
t380 = t390 * t502 + t398 * t464 - t413 * t431 + (-t414 - t461) * t473 + t542;
t379 = t389 * t473 + t397 * t431 - t398 * t432 + (-t390 - t430) * t474 + t547;
t1 = ((t537 * t483 + (t485 * t576 + t487 * t540) * t534) * t523 - (t537 * t482 + (t484 * t576 + t486 * t540) * t534) * t556 + (t537 * t505 + (t506 * t576 + t507 * t540) * t534) * t528) * t528 / 0.2e1 + t431 * ((t471 * t392 + t436 * t394 + t437 * t396) * t431 + (t391 * t471 + t393 * t436 + t395 * t437) * t432 + (t410 * t471 + t411 * t436 + t412 * t437) * t464) / 0.2e1 + t432 * ((t392 * t469 + t394 * t434 + t396 * t435) * t431 + (t469 * t391 + t434 * t393 + t435 * t395) * t432 + (t410 * t469 + t411 * t434 + t412 * t435) * t464) / 0.2e1 + t464 * ((t392 * t500 + t394 * t462 + t396 * t463) * t431 + (t391 * t500 + t393 * t462 + t395 * t463) * t432 + (t500 * t410 + t462 * t411 + t463 * t412) * t464) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t433 ^ 2 + t459 ^ 2 + t460 ^ 2) / 0.2e1 + m(4) * (t388 ^ 2 + t407 ^ 2 + t408 ^ 2) / 0.2e1 + m(5) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(6) * (t382 ^ 2 + t383 ^ 2 + t384 ^ 2) / 0.2e1 + m(7) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + ((t424 * t440 + t425 * t441 - t455 * t497 + t457 * t472 + t471 * t584) * t502 + (t401 * t440 + t403 * t441 - t415 * t497 + t419 * t472 + t471 * t586) * t474 + (t440 * t402 + t441 * t404 - t497 * t416 + t472 * t420 + t471 * t585) * t473) * t473 / 0.2e1 + ((t424 * t438 + t425 * t439 - t455 * t495 + t457 * t470 + t469 * t584) * t502 + (t438 * t401 + t439 * t403 - t495 * t415 + t470 * t419 + t586 * t469) * t474 + (t402 * t438 + t404 * t439 - t416 * t495 + t420 * t470 + t469 * t585) * t473) * t474 / 0.2e1 + ((t424 * t465 + t425 * t466 - t455 * t509 + t457 * t501 + t500 * t584) * t502 + (t401 * t465 + t403 * t466 - t415 * t509 + t419 * t501 + t500 * t586) * t474 + (t402 * t465 + t404 * t466 - t416 * t509 + t420 * t501 + t500 * t585) * t473) * t502 / 0.2e1 - ((t445 * t495 + t447 * t496 + t485 * t513 + t487 * t514) * t567 + (t479 * t495 + t480 * t496 + t506 * t513 + t507 * t514) * t537 + (-t444 * t495 - t446 * t496 - t484 * t513 - t486 * t514 + t587) * t566) * t580 * t566 / 0.2e1 + (t537 * ((t478 * t537 + t479 * t509 + t480 * t510) * t537 + ((t443 * t537 + t445 * t509 + t447 * t510) * t533 - (t442 * t537 + t444 * t509 + t446 * t510) * t536) * t534) + ((-t444 * t497 - t446 * t498 - t484 * t515 - t486 * t516) * t566 + (t479 * t497 + t480 * t498 + t506 * t515 + t507 * t516) * t537 + (t445 * t497 + t447 * t498 + t485 * t515 + t487 * t516 - t587) * t567) * t567) * t580 / 0.2e1;
T  = t1;
