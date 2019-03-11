% Calculate kinetic energy for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:47:46
% EndTime: 2019-03-09 22:47:49
% DurationCPUTime: 3.22s
% Computational Cost: add. (3413->384), mult. (5903->578), div. (0->0), fcn. (7108->14), ass. (0->173)
t580 = Icges(5,2) + Icges(6,3);
t531 = sin(qJ(2));
t532 = sin(qJ(1));
t534 = cos(qJ(2));
t535 = cos(qJ(1));
t568 = cos(pkin(6));
t549 = t535 * t568;
t504 = t531 * t532 - t534 * t549;
t505 = t531 * t549 + t532 * t534;
t527 = sin(pkin(6));
t561 = t527 * t535;
t460 = Icges(3,5) * t505 - Icges(3,6) * t504 - Icges(3,3) * t561;
t550 = t532 * t568;
t506 = t535 * t531 + t534 * t550;
t507 = -t531 * t550 + t535 * t534;
t564 = t527 * t532;
t461 = Icges(3,5) * t507 - Icges(3,6) * t506 + Icges(3,3) * t564;
t579 = (t460 * t535 - t461 * t532) * t527;
t560 = qJ(3) + qJ(4);
t524 = sin(t560);
t551 = cos(t560);
t481 = t505 * t551 - t524 * t561;
t526 = sin(pkin(12));
t528 = cos(pkin(12));
t444 = -t481 * t526 + t504 * t528;
t567 = t504 * t526;
t445 = t481 * t528 + t567;
t547 = t527 * t551;
t480 = t505 * t524 + t535 * t547;
t578 = -Icges(5,4) * t481 + Icges(6,5) * t445 - Icges(5,6) * t504 + Icges(6,6) * t444 + t580 * t480;
t483 = t507 * t551 + t524 * t564;
t446 = -t483 * t526 + t506 * t528;
t566 = t506 * t526;
t447 = t483 * t528 + t566;
t482 = t507 * t524 - t532 * t547;
t577 = -Icges(5,4) * t483 + Icges(6,5) * t447 - Icges(5,6) * t506 + Icges(6,6) * t446 + t580 * t482;
t497 = t524 * t568 + t531 * t547;
t562 = t527 * t534;
t478 = -t497 * t526 - t528 * t562;
t555 = t526 * t562;
t479 = t497 * t528 - t555;
t565 = t527 * t531;
t496 = t524 * t565 - t551 * t568;
t576 = -Icges(5,4) * t497 + Icges(6,5) * t479 + Icges(5,6) * t562 + Icges(6,6) * t478 + t580 * t496;
t533 = cos(qJ(3));
t571 = pkin(3) * t533;
t570 = pkin(5) * t528;
t563 = t527 * t533;
t474 = pkin(2) * t505 + pkin(9) * t504;
t475 = pkin(2) * t507 + pkin(9) * t506;
t556 = qJD(2) * t527;
t517 = t532 * t556;
t552 = t535 * t556;
t558 = t474 * t517 + t475 * t552;
t488 = qJD(3) * t506 + t517;
t557 = qJD(1) * (pkin(1) * t532 - pkin(8) * t561);
t518 = qJD(2) * t568 + qJD(1);
t530 = sin(qJ(3));
t554 = t530 * t564;
t553 = t530 * t561;
t456 = qJD(4) * t506 + t488;
t548 = t568 * t530;
t489 = qJD(3) * t504 - t552;
t425 = -pkin(3) * t553 + pkin(10) * t504 + t505 * t571;
t426 = pkin(3) * t554 + pkin(10) * t506 + t507 * t571;
t545 = t488 * t425 - t426 * t489 + t558;
t457 = qJD(4) * t504 + t489;
t508 = (pkin(2) * t531 - pkin(9) * t534) * t527;
t510 = qJD(1) * (pkin(1) * t535 + pkin(8) * t564);
t544 = t518 * t475 - t508 * t517 + t510;
t490 = (-qJD(3) - qJD(4)) * t562 + t518;
t438 = pkin(4) * t481 + qJ(5) * t480;
t543 = qJD(5) * t496 + t456 * t438 + t545;
t542 = -t474 * t518 - t508 * t552 - t557;
t466 = pkin(3) * t548 + (-pkin(10) * t534 + t531 * t571) * t527;
t509 = -qJD(3) * t562 + t518;
t541 = t509 * t426 - t466 * t488 + t544;
t439 = pkin(4) * t483 + qJ(5) * t482;
t540 = qJD(5) * t480 + t490 * t439 + t541;
t539 = -t425 * t509 + t489 * t466 + t542;
t458 = pkin(4) * t497 + qJ(5) * t496;
t538 = qJD(5) * t482 + t457 * t458 + t539;
t525 = pkin(12) + qJ(6);
t523 = cos(t525);
t522 = sin(t525);
t514 = rSges(2,1) * t535 - rSges(2,2) * t532;
t513 = rSges(2,1) * t532 + rSges(2,2) * t535;
t503 = t531 * t563 + t548;
t502 = -t530 * t565 + t533 * t568;
t495 = t568 * rSges(3,3) + (rSges(3,1) * t531 + rSges(3,2) * t534) * t527;
t494 = Icges(3,5) * t568 + (Icges(3,1) * t531 + Icges(3,4) * t534) * t527;
t493 = Icges(3,6) * t568 + (Icges(3,4) * t531 + Icges(3,2) * t534) * t527;
t492 = Icges(3,3) * t568 + (Icges(3,5) * t531 + Icges(3,6) * t534) * t527;
t487 = t507 * t533 + t554;
t486 = -t507 * t530 + t532 * t563;
t485 = t505 * t533 - t553;
t484 = -t505 * t530 - t533 * t561;
t473 = t497 * t523 - t522 * t562;
t472 = -t497 * t522 - t523 * t562;
t469 = rSges(3,1) * t507 - rSges(3,2) * t506 + rSges(3,3) * t564;
t468 = rSges(3,1) * t505 - rSges(3,2) * t504 - rSges(3,3) * t561;
t465 = Icges(3,1) * t507 - Icges(3,4) * t506 + Icges(3,5) * t564;
t464 = Icges(3,1) * t505 - Icges(3,4) * t504 - Icges(3,5) * t561;
t463 = Icges(3,4) * t507 - Icges(3,2) * t506 + Icges(3,6) * t564;
t462 = Icges(3,4) * t505 - Icges(3,2) * t504 - Icges(3,6) * t561;
t459 = rSges(4,1) * t503 + rSges(4,2) * t502 - rSges(4,3) * t562;
t455 = Icges(4,1) * t503 + Icges(4,4) * t502 - Icges(4,5) * t562;
t454 = Icges(4,4) * t503 + Icges(4,2) * t502 - Icges(4,6) * t562;
t453 = Icges(4,5) * t503 + Icges(4,6) * t502 - Icges(4,3) * t562;
t452 = qJD(6) * t496 + t490;
t451 = rSges(5,1) * t497 - rSges(5,2) * t496 - rSges(5,3) * t562;
t450 = Icges(5,1) * t497 - Icges(5,4) * t496 - Icges(5,5) * t562;
t448 = Icges(5,5) * t497 - Icges(5,6) * t496 - Icges(5,3) * t562;
t443 = t483 * t523 + t506 * t522;
t442 = -t483 * t522 + t506 * t523;
t441 = t481 * t523 + t504 * t522;
t440 = -t481 * t522 + t504 * t523;
t436 = qJD(6) * t480 + t457;
t435 = qJD(6) * t482 + t456;
t434 = rSges(4,1) * t487 + rSges(4,2) * t486 + rSges(4,3) * t506;
t433 = rSges(4,1) * t485 + rSges(4,2) * t484 + rSges(4,3) * t504;
t432 = Icges(4,1) * t487 + Icges(4,4) * t486 + Icges(4,5) * t506;
t431 = Icges(4,1) * t485 + Icges(4,4) * t484 + Icges(4,5) * t504;
t430 = Icges(4,4) * t487 + Icges(4,2) * t486 + Icges(4,6) * t506;
t429 = Icges(4,4) * t485 + Icges(4,2) * t484 + Icges(4,6) * t504;
t428 = Icges(4,5) * t487 + Icges(4,6) * t486 + Icges(4,3) * t506;
t427 = Icges(4,5) * t485 + Icges(4,6) * t484 + Icges(4,3) * t504;
t424 = rSges(5,1) * t483 - rSges(5,2) * t482 + rSges(5,3) * t506;
t423 = rSges(5,1) * t481 - rSges(5,2) * t480 + rSges(5,3) * t504;
t422 = Icges(5,1) * t483 - Icges(5,4) * t482 + Icges(5,5) * t506;
t421 = Icges(5,1) * t481 - Icges(5,4) * t480 + Icges(5,5) * t504;
t418 = Icges(5,5) * t483 - Icges(5,6) * t482 + Icges(5,3) * t506;
t417 = Icges(5,5) * t481 - Icges(5,6) * t480 + Icges(5,3) * t504;
t414 = rSges(6,1) * t479 + rSges(6,2) * t478 + rSges(6,3) * t496;
t413 = Icges(6,1) * t479 + Icges(6,4) * t478 + Icges(6,5) * t496;
t412 = Icges(6,4) * t479 + Icges(6,2) * t478 + Icges(6,6) * t496;
t410 = t469 * t518 - t495 * t517 + t510;
t409 = -t468 * t518 - t495 * t552 - t557;
t408 = rSges(7,1) * t473 + rSges(7,2) * t472 + rSges(7,3) * t496;
t407 = Icges(7,1) * t473 + Icges(7,4) * t472 + Icges(7,5) * t496;
t406 = Icges(7,4) * t473 + Icges(7,2) * t472 + Icges(7,6) * t496;
t405 = Icges(7,5) * t473 + Icges(7,6) * t472 + Icges(7,3) * t496;
t403 = (t468 * t532 + t469 * t535) * t556;
t402 = -pkin(5) * t555 + pkin(11) * t496 + t497 * t570;
t399 = rSges(6,1) * t447 + rSges(6,2) * t446 + rSges(6,3) * t482;
t398 = rSges(6,1) * t445 + rSges(6,2) * t444 + rSges(6,3) * t480;
t397 = Icges(6,1) * t447 + Icges(6,4) * t446 + Icges(6,5) * t482;
t396 = Icges(6,1) * t445 + Icges(6,4) * t444 + Icges(6,5) * t480;
t395 = Icges(6,4) * t447 + Icges(6,2) * t446 + Icges(6,6) * t482;
t394 = Icges(6,4) * t445 + Icges(6,2) * t444 + Icges(6,6) * t480;
t391 = rSges(7,1) * t443 + rSges(7,2) * t442 + rSges(7,3) * t482;
t390 = rSges(7,1) * t441 + rSges(7,2) * t440 + rSges(7,3) * t480;
t389 = Icges(7,1) * t443 + Icges(7,4) * t442 + Icges(7,5) * t482;
t388 = Icges(7,1) * t441 + Icges(7,4) * t440 + Icges(7,5) * t480;
t387 = Icges(7,4) * t443 + Icges(7,2) * t442 + Icges(7,6) * t482;
t386 = Icges(7,4) * t441 + Icges(7,2) * t440 + Icges(7,6) * t480;
t385 = Icges(7,5) * t443 + Icges(7,6) * t442 + Icges(7,3) * t482;
t384 = Icges(7,5) * t441 + Icges(7,6) * t440 + Icges(7,3) * t480;
t383 = pkin(5) * t566 + pkin(11) * t482 + t483 * t570;
t382 = pkin(5) * t567 + pkin(11) * t480 + t481 * t570;
t381 = t434 * t509 - t459 * t488 + t544;
t380 = -t433 * t509 + t459 * t489 + t542;
t379 = t433 * t488 - t434 * t489 + t558;
t378 = t424 * t490 - t451 * t456 + t541;
t377 = -t423 * t490 + t451 * t457 + t539;
t376 = t423 * t456 - t424 * t457 + t545;
t375 = t399 * t490 + (-t414 - t458) * t456 + t540;
t374 = t414 * t457 + (-t398 - t438) * t490 + t538;
t373 = t398 * t456 + (-t399 - t439) * t457 + t543;
t372 = t383 * t490 + t391 * t452 - t408 * t435 + (-t402 - t458) * t456 + t540;
t371 = -t390 * t452 + t402 * t457 + t408 * t436 + (-t382 - t438) * t490 + t538;
t370 = t382 * t456 + t390 * t435 - t391 * t436 + (-t383 - t439) * t457 + t543;
t1 = -((-t492 * t561 - t493 * t504 + t494 * t505) * t518 + ((-t463 * t504 + t465 * t505) * t532 + (t462 * t504 - t464 * t505 + t579) * t535) * t556) * t552 / 0.2e1 + ((t492 * t564 - t493 * t506 + t494 * t507) * t518 + (-(-t462 * t506 + t464 * t507) * t535 + (-t463 * t506 + t465 * t507 - t579) * t532) * t556) * t517 / 0.2e1 + m(6) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(7) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(5) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(4) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + m(3) * (t403 ^ 2 + t409 ^ 2 + t410 ^ 2) / 0.2e1 + t452 * ((t385 * t496 + t387 * t472 + t389 * t473) * t435 + (t384 * t496 + t386 * t472 + t388 * t473) * t436 + (t405 * t496 + t472 * t406 + t473 * t407) * t452) / 0.2e1 + t436 * ((t385 * t480 + t387 * t440 + t389 * t441) * t435 + (t384 * t480 + t386 * t440 + t388 * t441) * t436 + (t405 * t480 + t406 * t440 + t407 * t441) * t452) / 0.2e1 + t435 * ((t385 * t482 + t387 * t442 + t389 * t443) * t435 + (t384 * t482 + t386 * t442 + t388 * t443) * t436 + (t405 * t482 + t406 * t442 + t407 * t443) * t452) / 0.2e1 + t489 * ((t428 * t504 + t430 * t484 + t432 * t485) * t488 + (t427 * t504 + t429 * t484 + t431 * t485) * t489 + (t453 * t504 + t454 * t484 + t455 * t485) * t509) / 0.2e1 + t488 * ((t428 * t506 + t430 * t486 + t432 * t487) * t488 + (t427 * t506 + t429 * t486 + t431 * t487) * t489 + (t453 * t506 + t454 * t486 + t455 * t487) * t509) / 0.2e1 + t509 * ((-t428 * t562 + t430 * t502 + t432 * t503) * t488 + (-t427 * t562 + t429 * t502 + t431 * t503) * t489 + (-t453 * t562 + t454 * t502 + t455 * t503) * t509) / 0.2e1 + t518 * ((t568 * t461 + (t463 * t534 + t465 * t531) * t527) * t517 - (t568 * t460 + (t462 * t534 + t464 * t531) * t527) * t552 + (t568 * t492 + (t493 * t534 + t494 * t531) * t527) * t518) / 0.2e1 + ((t412 * t446 + t413 * t447 + t448 * t506 + t450 * t483 + t576 * t482) * t490 + (t394 * t446 + t396 * t447 + t417 * t506 + t421 * t483 + t578 * t482) * t457 + (t395 * t446 + t397 * t447 + t418 * t506 + t422 * t483 + t577 * t482) * t456) * t456 / 0.2e1 + ((t412 * t444 + t413 * t445 + t448 * t504 + t450 * t481 + t576 * t480) * t490 + (t394 * t444 + t396 * t445 + t417 * t504 + t481 * t421 + t578 * t480) * t457 + (t395 * t444 + t397 * t445 + t418 * t504 + t422 * t481 + t577 * t480) * t456) * t457 / 0.2e1 + ((t478 * t412 + t479 * t413 - t448 * t562 + t450 * t497 + t576 * t496) * t490 + (t394 * t478 + t396 * t479 - t417 * t562 + t421 * t497 + t578 * t496) * t457 + (t395 * t478 + t397 * t479 - t418 * t562 + t422 * t497 + t577 * t496) * t456) * t490 / 0.2e1 + (Icges(2,3) + m(2) * (t513 ^ 2 + t514 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
