% Calculate kinetic energy for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:18
% EndTime: 2019-03-09 14:08:21
% DurationCPUTime: 2.80s
% Computational Cost: add. (3363->387), mult. (5235->590), div. (0->0), fcn. (6180->14), ass. (0->177)
t531 = sin(pkin(6));
t533 = cos(pkin(6));
t539 = cos(qJ(2));
t540 = cos(qJ(1));
t569 = t539 * t540;
t536 = sin(qJ(2));
t537 = sin(qJ(1));
t572 = t536 * t537;
t505 = -t533 * t569 + t572;
t570 = t537 * t539;
t571 = t536 * t540;
t506 = t533 * t571 + t570;
t507 = t533 * t570 + t571;
t508 = -t533 * t572 + t569;
t573 = t531 * t540;
t575 = t531 * t537;
t549 = (Icges(3,5) * t506 - Icges(3,6) * t505 - Icges(3,3) * t573) * t540 - (Icges(3,5) * t508 - Icges(3,6) * t507 + Icges(3,3) * t575) * t537;
t585 = t531 * t549;
t530 = sin(pkin(12));
t532 = cos(pkin(12));
t482 = -t506 * t530 - t532 * t573;
t560 = t530 * t573;
t483 = t506 * t532 - t560;
t427 = Icges(4,5) * t483 + Icges(4,6) * t482 + Icges(4,3) * t505;
t461 = Icges(3,4) * t506 - Icges(3,2) * t505 - Icges(3,6) * t573;
t584 = -t427 + t461;
t484 = -t508 * t530 + t532 * t575;
t561 = t530 * t575;
t485 = t508 * t532 + t561;
t428 = Icges(4,5) * t485 + Icges(4,6) * t484 + Icges(4,3) * t507;
t462 = Icges(3,4) * t508 - Icges(3,2) * t507 + Icges(3,6) * t575;
t583 = t428 - t462;
t576 = t531 * t536;
t503 = -t530 * t576 + t532 * t533;
t577 = t530 * t533;
t504 = t532 * t576 + t577;
t574 = t531 * t539;
t452 = Icges(4,5) * t504 + Icges(4,6) * t503 - Icges(4,3) * t574;
t492 = Icges(3,6) * t533 + (Icges(3,4) * t536 + Icges(3,2) * t539) * t531;
t582 = t452 - t492;
t578 = t532 * pkin(3);
t476 = pkin(2) * t506 + qJ(3) * t505;
t477 = pkin(2) * t508 + qJ(3) * t507;
t563 = qJD(2) * t531;
t520 = t537 * t563;
t558 = t540 * t563;
t567 = t476 * t520 + t477 * t558;
t529 = pkin(12) + qJ(4);
t526 = cos(t529);
t566 = pkin(4) * t526;
t486 = qJD(4) * t507 + t520;
t564 = qJD(1) * (pkin(1) * t537 - pkin(8) * t573);
t562 = qJD(3) * t539;
t521 = qJD(2) * t533 + qJD(1);
t511 = qJD(1) * (pkin(1) * t540 + pkin(8) * t575);
t559 = qJD(3) * t505 + t521 * t477 + t511;
t457 = qJD(5) * t507 + t486;
t557 = qJ(5) + t529;
t525 = sin(t529);
t556 = pkin(4) * t525;
t555 = qJD(3) * t507 - t564;
t552 = cos(t557);
t509 = (pkin(2) * t536 - qJ(3) * t539) * t531;
t551 = (-rSges(4,1) * t504 - rSges(4,2) * t503 + rSges(4,3) * t574 - t509) * t563;
t550 = (-pkin(3) * t577 - (-pkin(9) * t539 + t536 * t578) * t531 - t509) * t563;
t487 = qJD(4) * t505 - t558;
t458 = qJD(5) * t505 + t487;
t548 = t531 * t552;
t490 = (-qJD(4) - qJD(5)) * t574 + t521;
t423 = -pkin(3) * t560 + pkin(9) * t505 + t506 * t578;
t424 = pkin(3) * t561 + pkin(9) * t507 + t508 * t578;
t547 = t423 * t520 + t424 * t558 - t531 * t562 + t567;
t546 = t521 * t424 + t537 * t550 + t559;
t394 = pkin(10) * t505 + t506 * t566 - t556 * t573;
t395 = pkin(10) * t507 + t508 * t566 + t556 * t575;
t545 = t486 * t394 - t395 * t487 + t547;
t437 = t556 * t533 + (-pkin(10) * t539 + t536 * t566) * t531;
t510 = -qJD(4) * t574 + t521;
t544 = t510 * t395 - t437 * t486 + t546;
t543 = (-t423 - t476) * t521 + t540 * t550 + t555;
t542 = -t394 * t510 + t487 * t437 + t543;
t538 = cos(qJ(6));
t535 = sin(qJ(6));
t522 = sin(t557);
t518 = rSges(2,1) * t540 - rSges(2,2) * t537;
t517 = rSges(2,1) * t537 + rSges(2,2) * t540;
t496 = rSges(3,3) * t533 + (rSges(3,1) * t536 + rSges(3,2) * t539) * t531;
t495 = t525 * t533 + t526 * t576;
t494 = -t525 * t576 + t526 * t533;
t493 = Icges(3,5) * t533 + (Icges(3,1) * t536 + Icges(3,4) * t539) * t531;
t491 = Icges(3,3) * t533 + (Icges(3,5) * t536 + Icges(3,6) * t539) * t531;
t489 = t533 * t522 + t536 * t548;
t488 = t522 * t576 - t533 * t552;
t481 = t508 * t526 + t525 * t575;
t480 = -t508 * t525 + t526 * t575;
t479 = t506 * t526 - t525 * t573;
t478 = -t506 * t525 - t526 * t573;
t475 = t508 * t552 + t522 * t575;
t474 = t508 * t522 - t537 * t548;
t473 = t506 * t552 - t522 * t573;
t472 = t506 * t522 + t540 * t548;
t471 = t489 * t538 - t535 * t574;
t470 = -t489 * t535 - t538 * t574;
t469 = rSges(3,1) * t508 - rSges(3,2) * t507 + rSges(3,3) * t575;
t468 = rSges(3,1) * t506 - rSges(3,2) * t505 - rSges(3,3) * t573;
t464 = Icges(3,1) * t508 - Icges(3,4) * t507 + Icges(3,5) * t575;
t463 = Icges(3,1) * t506 - Icges(3,4) * t505 - Icges(3,5) * t573;
t454 = Icges(4,1) * t504 + Icges(4,4) * t503 - Icges(4,5) * t574;
t453 = Icges(4,4) * t504 + Icges(4,2) * t503 - Icges(4,6) * t574;
t451 = qJD(6) * t488 + t490;
t450 = pkin(5) * t489 + pkin(11) * t488;
t449 = rSges(5,1) * t495 + rSges(5,2) * t494 - rSges(5,3) * t574;
t448 = Icges(5,1) * t495 + Icges(5,4) * t494 - Icges(5,5) * t574;
t447 = Icges(5,4) * t495 + Icges(5,2) * t494 - Icges(5,6) * t574;
t446 = Icges(5,5) * t495 + Icges(5,6) * t494 - Icges(5,3) * t574;
t445 = rSges(6,1) * t489 - rSges(6,2) * t488 - rSges(6,3) * t574;
t444 = Icges(6,1) * t489 - Icges(6,4) * t488 - Icges(6,5) * t574;
t443 = Icges(6,4) * t489 - Icges(6,2) * t488 - Icges(6,6) * t574;
t442 = Icges(6,5) * t489 - Icges(6,6) * t488 - Icges(6,3) * t574;
t441 = t475 * t538 + t507 * t535;
t440 = -t475 * t535 + t507 * t538;
t439 = t473 * t538 + t505 * t535;
t438 = -t473 * t535 + t505 * t538;
t436 = pkin(5) * t475 + pkin(11) * t474;
t435 = pkin(5) * t473 + pkin(11) * t472;
t434 = rSges(4,1) * t485 + rSges(4,2) * t484 + rSges(4,3) * t507;
t433 = rSges(4,1) * t483 + rSges(4,2) * t482 + rSges(4,3) * t505;
t432 = Icges(4,1) * t485 + Icges(4,4) * t484 + Icges(4,5) * t507;
t431 = Icges(4,1) * t483 + Icges(4,4) * t482 + Icges(4,5) * t505;
t430 = Icges(4,4) * t485 + Icges(4,2) * t484 + Icges(4,6) * t507;
t429 = Icges(4,4) * t483 + Icges(4,2) * t482 + Icges(4,6) * t505;
t426 = qJD(6) * t472 + t458;
t425 = qJD(6) * t474 + t457;
t422 = rSges(5,1) * t481 + rSges(5,2) * t480 + rSges(5,3) * t507;
t421 = rSges(5,1) * t479 + rSges(5,2) * t478 + rSges(5,3) * t505;
t420 = Icges(5,1) * t481 + Icges(5,4) * t480 + Icges(5,5) * t507;
t419 = Icges(5,1) * t479 + Icges(5,4) * t478 + Icges(5,5) * t505;
t418 = Icges(5,4) * t481 + Icges(5,2) * t480 + Icges(5,6) * t507;
t417 = Icges(5,4) * t479 + Icges(5,2) * t478 + Icges(5,6) * t505;
t416 = Icges(5,5) * t481 + Icges(5,6) * t480 + Icges(5,3) * t507;
t415 = Icges(5,5) * t479 + Icges(5,6) * t478 + Icges(5,3) * t505;
t414 = rSges(6,1) * t475 - rSges(6,2) * t474 + rSges(6,3) * t507;
t413 = rSges(6,1) * t473 - rSges(6,2) * t472 + rSges(6,3) * t505;
t412 = Icges(6,1) * t475 - Icges(6,4) * t474 + Icges(6,5) * t507;
t411 = Icges(6,1) * t473 - Icges(6,4) * t472 + Icges(6,5) * t505;
t410 = Icges(6,4) * t475 - Icges(6,2) * t474 + Icges(6,6) * t507;
t409 = Icges(6,4) * t473 - Icges(6,2) * t472 + Icges(6,6) * t505;
t408 = Icges(6,5) * t475 - Icges(6,6) * t474 + Icges(6,3) * t507;
t407 = Icges(6,5) * t473 - Icges(6,6) * t472 + Icges(6,3) * t505;
t406 = t469 * t521 - t496 * t520 + t511;
t405 = -t468 * t521 - t496 * t558 - t564;
t401 = rSges(7,1) * t471 + rSges(7,2) * t470 + rSges(7,3) * t488;
t399 = Icges(7,1) * t471 + Icges(7,4) * t470 + Icges(7,5) * t488;
t398 = Icges(7,4) * t471 + Icges(7,2) * t470 + Icges(7,6) * t488;
t397 = Icges(7,5) * t471 + Icges(7,6) * t470 + Icges(7,3) * t488;
t396 = (t468 * t537 + t469 * t540) * t563;
t391 = rSges(7,1) * t441 + rSges(7,2) * t440 + rSges(7,3) * t474;
t390 = rSges(7,1) * t439 + rSges(7,2) * t438 + rSges(7,3) * t472;
t389 = Icges(7,1) * t441 + Icges(7,4) * t440 + Icges(7,5) * t474;
t388 = Icges(7,1) * t439 + Icges(7,4) * t438 + Icges(7,5) * t472;
t387 = Icges(7,4) * t441 + Icges(7,2) * t440 + Icges(7,6) * t474;
t386 = Icges(7,4) * t439 + Icges(7,2) * t438 + Icges(7,6) * t472;
t385 = Icges(7,5) * t441 + Icges(7,6) * t440 + Icges(7,3) * t474;
t384 = Icges(7,5) * t439 + Icges(7,6) * t438 + Icges(7,3) * t472;
t383 = t434 * t521 + t537 * t551 + t559;
t382 = (-t433 - t476) * t521 + t540 * t551 + t555;
t381 = (-t562 + (t433 * t537 + t434 * t540) * qJD(2)) * t531 + t567;
t380 = t422 * t510 - t449 * t486 + t546;
t379 = -t421 * t510 + t449 * t487 + t543;
t378 = t421 * t486 - t422 * t487 + t547;
t377 = t414 * t490 - t445 * t457 + t544;
t376 = -t413 * t490 + t445 * t458 + t542;
t375 = t413 * t457 - t414 * t458 + t545;
t374 = t391 * t451 - t401 * t425 + t436 * t490 - t450 * t457 + t544;
t373 = -t390 * t451 + t401 * t426 - t435 * t490 + t450 * t458 + t542;
t372 = t390 * t425 - t391 * t426 + t435 * t457 - t436 * t458 + t545;
t1 = m(3) * (t396 ^ 2 + t405 ^ 2 + t406 ^ 2) / 0.2e1 + m(7) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + m(5) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + m(4) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + t426 * ((t385 * t472 + t387 * t438 + t389 * t439) * t425 + (t472 * t384 + t438 * t386 + t439 * t388) * t426 + (t397 * t472 + t398 * t438 + t399 * t439) * t451) / 0.2e1 + t451 * ((t385 * t488 + t387 * t470 + t389 * t471) * t425 + (t384 * t488 + t386 * t470 + t388 * t471) * t426 + (t488 * t397 + t470 * t398 + t471 * t399) * t451) / 0.2e1 + t425 * ((t474 * t385 + t440 * t387 + t441 * t389) * t425 + (t384 * t474 + t386 * t440 + t388 * t441) * t426 + (t397 * t474 + t398 * t440 + t399 * t441) * t451) / 0.2e1 + t457 * ((t408 * t507 - t410 * t474 + t412 * t475) * t457 + (t407 * t507 - t409 * t474 + t411 * t475) * t458 + (t442 * t507 - t443 * t474 + t444 * t475) * t490) / 0.2e1 + t458 * ((t408 * t505 - t410 * t472 + t412 * t473) * t457 + (t407 * t505 - t409 * t472 + t411 * t473) * t458 + (t442 * t505 - t443 * t472 + t444 * t473) * t490) / 0.2e1 + t490 * ((-t408 * t574 - t410 * t488 + t412 * t489) * t457 + (-t407 * t574 - t409 * t488 + t411 * t489) * t458 + (-t442 * t574 - t443 * t488 + t444 * t489) * t490) / 0.2e1 + t486 * ((t416 * t507 + t418 * t480 + t420 * t481) * t486 + (t415 * t507 + t417 * t480 + t419 * t481) * t487 + (t446 * t507 + t447 * t480 + t448 * t481) * t510) / 0.2e1 + t487 * ((t416 * t505 + t418 * t478 + t420 * t479) * t486 + (t415 * t505 + t417 * t478 + t419 * t479) * t487 + (t446 * t505 + t447 * t478 + t448 * t479) * t510) / 0.2e1 + t510 * ((-t416 * t574 + t418 * t494 + t420 * t495) * t486 + (-t415 * t574 + t417 * t494 + t419 * t495) * t487 + (-t446 * t574 + t447 * t494 + t448 * t495) * t510) / 0.2e1 + ((((t462 * t539 + t464 * t536) * t537 - (t461 * t539 + t463 * t536) * t540) * t531 - t549 * t533 + (t430 * t503 + t432 * t504) * t537 - (t429 * t503 + t431 * t504) * t540 + (t427 * t540 - t428 * t537) * t574) * t563 + (t533 * t491 + (t492 * t539 + t493 * t536) * t531 - t452 * t574 + t453 * t503 + t454 * t504) * t521) * t521 / 0.2e1 + (m(2) * (t517 ^ 2 + t518 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 - (((-t429 * t482 - t431 * t483 - t506 * t463 + t505 * t584 + t585) * t540 + (t430 * t482 + t432 * t483 + t464 * t506 + t505 * t583) * t537) * t563 + (t453 * t482 + t454 * t483 - t491 * t573 + t493 * t506 + t505 * t582) * t521) * t558 / 0.2e1 + (((-t429 * t484 - t431 * t485 - t463 * t508 + t507 * t584) * t540 + (t430 * t484 + t432 * t485 + t508 * t464 + t507 * t583 - t585) * t537) * t563 + (t453 * t484 + t454 * t485 + t491 * t575 + t493 * t508 + t507 * t582) * t521) * t520 / 0.2e1;
T  = t1;
