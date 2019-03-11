% Calculate kinetic energy for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:20
% EndTime: 2019-03-09 00:06:22
% DurationCPUTime: 2.74s
% Computational Cost: add. (3088->329), mult. (6807->508), div. (0->0), fcn. (8383->12), ass. (0->154)
t577 = Icges(6,1) + Icges(7,1);
t576 = Icges(6,4) + Icges(7,4);
t575 = Icges(6,5) + Icges(7,5);
t574 = Icges(6,2) + Icges(7,2);
t573 = Icges(6,6) + Icges(7,6);
t572 = Icges(6,3) + Icges(7,3);
t571 = rSges(7,3) + qJ(6);
t516 = sin(pkin(11));
t518 = cos(pkin(11));
t524 = cos(qJ(2));
t519 = cos(pkin(6));
t522 = sin(qJ(2));
t547 = t519 * t522;
t496 = t516 * t524 + t518 * t547;
t521 = sin(qJ(3));
t517 = sin(pkin(6));
t548 = t518 * t517;
t557 = cos(qJ(3));
t481 = t496 * t557 - t521 * t548;
t546 = t519 * t524;
t495 = t516 * t522 - t518 * t546;
t515 = qJ(4) + qJ(5);
t511 = sin(t515);
t512 = cos(t515);
t448 = -t481 * t511 + t495 * t512;
t449 = t481 * t512 + t495 * t511;
t537 = t517 * t557;
t480 = t496 * t521 + t518 * t537;
t570 = t573 * t448 + t575 * t449 + t572 * t480;
t498 = -t516 * t547 + t518 * t524;
t550 = t517 * t521;
t483 = t498 * t557 + t516 * t550;
t497 = t516 * t546 + t518 * t522;
t450 = -t483 * t511 + t497 * t512;
t451 = t483 * t512 + t497 * t511;
t482 = t498 * t521 - t516 * t537;
t569 = t573 * t450 + t575 * t451 + t572 * t482;
t568 = t574 * t448 + t576 * t449 + t573 * t480;
t567 = t574 * t450 + t576 * t451 + t573 * t482;
t566 = t576 * t448 + t577 * t449 + t575 * t480;
t565 = t576 * t450 + t577 * t451 + t575 * t482;
t500 = t519 * t521 + t522 * t537;
t549 = t517 * t524;
t475 = -t500 * t511 - t512 * t549;
t476 = t500 * t512 - t511 * t549;
t499 = -t519 * t557 + t522 * t550;
t564 = t573 * t475 + t575 * t476 + t572 * t499;
t563 = t574 * t475 + t576 * t476 + t573 * t499;
t562 = t576 * t475 + t577 * t476 + t575 * t499;
t561 = qJD(2) ^ 2;
t523 = cos(qJ(4));
t555 = t523 * pkin(4);
t520 = sin(qJ(4));
t553 = t495 * t520;
t552 = t497 * t520;
t551 = t516 * t517;
t534 = pkin(5) * t511;
t542 = pkin(5) * t512;
t545 = rSges(7,1) * t449 + rSges(7,2) * t448 + t571 * t480 + t481 * t542 + t495 * t534;
t544 = rSges(7,1) * t451 + rSges(7,2) * t450 + t571 * t482 + t483 * t542 + t497 * t534;
t543 = rSges(7,1) * t476 + rSges(7,2) * t475 + t571 * t499 + t500 * t542 - t534 * t549;
t540 = qJD(2) * t517;
t507 = t516 * t540;
t486 = qJD(3) * t497 + t507;
t510 = qJD(2) * t519;
t538 = t520 * t549;
t446 = qJD(4) * t482 + t486;
t536 = t518 * t540;
t472 = pkin(2) * t496 + pkin(8) * t495;
t473 = pkin(2) * t498 + pkin(8) * t497;
t535 = t472 * t507 + t473 * t536 + qJD(1);
t487 = qJD(3) * t495 - t536;
t502 = -qJD(3) * t549 + t510;
t501 = (pkin(2) * t522 - pkin(8) * t524) * t517;
t533 = t473 * t510 - t501 * t507;
t447 = qJD(4) * t480 + t487;
t479 = qJD(4) * t499 + t502;
t444 = t481 * pkin(3) + t480 * pkin(9);
t445 = t483 * pkin(3) + t482 * pkin(9);
t532 = t486 * t444 - t445 * t487 + t535;
t531 = (-t472 * t519 - t501 * t548) * qJD(2);
t474 = t500 * pkin(3) + t499 * pkin(9);
t530 = t502 * t445 - t474 * t486 + t533;
t387 = pkin(4) * t553 + pkin(10) * t480 + t481 * t555;
t388 = pkin(4) * t552 + pkin(10) * t482 + t483 * t555;
t529 = t446 * t387 - t388 * t447 + t532;
t528 = -t444 * t502 + t487 * t474 + t531;
t427 = -pkin(4) * t538 + pkin(10) * t499 + t500 * t555;
t527 = t479 * t388 - t427 * t446 + t530;
t526 = -t387 * t479 + t447 * t427 + t528;
t491 = rSges(3,3) * t519 + (rSges(3,1) * t522 + rSges(3,2) * t524) * t517;
t490 = Icges(3,5) * t519 + (Icges(3,1) * t522 + Icges(3,4) * t524) * t517;
t489 = Icges(3,6) * t519 + (Icges(3,4) * t522 + Icges(3,2) * t524) * t517;
t488 = Icges(3,3) * t519 + (Icges(3,5) * t522 + Icges(3,6) * t524) * t517;
t485 = t500 * t523 - t538;
t484 = -t500 * t520 - t523 * t549;
t470 = rSges(4,1) * t500 - rSges(4,2) * t499 - rSges(4,3) * t549;
t469 = Icges(4,1) * t500 - Icges(4,4) * t499 - Icges(4,5) * t549;
t468 = Icges(4,4) * t500 - Icges(4,2) * t499 - Icges(4,6) * t549;
t467 = Icges(4,5) * t500 - Icges(4,6) * t499 - Icges(4,3) * t549;
t464 = rSges(3,1) * t498 - rSges(3,2) * t497 + rSges(3,3) * t551;
t463 = rSges(3,1) * t496 - rSges(3,2) * t495 - rSges(3,3) * t548;
t462 = Icges(3,1) * t498 - Icges(3,4) * t497 + Icges(3,5) * t551;
t461 = Icges(3,1) * t496 - Icges(3,4) * t495 - Icges(3,5) * t548;
t460 = Icges(3,4) * t498 - Icges(3,2) * t497 + Icges(3,6) * t551;
t459 = Icges(3,4) * t496 - Icges(3,2) * t495 - Icges(3,6) * t548;
t458 = Icges(3,5) * t498 - Icges(3,6) * t497 + Icges(3,3) * t551;
t457 = Icges(3,5) * t496 - Icges(3,6) * t495 - Icges(3,3) * t548;
t456 = qJD(5) * t499 + t479;
t455 = t483 * t523 + t552;
t454 = -t483 * t520 + t497 * t523;
t453 = t481 * t523 + t553;
t452 = -t481 * t520 + t495 * t523;
t441 = (-t463 * t519 - t491 * t548) * qJD(2);
t440 = (t464 * t519 - t491 * t551) * qJD(2);
t439 = rSges(5,1) * t485 + rSges(5,2) * t484 + rSges(5,3) * t499;
t438 = Icges(5,1) * t485 + Icges(5,4) * t484 + Icges(5,5) * t499;
t437 = Icges(5,4) * t485 + Icges(5,2) * t484 + Icges(5,6) * t499;
t436 = Icges(5,5) * t485 + Icges(5,6) * t484 + Icges(5,3) * t499;
t435 = rSges(4,1) * t483 - rSges(4,2) * t482 + rSges(4,3) * t497;
t434 = rSges(4,1) * t481 - rSges(4,2) * t480 + rSges(4,3) * t495;
t433 = Icges(4,1) * t483 - Icges(4,4) * t482 + Icges(4,5) * t497;
t432 = Icges(4,1) * t481 - Icges(4,4) * t480 + Icges(4,5) * t495;
t431 = Icges(4,4) * t483 - Icges(4,2) * t482 + Icges(4,6) * t497;
t430 = Icges(4,4) * t481 - Icges(4,2) * t480 + Icges(4,6) * t495;
t429 = Icges(4,5) * t483 - Icges(4,6) * t482 + Icges(4,3) * t497;
t428 = Icges(4,5) * t481 - Icges(4,6) * t480 + Icges(4,3) * t495;
t426 = rSges(6,1) * t476 + rSges(6,2) * t475 + rSges(6,3) * t499;
t418 = qJD(5) * t480 + t447;
t417 = qJD(5) * t482 + t446;
t415 = qJD(1) + (t463 * t516 + t464 * t518) * t540;
t413 = rSges(5,1) * t455 + rSges(5,2) * t454 + rSges(5,3) * t482;
t412 = rSges(5,1) * t453 + rSges(5,2) * t452 + rSges(5,3) * t480;
t411 = Icges(5,1) * t455 + Icges(5,4) * t454 + Icges(5,5) * t482;
t410 = Icges(5,1) * t453 + Icges(5,4) * t452 + Icges(5,5) * t480;
t409 = Icges(5,4) * t455 + Icges(5,2) * t454 + Icges(5,6) * t482;
t408 = Icges(5,4) * t453 + Icges(5,2) * t452 + Icges(5,6) * t480;
t407 = Icges(5,5) * t455 + Icges(5,6) * t454 + Icges(5,3) * t482;
t406 = Icges(5,5) * t453 + Icges(5,6) * t452 + Icges(5,3) * t480;
t404 = rSges(6,1) * t451 + rSges(6,2) * t450 + rSges(6,3) * t482;
t402 = rSges(6,1) * t449 + rSges(6,2) * t448 + rSges(6,3) * t480;
t382 = -t434 * t502 + t470 * t487 + t531;
t381 = t435 * t502 - t470 * t486 + t533;
t380 = t434 * t486 - t435 * t487 + t535;
t379 = -t412 * t479 + t439 * t447 + t528;
t378 = t413 * t479 - t439 * t446 + t530;
t377 = t412 * t446 - t413 * t447 + t532;
t376 = -t402 * t456 + t418 * t426 + t526;
t375 = t404 * t456 - t417 * t426 + t527;
t374 = t402 * t417 - t404 * t418 + t529;
t373 = qJD(6) * t482 + t418 * t543 - t456 * t545 + t526;
t372 = qJD(6) * t480 - t417 * t543 + t456 * t544 + t527;
t371 = qJD(6) * t499 + t417 * t545 - t418 * t544 + t529;
t1 = t487 * ((t495 * t429 - t480 * t431 + t481 * t433) * t486 + (t495 * t428 - t480 * t430 + t481 * t432) * t487 + (t495 * t467 - t480 * t468 + t481 * t469) * t502) / 0.2e1 + t502 * ((-t429 * t549 - t499 * t431 + t500 * t433) * t486 + (-t428 * t549 - t499 * t430 + t500 * t432) * t487 + (-t467 * t549 - t499 * t468 + t500 * t469) * t502) / 0.2e1 + t446 * ((t482 * t407 + t454 * t409 + t455 * t411) * t446 + (t482 * t406 + t408 * t454 + t455 * t410) * t447 + (t482 * t436 + t454 * t437 + t455 * t438) * t479) / 0.2e1 + t447 * ((t480 * t407 + t452 * t409 + t453 * t411) * t446 + (t480 * t406 + t452 * t408 + t453 * t410) * t447 + (t480 * t436 + t452 * t437 + t453 * t438) * t479) / 0.2e1 + t479 * ((t499 * t407 + t484 * t409 + t485 * t411) * t446 + (t499 * t406 + t484 * t408 + t485 * t410) * t447 + (t499 * t436 + t484 * t437 + t485 * t438) * t479) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t415 ^ 2 + t440 ^ 2 + t441 ^ 2) / 0.2e1 + m(4) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(5) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(6) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(7) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + t486 * ((t497 * t429 - t482 * t431 + t483 * t433) * t486 + (t497 * t428 - t482 * t430 + t483 * t432) * t487 + (t497 * t467 - t482 * t468 + t483 * t469) * t502) / 0.2e1 - t561 * ((-t458 * t548 - t460 * t495 + t462 * t496) * t551 - (-t457 * t548 - t459 * t495 + t461 * t496) * t548 + (-t488 * t548 - t489 * t495 + t490 * t496) * t519) * t548 / 0.2e1 + ((t450 * t563 + t451 * t562 + t482 * t564) * t456 + (t568 * t450 + t566 * t451 + t482 * t570) * t418 + (t567 * t450 + t565 * t451 + t569 * t482) * t417) * t417 / 0.2e1 + ((t448 * t563 + t449 * t562 + t480 * t564) * t456 + (t568 * t448 + t566 * t449 + t570 * t480) * t418 + (t448 * t567 + t449 * t565 + t480 * t569) * t417) * t418 / 0.2e1 + ((t563 * t475 + t562 * t476 + t564 * t499) * t456 + (t568 * t475 + t566 * t476 + t499 * t570) * t418 + (t475 * t567 + t476 * t565 + t499 * t569) * t417) * t456 / 0.2e1 + (t519 * (t519 ^ 2 * t488 + (((t460 * t524 + t462 * t522) * t516 - (t459 * t524 + t461 * t522) * t518) * t517 + (-t457 * t518 + t458 * t516 + t489 * t524 + t490 * t522) * t519) * t517) + ((t458 * t551 - t460 * t497 + t462 * t498) * t551 - (t457 * t551 - t459 * t497 + t461 * t498) * t548 + (t488 * t551 - t489 * t497 + t490 * t498) * t519) * t551) * t561 / 0.2e1;
T  = t1;
