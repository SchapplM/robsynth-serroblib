% Calculate kinetic energy for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:34:10
% EndTime: 2019-03-08 21:34:13
% DurationCPUTime: 2.65s
% Computational Cost: add. (2869->319), mult. (6346->470), div. (0->0), fcn. (7814->12), ass. (0->148)
t595 = Icges(6,1) + Icges(7,1);
t594 = -Icges(6,4) + Icges(7,5);
t593 = Icges(7,4) + Icges(6,5);
t592 = Icges(4,2) + Icges(5,3);
t591 = Icges(6,2) + Icges(7,3);
t590 = Icges(7,2) + Icges(6,3);
t589 = -Icges(6,6) + Icges(7,6);
t588 = rSges(7,1) + pkin(5);
t587 = rSges(7,3) + qJ(6);
t529 = sin(pkin(10));
t532 = cos(pkin(10));
t537 = cos(qJ(2));
t533 = cos(pkin(6));
t536 = sin(qJ(2));
t559 = t533 * t536;
t514 = t529 * t537 + t532 * t559;
t530 = sin(pkin(6));
t535 = sin(qJ(3));
t561 = t530 * t535;
t567 = cos(qJ(3));
t499 = t514 * t567 - t532 * t561;
t558 = t533 * t537;
t513 = t529 * t536 - t532 * t558;
t551 = pkin(11) + qJ(5);
t527 = sin(t551);
t546 = cos(t551);
t465 = t499 * t527 - t513 * t546;
t466 = t499 * t546 + t513 * t527;
t549 = t530 * t567;
t498 = t514 * t535 + t532 * t549;
t586 = t465 * t591 + t466 * t594 + t498 * t589;
t516 = -t529 * t559 + t532 * t537;
t501 = t516 * t567 + t529 * t561;
t515 = t529 * t558 + t532 * t536;
t467 = t501 * t527 - t515 * t546;
t468 = t501 * t546 + t515 * t527;
t500 = t516 * t535 - t529 * t549;
t585 = t467 * t591 + t468 * t594 + t500 * t589;
t584 = t465 * t589 + t466 * t593 + t498 * t590;
t583 = t467 * t589 + t468 * t593 + t500 * t590;
t582 = t594 * t465 + t466 * t595 + t593 * t498;
t581 = t594 * t467 + t468 * t595 + t593 * t500;
t528 = sin(pkin(11));
t531 = cos(pkin(11));
t469 = -t499 * t528 + t513 * t531;
t565 = t513 * t528;
t470 = t499 * t531 + t565;
t580 = -Icges(4,4) * t499 + Icges(5,5) * t470 - Icges(4,6) * t513 + Icges(5,6) * t469 + t498 * t592;
t471 = -t501 * t528 + t515 * t531;
t564 = t515 * t528;
t472 = t501 * t531 + t564;
t579 = -Icges(4,4) * t501 + Icges(5,5) * t472 - Icges(4,6) * t515 + Icges(5,6) * t471 + t500 * t592;
t518 = t533 * t535 + t536 * t549;
t560 = t530 * t537;
t491 = t518 * t527 + t546 * t560;
t492 = t518 * t546 - t527 * t560;
t517 = -t533 * t567 + t536 * t561;
t578 = t491 * t591 + t492 * t594 + t517 * t589;
t577 = t491 * t589 + t492 * t593 + t517 * t590;
t576 = t594 * t491 + t492 * t595 + t593 * t517;
t496 = -t518 * t528 - t531 * t560;
t550 = t528 * t560;
t497 = t518 * t531 - t550;
t575 = -Icges(4,4) * t518 + Icges(5,5) * t497 + Icges(4,6) * t560 + Icges(5,6) * t496 + t517 * t592;
t574 = qJD(2) ^ 2;
t566 = pkin(4) * t531;
t563 = t529 * t530;
t562 = t530 * t532;
t556 = rSges(7,2) * t498 + t587 * t465 + t588 * t466;
t555 = rSges(7,2) * t500 + t587 * t467 + t588 * t468;
t554 = rSges(7,2) * t517 + t587 * t491 + t588 * t492;
t553 = qJD(2) * t530;
t523 = t529 * t553;
t502 = qJD(3) * t515 + t523;
t526 = qJD(2) * t533;
t548 = t532 * t553;
t488 = pkin(2) * t514 + pkin(8) * t513;
t489 = pkin(2) * t516 + pkin(8) * t515;
t547 = t488 * t523 + t489 * t548 + qJD(1);
t503 = qJD(3) * t513 - t548;
t520 = -qJD(3) * t560 + t526;
t461 = pkin(3) * t499 + qJ(4) * t498;
t545 = qJD(4) * t517 + t502 * t461 + t547;
t519 = (pkin(2) * t536 - pkin(8) * t537) * t530;
t544 = t489 * t526 - t519 * t523;
t462 = pkin(3) * t501 + qJ(4) * t500;
t543 = qJD(4) * t498 + t520 * t462 + t544;
t542 = (-t488 * t533 - t519 * t562) * qJD(2);
t404 = pkin(4) * t565 + pkin(9) * t498 + t499 * t566;
t405 = pkin(4) * t564 + pkin(9) * t500 + t501 * t566;
t541 = t502 * t404 + (-t405 - t462) * t503 + t545;
t490 = t518 * pkin(3) + t517 * qJ(4);
t540 = qJD(4) * t500 + t503 * t490 + t542;
t443 = -pkin(4) * t550 + pkin(9) * t517 + t518 * t566;
t539 = t520 * t405 + (-t443 - t490) * t502 + t543;
t538 = t503 * t443 + (-t404 - t461) * t520 + t540;
t509 = t533 * rSges(3,3) + (rSges(3,1) * t536 + rSges(3,2) * t537) * t530;
t508 = Icges(3,5) * t533 + (Icges(3,1) * t536 + Icges(3,4) * t537) * t530;
t507 = Icges(3,6) * t533 + (Icges(3,4) * t536 + Icges(3,2) * t537) * t530;
t506 = Icges(3,3) * t533 + (Icges(3,5) * t536 + Icges(3,6) * t537) * t530;
t495 = qJD(5) * t517 + t520;
t486 = t518 * rSges(4,1) - t517 * rSges(4,2) - rSges(4,3) * t560;
t485 = Icges(4,1) * t518 - Icges(4,4) * t517 - Icges(4,5) * t560;
t483 = Icges(4,5) * t518 - Icges(4,6) * t517 - Icges(4,3) * t560;
t480 = rSges(3,1) * t516 - rSges(3,2) * t515 + rSges(3,3) * t563;
t479 = rSges(3,1) * t514 - rSges(3,2) * t513 - rSges(3,3) * t562;
t478 = Icges(3,1) * t516 - Icges(3,4) * t515 + Icges(3,5) * t563;
t477 = Icges(3,1) * t514 - Icges(3,4) * t513 - Icges(3,5) * t562;
t476 = Icges(3,4) * t516 - Icges(3,2) * t515 + Icges(3,6) * t563;
t475 = Icges(3,4) * t514 - Icges(3,2) * t513 - Icges(3,6) * t562;
t474 = Icges(3,5) * t516 - Icges(3,6) * t515 + Icges(3,3) * t563;
t473 = Icges(3,5) * t514 - Icges(3,6) * t513 - Icges(3,3) * t562;
t464 = qJD(5) * t498 + t503;
t463 = qJD(5) * t500 + t502;
t457 = (-t479 * t533 - t509 * t562) * qJD(2);
t456 = (t480 * t533 - t509 * t563) * qJD(2);
t455 = rSges(5,1) * t497 + rSges(5,2) * t496 + rSges(5,3) * t517;
t454 = rSges(4,1) * t501 - rSges(4,2) * t500 + rSges(4,3) * t515;
t453 = rSges(4,1) * t499 - rSges(4,2) * t498 + rSges(4,3) * t513;
t452 = Icges(5,1) * t497 + Icges(5,4) * t496 + Icges(5,5) * t517;
t451 = Icges(5,4) * t497 + Icges(5,2) * t496 + Icges(5,6) * t517;
t449 = Icges(4,1) * t501 - Icges(4,4) * t500 + Icges(4,5) * t515;
t448 = Icges(4,1) * t499 - Icges(4,4) * t498 + Icges(4,5) * t513;
t445 = Icges(4,5) * t501 - Icges(4,6) * t500 + Icges(4,3) * t515;
t444 = Icges(4,5) * t499 - Icges(4,6) * t498 + Icges(4,3) * t513;
t442 = rSges(6,1) * t492 - rSges(6,2) * t491 + rSges(6,3) * t517;
t433 = qJD(1) + (t479 * t529 + t480 * t532) * t553;
t429 = rSges(5,1) * t472 + rSges(5,2) * t471 + rSges(5,3) * t500;
t428 = rSges(5,1) * t470 + rSges(5,2) * t469 + rSges(5,3) * t498;
t427 = Icges(5,1) * t472 + Icges(5,4) * t471 + Icges(5,5) * t500;
t426 = Icges(5,1) * t470 + Icges(5,4) * t469 + Icges(5,5) * t498;
t425 = Icges(5,4) * t472 + Icges(5,2) * t471 + Icges(5,6) * t500;
t424 = Icges(5,4) * t470 + Icges(5,2) * t469 + Icges(5,6) * t498;
t421 = rSges(6,1) * t468 - rSges(6,2) * t467 + rSges(6,3) * t500;
t419 = rSges(6,1) * t466 - rSges(6,2) * t465 + rSges(6,3) * t498;
t401 = -t453 * t520 + t486 * t503 + t542;
t400 = t454 * t520 - t486 * t502 + t544;
t399 = t453 * t502 - t454 * t503 + t547;
t398 = t455 * t503 + (-t428 - t461) * t520 + t540;
t397 = t429 * t520 + (-t455 - t490) * t502 + t543;
t396 = t428 * t502 + (-t429 - t462) * t503 + t545;
t395 = -t419 * t495 + t442 * t464 + t538;
t394 = t421 * t495 - t442 * t463 + t539;
t393 = t419 * t463 - t421 * t464 + t541;
t392 = qJD(6) * t467 + t464 * t554 - t495 * t556 + t538;
t391 = qJD(6) * t465 - t463 * t554 + t495 * t555 + t539;
t390 = qJD(6) * t491 + t463 * t556 - t464 * t555 + t541;
t1 = m(7) * (t390 ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(6) * (t393 ^ 2 + t394 ^ 2 + t395 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(5) * (t396 ^ 2 + t397 ^ 2 + t398 ^ 2) / 0.2e1 + m(4) * (t399 ^ 2 + t400 ^ 2 + t401 ^ 2) / 0.2e1 + m(3) * (t433 ^ 2 + t456 ^ 2 + t457 ^ 2) / 0.2e1 - t574 * ((-t474 * t562 - t476 * t513 + t478 * t514) * t563 - (-t473 * t562 - t475 * t513 + t477 * t514) * t562 + (-t506 * t562 - t507 * t513 + t508 * t514) * t533) * t562 / 0.2e1 + ((t578 * t467 + t576 * t468 + t577 * t500) * t495 + (t586 * t467 + t582 * t468 + t584 * t500) * t464 + (t585 * t467 + t581 * t468 + t583 * t500) * t463) * t463 / 0.2e1 + ((t578 * t465 + t576 * t466 + t577 * t498) * t495 + (t586 * t465 + t582 * t466 + t584 * t498) * t464 + (t585 * t465 + t581 * t466 + t583 * t498) * t463) * t464 / 0.2e1 + ((t578 * t491 + t576 * t492 + t577 * t517) * t495 + (t586 * t491 + t582 * t492 + t584 * t517) * t464 + (t585 * t491 + t581 * t492 + t583 * t517) * t463) * t495 / 0.2e1 + ((t451 * t471 + t452 * t472 + t483 * t515 + t485 * t501 + t575 * t500) * t520 + (t424 * t471 + t426 * t472 + t444 * t515 + t448 * t501 + t580 * t500) * t503 + (t425 * t471 + t427 * t472 + t445 * t515 + t449 * t501 + t579 * t500) * t502) * t502 / 0.2e1 + ((t451 * t469 + t452 * t470 + t483 * t513 + t485 * t499 + t575 * t498) * t520 + (t424 * t469 + t426 * t470 + t444 * t513 + t448 * t499 + t580 * t498) * t503 + (t425 * t469 + t427 * t470 + t445 * t513 + t449 * t499 + t579 * t498) * t502) * t503 / 0.2e1 + ((t451 * t496 + t452 * t497 - t483 * t560 + t518 * t485 + t575 * t517) * t520 + (t424 * t496 + t426 * t497 - t444 * t560 + t518 * t448 + t580 * t517) * t503 + (t425 * t496 + t427 * t497 - t445 * t560 + t518 * t449 + t579 * t517) * t502) * t520 / 0.2e1 + (t533 * (t533 ^ 2 * t506 + (((t476 * t537 + t478 * t536) * t529 - (t475 * t537 + t477 * t536) * t532) * t530 + (-t473 * t532 + t474 * t529 + t507 * t537 + t508 * t536) * t533) * t530) + ((t474 * t563 - t476 * t515 + t478 * t516) * t563 - (t473 * t563 - t475 * t515 + t477 * t516) * t562 + (t506 * t563 - t507 * t515 + t508 * t516) * t533) * t563) * t574 / 0.2e1;
T  = t1;
