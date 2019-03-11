% Calculate kinetic energy for
% S6PRRRRP4
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:06
% EndTime: 2019-03-09 00:13:08
% DurationCPUTime: 2.69s
% Computational Cost: add. (3037->324), mult. (6724->500), div. (0->0), fcn. (8297->12), ass. (0->153)
t587 = Icges(6,1) + Icges(7,1);
t586 = -Icges(6,4) + Icges(7,5);
t585 = Icges(7,4) + Icges(6,5);
t584 = Icges(6,2) + Icges(7,3);
t583 = Icges(7,2) + Icges(6,3);
t582 = -Icges(6,6) + Icges(7,6);
t581 = rSges(7,1) + pkin(5);
t580 = rSges(7,3) + qJ(6);
t527 = sin(pkin(11));
t529 = cos(pkin(11));
t535 = cos(qJ(2));
t530 = cos(pkin(6));
t533 = sin(qJ(2));
t557 = t530 * t533;
t513 = t527 * t535 + t529 * t557;
t532 = sin(qJ(3));
t528 = sin(pkin(6));
t558 = t529 * t528;
t566 = cos(qJ(3));
t496 = t513 * t566 - t532 * t558;
t556 = t530 * t535;
t512 = t527 * t533 - t529 * t556;
t555 = qJ(4) + qJ(5);
t526 = sin(t555);
t545 = cos(t555);
t463 = t496 * t526 - t512 * t545;
t464 = t496 * t545 + t512 * t526;
t548 = t528 * t566;
t495 = t513 * t532 + t529 * t548;
t579 = t584 * t463 + t586 * t464 + t582 * t495;
t515 = -t527 * t557 + t529 * t535;
t560 = t528 * t532;
t498 = t515 * t566 + t527 * t560;
t514 = t527 * t556 + t529 * t533;
t465 = t498 * t526 - t514 * t545;
t466 = t498 * t545 + t514 * t526;
t497 = t515 * t532 - t527 * t548;
t578 = t584 * t465 + t586 * t466 + t582 * t497;
t577 = t582 * t463 + t585 * t464 + t583 * t495;
t576 = t582 * t465 + t585 * t466 + t583 * t497;
t575 = t586 * t463 + t587 * t464 + t585 * t495;
t574 = t586 * t465 + t587 * t466 + t585 * t497;
t517 = t530 * t532 + t533 * t548;
t559 = t528 * t535;
t490 = t517 * t526 + t545 * t559;
t491 = t517 * t545 - t526 * t559;
t516 = -t530 * t566 + t533 * t560;
t573 = t584 * t490 + t586 * t491 + t582 * t516;
t572 = t582 * t490 + t585 * t491 + t583 * t516;
t571 = t586 * t490 + t587 * t491 + t585 * t516;
t570 = qJD(2) ^ 2;
t534 = cos(qJ(4));
t565 = pkin(4) * t534;
t531 = sin(qJ(4));
t563 = t512 * t531;
t562 = t514 * t531;
t561 = t527 * t528;
t554 = rSges(7,2) * t495 + t580 * t463 + t581 * t464;
t553 = rSges(7,2) * t497 + t580 * t465 + t581 * t466;
t552 = rSges(7,2) * t516 + t580 * t490 + t581 * t491;
t551 = qJD(2) * t528;
t522 = t527 * t551;
t501 = qJD(3) * t514 + t522;
t525 = qJD(2) * t530;
t549 = t531 * t559;
t461 = qJD(4) * t497 + t501;
t547 = t529 * t551;
t487 = pkin(2) * t513 + pkin(8) * t512;
t488 = pkin(2) * t515 + pkin(8) * t514;
t546 = t487 * t522 + t488 * t547 + qJD(1);
t502 = qJD(3) * t512 - t547;
t519 = -qJD(3) * t559 + t525;
t518 = (pkin(2) * t533 - pkin(8) * t535) * t528;
t544 = t488 * t525 - t518 * t522;
t462 = qJD(4) * t495 + t502;
t494 = qJD(4) * t516 + t519;
t459 = t496 * pkin(3) + t495 * pkin(9);
t460 = t498 * pkin(3) + t497 * pkin(9);
t543 = t501 * t459 - t460 * t502 + t546;
t542 = (-t487 * t530 - t518 * t558) * qJD(2);
t489 = t517 * pkin(3) + t516 * pkin(9);
t541 = t519 * t460 - t489 * t501 + t544;
t400 = pkin(4) * t563 + pkin(10) * t495 + t496 * t565;
t401 = pkin(4) * t562 + pkin(10) * t497 + t498 * t565;
t540 = t461 * t400 - t401 * t462 + t543;
t539 = -t459 * t519 + t502 * t489 + t542;
t441 = -pkin(4) * t549 + pkin(10) * t516 + t517 * t565;
t538 = t494 * t401 - t441 * t461 + t541;
t537 = -t400 * t494 + t462 * t441 + t539;
t508 = rSges(3,3) * t530 + (rSges(3,1) * t533 + rSges(3,2) * t535) * t528;
t507 = Icges(3,5) * t530 + (Icges(3,1) * t533 + Icges(3,4) * t535) * t528;
t506 = Icges(3,6) * t530 + (Icges(3,4) * t533 + Icges(3,2) * t535) * t528;
t505 = Icges(3,3) * t530 + (Icges(3,5) * t533 + Icges(3,6) * t535) * t528;
t500 = t517 * t534 - t549;
t499 = -t517 * t531 - t534 * t559;
t485 = rSges(4,1) * t517 - rSges(4,2) * t516 - rSges(4,3) * t559;
t484 = Icges(4,1) * t517 - Icges(4,4) * t516 - Icges(4,5) * t559;
t483 = Icges(4,4) * t517 - Icges(4,2) * t516 - Icges(4,6) * t559;
t482 = Icges(4,5) * t517 - Icges(4,6) * t516 - Icges(4,3) * t559;
t479 = rSges(3,1) * t515 - rSges(3,2) * t514 + rSges(3,3) * t561;
t478 = rSges(3,1) * t513 - rSges(3,2) * t512 - rSges(3,3) * t558;
t477 = Icges(3,1) * t515 - Icges(3,4) * t514 + Icges(3,5) * t561;
t476 = Icges(3,1) * t513 - Icges(3,4) * t512 - Icges(3,5) * t558;
t475 = Icges(3,4) * t515 - Icges(3,2) * t514 + Icges(3,6) * t561;
t474 = Icges(3,4) * t513 - Icges(3,2) * t512 - Icges(3,6) * t558;
t473 = Icges(3,5) * t515 - Icges(3,6) * t514 + Icges(3,3) * t561;
t472 = Icges(3,5) * t513 - Icges(3,6) * t512 - Icges(3,3) * t558;
t471 = qJD(5) * t516 + t494;
t470 = t498 * t534 + t562;
t469 = -t498 * t531 + t514 * t534;
t468 = t496 * t534 + t563;
t467 = -t496 * t531 + t512 * t534;
t455 = (-t478 * t530 - t508 * t558) * qJD(2);
t454 = (t479 * t530 - t508 * t561) * qJD(2);
t453 = rSges(5,1) * t500 + rSges(5,2) * t499 + rSges(5,3) * t516;
t452 = Icges(5,1) * t500 + Icges(5,4) * t499 + Icges(5,5) * t516;
t451 = Icges(5,4) * t500 + Icges(5,2) * t499 + Icges(5,6) * t516;
t450 = Icges(5,5) * t500 + Icges(5,6) * t499 + Icges(5,3) * t516;
t449 = rSges(4,1) * t498 - rSges(4,2) * t497 + rSges(4,3) * t514;
t448 = rSges(4,1) * t496 - rSges(4,2) * t495 + rSges(4,3) * t512;
t447 = Icges(4,1) * t498 - Icges(4,4) * t497 + Icges(4,5) * t514;
t446 = Icges(4,1) * t496 - Icges(4,4) * t495 + Icges(4,5) * t512;
t445 = Icges(4,4) * t498 - Icges(4,2) * t497 + Icges(4,6) * t514;
t444 = Icges(4,4) * t496 - Icges(4,2) * t495 + Icges(4,6) * t512;
t443 = Icges(4,5) * t498 - Icges(4,6) * t497 + Icges(4,3) * t514;
t442 = Icges(4,5) * t496 - Icges(4,6) * t495 + Icges(4,3) * t512;
t440 = rSges(6,1) * t491 - rSges(6,2) * t490 + rSges(6,3) * t516;
t432 = qJD(5) * t495 + t462;
t431 = qJD(5) * t497 + t461;
t429 = qJD(1) + (t478 * t527 + t479 * t529) * t551;
t426 = rSges(5,1) * t470 + rSges(5,2) * t469 + rSges(5,3) * t497;
t425 = rSges(5,1) * t468 + rSges(5,2) * t467 + rSges(5,3) * t495;
t424 = Icges(5,1) * t470 + Icges(5,4) * t469 + Icges(5,5) * t497;
t423 = Icges(5,1) * t468 + Icges(5,4) * t467 + Icges(5,5) * t495;
t422 = Icges(5,4) * t470 + Icges(5,2) * t469 + Icges(5,6) * t497;
t421 = Icges(5,4) * t468 + Icges(5,2) * t467 + Icges(5,6) * t495;
t420 = Icges(5,5) * t470 + Icges(5,6) * t469 + Icges(5,3) * t497;
t419 = Icges(5,5) * t468 + Icges(5,6) * t467 + Icges(5,3) * t495;
t417 = rSges(6,1) * t466 - rSges(6,2) * t465 + rSges(6,3) * t497;
t415 = rSges(6,1) * t464 - rSges(6,2) * t463 + rSges(6,3) * t495;
t397 = -t448 * t519 + t485 * t502 + t542;
t396 = t449 * t519 - t485 * t501 + t544;
t395 = t448 * t501 - t449 * t502 + t546;
t394 = -t425 * t494 + t453 * t462 + t539;
t393 = t426 * t494 - t453 * t461 + t541;
t392 = t425 * t461 - t426 * t462 + t543;
t391 = -t415 * t471 + t432 * t440 + t537;
t390 = t417 * t471 - t431 * t440 + t538;
t389 = t415 * t431 - t417 * t432 + t540;
t388 = qJD(6) * t465 + t432 * t552 - t471 * t554 + t537;
t387 = qJD(6) * t463 - t431 * t552 + t471 * t553 + t538;
t386 = qJD(6) * t490 + t431 * t554 - t432 * t553 + t540;
t1 = -t570 * ((-t473 * t558 - t475 * t512 + t477 * t513) * t561 - (-t472 * t558 - t474 * t512 + t476 * t513) * t558 + (-t505 * t558 - t506 * t512 + t507 * t513) * t530) * t558 / 0.2e1 + t462 * ((t420 * t495 + t422 * t467 + t424 * t468) * t461 + (t495 * t419 + t467 * t421 + t468 * t423) * t462 + (t450 * t495 + t451 * t467 + t452 * t468) * t494) / 0.2e1 + t494 * ((t420 * t516 + t422 * t499 + t424 * t500) * t461 + (t419 * t516 + t421 * t499 + t423 * t500) * t462 + (t516 * t450 + t499 * t451 + t500 * t452) * t494) / 0.2e1 + t461 * ((t497 * t420 + t469 * t422 + t470 * t424) * t461 + (t419 * t497 + t421 * t469 + t423 * t470) * t462 + (t450 * t497 + t451 * t469 + t452 * t470) * t494) / 0.2e1 + t501 * ((t514 * t443 - t497 * t445 + t498 * t447) * t501 + (t442 * t514 - t444 * t497 + t446 * t498) * t502 + (t482 * t514 - t483 * t497 + t484 * t498) * t519) / 0.2e1 + t502 * ((t443 * t512 - t445 * t495 + t447 * t496) * t501 + (t512 * t442 - t495 * t444 + t496 * t446) * t502 + (t482 * t512 - t483 * t495 + t484 * t496) * t519) / 0.2e1 + t519 * ((-t443 * t559 - t445 * t516 + t447 * t517) * t501 + (-t442 * t559 - t444 * t516 + t446 * t517) * t502 + (-t482 * t559 - t516 * t483 + t517 * t484) * t519) / 0.2e1 + m(7) * (t386 ^ 2 + t387 ^ 2 + t388 ^ 2) / 0.2e1 + m(6) * (t389 ^ 2 + t390 ^ 2 + t391 ^ 2) / 0.2e1 + m(5) * (t392 ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(4) * (t395 ^ 2 + t396 ^ 2 + t397 ^ 2) / 0.2e1 + m(3) * (t429 ^ 2 + t454 ^ 2 + t455 ^ 2) / 0.2e1 + ((t465 * t573 + t466 * t571 + t497 * t572) * t471 + (t465 * t579 + t575 * t466 + t577 * t497) * t432 + (t578 * t465 + t574 * t466 + t576 * t497) * t431) * t431 / 0.2e1 + ((t463 * t573 + t464 * t571 + t495 * t572) * t471 + (t579 * t463 + t575 * t464 + t577 * t495) * t432 + (t463 * t578 + t464 * t574 + t495 * t576) * t431) * t432 / 0.2e1 + ((t573 * t490 + t571 * t491 + t572 * t516) * t471 + (t490 * t579 + t575 * t491 + t577 * t516) * t432 + (t490 * t578 + t491 * t574 + t516 * t576) * t431) * t471 / 0.2e1 + (((t473 * t561 - t475 * t514 + t477 * t515) * t561 - (t472 * t561 - t474 * t514 + t476 * t515) * t558 + (t505 * t561 - t506 * t514 + t507 * t515) * t530) * t561 + t530 * (t530 ^ 2 * t505 + (((t475 * t535 + t477 * t533) * t527 - (t474 * t535 + t476 * t533) * t529) * t528 + (-t472 * t529 + t473 * t527 + t506 * t535 + t507 * t533) * t530) * t528)) * t570 / 0.2e1;
T  = t1;
