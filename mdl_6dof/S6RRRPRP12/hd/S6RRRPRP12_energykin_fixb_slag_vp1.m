% Calculate kinetic energy for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP12_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:36
% EndTime: 2019-03-09 17:51:39
% DurationCPUTime: 2.88s
% Computational Cost: add. (2270->286), mult. (5682->421), div. (0->0), fcn. (6890->10), ass. (0->138)
t586 = Icges(4,1) + Icges(5,2);
t585 = Icges(5,1) + Icges(4,3);
t584 = Icges(6,1) + Icges(7,1);
t583 = -Icges(4,4) - Icges(5,6);
t582 = Icges(5,4) - Icges(4,5);
t581 = Icges(6,4) - Icges(7,5);
t580 = Icges(7,4) + Icges(6,5);
t579 = Icges(5,5) - Icges(4,6);
t578 = Icges(4,2) + Icges(5,3);
t577 = Icges(6,2) + Icges(7,3);
t576 = Icges(7,2) + Icges(6,3);
t575 = Icges(6,6) - Icges(7,6);
t574 = rSges(7,1) + pkin(5);
t573 = rSges(7,3) + qJ(6);
t515 = sin(qJ(2));
t516 = sin(qJ(1));
t517 = cos(qJ(2));
t518 = cos(qJ(1));
t543 = cos(pkin(6));
t529 = t518 * t543;
t496 = t515 * t516 - t517 * t529;
t497 = t515 * t529 + t516 * t517;
t513 = sin(pkin(6));
t540 = t513 * t518;
t459 = Icges(3,5) * t497 - Icges(3,6) * t496 - Icges(3,3) * t540;
t530 = t516 * t543;
t498 = t518 * t515 + t517 * t530;
t499 = -t515 * t530 + t518 * t517;
t542 = t513 * t516;
t460 = Icges(3,5) * t499 - Icges(3,6) * t498 + Icges(3,3) * t542;
t572 = t513 * (t459 * t518 - t460 * t516);
t546 = cos(qJ(3));
t533 = t513 * t546;
t544 = sin(qJ(3));
t480 = t497 * t544 + t518 * t533;
t514 = sin(qJ(5));
t545 = cos(qJ(5));
t445 = -t480 * t545 + t496 * t514;
t446 = t480 * t514 + t496 * t545;
t532 = t513 * t544;
t481 = t497 * t546 - t518 * t532;
t571 = t445 * t577 - t446 * t581 - t481 * t575;
t482 = t499 * t544 - t516 * t533;
t447 = -t482 * t545 + t498 * t514;
t448 = t482 * t514 + t498 * t545;
t483 = t499 * t546 + t516 * t532;
t570 = t447 * t577 - t448 * t581 - t483 * t575;
t569 = -t445 * t575 + t446 * t580 + t481 * t576;
t568 = -t447 * t575 + t448 * t580 + t483 * t576;
t567 = -t445 * t581 + t446 * t584 + t481 * t580;
t566 = -t447 * t581 + t448 * t584 + t483 * t580;
t494 = t515 * t532 - t543 * t546;
t541 = t513 * t517;
t478 = t494 * t545 + t514 * t541;
t479 = t494 * t514 - t541 * t545;
t495 = t515 * t533 + t543 * t544;
t565 = -t478 * t577 - t479 * t581 - t495 * t575;
t564 = t478 * t575 + t479 * t580 + t495 * t576;
t563 = t478 * t581 + t479 * t584 + t495 * t580;
t562 = t480 * t578 + t481 * t583 + t496 * t579;
t561 = t482 * t578 + t483 * t583 + t498 * t579;
t560 = t480 * t579 - t481 * t582 + t496 * t585;
t559 = t482 * t579 - t483 * t582 + t498 * t585;
t558 = t583 * t480 + t481 * t586 - t582 * t496;
t557 = t583 * t482 + t483 * t586 - t582 * t498;
t556 = t494 * t578 + t495 * t583 - t541 * t579;
t555 = t583 * t494 + t495 * t586 + t582 * t541;
t554 = t494 * t579 - t495 * t582 - t541 * t585;
t539 = rSges(7,2) * t481 + t573 * t445 + t574 * t446;
t538 = rSges(7,2) * t483 + t573 * t447 + t574 * t448;
t537 = rSges(7,2) * t495 - t573 * t478 + t574 * t479;
t471 = pkin(2) * t497 + pkin(9) * t496;
t472 = pkin(2) * t499 + pkin(9) * t498;
t534 = qJD(2) * t513;
t509 = t516 * t534;
t531 = t518 * t534;
t536 = t471 * t509 + t472 * t531;
t484 = qJD(3) * t498 + t509;
t535 = qJD(1) * (pkin(1) * t516 - pkin(8) * t540);
t510 = qJD(2) * t543 + qJD(1);
t440 = pkin(3) * t481 + qJ(4) * t480;
t528 = qJD(4) * t494 + t484 * t440 + t536;
t485 = qJD(3) * t496 - t531;
t501 = -qJD(3) * t541 + t510;
t500 = (pkin(2) * t515 - pkin(9) * t517) * t513;
t502 = qJD(1) * (pkin(1) * t518 + pkin(8) * t542);
t526 = t510 * t472 - t500 * t509 + t502;
t441 = pkin(3) * t483 + qJ(4) * t482;
t525 = qJD(4) * t480 + t501 * t441 + t526;
t449 = pkin(4) * t496 + pkin(10) * t481;
t450 = pkin(4) * t498 + pkin(10) * t483;
t524 = t484 * t449 + (-t441 - t450) * t485 + t528;
t523 = -t471 * t510 - t500 * t531 - t535;
t470 = pkin(3) * t495 + qJ(4) * t494;
t522 = qJD(4) * t482 + t485 * t470 + t523;
t486 = -pkin(4) * t541 + pkin(10) * t495;
t521 = t501 * t450 + (-t470 - t486) * t484 + t525;
t520 = t485 * t486 + (-t440 - t449) * t501 + t522;
t505 = rSges(2,1) * t518 - rSges(2,2) * t516;
t504 = rSges(2,1) * t516 + rSges(2,2) * t518;
t490 = t543 * rSges(3,3) + (rSges(3,1) * t515 + rSges(3,2) * t517) * t513;
t489 = Icges(3,5) * t543 + (Icges(3,1) * t515 + Icges(3,4) * t517) * t513;
t488 = Icges(3,6) * t543 + (Icges(3,4) * t515 + Icges(3,2) * t517) * t513;
t487 = Icges(3,3) * t543 + (Icges(3,5) * t515 + Icges(3,6) * t517) * t513;
t473 = qJD(5) * t495 + t501;
t467 = rSges(3,1) * t499 - rSges(3,2) * t498 + rSges(3,3) * t542;
t466 = rSges(3,1) * t497 - rSges(3,2) * t496 - rSges(3,3) * t540;
t464 = Icges(3,1) * t499 - Icges(3,4) * t498 + Icges(3,5) * t542;
t463 = Icges(3,1) * t497 - Icges(3,4) * t496 - Icges(3,5) * t540;
t462 = Icges(3,4) * t499 - Icges(3,2) * t498 + Icges(3,6) * t542;
t461 = Icges(3,4) * t497 - Icges(3,2) * t496 - Icges(3,6) * t540;
t458 = rSges(4,1) * t495 - rSges(4,2) * t494 - rSges(4,3) * t541;
t457 = -rSges(5,1) * t541 - rSges(5,2) * t495 + rSges(5,3) * t494;
t443 = qJD(5) * t481 + t485;
t442 = qJD(5) * t483 + t484;
t434 = rSges(4,1) * t483 - rSges(4,2) * t482 + rSges(4,3) * t498;
t433 = rSges(4,1) * t481 - rSges(4,2) * t480 + rSges(4,3) * t496;
t432 = rSges(5,1) * t498 - rSges(5,2) * t483 + rSges(5,3) * t482;
t431 = rSges(5,1) * t496 - rSges(5,2) * t481 + rSges(5,3) * t480;
t418 = rSges(6,1) * t479 + rSges(6,2) * t478 + rSges(6,3) * t495;
t409 = t467 * t510 - t490 * t509 + t502;
t408 = -t466 * t510 - t490 * t531 - t535;
t407 = (t466 * t516 + t467 * t518) * t534;
t404 = rSges(6,1) * t448 - rSges(6,2) * t447 + rSges(6,3) * t483;
t402 = rSges(6,1) * t446 - rSges(6,2) * t445 + rSges(6,3) * t481;
t388 = t434 * t501 - t458 * t484 + t526;
t387 = -t433 * t501 + t458 * t485 + t523;
t386 = t433 * t484 - t434 * t485 + t536;
t385 = t432 * t501 + (-t457 - t470) * t484 + t525;
t384 = t457 * t485 + (-t431 - t440) * t501 + t522;
t383 = t431 * t484 + (-t432 - t441) * t485 + t528;
t382 = t404 * t473 - t418 * t442 + t521;
t381 = -t402 * t473 + t418 * t443 + t520;
t380 = t402 * t442 - t404 * t443 + t524;
t379 = qJD(6) * t445 - t442 * t537 + t473 * t538 + t521;
t378 = qJD(6) * t447 + t443 * t537 - t473 * t539 + t520;
t377 = -qJD(6) * t478 + t442 * t539 - t443 * t538 + t524;
t1 = t510 * ((t543 * t460 + (t462 * t517 + t464 * t515) * t513) * t509 - (t543 * t459 + (t461 * t517 + t463 * t515) * t513) * t531 + (t543 * t487 + (t488 * t517 + t489 * t515) * t513) * t510) / 0.2e1 + m(7) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(6) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(5) * (t383 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(4) * (t386 ^ 2 + t387 ^ 2 + t388 ^ 2) / 0.2e1 + m(3) * (t407 ^ 2 + t408 ^ 2 + t409 ^ 2) / 0.2e1 - ((-t487 * t540 - t488 * t496 + t489 * t497) * t510 + ((-t462 * t496 + t464 * t497) * t516 + (t496 * t461 - t497 * t463 + t572) * t518) * t534) * t531 / 0.2e1 + ((t487 * t542 - t488 * t498 + t489 * t499) * t510 + (-(-t461 * t498 + t463 * t499) * t518 + (-t498 * t462 + t499 * t464 - t572) * t516) * t534) * t509 / 0.2e1 + ((t565 * t447 + t563 * t448 + t564 * t483) * t473 + (t571 * t447 + t567 * t448 + t569 * t483) * t443 + (t570 * t447 + t566 * t448 + t568 * t483) * t442) * t442 / 0.2e1 + ((t565 * t445 + t563 * t446 + t564 * t481) * t473 + (t571 * t445 + t567 * t446 + t569 * t481) * t443 + (t570 * t445 + t566 * t446 + t568 * t481) * t442) * t443 / 0.2e1 + ((-t565 * t478 + t563 * t479 + t564 * t495) * t473 + (-t571 * t478 + t567 * t479 + t569 * t495) * t443 + (-t570 * t478 + t566 * t479 + t568 * t495) * t442) * t473 / 0.2e1 + ((t556 * t482 + t555 * t483 + t554 * t498) * t501 + (t562 * t482 + t558 * t483 + t560 * t498) * t485 + (t561 * t482 + t557 * t483 + t559 * t498) * t484) * t484 / 0.2e1 + ((t556 * t480 + t555 * t481 + t554 * t496) * t501 + (t562 * t480 + t558 * t481 + t560 * t496) * t485 + (t561 * t480 + t557 * t481 + t559 * t496) * t484) * t485 / 0.2e1 + ((t556 * t494 + t555 * t495 - t554 * t541) * t501 + (t562 * t494 + t558 * t495 - t560 * t541) * t485 + (t561 * t494 + t557 * t495 - t559 * t541) * t484) * t501 / 0.2e1 + (Icges(2,3) + m(2) * (t504 ^ 2 + t505 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
