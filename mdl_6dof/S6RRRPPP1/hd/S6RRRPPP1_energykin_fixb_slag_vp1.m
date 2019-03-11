% Calculate kinetic energy for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
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
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:14:26
% EndTime: 2019-03-09 15:14:29
% DurationCPUTime: 2.37s
% Computational Cost: add. (2049->271), mult. (5431->398), div. (0->0), fcn. (6505->10), ass. (0->143)
t573 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t572 = Icges(6,1) + Icges(7,1) + Icges(5,3);
t571 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t570 = -Icges(6,4) + Icges(5,5) + Icges(7,5);
t569 = Icges(7,4) + Icges(6,5) - Icges(5,6);
t568 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t567 = rSges(7,1) + pkin(5);
t566 = rSges(7,3) + qJ(6);
t501 = sin(qJ(3));
t504 = cos(qJ(3));
t506 = cos(qJ(1));
t503 = sin(qJ(1));
t505 = cos(qJ(2));
t542 = t503 * t505;
t479 = -t501 * t542 - t504 * t506;
t480 = -t501 * t506 + t504 * t542;
t502 = sin(qJ(2));
t550 = sin(pkin(6));
t551 = cos(pkin(10));
t525 = t550 * t551;
t515 = t502 * t525;
t552 = cos(pkin(6));
t527 = t552 * t551;
t549 = sin(pkin(10));
t433 = -t479 * t527 + t480 * t549 - t503 * t515;
t524 = t550 * t549;
t514 = t502 * t524;
t526 = t552 * t549;
t434 = t479 * t526 + t480 * t551 + t503 * t514;
t530 = t502 * t552;
t457 = -t479 * t550 + t503 * t530;
t565 = t571 * t433 + t573 * t434 + t570 * t457;
t541 = t505 * t506;
t481 = -t501 * t541 + t503 * t504;
t482 = t501 * t503 + t504 * t541;
t435 = -t481 * t527 + t482 * t549 - t506 * t515;
t436 = t481 * t526 + t482 * t551 + t506 * t514;
t458 = -t481 * t550 + t506 * t530;
t564 = t571 * t435 + t573 * t436 + t570 * t458;
t563 = t568 * t433 + t571 * t434 + t569 * t457;
t562 = t568 * t435 + t571 * t436 + t569 * t458;
t561 = t569 * t433 + t570 * t434 + t572 * t457;
t560 = t569 * t435 + t570 * t436 + t572 * t458;
t453 = t505 * t525 + (t501 * t527 + t504 * t549) * t502;
t544 = t502 * t504;
t546 = t501 * t502;
t454 = -t505 * t524 - t526 * t546 + t544 * t551;
t478 = -t505 * t552 + t546 * t550;
t559 = t571 * t453 + t573 * t454 + t570 * t478;
t558 = t568 * t453 + t571 * t454 + t569 * t478;
t557 = t569 * t453 + t570 * t454 + t572 * t478;
t548 = Icges(3,4) * t502;
t547 = Icges(3,4) * t505;
t545 = t502 * t503;
t543 = t502 * t506;
t540 = rSges(7,2) * t433 + t566 * t434 + t457 * t567;
t539 = rSges(7,2) * t435 + t566 * t436 + t458 * t567;
t411 = pkin(4) * t434 + qJ(5) * t433;
t439 = t480 * pkin(3) + qJ(4) * t457;
t538 = -t411 - t439;
t412 = pkin(4) * t436 + qJ(5) * t435;
t440 = t482 * pkin(3) + qJ(4) * t458;
t537 = -t412 - t440;
t536 = rSges(7,2) * t453 + t566 * t454 + t478 * t567;
t428 = pkin(4) * t454 + qJ(5) * t453;
t459 = pkin(3) * t544 + qJ(4) * t478;
t535 = -t428 - t459;
t528 = pkin(2) * t505 + pkin(9) * t502;
t483 = t528 * t503;
t484 = t528 * t506;
t532 = qJD(2) * t506;
t533 = qJD(2) * t503;
t534 = t483 * t533 + t484 * t532;
t531 = qJD(3) * t502;
t485 = t506 * t531 + t533;
t529 = qJD(4) * t478 + t485 * t439 + t534;
t523 = rSges(3,1) * t505 - rSges(3,2) * t502;
t522 = Icges(3,1) * t505 - t548;
t521 = -Icges(3,2) * t502 + t547;
t520 = Icges(3,5) * t505 - Icges(3,6) * t502;
t466 = -Icges(3,6) * t506 + t503 * t521;
t469 = -Icges(3,5) * t506 + t503 * t522;
t519 = t466 * t502 - t469 * t505;
t467 = Icges(3,6) * t503 + t506 * t521;
t470 = Icges(3,5) * t503 + t506 * t522;
t518 = -t467 * t502 + t470 * t505;
t492 = Icges(3,2) * t505 + t548;
t493 = Icges(3,1) * t502 + t547;
t517 = -t492 * t502 + t493 * t505;
t490 = qJD(1) * (pkin(1) * t506 + pkin(8) * t503);
t498 = pkin(2) * t502 - pkin(9) * t505;
t516 = qJD(1) * t484 - t498 * t533 + t490;
t513 = qJD(5) * t453 + t485 * t411 + t529;
t500 = -qJD(3) * t505 + qJD(1);
t512 = qJD(4) * t457 + t500 * t440 + t516;
t499 = pkin(1) * t503 - pkin(8) * t506;
t511 = (-t483 - t499) * qJD(1) - t498 * t532;
t510 = qJD(5) * t433 + t500 * t412 + t512;
t486 = t503 * t531 - t532;
t509 = qJD(4) * t458 + t486 * t459 + t511;
t508 = qJD(5) * t435 + t486 * t428 + t509;
t496 = rSges(2,1) * t506 - rSges(2,2) * t503;
t495 = rSges(2,1) * t503 + rSges(2,2) * t506;
t494 = rSges(3,1) * t502 + rSges(3,2) * t505;
t491 = Icges(3,5) * t502 + Icges(3,6) * t505;
t473 = rSges(3,3) * t503 + t506 * t523;
t472 = -rSges(3,3) * t506 + t503 * t523;
t471 = -rSges(4,3) * t505 + (rSges(4,1) * t504 - rSges(4,2) * t501) * t502;
t468 = -Icges(4,5) * t505 + (Icges(4,1) * t504 - Icges(4,4) * t501) * t502;
t465 = -Icges(4,6) * t505 + (Icges(4,4) * t504 - Icges(4,2) * t501) * t502;
t464 = Icges(3,3) * t503 + t506 * t520;
t463 = -Icges(3,3) * t506 + t503 * t520;
t462 = -Icges(4,3) * t505 + (Icges(4,5) * t504 - Icges(4,6) * t501) * t502;
t450 = rSges(4,1) * t482 + rSges(4,2) * t481 + rSges(4,3) * t543;
t449 = rSges(4,1) * t480 + rSges(4,2) * t479 + rSges(4,3) * t545;
t448 = Icges(4,1) * t482 + Icges(4,4) * t481 + Icges(4,5) * t543;
t447 = Icges(4,1) * t480 + Icges(4,4) * t479 + Icges(4,5) * t545;
t446 = Icges(4,4) * t482 + Icges(4,2) * t481 + Icges(4,6) * t543;
t445 = Icges(4,4) * t480 + Icges(4,2) * t479 + Icges(4,6) * t545;
t444 = Icges(4,5) * t482 + Icges(4,6) * t481 + Icges(4,3) * t543;
t443 = Icges(4,5) * t480 + Icges(4,6) * t479 + Icges(4,3) * t545;
t442 = qJD(1) * t473 - t494 * t533 + t490;
t441 = -t494 * t532 + (-t472 - t499) * qJD(1);
t438 = (t472 * t503 + t473 * t506) * qJD(2);
t426 = rSges(6,1) * t478 - rSges(6,2) * t454 + rSges(6,3) * t453;
t424 = rSges(5,1) * t454 - rSges(5,2) * t453 + rSges(5,3) * t478;
t409 = t450 * t500 - t471 * t485 + t516;
t408 = -t449 * t500 + t471 * t486 + t511;
t406 = rSges(6,1) * t458 - rSges(6,2) * t436 + rSges(6,3) * t435;
t404 = rSges(6,1) * t457 - rSges(6,2) * t434 + rSges(6,3) * t433;
t402 = rSges(5,1) * t436 - rSges(5,2) * t435 + rSges(5,3) * t458;
t401 = rSges(5,1) * t434 - rSges(5,2) * t433 + rSges(5,3) * t457;
t382 = t449 * t485 - t450 * t486 + t534;
t381 = t402 * t500 + (-t424 - t459) * t485 + t512;
t380 = t424 * t486 + (-t401 - t439) * t500 + t509;
t379 = t401 * t485 + (-t402 - t440) * t486 + t529;
t378 = t406 * t500 + (-t426 + t535) * t485 + t510;
t377 = t426 * t486 + (-t404 + t538) * t500 + t508;
t376 = t404 * t485 + (-t406 + t537) * t486 + t513;
t375 = qJD(6) * t434 + t539 * t500 + (t535 - t536) * t485 + t510;
t374 = qJD(6) * t436 + t536 * t486 + (t538 - t540) * t500 + t508;
t373 = qJD(6) * t454 + t540 * t485 + (t537 - t539) * t486 + t513;
t1 = qJD(1) * ((t505 * t492 + t502 * t493) * qJD(1) + ((t467 * t505 + t470 * t502) * t503 - (t466 * t505 + t469 * t502) * t506) * qJD(2)) / 0.2e1 + m(5) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + m(6) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(7) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(4) * (t382 ^ 2 + t408 ^ 2 + t409 ^ 2) / 0.2e1 + m(3) * (t438 ^ 2 + t441 ^ 2 + t442 ^ 2) / 0.2e1 - ((-t506 * t491 + t503 * t517) * qJD(1) + (t506 ^ 2 * t463 + (t518 * t503 + (-t464 + t519) * t506) * t503) * qJD(2)) * t532 / 0.2e1 + ((t503 * t491 + t506 * t517) * qJD(1) + (t503 ^ 2 * t464 + (t519 * t506 + (-t463 + t518) * t503) * t506) * qJD(2)) * t533 / 0.2e1 + (m(2) * (t495 ^ 2 + t496 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t558 * t435 + t559 * t436 + t557 * t458 + t462 * t543 + t465 * t481 + t468 * t482) * t500 + (t563 * t435 + t565 * t436 + t443 * t543 + t445 * t481 + t447 * t482 + t561 * t458) * t486 + (t562 * t435 + t564 * t436 + t444 * t543 + t481 * t446 + t482 * t448 + t560 * t458) * t485) * t485 / 0.2e1 + ((t558 * t433 + t559 * t434 + t557 * t457 + t462 * t545 + t465 * t479 + t468 * t480) * t500 + (t563 * t433 + t565 * t434 + t443 * t545 + t479 * t445 + t480 * t447 + t561 * t457) * t486 + (t562 * t433 + t564 * t434 + t444 * t545 + t446 * t479 + t448 * t480 + t560 * t457) * t485) * t486 / 0.2e1 + ((-t505 * t462 + (-t465 * t501 + t468 * t504) * t502 + t557 * t478 + t559 * t454 + t558 * t453) * t500 + (-t443 * t505 + (-t445 * t501 + t447 * t504) * t502 + t561 * t478 + t565 * t454 + t563 * t453) * t486 + (-t444 * t505 + (-t446 * t501 + t448 * t504) * t502 + t560 * t478 + t564 * t454 + t562 * t453) * t485) * t500 / 0.2e1;
T  = t1;
