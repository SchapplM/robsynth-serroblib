% Calculate kinetic energy for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR14_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR14_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR14_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:32
% EndTime: 2019-03-09 11:32:35
% DurationCPUTime: 2.87s
% Computational Cost: add. (1789->294), mult. (4477->433), div. (0->0), fcn. (5258->10), ass. (0->141)
t575 = Icges(3,1) + Icges(4,2);
t574 = Icges(5,1) + Icges(6,2);
t573 = Icges(6,1) + Icges(5,3);
t572 = Icges(3,4) + Icges(4,6);
t571 = Icges(5,4) + Icges(6,6);
t570 = Icges(6,4) - Icges(5,5);
t569 = Icges(3,5) - Icges(4,4);
t568 = Icges(6,5) - Icges(5,6);
t567 = Icges(3,2) + Icges(4,3);
t566 = Icges(5,2) + Icges(6,3);
t565 = Icges(3,6) - Icges(4,5);
t564 = Icges(3,3) + Icges(4,1);
t502 = sin(pkin(6));
t503 = cos(pkin(6));
t508 = cos(qJ(2));
t509 = cos(qJ(1));
t532 = t508 * t509;
t505 = sin(qJ(2));
t506 = sin(qJ(1));
t535 = t505 * t506;
t484 = -t503 * t532 + t535;
t533 = t506 * t508;
t534 = t505 * t509;
t485 = t503 * t534 + t533;
t486 = t503 * t533 + t534;
t487 = -t503 * t535 + t532;
t536 = t502 * t509;
t537 = t502 * t506;
t546 = (t565 * t484 - t569 * t485 + t564 * t536) * t509 + (-t565 * t486 + t569 * t487 + t564 * t537) * t506;
t563 = t546 * t502;
t539 = sin(qJ(4));
t525 = t502 * t539;
t540 = cos(qJ(4));
t460 = t484 * t540 + t509 * t525;
t526 = t502 * t540;
t461 = -t484 * t539 + t509 * t526;
t562 = -t566 * t460 + t571 * t461 + t568 * t485;
t458 = -t486 * t540 + t506 * t525;
t459 = t486 * t539 + t506 * t526;
t561 = t566 * t458 - t571 * t459 + t568 * t487;
t560 = t568 * t458 - t570 * t459 + t573 * t487;
t559 = -t568 * t460 + t570 * t461 + t573 * t485;
t558 = -t571 * t460 + t574 * t461 + t570 * t485;
t557 = t571 * t458 - t574 * t459 + t570 * t487;
t482 = t503 * t539 + t508 * t526;
t483 = t503 * t540 - t508 * t525;
t538 = t502 * t505;
t556 = t566 * t482 - t571 * t483 + t568 * t538;
t555 = t571 * t482 - t574 * t483 + t570 * t538;
t554 = t568 * t482 - t570 * t483 + t573 * t538;
t553 = t567 * t486 - t572 * t487 - t565 * t537;
t552 = t567 * t484 - t572 * t485 + t565 * t536;
t551 = -t572 * t486 + t575 * t487 + t569 * t537;
t550 = t572 * t484 - t575 * t485 + t569 * t536;
t549 = t565 * t503 + (t572 * t505 + t567 * t508) * t502;
t548 = t569 * t503 + (t575 * t505 + t572 * t508) * t502;
t547 = t564 * t503 + (t569 * t505 + t565 * t508) * t502;
t448 = pkin(2) * t485 + qJ(3) * t484;
t449 = pkin(2) * t487 + qJ(3) * t486;
t529 = qJD(2) * t502;
t499 = t506 * t529;
t524 = t509 * t529;
t531 = t448 * t499 + t449 * t524;
t462 = qJD(4) * t487 + t499;
t530 = qJD(1) * (pkin(1) * t506 - pkin(8) * t536);
t528 = qJD(3) * t508;
t500 = qJD(2) * t503 + qJD(1);
t490 = qJD(1) * (pkin(1) * t509 + pkin(8) * t537);
t527 = qJD(3) * t484 + t500 * t449 + t490;
t489 = qJD(4) * t538 + t500;
t523 = qJD(3) * t486 - t530;
t488 = (pkin(2) * t505 - qJ(3) * t508) * t502;
t520 = (-rSges(4,1) * t503 - (-rSges(4,2) * t505 - rSges(4,3) * t508) * t502 - t488) * t529;
t519 = (-pkin(3) * t503 - pkin(9) * t538 - t488) * t529;
t463 = qJD(4) * t485 - t524;
t465 = pkin(3) * t537 + pkin(9) * t487;
t466 = -pkin(3) * t536 + pkin(9) * t485;
t516 = t465 * t524 + t466 * t499 - t502 * t528 + t531;
t411 = -pkin(4) * t461 - qJ(5) * t460;
t515 = qJD(5) * t482 + t462 * t411 + t516;
t514 = t500 * t465 + t506 * t519 + t527;
t410 = pkin(4) * t459 + qJ(5) * t458;
t513 = -qJD(5) * t460 + t489 * t410 + t514;
t512 = (-t448 - t466) * t500 + t509 * t519 + t523;
t447 = pkin(4) * t483 + qJ(5) * t482;
t511 = qJD(5) * t458 + t463 * t447 + t512;
t507 = cos(qJ(6));
t504 = sin(qJ(6));
t494 = rSges(2,1) * t509 - rSges(2,2) * t506;
t493 = rSges(2,1) * t506 + rSges(2,2) * t509;
t473 = rSges(3,3) * t503 + (rSges(3,1) * t505 + rSges(3,2) * t508) * t502;
t464 = pkin(5) * t538 + pkin(10) * t483;
t457 = t482 * t504 + t507 * t538;
t456 = t482 * t507 - t504 * t538;
t450 = qJD(6) * t483 + t489;
t446 = rSges(3,1) * t487 - rSges(3,2) * t486 + rSges(3,3) * t537;
t445 = rSges(3,1) * t485 - rSges(3,2) * t484 - rSges(3,3) * t536;
t444 = -rSges(4,1) * t536 - rSges(4,2) * t485 + rSges(4,3) * t484;
t443 = rSges(4,1) * t537 - rSges(4,2) * t487 + rSges(4,3) * t486;
t427 = rSges(5,1) * t483 - rSges(5,2) * t482 + rSges(5,3) * t538;
t426 = rSges(6,1) * t538 - rSges(6,2) * t483 + rSges(6,3) * t482;
t419 = pkin(5) * t487 + pkin(10) * t459;
t418 = pkin(5) * t485 - pkin(10) * t461;
t417 = t458 * t504 + t487 * t507;
t416 = t458 * t507 - t487 * t504;
t415 = -t460 * t504 + t485 * t507;
t414 = -t460 * t507 - t485 * t504;
t413 = -qJD(6) * t461 + t463;
t412 = qJD(6) * t459 + t462;
t407 = rSges(6,1) * t487 - rSges(6,2) * t459 + rSges(6,3) * t458;
t406 = rSges(6,1) * t485 + rSges(6,2) * t461 - rSges(6,3) * t460;
t405 = -rSges(5,1) * t461 + rSges(5,2) * t460 + rSges(5,3) * t485;
t404 = rSges(5,1) * t459 - rSges(5,2) * t458 + rSges(5,3) * t487;
t391 = rSges(7,1) * t457 + rSges(7,2) * t456 + rSges(7,3) * t483;
t390 = Icges(7,1) * t457 + Icges(7,4) * t456 + Icges(7,5) * t483;
t389 = Icges(7,4) * t457 + Icges(7,2) * t456 + Icges(7,6) * t483;
t388 = Icges(7,5) * t457 + Icges(7,6) * t456 + Icges(7,3) * t483;
t386 = t446 * t500 - t473 * t499 + t490;
t385 = -t445 * t500 - t473 * t524 - t530;
t384 = (t445 * t506 + t446 * t509) * t529;
t383 = rSges(7,1) * t417 + rSges(7,2) * t416 + rSges(7,3) * t459;
t382 = rSges(7,1) * t415 + rSges(7,2) * t414 - rSges(7,3) * t461;
t381 = Icges(7,1) * t417 + Icges(7,4) * t416 + Icges(7,5) * t459;
t380 = Icges(7,1) * t415 + Icges(7,4) * t414 - Icges(7,5) * t461;
t379 = Icges(7,4) * t417 + Icges(7,2) * t416 + Icges(7,6) * t459;
t378 = Icges(7,4) * t415 + Icges(7,2) * t414 - Icges(7,6) * t461;
t377 = Icges(7,5) * t417 + Icges(7,6) * t416 + Icges(7,3) * t459;
t376 = Icges(7,5) * t415 + Icges(7,6) * t414 - Icges(7,3) * t461;
t375 = t443 * t500 + t506 * t520 + t527;
t374 = (-t444 - t448) * t500 + t509 * t520 + t523;
t373 = (-t528 + (t443 * t509 + t444 * t506) * qJD(2)) * t502 + t531;
t372 = t404 * t489 - t427 * t462 + t514;
t371 = -t405 * t489 + t427 * t463 + t512;
t370 = -t404 * t463 + t405 * t462 + t516;
t369 = t407 * t489 + (-t426 - t447) * t462 + t513;
t368 = t426 * t463 + (-t406 - t411) * t489 + t511;
t367 = t406 * t462 + (-t407 - t410) * t463 + t515;
t366 = t383 * t450 - t391 * t412 + t419 * t489 + (-t447 - t464) * t462 + t513;
t365 = -t382 * t450 + t391 * t413 + t463 * t464 + (-t411 - t418) * t489 + t511;
t364 = t382 * t412 - t383 * t413 + t418 * t462 + (-t410 - t419) * t463 + t515;
t1 = m(4) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(3) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + m(7) * (t364 ^ 2 + t365 ^ 2 + t366 ^ 2) / 0.2e1 + m(6) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(5) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + t450 * ((t377 * t483 + t379 * t456 + t381 * t457) * t412 + (t376 * t483 + t378 * t456 + t380 * t457) * t413 + (t483 * t388 + t456 * t389 + t457 * t390) * t450) / 0.2e1 + t412 * ((t459 * t377 + t416 * t379 + t417 * t381) * t412 + (t376 * t459 + t378 * t416 + t380 * t417) * t413 + (t388 * t459 + t416 * t389 + t417 * t390) * t450) / 0.2e1 + t413 * ((-t377 * t461 + t379 * t414 + t381 * t415) * t412 + (-t461 * t376 + t414 * t378 + t415 * t380) * t413 + (-t388 * t461 + t389 * t414 + t390 * t415) * t450) / 0.2e1 + ((t556 * t458 - t555 * t459 + t554 * t487) * t489 + (t562 * t458 - t558 * t459 + t559 * t487) * t463 + (t561 * t458 - t557 * t459 + t560 * t487) * t462) * t462 / 0.2e1 + ((-t556 * t460 + t555 * t461 + t554 * t485) * t489 + (-t562 * t460 + t558 * t461 + t559 * t485) * t463 + (-t561 * t460 + t557 * t461 + t560 * t485) * t462) * t463 / 0.2e1 + ((t556 * t482 - t555 * t483 + t554 * t538) * t489 + (t562 * t482 - t558 * t483 + t559 * t538) * t463 + (t561 * t482 - t557 * t483 + t560 * t538) * t462) * t489 / 0.2e1 + ((t546 * t503 + ((t550 * t505 + t552 * t508) * t509 + (t551 * t505 - t553 * t508) * t506) * t502) * t529 + (t547 * t503 + (t548 * t505 + t549 * t508) * t502) * t500) * t500 / 0.2e1 + (Icges(2,3) + m(2) * (t493 ^ 2 + t494 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t552 * t486 + t550 * t487) * t509 + (t553 * t486 + t551 * t487 + t563) * t506) * t529 + (-t549 * t486 + t548 * t487 + t547 * t537) * t500) * t499 / 0.2e1 - (((-t552 * t484 + t550 * t485 - t563) * t509 + (t553 * t484 + t551 * t485) * t506) * t529 + (-t549 * t484 + t548 * t485 - t547 * t536) * t500) * t524 / 0.2e1;
T  = t1;
