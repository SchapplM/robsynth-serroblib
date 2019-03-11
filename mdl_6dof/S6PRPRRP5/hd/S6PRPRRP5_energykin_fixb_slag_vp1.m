% Calculate kinetic energy for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:14:01
% EndTime: 2019-03-08 20:14:03
% DurationCPUTime: 2.70s
% Computational Cost: add. (2018->280), mult. (5147->425), div. (0->0), fcn. (6193->10), ass. (0->139)
t580 = Icges(3,1) + Icges(4,2);
t579 = Icges(6,1) + Icges(7,1);
t578 = Icges(3,4) + Icges(4,6);
t577 = Icges(6,4) + Icges(7,4);
t576 = Icges(3,5) - Icges(4,4);
t575 = Icges(6,5) + Icges(7,5);
t574 = Icges(3,2) + Icges(4,3);
t573 = Icges(6,2) + Icges(7,2);
t572 = Icges(3,6) - Icges(4,5);
t571 = Icges(6,6) + Icges(7,6);
t570 = Icges(3,3) + Icges(4,1);
t569 = Icges(6,3) + Icges(7,3);
t568 = rSges(7,3) + qJ(6);
t504 = cos(pkin(6));
t502 = sin(pkin(6));
t503 = cos(pkin(10));
t535 = t503 * t502;
t501 = sin(pkin(10));
t538 = t501 * t502;
t508 = sin(qJ(2));
t510 = cos(qJ(2));
t551 = t570 * t504 + (t576 * t508 + t572 * t510) * t502;
t533 = t504 * t510;
t486 = t501 * t508 - t503 * t533;
t534 = t504 * t508;
t487 = t501 * t510 + t503 * t534;
t552 = -t572 * t486 + t576 * t487 - t570 * t535;
t488 = t501 * t533 + t503 * t508;
t489 = -t501 * t534 + t503 * t510;
t553 = -t572 * t488 + t576 * t489 + t570 * t538;
t567 = -t551 * t504 + t535 * t552 - t538 * t553;
t507 = sin(qJ(4));
t543 = cos(qJ(4));
t524 = t502 * t543;
t464 = t488 * t507 + t501 * t524;
t506 = sin(qJ(5));
t509 = cos(qJ(5));
t429 = -t464 * t506 + t489 * t509;
t539 = t489 * t506;
t430 = t464 * t509 + t539;
t537 = t502 * t507;
t463 = -t488 * t543 + t501 * t537;
t566 = t571 * t429 + t575 * t430 + t569 * t463;
t466 = t486 * t507 - t503 * t524;
t431 = -t466 * t506 + t487 * t509;
t540 = t487 * t506;
t432 = t466 * t509 + t540;
t465 = t486 * t543 + t507 * t535;
t565 = t571 * t431 + t575 * t432 - t569 * t465;
t564 = t573 * t429 + t577 * t430 + t571 * t463;
t563 = t573 * t431 + t577 * t432 - t571 * t465;
t562 = t577 * t429 + t579 * t430 + t575 * t463;
t561 = t577 * t431 + t579 * t432 - t575 * t465;
t491 = t504 * t543 - t510 * t537;
t536 = t502 * t508;
t467 = -t491 * t506 + t509 * t536;
t525 = t506 * t536;
t468 = t491 * t509 + t525;
t490 = t504 * t507 + t510 * t524;
t560 = t571 * t467 + t575 * t468 + t569 * t490;
t559 = t573 * t467 + t577 * t468 + t571 * t490;
t558 = t577 * t467 + t579 * t468 + t575 * t490;
t557 = t574 * t488 - t578 * t489 - t572 * t538;
t556 = t574 * t486 - t578 * t487 + t572 * t535;
t555 = -t578 * t488 + t580 * t489 + t576 * t538;
t554 = t578 * t486 - t580 * t487 + t576 * t535;
t550 = t572 * t504 + (t578 * t508 + t574 * t510) * t502;
t549 = t576 * t504 + (t580 * t508 + t578 * t510) * t502;
t547 = qJD(2) ^ 2;
t542 = pkin(5) * t509;
t532 = rSges(7,1) * t430 + rSges(7,2) * t429 + pkin(5) * t539 + t568 * t463 + t464 * t542;
t531 = rSges(7,1) * t432 + rSges(7,2) * t431 + pkin(5) * t540 - t568 * t465 + t466 * t542;
t530 = rSges(7,1) * t468 + rSges(7,2) * t467 + pkin(5) * t525 + t568 * t490 + t491 * t542;
t457 = pkin(2) * t489 + qJ(3) * t488;
t500 = qJD(2) * t504;
t529 = qJD(3) * t486 + t457 * t500;
t528 = qJD(2) * t502;
t497 = t501 * t528;
t469 = qJD(4) * t489 + t497;
t493 = qJD(4) * t536 + t500;
t527 = qJD(3) * t510;
t523 = t503 * t528;
t456 = pkin(2) * t487 + qJ(3) * t486;
t522 = t456 * t497 + t457 * t523 + qJD(1);
t492 = (pkin(2) * t508 - qJ(3) * t510) * t502;
t520 = (-t504 * rSges(4,1) - (-rSges(4,2) * t508 - rSges(4,3) * t510) * t502 - t492) * t502;
t519 = (-pkin(3) * t504 - pkin(8) * t536 - t492) * t502;
t470 = qJD(4) * t487 - t523;
t471 = pkin(3) * t538 + pkin(8) * t489;
t472 = -pkin(3) * t535 + pkin(8) * t487;
t516 = t471 * t523 + t472 * t497 - t502 * t527 + t522;
t515 = qJD(2) * t501 * t519 + t471 * t500 + t529;
t425 = pkin(4) * t464 + pkin(9) * t463;
t426 = pkin(4) * t466 - pkin(9) * t465;
t514 = -t470 * t425 + t469 * t426 + t516;
t485 = qJD(3) * t488;
t513 = t485 + ((-t456 - t472) * t504 + t503 * t519) * qJD(2);
t458 = pkin(4) * t491 + pkin(9) * t490;
t512 = t493 * t425 - t458 * t469 + t515;
t511 = -t426 * t493 + t470 * t458 + t513;
t479 = t504 * rSges(3,3) + (rSges(3,1) * t508 + rSges(3,2) * t510) * t502;
t461 = qJD(5) * t490 + t493;
t454 = rSges(5,1) * t491 - rSges(5,2) * t490 + rSges(5,3) * t536;
t453 = Icges(5,1) * t491 - Icges(5,4) * t490 + Icges(5,5) * t536;
t452 = Icges(5,4) * t491 - Icges(5,2) * t490 + Icges(5,6) * t536;
t451 = Icges(5,5) * t491 - Icges(5,6) * t490 + Icges(5,3) * t536;
t450 = rSges(3,1) * t489 - rSges(3,2) * t488 + rSges(3,3) * t538;
t449 = rSges(3,1) * t487 - rSges(3,2) * t486 - rSges(3,3) * t535;
t448 = -rSges(4,1) * t535 - rSges(4,2) * t487 + rSges(4,3) * t486;
t447 = rSges(4,1) * t538 - rSges(4,2) * t489 + rSges(4,3) * t488;
t428 = -qJD(5) * t465 + t470;
t427 = qJD(5) * t463 + t469;
t422 = (-t449 * t504 - t479 * t535) * qJD(2);
t421 = (t450 * t504 - t479 * t538) * qJD(2);
t420 = rSges(6,1) * t468 + rSges(6,2) * t467 + rSges(6,3) * t490;
t412 = rSges(5,1) * t466 + rSges(5,2) * t465 + rSges(5,3) * t487;
t411 = rSges(5,1) * t464 - rSges(5,2) * t463 + rSges(5,3) * t489;
t410 = Icges(5,1) * t466 + Icges(5,4) * t465 + Icges(5,5) * t487;
t409 = Icges(5,1) * t464 - Icges(5,4) * t463 + Icges(5,5) * t489;
t408 = Icges(5,4) * t466 + Icges(5,2) * t465 + Icges(5,6) * t487;
t407 = Icges(5,4) * t464 - Icges(5,2) * t463 + Icges(5,6) * t489;
t406 = Icges(5,5) * t466 + Icges(5,6) * t465 + Icges(5,3) * t487;
t405 = Icges(5,5) * t464 - Icges(5,6) * t463 + Icges(5,3) * t489;
t402 = qJD(1) + (t449 * t501 + t450 * t503) * t528;
t401 = rSges(6,1) * t432 + rSges(6,2) * t431 - rSges(6,3) * t465;
t399 = rSges(6,1) * t430 + rSges(6,2) * t429 + rSges(6,3) * t463;
t383 = t485 + ((-t448 - t456) * t504 + t503 * t520) * qJD(2);
t382 = (t447 * t504 + t501 * t520) * qJD(2) + t529;
t381 = (-t527 + (t447 * t503 + t448 * t501) * qJD(2)) * t502 + t522;
t380 = -t412 * t493 + t454 * t470 + t513;
t379 = t411 * t493 - t454 * t469 + t515;
t378 = -t470 * t411 + t469 * t412 + t516;
t377 = -t401 * t461 + t420 * t428 + t511;
t376 = t399 * t461 - t420 * t427 + t512;
t375 = -t428 * t399 + t427 * t401 + t514;
t374 = qJD(6) * t463 + t428 * t530 - t461 * t531 + t511;
t373 = -qJD(6) * t465 - t427 * t530 + t461 * t532 + t512;
t372 = qJD(6) * t490 + t427 * t531 - t428 * t532 + t514;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t402 ^ 2 + t421 ^ 2 + t422 ^ 2) / 0.2e1 + m(4) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(7) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + m(5) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + t469 * ((t489 * t405 - t463 * t407 + t464 * t409) * t469 + (t406 * t489 - t408 * t463 + t410 * t464) * t470 + (t451 * t489 - t452 * t463 + t453 * t464) * t493) / 0.2e1 + t470 * ((t405 * t487 + t407 * t465 + t409 * t466) * t469 + (t487 * t406 + t465 * t408 + t466 * t410) * t470 + (t451 * t487 + t452 * t465 + t453 * t466) * t493) / 0.2e1 + t493 * ((t405 * t536 - t407 * t490 + t409 * t491) * t469 + (t406 * t536 - t490 * t408 + t491 * t410) * t470 + (t451 * t536 - t490 * t452 + t491 * t453) * t493) / 0.2e1 + ((t429 * t559 + t430 * t558 + t463 * t560) * t461 + (t429 * t563 + t430 * t561 + t463 * t565) * t428 + (t564 * t429 + t562 * t430 + t566 * t463) * t427) * t427 / 0.2e1 + ((t431 * t559 + t432 * t558 - t465 * t560) * t461 + (t563 * t431 + t561 * t432 - t565 * t465) * t428 + (t564 * t431 + t562 * t432 - t465 * t566) * t427) * t428 / 0.2e1 + ((t559 * t467 + t558 * t468 + t560 * t490) * t461 + (t467 * t563 + t468 * t561 + t490 * t565) * t428 + (t564 * t467 + t562 * t468 + t490 * t566) * t427) * t461 / 0.2e1 - ((t486 * t557 + t487 * t555) * t538 + (-t486 * t550 + t487 * t549) * t504 + (-t486 * t556 + t487 * t554 + t567) * t535) * t547 * t535 / 0.2e1 + ((t551 * t504 ^ 2 + (((t508 * t554 + t510 * t556) * t503 + (t508 * t555 - t510 * t557) * t501) * t502 + (t501 * t553 - t503 * t552 + t508 * t549 + t510 * t550) * t504) * t502) * t504 + ((-t488 * t556 + t489 * t554) * t535 + (-t488 * t550 + t489 * t549) * t504 + (t488 * t557 + t489 * t555 - t567) * t538) * t538) * t547 / 0.2e1;
T  = t1;
