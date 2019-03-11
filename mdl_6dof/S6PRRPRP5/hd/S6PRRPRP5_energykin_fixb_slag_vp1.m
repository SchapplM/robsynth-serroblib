% Calculate kinetic energy for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:26
% EndTime: 2019-03-08 21:45:28
% DurationCPUTime: 2.57s
% Computational Cost: add. (2201->278), mult. (5630->412), div. (0->0), fcn. (6856->10), ass. (0->134)
t579 = Icges(4,1) + Icges(5,2);
t578 = Icges(5,1) + Icges(4,3);
t577 = Icges(6,1) + Icges(7,1);
t576 = -Icges(4,4) - Icges(5,6);
t575 = Icges(5,4) - Icges(4,5);
t574 = Icges(6,4) - Icges(7,5);
t573 = Icges(7,4) + Icges(6,5);
t572 = Icges(5,5) - Icges(4,6);
t571 = Icges(4,2) + Icges(5,3);
t570 = Icges(6,2) + Icges(7,3);
t569 = Icges(7,2) + Icges(6,3);
t568 = Icges(6,6) - Icges(7,6);
t567 = rSges(7,1) + pkin(5);
t566 = rSges(7,3) + qJ(6);
t509 = sin(pkin(10));
t511 = cos(pkin(10));
t515 = cos(qJ(2));
t512 = cos(pkin(6));
t514 = sin(qJ(2));
t534 = t512 * t514;
t496 = t509 * t515 + t511 * t534;
t510 = sin(pkin(6));
t540 = cos(qJ(3));
t527 = t510 * t540;
t538 = sin(qJ(3));
t479 = t496 * t538 + t511 * t527;
t533 = t512 * t515;
t495 = t509 * t514 - t511 * t533;
t513 = sin(qJ(5));
t539 = cos(qJ(5));
t445 = -t479 * t539 + t495 * t513;
t446 = t479 * t513 + t495 * t539;
t526 = t510 * t538;
t480 = t496 * t540 - t511 * t526;
t565 = t570 * t445 - t574 * t446 - t568 * t480;
t498 = -t509 * t534 + t511 * t515;
t481 = t498 * t538 - t509 * t527;
t497 = t509 * t533 + t511 * t514;
t447 = -t481 * t539 + t497 * t513;
t448 = t481 * t513 + t497 * t539;
t482 = t498 * t540 + t509 * t526;
t564 = t570 * t447 - t574 * t448 - t568 * t482;
t563 = -t568 * t445 + t573 * t446 + t569 * t480;
t562 = -t568 * t447 + t573 * t448 + t569 * t482;
t561 = -t574 * t445 + t577 * t446 + t573 * t480;
t560 = -t574 * t447 + t577 * t448 + t573 * t482;
t559 = t571 * t479 + t576 * t480 + t572 * t495;
t558 = t571 * t481 + t576 * t482 + t572 * t497;
t557 = t572 * t479 - t575 * t480 + t578 * t495;
t556 = t572 * t481 - t575 * t482 + t578 * t497;
t555 = t576 * t479 + t579 * t480 - t575 * t495;
t554 = t576 * t481 + t579 * t482 - t575 * t497;
t499 = -t512 * t540 + t514 * t526;
t535 = t510 * t515;
t483 = t499 * t539 + t513 * t535;
t484 = t499 * t513 - t535 * t539;
t500 = t512 * t538 + t514 * t527;
t553 = -t570 * t483 - t574 * t484 - t568 * t500;
t552 = t568 * t483 + t573 * t484 + t569 * t500;
t551 = t574 * t483 + t577 * t484 + t573 * t500;
t550 = t571 * t499 + t576 * t500 - t572 * t535;
t549 = t576 * t499 + t579 * t500 + t575 * t535;
t548 = t572 * t499 - t575 * t500 - t578 * t535;
t547 = qJD(2) ^ 2;
t537 = t509 * t510;
t536 = t510 * t511;
t532 = rSges(7,2) * t480 + t566 * t445 + t567 * t446;
t531 = rSges(7,2) * t482 + t566 * t447 + t567 * t448;
t530 = rSges(7,2) * t500 - t566 * t483 + t567 * t484;
t529 = qJD(2) * t510;
t506 = t509 * t529;
t485 = qJD(3) * t497 + t506;
t508 = qJD(2) * t512;
t525 = t511 * t529;
t471 = pkin(2) * t496 + pkin(8) * t495;
t472 = pkin(2) * t498 + pkin(8) * t497;
t524 = t471 * t506 + t472 * t525 + qJD(1);
t486 = qJD(3) * t495 - t525;
t502 = -qJD(3) * t535 + t508;
t440 = pkin(3) * t480 + qJ(4) * t479;
t523 = qJD(4) * t499 + t485 * t440 + t524;
t501 = (pkin(2) * t514 - pkin(8) * t515) * t510;
t522 = t472 * t508 - t501 * t506;
t441 = pkin(3) * t482 + qJ(4) * t481;
t521 = qJD(4) * t479 + t502 * t441 + t522;
t520 = (-t471 * t512 - t501 * t536) * qJD(2);
t450 = pkin(4) * t495 + pkin(9) * t480;
t451 = pkin(4) * t497 + pkin(9) * t482;
t519 = t485 * t450 + (-t441 - t451) * t486 + t523;
t473 = pkin(3) * t500 + qJ(4) * t499;
t518 = qJD(4) * t481 + t486 * t473 + t520;
t487 = -pkin(4) * t535 + t500 * pkin(9);
t517 = t502 * t451 + (-t473 - t487) * t485 + t521;
t516 = t486 * t487 + (-t440 - t450) * t502 + t518;
t491 = t512 * rSges(3,3) + (rSges(3,1) * t514 + rSges(3,2) * t515) * t510;
t490 = Icges(3,5) * t512 + (Icges(3,1) * t514 + Icges(3,4) * t515) * t510;
t489 = Icges(3,6) * t512 + (Icges(3,4) * t514 + Icges(3,2) * t515) * t510;
t488 = Icges(3,3) * t512 + (Icges(3,5) * t514 + Icges(3,6) * t515) * t510;
t478 = qJD(5) * t500 + t502;
t469 = t500 * rSges(4,1) - t499 * rSges(4,2) - rSges(4,3) * t535;
t468 = -rSges(5,1) * t535 - t500 * rSges(5,2) + t499 * rSges(5,3);
t459 = rSges(3,1) * t498 - rSges(3,2) * t497 + rSges(3,3) * t537;
t458 = rSges(3,1) * t496 - rSges(3,2) * t495 - rSges(3,3) * t536;
t457 = Icges(3,1) * t498 - Icges(3,4) * t497 + Icges(3,5) * t537;
t456 = Icges(3,1) * t496 - Icges(3,4) * t495 - Icges(3,5) * t536;
t455 = Icges(3,4) * t498 - Icges(3,2) * t497 + Icges(3,6) * t537;
t454 = Icges(3,4) * t496 - Icges(3,2) * t495 - Icges(3,6) * t536;
t453 = Icges(3,5) * t498 - Icges(3,6) * t497 + Icges(3,3) * t537;
t452 = Icges(3,5) * t496 - Icges(3,6) * t495 - Icges(3,3) * t536;
t444 = qJD(5) * t480 + t486;
t443 = qJD(5) * t482 + t485;
t435 = (-t458 * t512 - t491 * t536) * qJD(2);
t434 = (t459 * t512 - t491 * t537) * qJD(2);
t433 = rSges(6,1) * t484 + rSges(6,2) * t483 + rSges(6,3) * t500;
t425 = rSges(4,1) * t482 - rSges(4,2) * t481 + rSges(4,3) * t497;
t424 = rSges(4,1) * t480 - rSges(4,2) * t479 + rSges(4,3) * t495;
t423 = rSges(5,1) * t497 - rSges(5,2) * t482 + rSges(5,3) * t481;
t422 = rSges(5,1) * t495 - rSges(5,2) * t480 + rSges(5,3) * t479;
t408 = qJD(1) + (t458 * t509 + t459 * t511) * t529;
t405 = rSges(6,1) * t448 - rSges(6,2) * t447 + rSges(6,3) * t482;
t403 = rSges(6,1) * t446 - rSges(6,2) * t445 + rSges(6,3) * t480;
t389 = -t424 * t502 + t469 * t486 + t520;
t388 = t425 * t502 - t469 * t485 + t522;
t387 = t424 * t485 - t425 * t486 + t524;
t386 = t468 * t486 + (-t422 - t440) * t502 + t518;
t385 = t423 * t502 + (-t468 - t473) * t485 + t521;
t384 = t422 * t485 + (-t423 - t441) * t486 + t523;
t383 = -t403 * t478 + t433 * t444 + t516;
t382 = t405 * t478 - t433 * t443 + t517;
t381 = t403 * t443 - t405 * t444 + t519;
t380 = qJD(6) * t447 + t444 * t530 - t478 * t532 + t516;
t379 = qJD(6) * t445 - t443 * t530 + t478 * t531 + t517;
t378 = -qJD(6) * t483 + t443 * t532 - t444 * t531 + t519;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t408 ^ 2 + t434 ^ 2 + t435 ^ 2) / 0.2e1 + m(4) * (t387 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 - t547 * ((-t453 * t536 - t455 * t495 + t457 * t496) * t537 - (-t452 * t536 - t454 * t495 + t456 * t496) * t536 + (-t488 * t536 - t489 * t495 + t490 * t496) * t512) * t536 / 0.2e1 + m(5) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + m(6) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(7) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + ((t553 * t447 + t551 * t448 + t552 * t482) * t478 + (t565 * t447 + t561 * t448 + t563 * t482) * t444 + (t564 * t447 + t560 * t448 + t562 * t482) * t443) * t443 / 0.2e1 + ((t553 * t445 + t551 * t446 + t552 * t480) * t478 + (t565 * t445 + t561 * t446 + t563 * t480) * t444 + (t564 * t445 + t560 * t446 + t562 * t480) * t443) * t444 / 0.2e1 + ((-t553 * t483 + t551 * t484 + t552 * t500) * t478 + (-t565 * t483 + t561 * t484 + t563 * t500) * t444 + (-t564 * t483 + t560 * t484 + t562 * t500) * t443) * t478 / 0.2e1 + ((t550 * t481 + t549 * t482 + t548 * t497) * t502 + (t559 * t481 + t555 * t482 + t557 * t497) * t486 + (t558 * t481 + t554 * t482 + t556 * t497) * t485) * t485 / 0.2e1 + ((t550 * t479 + t549 * t480 + t548 * t495) * t502 + (t559 * t479 + t555 * t480 + t557 * t495) * t486 + (t558 * t479 + t554 * t480 + t556 * t495) * t485) * t486 / 0.2e1 + ((t550 * t499 + t549 * t500 - t548 * t535) * t502 + (t559 * t499 + t555 * t500 - t557 * t535) * t486 + (t558 * t499 + t554 * t500 - t556 * t535) * t485) * t502 / 0.2e1 + (t512 * (t512 ^ 2 * t488 + (((t455 * t515 + t457 * t514) * t509 - (t454 * t515 + t456 * t514) * t511) * t510 + (-t452 * t511 + t453 * t509 + t489 * t515 + t490 * t514) * t512) * t510) + ((t453 * t537 - t455 * t497 + t457 * t498) * t537 - (t452 * t537 - t454 * t497 + t456 * t498) * t536 + (t488 * t537 - t489 * t497 + t490 * t498) * t512) * t537) * t547 / 0.2e1;
T  = t1;
