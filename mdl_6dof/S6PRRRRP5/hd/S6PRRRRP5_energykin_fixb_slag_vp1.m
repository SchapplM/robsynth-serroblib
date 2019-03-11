% Calculate kinetic energy for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:11
% EndTime: 2019-03-09 00:20:13
% DurationCPUTime: 3.06s
% Computational Cost: add. (5339->334), mult. (14610->523), div. (0->0), fcn. (18857->14), ass. (0->159)
t633 = Icges(6,1) + Icges(7,1);
t632 = Icges(6,4) + Icges(7,4);
t631 = Icges(6,5) + Icges(7,5);
t630 = Icges(6,2) + Icges(7,2);
t629 = Icges(6,6) + Icges(7,6);
t628 = Icges(6,3) + Icges(7,3);
t627 = rSges(7,3) + qJ(6);
t566 = sin(pkin(12));
t568 = cos(pkin(12));
t574 = sin(qJ(2));
t569 = cos(pkin(6));
t576 = cos(qJ(2));
t600 = t569 * t576;
t556 = -t566 * t574 + t568 * t600;
t601 = t569 * t574;
t557 = t566 * t576 + t568 * t601;
t573 = sin(qJ(3));
t567 = sin(pkin(6));
t608 = sin(pkin(7));
t591 = t567 * t608;
t609 = cos(pkin(7));
t613 = cos(qJ(3));
t521 = t557 * t613 + (t556 * t609 - t568 * t591) * t573;
t572 = sin(qJ(4));
t592 = t567 * t609;
t584 = -t556 * t608 - t568 * t592;
t612 = cos(qJ(4));
t505 = t521 * t612 + t572 * t584;
t589 = t613 * t608;
t588 = t567 * t589;
t590 = t609 * t613;
t520 = -t556 * t590 + t557 * t573 + t568 * t588;
t571 = sin(qJ(5));
t575 = cos(qJ(5));
t477 = -t505 * t571 + t520 * t575;
t607 = t520 * t571;
t478 = t505 * t575 + t607;
t504 = t521 * t572 - t584 * t612;
t626 = t629 * t477 + t631 * t478 + t628 * t504;
t558 = -t566 * t600 - t568 * t574;
t559 = -t566 * t601 + t568 * t576;
t523 = t559 * t613 + (t558 * t609 + t566 * t591) * t573;
t583 = -t558 * t608 + t566 * t592;
t507 = t523 * t612 + t572 * t583;
t522 = -t558 * t590 + t559 * t573 - t566 * t588;
t479 = -t507 * t571 + t522 * t575;
t606 = t522 * t571;
t480 = t507 * t575 + t606;
t506 = t523 * t572 - t583 * t612;
t625 = t629 * t479 + t631 * t480 + t628 * t506;
t624 = t630 * t477 + t632 * t478 + t629 * t504;
t623 = t630 * t479 + t632 * t480 + t629 * t506;
t622 = t632 * t477 + t633 * t478 + t631 * t504;
t621 = t632 * t479 + t633 * t480 + t631 * t506;
t545 = t569 * t608 * t573 + (t573 * t576 * t609 + t574 * t613) * t567;
t582 = t569 * t609 - t576 * t591;
t526 = t545 * t612 + t572 * t582;
t603 = t567 * t574;
t544 = -t567 * t576 * t590 - t569 * t589 + t573 * t603;
t502 = -t526 * t571 + t544 * t575;
t605 = t544 * t571;
t503 = t526 * t575 + t605;
t525 = t545 * t572 - t582 * t612;
t620 = t629 * t502 + t631 * t503 + t628 * t525;
t619 = t630 * t502 + t632 * t503 + t629 * t525;
t618 = t632 * t502 + t633 * t503 + t631 * t525;
t617 = qJD(2) ^ 2;
t611 = pkin(5) * t575;
t604 = t566 * t567;
t602 = t568 * t567;
t599 = rSges(7,1) * t478 + rSges(7,2) * t477 + pkin(5) * t607 + t627 * t504 + t505 * t611;
t598 = rSges(7,1) * t480 + rSges(7,2) * t479 + pkin(5) * t606 + t627 * t506 + t507 * t611;
t597 = rSges(7,1) * t503 + rSges(7,2) * t502 + pkin(5) * t605 + t627 * t525 + t526 * t611;
t596 = qJD(2) * t567;
t563 = t566 * t596;
t537 = qJD(3) * t583 + t563;
t565 = qJD(2) * t569;
t547 = qJD(3) * t582 + t565;
t500 = qJD(4) * t522 + t537;
t517 = qJD(4) * t544 + t547;
t594 = t568 * t596;
t527 = t557 * pkin(2) + pkin(9) * t584;
t528 = t559 * pkin(2) + pkin(9) * t583;
t593 = t527 * t563 + t528 * t594 + qJD(1);
t538 = qJD(3) * t584 - t594;
t546 = pkin(2) * t603 + pkin(9) * t582;
t587 = t528 * t565 - t546 * t563;
t501 = qJD(4) * t520 + t538;
t495 = pkin(3) * t521 + pkin(10) * t520;
t496 = pkin(3) * t523 + pkin(10) * t522;
t586 = t537 * t495 - t496 * t538 + t593;
t585 = (-t527 * t569 - t546 * t602) * qJD(2);
t514 = pkin(3) * t545 + pkin(10) * t544;
t581 = t547 * t496 - t514 * t537 + t587;
t474 = pkin(4) * t505 + pkin(11) * t504;
t475 = pkin(4) * t507 + pkin(11) * t506;
t580 = t500 * t474 - t475 * t501 + t586;
t579 = -t495 * t547 + t538 * t514 + t585;
t497 = pkin(4) * t526 + pkin(11) * t525;
t578 = t517 * t475 - t497 * t500 + t581;
t577 = -t474 * t517 + t501 * t497 + t579;
t553 = t569 * rSges(3,3) + (rSges(3,1) * t574 + rSges(3,2) * t576) * t567;
t552 = Icges(3,5) * t569 + (Icges(3,1) * t574 + Icges(3,4) * t576) * t567;
t551 = Icges(3,6) * t569 + (Icges(3,4) * t574 + Icges(3,2) * t576) * t567;
t550 = Icges(3,3) * t569 + (Icges(3,5) * t574 + Icges(3,6) * t576) * t567;
t536 = rSges(3,1) * t559 + rSges(3,2) * t558 + rSges(3,3) * t604;
t535 = rSges(3,1) * t557 + rSges(3,2) * t556 - rSges(3,3) * t602;
t534 = Icges(3,1) * t559 + Icges(3,4) * t558 + Icges(3,5) * t604;
t533 = Icges(3,1) * t557 + Icges(3,4) * t556 - Icges(3,5) * t602;
t532 = Icges(3,4) * t559 + Icges(3,2) * t558 + Icges(3,6) * t604;
t531 = Icges(3,4) * t557 + Icges(3,2) * t556 - Icges(3,6) * t602;
t530 = Icges(3,5) * t559 + Icges(3,6) * t558 + Icges(3,3) * t604;
t529 = Icges(3,5) * t557 + Icges(3,6) * t556 - Icges(3,3) * t602;
t513 = (-t535 * t569 - t553 * t602) * qJD(2);
t512 = (t536 * t569 - t553 * t604) * qJD(2);
t511 = t545 * rSges(4,1) - t544 * rSges(4,2) + rSges(4,3) * t582;
t510 = Icges(4,1) * t545 - Icges(4,4) * t544 + Icges(4,5) * t582;
t509 = Icges(4,4) * t545 - Icges(4,2) * t544 + Icges(4,6) * t582;
t508 = Icges(4,5) * t545 - Icges(4,6) * t544 + Icges(4,3) * t582;
t499 = qJD(1) + (t535 * t566 + t536 * t568) * t596;
t494 = qJD(5) * t525 + t517;
t492 = rSges(5,1) * t526 - rSges(5,2) * t525 + rSges(5,3) * t544;
t491 = t523 * rSges(4,1) - t522 * rSges(4,2) + rSges(4,3) * t583;
t490 = t521 * rSges(4,1) - t520 * rSges(4,2) + rSges(4,3) * t584;
t489 = Icges(5,1) * t526 - Icges(5,4) * t525 + Icges(5,5) * t544;
t488 = Icges(5,4) * t526 - Icges(5,2) * t525 + Icges(5,6) * t544;
t487 = Icges(5,5) * t526 - Icges(5,6) * t525 + Icges(5,3) * t544;
t486 = Icges(4,1) * t523 - Icges(4,4) * t522 + Icges(4,5) * t583;
t485 = Icges(4,1) * t521 - Icges(4,4) * t520 + Icges(4,5) * t584;
t484 = Icges(4,4) * t523 - Icges(4,2) * t522 + Icges(4,6) * t583;
t483 = Icges(4,4) * t521 - Icges(4,2) * t520 + Icges(4,6) * t584;
t482 = Icges(4,5) * t523 - Icges(4,6) * t522 + Icges(4,3) * t583;
t481 = Icges(4,5) * t521 - Icges(4,6) * t520 + Icges(4,3) * t584;
t473 = qJD(5) * t504 + t501;
t472 = qJD(5) * t506 + t500;
t470 = rSges(6,1) * t503 + rSges(6,2) * t502 + rSges(6,3) * t525;
t468 = rSges(5,1) * t507 - rSges(5,2) * t506 + rSges(5,3) * t522;
t467 = rSges(5,1) * t505 - rSges(5,2) * t504 + rSges(5,3) * t520;
t466 = Icges(5,1) * t507 - Icges(5,4) * t506 + Icges(5,5) * t522;
t465 = Icges(5,1) * t505 - Icges(5,4) * t504 + Icges(5,5) * t520;
t462 = Icges(5,4) * t507 - Icges(5,2) * t506 + Icges(5,6) * t522;
t461 = Icges(5,4) * t505 - Icges(5,2) * t504 + Icges(5,6) * t520;
t458 = Icges(5,5) * t507 - Icges(5,6) * t506 + Icges(5,3) * t522;
t457 = Icges(5,5) * t505 - Icges(5,6) * t504 + Icges(5,3) * t520;
t451 = -t490 * t547 + t511 * t538 + t585;
t450 = t491 * t547 - t511 * t537 + t587;
t449 = rSges(6,1) * t480 + rSges(6,2) * t479 + rSges(6,3) * t506;
t447 = rSges(6,1) * t478 + rSges(6,2) * t477 + rSges(6,3) * t504;
t431 = t490 * t537 - t491 * t538 + t593;
t430 = -t467 * t517 + t492 * t501 + t579;
t429 = t468 * t517 - t492 * t500 + t581;
t428 = t467 * t500 - t468 * t501 + t586;
t427 = -t447 * t494 + t470 * t473 + t577;
t426 = t449 * t494 - t470 * t472 + t578;
t425 = t447 * t472 - t449 * t473 + t580;
t424 = qJD(6) * t506 + t473 * t597 - t494 * t599 + t577;
t423 = qJD(6) * t504 - t472 * t597 + t494 * t598 + t578;
t422 = qJD(6) * t525 + t472 * t599 - t473 * t598 + t580;
t1 = m(3) * (t499 ^ 2 + t512 ^ 2 + t513 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(6) * (t425 ^ 2 + t426 ^ 2 + t427 ^ 2) / 0.2e1 + m(7) * (t422 ^ 2 + t423 ^ 2 + t424 ^ 2) / 0.2e1 + m(5) * (t428 ^ 2 + t429 ^ 2 + t430 ^ 2) / 0.2e1 + m(4) * (t431 ^ 2 + t450 ^ 2 + t451 ^ 2) / 0.2e1 + t500 * ((t522 * t458 - t506 * t462 + t507 * t466) * t500 + (t457 * t522 - t461 * t506 + t465 * t507) * t501 + (t487 * t522 - t488 * t506 + t489 * t507) * t517) / 0.2e1 + t517 * ((t458 * t544 - t462 * t525 + t466 * t526) * t500 + (t457 * t544 - t525 * t461 + t526 * t465) * t501 + (t544 * t487 - t525 * t488 + t526 * t489) * t517) / 0.2e1 + t501 * ((t458 * t520 - t462 * t504 + t466 * t505) * t500 + (t520 * t457 - t504 * t461 + t505 * t465) * t501 + (t487 * t520 - t488 * t504 + t489 * t505) * t517) / 0.2e1 + t537 * ((t583 * t482 - t522 * t484 + t523 * t486) * t537 + (t481 * t583 - t522 * t483 + t523 * t485) * t538 + (t508 * t583 - t522 * t509 + t523 * t510) * t547) / 0.2e1 + t538 * ((t482 * t584 - t520 * t484 + t521 * t486) * t537 + (t584 * t481 - t520 * t483 + t521 * t485) * t538 + (t508 * t584 - t520 * t509 + t521 * t510) * t547) / 0.2e1 + t547 * ((t482 * t582 - t544 * t484 + t545 * t486) * t537 + (t481 * t582 - t544 * t483 + t545 * t485) * t538 + (t582 * t508 - t544 * t509 + t545 * t510) * t547) / 0.2e1 - t617 * ((-t530 * t602 + t532 * t556 + t534 * t557) * t604 - (-t529 * t602 + t531 * t556 + t533 * t557) * t602 + (-t550 * t602 + t551 * t556 + t552 * t557) * t569) * t602 / 0.2e1 + ((t479 * t619 + t480 * t618 + t506 * t620) * t494 + (t624 * t479 + t622 * t480 + t506 * t626) * t473 + (t623 * t479 + t621 * t480 + t625 * t506) * t472) * t472 / 0.2e1 + ((t477 * t619 + t478 * t618 + t504 * t620) * t494 + (t624 * t477 + t622 * t478 + t626 * t504) * t473 + (t477 * t623 + t478 * t621 + t504 * t625) * t472) * t473 / 0.2e1 + ((t619 * t502 + t618 * t503 + t620 * t525) * t494 + (t624 * t502 + t622 * t503 + t525 * t626) * t473 + (t502 * t623 + t503 * t621 + t525 * t625) * t472) * t494 / 0.2e1 + (t569 * (t569 ^ 2 * t550 + (((t532 * t576 + t534 * t574) * t566 - (t531 * t576 + t533 * t574) * t568) * t567 + (-t529 * t568 + t530 * t566 + t551 * t576 + t552 * t574) * t569) * t567) + ((t530 * t604 + t532 * t558 + t534 * t559) * t604 - (t529 * t604 + t531 * t558 + t533 * t559) * t602 + (t550 * t604 + t551 * t558 + t552 * t559) * t569) * t604) * t617 / 0.2e1;
T  = t1;
