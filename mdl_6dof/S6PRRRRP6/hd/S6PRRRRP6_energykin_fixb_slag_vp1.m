% Calculate kinetic energy for
% S6PRRRRP6
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:33
% EndTime: 2019-03-09 00:28:36
% DurationCPUTime: 3.14s
% Computational Cost: add. (5282->329), mult. (14495->516), div. (0->0), fcn. (18715->14), ass. (0->156)
t643 = Icges(6,1) + Icges(7,1);
t642 = -Icges(6,4) + Icges(7,5);
t641 = Icges(7,4) + Icges(6,5);
t640 = Icges(6,2) + Icges(7,3);
t639 = Icges(7,2) + Icges(6,3);
t638 = -Icges(6,6) + Icges(7,6);
t637 = rSges(7,1) + pkin(5);
t636 = rSges(7,3) + qJ(6);
t584 = sin(pkin(12));
t586 = cos(pkin(12));
t591 = sin(qJ(2));
t587 = cos(pkin(6));
t592 = cos(qJ(2));
t613 = t587 * t592;
t575 = -t584 * t591 + t586 * t613;
t614 = t587 * t591;
t576 = t584 * t592 + t586 * t614;
t590 = sin(qJ(3));
t585 = sin(pkin(6));
t618 = sin(pkin(7));
t604 = t585 * t618;
t619 = cos(pkin(7));
t622 = cos(qJ(3));
t536 = t576 * t622 + (t575 * t619 - t586 * t604) * t590;
t605 = t585 * t619;
t562 = -t575 * t618 - t586 * t605;
t589 = sin(qJ(4));
t621 = cos(qJ(4));
t518 = t536 * t621 + t562 * t589;
t602 = t622 * t618;
t601 = t585 * t602;
t603 = t619 * t622;
t535 = -t575 * t603 + t576 * t590 + t586 * t601;
t588 = sin(qJ(5));
t620 = cos(qJ(5));
t490 = t518 * t588 - t535 * t620;
t491 = t518 * t620 + t535 * t588;
t517 = t536 * t589 - t562 * t621;
t635 = t490 * t640 + t491 * t642 + t517 * t638;
t577 = -t584 * t613 - t586 * t591;
t578 = -t584 * t614 + t586 * t592;
t538 = t578 * t622 + (t577 * t619 + t584 * t604) * t590;
t563 = -t577 * t618 + t584 * t605;
t520 = t538 * t621 + t563 * t589;
t537 = -t577 * t603 + t578 * t590 - t584 * t601;
t492 = t520 * t588 - t537 * t620;
t493 = t520 * t620 + t537 * t588;
t519 = t538 * t589 - t563 * t621;
t634 = t492 * t640 + t493 * t642 + t519 * t638;
t633 = t490 * t638 + t491 * t641 + t517 * t639;
t632 = t492 * t638 + t493 * t641 + t519 * t639;
t631 = t642 * t490 + t491 * t643 + t641 * t517;
t630 = t642 * t492 + t493 * t643 + t641 * t519;
t561 = t587 * t618 * t590 + (t590 * t592 * t619 + t591 * t622) * t585;
t574 = t587 * t619 - t592 * t604;
t541 = t561 * t621 + t574 * t589;
t615 = t585 * t591;
t560 = -t585 * t592 * t603 - t587 * t602 + t590 * t615;
t515 = t541 * t588 - t560 * t620;
t516 = t541 * t620 + t560 * t588;
t540 = t561 * t589 - t574 * t621;
t629 = t515 * t640 + t516 * t642 + t540 * t638;
t628 = t515 * t638 + t516 * t641 + t540 * t639;
t627 = t642 * t515 + t516 * t643 + t641 * t540;
t626 = qJD(2) ^ 2;
t617 = t584 * t585;
t616 = t585 * t586;
t612 = rSges(7,2) * t517 + t636 * t490 + t637 * t491;
t611 = rSges(7,2) * t519 + t636 * t492 + t637 * t493;
t610 = rSges(7,2) * t540 + t636 * t515 + t637 * t516;
t609 = qJD(2) * t585;
t582 = t584 * t609;
t552 = qJD(3) * t563 + t582;
t583 = qJD(2) * t587;
t565 = qJD(3) * t574 + t583;
t513 = qJD(4) * t537 + t552;
t532 = qJD(4) * t560 + t565;
t607 = t586 * t609;
t542 = t576 * pkin(2) + pkin(9) * t562;
t543 = t578 * pkin(2) + pkin(9) * t563;
t606 = t542 * t582 + t543 * t607 + qJD(1);
t553 = qJD(3) * t562 - t607;
t564 = pkin(2) * t615 + pkin(9) * t574;
t600 = t543 * t583 - t564 * t582;
t514 = qJD(4) * t535 + t553;
t508 = pkin(3) * t536 + pkin(10) * t535;
t509 = pkin(3) * t538 + pkin(10) * t537;
t599 = t552 * t508 - t509 * t553 + t606;
t598 = (-t542 * t587 - t564 * t616) * qJD(2);
t527 = pkin(3) * t561 + pkin(10) * t560;
t597 = t565 * t509 - t527 * t552 + t600;
t487 = pkin(4) * t518 + pkin(11) * t517;
t488 = pkin(4) * t520 + pkin(11) * t519;
t596 = t513 * t487 - t488 * t514 + t599;
t595 = -t508 * t565 + t553 * t527 + t598;
t510 = pkin(4) * t541 + pkin(11) * t540;
t594 = t532 * t488 - t510 * t513 + t597;
t593 = -t487 * t532 + t514 * t510 + t595;
t571 = t587 * rSges(3,3) + (rSges(3,1) * t591 + rSges(3,2) * t592) * t585;
t570 = Icges(3,5) * t587 + (Icges(3,1) * t591 + Icges(3,4) * t592) * t585;
t569 = Icges(3,6) * t587 + (Icges(3,4) * t591 + Icges(3,2) * t592) * t585;
t568 = Icges(3,3) * t587 + (Icges(3,5) * t591 + Icges(3,6) * t592) * t585;
t551 = rSges(3,1) * t578 + rSges(3,2) * t577 + rSges(3,3) * t617;
t550 = rSges(3,1) * t576 + rSges(3,2) * t575 - rSges(3,3) * t616;
t549 = Icges(3,1) * t578 + Icges(3,4) * t577 + Icges(3,5) * t617;
t548 = Icges(3,1) * t576 + Icges(3,4) * t575 - Icges(3,5) * t616;
t547 = Icges(3,4) * t578 + Icges(3,2) * t577 + Icges(3,6) * t617;
t546 = Icges(3,4) * t576 + Icges(3,2) * t575 - Icges(3,6) * t616;
t545 = Icges(3,5) * t578 + Icges(3,6) * t577 + Icges(3,3) * t617;
t544 = Icges(3,5) * t576 + Icges(3,6) * t575 - Icges(3,3) * t616;
t526 = (-t550 * t587 - t571 * t616) * qJD(2);
t525 = (t551 * t587 - t571 * t617) * qJD(2);
t524 = rSges(4,1) * t561 - rSges(4,2) * t560 + rSges(4,3) * t574;
t523 = Icges(4,1) * t561 - Icges(4,4) * t560 + Icges(4,5) * t574;
t522 = Icges(4,4) * t561 - Icges(4,2) * t560 + Icges(4,6) * t574;
t521 = Icges(4,5) * t561 - Icges(4,6) * t560 + Icges(4,3) * t574;
t512 = qJD(1) + (t550 * t584 + t551 * t586) * t609;
t507 = qJD(5) * t540 + t532;
t505 = rSges(5,1) * t541 - rSges(5,2) * t540 + rSges(5,3) * t560;
t504 = rSges(4,1) * t538 - rSges(4,2) * t537 + rSges(4,3) * t563;
t503 = rSges(4,1) * t536 - rSges(4,2) * t535 + rSges(4,3) * t562;
t502 = Icges(5,1) * t541 - Icges(5,4) * t540 + Icges(5,5) * t560;
t501 = Icges(5,4) * t541 - Icges(5,2) * t540 + Icges(5,6) * t560;
t500 = Icges(5,5) * t541 - Icges(5,6) * t540 + Icges(5,3) * t560;
t499 = Icges(4,1) * t538 - Icges(4,4) * t537 + Icges(4,5) * t563;
t498 = Icges(4,1) * t536 - Icges(4,4) * t535 + Icges(4,5) * t562;
t497 = Icges(4,4) * t538 - Icges(4,2) * t537 + Icges(4,6) * t563;
t496 = Icges(4,4) * t536 - Icges(4,2) * t535 + Icges(4,6) * t562;
t495 = Icges(4,5) * t538 - Icges(4,6) * t537 + Icges(4,3) * t563;
t494 = Icges(4,5) * t536 - Icges(4,6) * t535 + Icges(4,3) * t562;
t485 = qJD(5) * t517 + t514;
t484 = qJD(5) * t519 + t513;
t482 = rSges(6,1) * t516 - rSges(6,2) * t515 + rSges(6,3) * t540;
t480 = rSges(5,1) * t520 - rSges(5,2) * t519 + rSges(5,3) * t537;
t479 = rSges(5,1) * t518 - rSges(5,2) * t517 + rSges(5,3) * t535;
t478 = Icges(5,1) * t520 - Icges(5,4) * t519 + Icges(5,5) * t537;
t477 = Icges(5,1) * t518 - Icges(5,4) * t517 + Icges(5,5) * t535;
t474 = Icges(5,4) * t520 - Icges(5,2) * t519 + Icges(5,6) * t537;
t473 = Icges(5,4) * t518 - Icges(5,2) * t517 + Icges(5,6) * t535;
t470 = Icges(5,5) * t520 - Icges(5,6) * t519 + Icges(5,3) * t537;
t469 = Icges(5,5) * t518 - Icges(5,6) * t517 + Icges(5,3) * t535;
t462 = -t503 * t565 + t524 * t553 + t598;
t461 = t504 * t565 - t524 * t552 + t600;
t460 = rSges(6,1) * t493 - rSges(6,2) * t492 + rSges(6,3) * t519;
t458 = rSges(6,1) * t491 - rSges(6,2) * t490 + rSges(6,3) * t517;
t444 = t503 * t552 - t504 * t553 + t606;
t443 = -t479 * t532 + t505 * t514 + t595;
t442 = t480 * t532 - t505 * t513 + t597;
t441 = t479 * t513 - t480 * t514 + t599;
t440 = -t458 * t507 + t482 * t485 + t593;
t439 = t460 * t507 - t482 * t484 + t594;
t438 = t458 * t484 - t460 * t485 + t596;
t437 = qJD(6) * t492 + t485 * t610 - t507 * t612 + t593;
t436 = qJD(6) * t490 - t484 * t610 + t507 * t611 + t594;
t435 = qJD(6) * t515 + t484 * t612 - t485 * t611 + t596;
t1 = t565 * ((t495 * t574 - t497 * t560 + t499 * t561) * t552 + (t494 * t574 - t496 * t560 + t498 * t561) * t553 + (t521 * t574 - t522 * t560 + t523 * t561) * t565) / 0.2e1 + m(3) * (t512 ^ 2 + t525 ^ 2 + t526 ^ 2) / 0.2e1 + m(6) * (t438 ^ 2 + t439 ^ 2 + t440 ^ 2) / 0.2e1 + m(7) * (t435 ^ 2 + t436 ^ 2 + t437 ^ 2) / 0.2e1 + m(5) * (t441 ^ 2 + t442 ^ 2 + t443 ^ 2) / 0.2e1 + m(4) * (t444 ^ 2 + t461 ^ 2 + t462 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t513 * ((t537 * t470 - t519 * t474 + t520 * t478) * t513 + (t469 * t537 - t473 * t519 + t477 * t520) * t514 + (t500 * t537 - t501 * t519 + t502 * t520) * t532) / 0.2e1 + t514 * ((t470 * t535 - t474 * t517 + t478 * t518) * t513 + (t535 * t469 - t517 * t473 + t518 * t477) * t514 + (t500 * t535 - t501 * t517 + t502 * t518) * t532) / 0.2e1 + t532 * ((t470 * t560 - t474 * t540 + t478 * t541) * t513 + (t469 * t560 - t473 * t540 + t477 * t541) * t514 + (t500 * t560 - t501 * t540 + t502 * t541) * t532) / 0.2e1 + t552 * ((t495 * t563 - t497 * t537 + t499 * t538) * t552 + (t494 * t563 - t496 * t537 + t498 * t538) * t553 + (t521 * t563 - t522 * t537 + t523 * t538) * t565) / 0.2e1 + t553 * ((t495 * t562 - t497 * t535 + t499 * t536) * t552 + (t494 * t562 - t496 * t535 + t498 * t536) * t553 + (t521 * t562 - t522 * t535 + t523 * t536) * t565) / 0.2e1 - t626 * ((-t545 * t616 + t547 * t575 + t549 * t576) * t617 - (-t544 * t616 + t546 * t575 + t548 * t576) * t616 + (-t568 * t616 + t569 * t575 + t570 * t576) * t587) * t616 / 0.2e1 + ((t492 * t629 + t493 * t627 + t519 * t628) * t507 + (t492 * t635 + t631 * t493 + t633 * t519) * t485 + (t634 * t492 + t630 * t493 + t632 * t519) * t484) * t484 / 0.2e1 + ((t490 * t629 + t491 * t627 + t517 * t628) * t507 + (t635 * t490 + t631 * t491 + t633 * t517) * t485 + (t490 * t634 + t491 * t630 + t517 * t632) * t484) * t485 / 0.2e1 + ((t515 * t629 + t516 * t627 + t540 * t628) * t507 + (t515 * t635 + t631 * t516 + t633 * t540) * t485 + (t515 * t634 + t516 * t630 + t540 * t632) * t484) * t507 / 0.2e1 + (t587 * (t587 ^ 2 * t568 + (((t547 * t592 + t549 * t591) * t584 - (t546 * t592 + t548 * t591) * t586) * t585 + (-t544 * t586 + t545 * t584 + t569 * t592 + t570 * t591) * t587) * t585) + ((t545 * t617 + t547 * t577 + t549 * t578) * t617 - (t544 * t617 + t546 * t577 + t548 * t578) * t616 + (t568 * t617 + t569 * t577 + t570 * t578) * t587) * t617) * t626 / 0.2e1;
T  = t1;
