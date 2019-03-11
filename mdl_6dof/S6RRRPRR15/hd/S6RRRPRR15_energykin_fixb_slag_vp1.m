% Calculate kinetic energy for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR15_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR15_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:31
% EndTime: 2019-03-09 20:24:34
% DurationCPUTime: 3.17s
% Computational Cost: add. (4152->346), mult. (11255->530), div. (0->0), fcn. (14306->14), ass. (0->164)
t634 = Icges(4,1) + Icges(5,2);
t633 = Icges(5,1) + Icges(4,3);
t632 = -Icges(4,4) - Icges(5,6);
t631 = -Icges(5,4) + Icges(4,5);
t630 = Icges(5,5) - Icges(4,6);
t629 = Icges(4,2) + Icges(5,3);
t573 = sin(pkin(6));
t574 = cos(pkin(6));
t581 = cos(qJ(2));
t582 = cos(qJ(1));
t604 = t581 * t582;
t578 = sin(qJ(2));
t579 = sin(qJ(1));
t607 = t578 * t579;
t558 = t574 * t604 - t607;
t605 = t579 * t581;
t606 = t578 * t582;
t559 = t574 * t606 + t605;
t560 = -t574 * t605 - t606;
t561 = -t574 * t607 + t604;
t608 = t573 * t582;
t609 = t573 * t579;
t592 = (Icges(3,5) * t559 + Icges(3,6) * t558 - Icges(3,3) * t608) * t582 - (Icges(3,5) * t561 + Icges(3,6) * t560 + Icges(3,3) * t609) * t579;
t628 = t573 * t592;
t577 = sin(qJ(3));
t611 = sin(pkin(7));
t614 = cos(qJ(3));
t593 = t614 * t611;
t591 = t573 * t593;
t612 = cos(pkin(7));
t594 = t612 * t614;
t522 = -t558 * t594 + t559 * t577 + t582 * t591;
t596 = t577 * t611;
t597 = t577 * t612;
t523 = t558 * t597 + t559 * t614 - t596 * t608;
t599 = t573 * t612;
t546 = -t558 * t611 - t582 * t599;
t627 = t629 * t522 + t632 * t523 + t630 * t546;
t524 = -t560 * t594 + t561 * t577 - t579 * t591;
t598 = t573 * t611;
t525 = t561 * t614 + (t560 * t612 + t579 * t598) * t577;
t547 = -t560 * t611 + t579 * t599;
t626 = t629 * t524 + t632 * t525 + t630 * t547;
t625 = t630 * t522 + t631 * t523 + t633 * t546;
t624 = t630 * t524 + t631 * t525 + t633 * t547;
t623 = t632 * t522 + t634 * t523 + t631 * t546;
t622 = t632 * t524 + t634 * t525 + t631 * t547;
t610 = t573 * t578;
t544 = -t573 * t581 * t594 - t574 * t593 + t577 * t610;
t545 = t574 * t596 + (t578 * t614 + t581 * t597) * t573;
t557 = t574 * t612 - t581 * t598;
t621 = t629 * t544 + t632 * t545 + t630 * t557;
t620 = t630 * t544 + t631 * t545 + t633 * t557;
t619 = t632 * t544 + t634 * t545 + t631 * t557;
t613 = cos(qJ(5));
t527 = t559 * pkin(2) + pkin(10) * t546;
t528 = t561 * pkin(2) + pkin(10) * t547;
t601 = qJD(2) * t573;
t570 = t579 * t601;
t600 = t582 * t601;
t603 = t527 * t570 + t528 * t600;
t537 = qJD(3) * t547 + t570;
t602 = qJD(1) * (pkin(1) * t579 - pkin(9) * t608);
t571 = qJD(2) * t574 + qJD(1);
t488 = qJD(5) * t525 + t537;
t548 = qJD(3) * t557 + t571;
t484 = pkin(3) * t523 + qJ(4) * t522;
t595 = qJD(4) * t544 + t537 * t484 + t603;
t510 = qJD(5) * t545 + t548;
t538 = qJD(3) * t546 - t600;
t489 = qJD(5) * t523 + t538;
t549 = pkin(2) * t610 + pkin(10) * t557;
t562 = qJD(1) * (pkin(1) * t582 + pkin(9) * t609);
t590 = t571 * t528 - t549 * t570 + t562;
t485 = pkin(3) * t525 + qJ(4) * t524;
t589 = qJD(4) * t522 + t548 * t485 + t590;
t499 = pkin(4) * t546 + pkin(11) * t523;
t500 = pkin(4) * t547 + pkin(11) * t525;
t588 = t537 * t499 + (-t485 - t500) * t538 + t595;
t587 = -t527 * t571 - t549 * t600 - t602;
t509 = pkin(3) * t545 + qJ(4) * t544;
t586 = qJD(4) * t524 + t538 * t509 + t587;
t526 = pkin(4) * t557 + pkin(11) * t545;
t585 = t548 * t500 + (-t509 - t526) * t537 + t589;
t584 = t538 * t526 + (-t484 - t499) * t548 + t586;
t580 = cos(qJ(6));
t576 = sin(qJ(5));
t575 = sin(qJ(6));
t568 = rSges(2,1) * t582 - rSges(2,2) * t579;
t567 = rSges(2,1) * t579 + rSges(2,2) * t582;
t555 = rSges(3,3) * t574 + (rSges(3,1) * t578 + rSges(3,2) * t581) * t573;
t554 = Icges(3,5) * t574 + (Icges(3,1) * t578 + Icges(3,4) * t581) * t573;
t553 = Icges(3,6) * t574 + (Icges(3,4) * t578 + Icges(3,2) * t581) * t573;
t552 = Icges(3,3) * t574 + (Icges(3,5) * t578 + Icges(3,6) * t581) * t573;
t536 = rSges(3,1) * t561 + rSges(3,2) * t560 + rSges(3,3) * t609;
t535 = rSges(3,1) * t559 + rSges(3,2) * t558 - rSges(3,3) * t608;
t534 = Icges(3,1) * t561 + Icges(3,4) * t560 + Icges(3,5) * t609;
t533 = Icges(3,1) * t559 + Icges(3,4) * t558 - Icges(3,5) * t608;
t532 = Icges(3,4) * t561 + Icges(3,2) * t560 + Icges(3,6) * t609;
t531 = Icges(3,4) * t559 + Icges(3,2) * t558 - Icges(3,6) * t608;
t521 = t544 * t576 + t557 * t613;
t520 = -t544 * t613 + t557 * t576;
t508 = rSges(5,1) * t557 - rSges(5,2) * t545 + rSges(5,3) * t544;
t507 = rSges(4,1) * t545 - rSges(4,2) * t544 + rSges(4,3) * t557;
t498 = t524 * t576 + t547 * t613;
t497 = -t524 * t613 + t547 * t576;
t496 = t522 * t576 + t546 * t613;
t495 = -t522 * t613 + t546 * t576;
t493 = t521 * t580 + t545 * t575;
t492 = -t521 * t575 + t545 * t580;
t491 = t536 * t571 - t555 * t570 + t562;
t490 = -t535 * t571 - t555 * t600 - t602;
t487 = (t535 * t579 + t536 * t582) * t601;
t482 = pkin(5) * t521 + pkin(12) * t520;
t481 = qJD(6) * t520 + t510;
t478 = rSges(5,1) * t547 - rSges(5,2) * t525 + rSges(5,3) * t524;
t477 = rSges(5,1) * t546 - rSges(5,2) * t523 + rSges(5,3) * t522;
t476 = rSges(4,1) * t525 - rSges(4,2) * t524 + rSges(4,3) * t547;
t475 = rSges(4,1) * t523 - rSges(4,2) * t522 + rSges(4,3) * t546;
t462 = rSges(6,1) * t521 - rSges(6,2) * t520 + rSges(6,3) * t545;
t461 = Icges(6,1) * t521 - Icges(6,4) * t520 + Icges(6,5) * t545;
t460 = Icges(6,4) * t521 - Icges(6,2) * t520 + Icges(6,6) * t545;
t459 = Icges(6,5) * t521 - Icges(6,6) * t520 + Icges(6,3) * t545;
t458 = t498 * t580 + t525 * t575;
t457 = -t498 * t575 + t525 * t580;
t456 = t496 * t580 + t523 * t575;
t455 = -t496 * t575 + t523 * t580;
t453 = pkin(5) * t498 + pkin(12) * t497;
t452 = pkin(5) * t496 + pkin(12) * t495;
t451 = qJD(6) * t495 + t489;
t450 = qJD(6) * t497 + t488;
t449 = rSges(6,1) * t498 - rSges(6,2) * t497 + rSges(6,3) * t525;
t448 = rSges(6,1) * t496 - rSges(6,2) * t495 + rSges(6,3) * t523;
t447 = Icges(6,1) * t498 - Icges(6,4) * t497 + Icges(6,5) * t525;
t446 = Icges(6,1) * t496 - Icges(6,4) * t495 + Icges(6,5) * t523;
t445 = Icges(6,4) * t498 - Icges(6,2) * t497 + Icges(6,6) * t525;
t444 = Icges(6,4) * t496 - Icges(6,2) * t495 + Icges(6,6) * t523;
t443 = Icges(6,5) * t498 - Icges(6,6) * t497 + Icges(6,3) * t525;
t442 = Icges(6,5) * t496 - Icges(6,6) * t495 + Icges(6,3) * t523;
t441 = rSges(7,1) * t493 + rSges(7,2) * t492 + rSges(7,3) * t520;
t440 = Icges(7,1) * t493 + Icges(7,4) * t492 + Icges(7,5) * t520;
t439 = Icges(7,4) * t493 + Icges(7,2) * t492 + Icges(7,6) * t520;
t438 = Icges(7,5) * t493 + Icges(7,6) * t492 + Icges(7,3) * t520;
t437 = rSges(7,1) * t458 + rSges(7,2) * t457 + rSges(7,3) * t497;
t436 = rSges(7,1) * t456 + rSges(7,2) * t455 + rSges(7,3) * t495;
t435 = t476 * t548 - t507 * t537 + t590;
t434 = -t475 * t548 + t507 * t538 + t587;
t433 = Icges(7,1) * t458 + Icges(7,4) * t457 + Icges(7,5) * t497;
t432 = Icges(7,1) * t456 + Icges(7,4) * t455 + Icges(7,5) * t495;
t431 = Icges(7,4) * t458 + Icges(7,2) * t457 + Icges(7,6) * t497;
t430 = Icges(7,4) * t456 + Icges(7,2) * t455 + Icges(7,6) * t495;
t429 = Icges(7,5) * t458 + Icges(7,6) * t457 + Icges(7,3) * t497;
t428 = Icges(7,5) * t456 + Icges(7,6) * t455 + Icges(7,3) * t495;
t427 = t475 * t537 - t476 * t538 + t603;
t426 = t478 * t548 + (-t508 - t509) * t537 + t589;
t425 = t508 * t538 + (-t477 - t484) * t548 + t586;
t424 = t477 * t537 + (-t478 - t485) * t538 + t595;
t423 = t449 * t510 - t462 * t488 + t585;
t422 = -t448 * t510 + t462 * t489 + t584;
t421 = t448 * t488 - t449 * t489 + t588;
t420 = t437 * t481 - t441 * t450 + t453 * t510 - t482 * t488 + t585;
t419 = -t436 * t481 + t441 * t451 - t452 * t510 + t482 * t489 + t584;
t418 = t436 * t450 - t437 * t451 + t452 * t488 - t453 * t489 + t588;
t1 = ((t552 * t609 + t553 * t560 + t554 * t561) * t571 + (-(t531 * t560 + t533 * t561) * t582 + (t560 * t532 + t561 * t534 - t628) * t579) * t601) * t570 / 0.2e1 - ((-t552 * t608 + t553 * t558 + t554 * t559) * t571 + ((t532 * t558 + t534 * t559) * t579 + (-t558 * t531 - t559 * t533 + t628) * t582) * t601) * t600 / 0.2e1 + t451 * ((t429 * t495 + t431 * t455 + t433 * t456) * t450 + (t428 * t495 + t430 * t455 + t432 * t456) * t451 + (t438 * t495 + t439 * t455 + t440 * t456) * t481) / 0.2e1 + t481 * ((t429 * t520 + t431 * t492 + t433 * t493) * t450 + (t428 * t520 + t430 * t492 + t432 * t493) * t451 + (t438 * t520 + t439 * t492 + t440 * t493) * t481) / 0.2e1 + t450 * ((t429 * t497 + t431 * t457 + t433 * t458) * t450 + (t428 * t497 + t430 * t457 + t432 * t458) * t451 + (t438 * t497 + t439 * t457 + t440 * t458) * t481) / 0.2e1 + t510 * ((t443 * t545 - t445 * t520 + t447 * t521) * t488 + (t442 * t545 - t444 * t520 + t446 * t521) * t489 + (t459 * t545 - t460 * t520 + t461 * t521) * t510) / 0.2e1 + t489 * ((t443 * t523 - t445 * t495 + t447 * t496) * t488 + (t442 * t523 - t444 * t495 + t446 * t496) * t489 + (t459 * t523 - t460 * t495 + t461 * t496) * t510) / 0.2e1 + t488 * ((t443 * t525 - t445 * t497 + t447 * t498) * t488 + (t442 * t525 - t444 * t497 + t446 * t498) * t489 + (t459 * t525 - t460 * t497 + t461 * t498) * t510) / 0.2e1 + t571 * ((t574 * t552 + (t553 * t581 + t554 * t578) * t573) * t571 + (((t532 * t581 + t534 * t578) * t579 - (t531 * t581 + t533 * t578) * t582) * t573 - t592 * t574) * t601) / 0.2e1 + m(7) * (t418 ^ 2 + t419 ^ 2 + t420 ^ 2) / 0.2e1 + m(6) * (t421 ^ 2 + t422 ^ 2 + t423 ^ 2) / 0.2e1 + m(5) * (t424 ^ 2 + t425 ^ 2 + t426 ^ 2) / 0.2e1 + m(4) * (t427 ^ 2 + t434 ^ 2 + t435 ^ 2) / 0.2e1 + m(3) * (t487 ^ 2 + t490 ^ 2 + t491 ^ 2) / 0.2e1 + ((t524 * t621 + t525 * t619 + t547 * t620) * t548 + (t524 * t627 + t525 * t623 + t547 * t625) * t538 + (t524 * t626 + t525 * t622 + t547 * t624) * t537) * t537 / 0.2e1 + ((t522 * t621 + t523 * t619 + t546 * t620) * t548 + (t522 * t627 + t523 * t623 + t546 * t625) * t538 + (t522 * t626 + t523 * t622 + t546 * t624) * t537) * t538 / 0.2e1 + ((t544 * t621 + t545 * t619 + t557 * t620) * t548 + (t544 * t627 + t545 * t623 + t557 * t625) * t538 + (t544 * t626 + t545 * t622 + t557 * t624) * t537) * t548 / 0.2e1 + (Icges(2,3) + m(2) * (t567 ^ 2 + t568 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
