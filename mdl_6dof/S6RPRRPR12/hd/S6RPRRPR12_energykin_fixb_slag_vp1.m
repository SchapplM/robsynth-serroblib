% Calculate kinetic energy for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR12_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR12_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:48:03
% EndTime: 2019-03-09 05:48:06
% DurationCPUTime: 2.50s
% Computational Cost: add. (4340->302), mult. (11887->458), div. (0->0), fcn. (15369->14), ass. (0->144)
t604 = Icges(5,1) + Icges(6,2);
t603 = Icges(6,1) + Icges(5,3);
t602 = -Icges(5,4) - Icges(6,6);
t601 = -Icges(6,4) + Icges(5,5);
t600 = Icges(6,5) - Icges(5,6);
t599 = Icges(5,2) + Icges(6,3);
t547 = cos(pkin(6));
t580 = cos(pkin(12));
t582 = sin(qJ(1));
t569 = t582 * t580;
t545 = sin(pkin(12));
t552 = cos(qJ(1));
t576 = t552 * t545;
t534 = t547 * t576 + t569;
t550 = sin(qJ(3));
t571 = t552 * t580;
t574 = t582 * t545;
t560 = -t547 * t571 + t574;
t581 = cos(pkin(7));
t555 = t560 * t581;
t546 = sin(pkin(6));
t579 = sin(pkin(7));
t572 = t546 * t579;
t584 = cos(qJ(3));
t513 = t534 * t584 + (-t552 * t572 - t555) * t550;
t573 = t546 * t581;
t526 = -t552 * t573 + t560 * t579;
t549 = sin(qJ(4));
t583 = cos(qJ(4));
t493 = t513 * t549 - t526 * t583;
t494 = t513 * t583 + t526 * t549;
t570 = t584 * t579;
t577 = t546 * t552;
t512 = t534 * t550 + t555 * t584 + t570 * t577;
t598 = t599 * t493 + t602 * t494 + t600 * t512;
t535 = -t547 * t574 + t571;
t561 = t547 * t569 + t576;
t589 = t561 * t581 - t582 * t572;
t515 = t535 * t584 - t589 * t550;
t527 = t561 * t579 + t573 * t582;
t495 = t515 * t549 - t527 * t583;
t496 = t515 * t583 + t527 * t549;
t514 = t535 * t550 + t589 * t584;
t597 = t599 * t495 + t602 * t496 + t600 * t514;
t596 = t600 * t493 + t601 * t494 + t603 * t512;
t595 = t600 * t495 + t601 * t496 + t603 * t514;
t594 = t602 * t493 + t604 * t494 + t601 * t512;
t593 = t602 * t495 + t604 * t496 + t601 * t514;
t567 = t581 * t580;
t525 = t547 * t579 * t550 + (t545 * t584 + t550 * t567) * t546;
t533 = t547 * t581 - t572 * t580;
t510 = t525 * t549 - t533 * t583;
t511 = t525 * t583 + t533 * t549;
t578 = t546 * t545;
t524 = -t546 * t567 * t584 - t547 * t570 + t550 * t578;
t592 = t599 * t510 + t602 * t511 + t600 * t524;
t591 = t600 * t510 + t601 * t511 + t603 * t524;
t590 = t602 * t510 + t604 * t511 + t601 * t524;
t522 = qJD(3) * t526;
t498 = qJD(4) * t512 + t522;
t523 = qJD(3) * t527;
t499 = qJD(4) * t514 + t523;
t530 = qJD(3) * t533 + qJD(1);
t575 = t546 * t582;
t516 = qJD(4) * t524 + t530;
t568 = -qJD(2) * t577 + qJD(1) * (t552 * pkin(1) + qJ(2) * t575);
t537 = pkin(1) * t582 - qJ(2) * t577;
t543 = qJD(2) * t575;
t566 = t543 + (-t534 * pkin(2) - pkin(9) * t526 - t537) * qJD(1);
t564 = qJD(1) * (t535 * pkin(2) + pkin(9) * t527) + t568;
t485 = pkin(3) * t513 + pkin(10) * t512;
t486 = pkin(3) * t515 + pkin(10) * t514;
t544 = qJD(2) * t547;
t563 = t485 * t523 - t486 * t522 + t544;
t456 = pkin(4) * t494 + qJ(5) * t493;
t562 = qJD(5) * t510 + t499 * t456 + t563;
t506 = pkin(3) * t525 + pkin(10) * t524;
t559 = -t485 * t530 + t506 * t522 + t566;
t558 = t530 * t486 - t506 * t523 + t564;
t484 = pkin(4) * t511 + qJ(5) * t510;
t557 = qJD(5) * t495 + t498 * t484 + t559;
t457 = pkin(4) * t496 + qJ(5) * t495;
t554 = qJD(5) * t493 + t516 * t457 + t558;
t551 = cos(qJ(6));
t548 = sin(qJ(6));
t541 = t552 * rSges(2,1) - rSges(2,2) * t582;
t540 = rSges(2,1) * t582 + t552 * rSges(2,2);
t505 = qJD(1) * (t535 * rSges(3,1) - rSges(3,2) * t561 + rSges(3,3) * t575) + t568;
t504 = t543 + (-t534 * rSges(3,1) + rSges(3,2) * t560 + rSges(3,3) * t577 - t537) * qJD(1);
t503 = rSges(4,1) * t525 - rSges(4,2) * t524 + rSges(4,3) * t533;
t502 = Icges(4,1) * t525 - Icges(4,4) * t524 + Icges(4,5) * t533;
t501 = Icges(4,4) * t525 - Icges(4,2) * t524 + Icges(4,6) * t533;
t500 = Icges(4,5) * t525 - Icges(4,6) * t524 + Icges(4,3) * t533;
t497 = pkin(5) * t524 + pkin(11) * t511;
t490 = t510 * t548 + t524 * t551;
t489 = t510 * t551 - t524 * t548;
t487 = qJD(6) * t511 + t516;
t481 = rSges(4,1) * t515 - rSges(4,2) * t514 + rSges(4,3) * t527;
t480 = rSges(4,1) * t513 - rSges(4,2) * t512 + rSges(4,3) * t526;
t479 = Icges(4,1) * t515 - Icges(4,4) * t514 + Icges(4,5) * t527;
t478 = Icges(4,1) * t513 - Icges(4,4) * t512 + Icges(4,5) * t526;
t477 = Icges(4,4) * t515 - Icges(4,2) * t514 + Icges(4,6) * t527;
t476 = Icges(4,4) * t513 - Icges(4,2) * t512 + Icges(4,6) * t526;
t475 = Icges(4,5) * t515 - Icges(4,6) * t514 + Icges(4,3) * t527;
t474 = Icges(4,5) * t513 - Icges(4,6) * t512 + Icges(4,3) * t526;
t473 = rSges(5,1) * t511 - rSges(5,2) * t510 + rSges(5,3) * t524;
t472 = rSges(6,1) * t524 - rSges(6,2) * t511 + rSges(6,3) * t510;
t465 = pkin(5) * t514 + pkin(11) * t496;
t464 = pkin(5) * t512 + pkin(11) * t494;
t463 = t495 * t548 + t514 * t551;
t462 = t495 * t551 - t514 * t548;
t461 = t493 * t548 + t512 * t551;
t460 = t493 * t551 - t512 * t548;
t459 = qJD(6) * t496 + t499;
t458 = qJD(6) * t494 + t498;
t453 = rSges(5,1) * t496 - rSges(5,2) * t495 + rSges(5,3) * t514;
t452 = rSges(5,1) * t494 - rSges(5,2) * t493 + rSges(5,3) * t512;
t451 = rSges(6,1) * t514 - rSges(6,2) * t496 + rSges(6,3) * t495;
t450 = rSges(6,1) * t512 - rSges(6,2) * t494 + rSges(6,3) * t493;
t437 = rSges(7,1) * t490 + rSges(7,2) * t489 + rSges(7,3) * t511;
t436 = Icges(7,1) * t490 + Icges(7,4) * t489 + Icges(7,5) * t511;
t435 = Icges(7,4) * t490 + Icges(7,2) * t489 + Icges(7,6) * t511;
t434 = Icges(7,5) * t490 + Icges(7,6) * t489 + Icges(7,3) * t511;
t432 = t481 * t530 - t503 * t523 + t564;
t431 = -t480 * t530 + t503 * t522 + t566;
t430 = t544 + (t480 * t527 - t481 * t526) * qJD(3);
t429 = rSges(7,1) * t463 + rSges(7,2) * t462 + rSges(7,3) * t496;
t428 = rSges(7,1) * t461 + rSges(7,2) * t460 + rSges(7,3) * t494;
t427 = Icges(7,1) * t463 + Icges(7,4) * t462 + Icges(7,5) * t496;
t426 = Icges(7,1) * t461 + Icges(7,4) * t460 + Icges(7,5) * t494;
t425 = Icges(7,4) * t463 + Icges(7,2) * t462 + Icges(7,6) * t496;
t424 = Icges(7,4) * t461 + Icges(7,2) * t460 + Icges(7,6) * t494;
t423 = Icges(7,5) * t463 + Icges(7,6) * t462 + Icges(7,3) * t496;
t422 = Icges(7,5) * t461 + Icges(7,6) * t460 + Icges(7,3) * t494;
t421 = t453 * t516 - t473 * t499 + t558;
t420 = -t452 * t516 + t473 * t498 + t559;
t419 = t452 * t499 - t453 * t498 + t563;
t418 = t451 * t516 + (-t472 - t484) * t499 + t554;
t417 = t472 * t498 + (-t450 - t456) * t516 + t557;
t416 = t450 * t499 + (-t451 - t457) * t498 + t562;
t415 = t429 * t487 - t437 * t459 + t465 * t516 + (-t484 - t497) * t499 + t554;
t414 = -t428 * t487 + t437 * t458 + t497 * t498 + (-t456 - t464) * t516 + t557;
t413 = t428 * t459 - t429 * t458 + t464 * t499 + (-t457 - t465) * t498 + t562;
t1 = ((t500 * t527 - t501 * t514 + t502 * t515) * t530 + ((t475 * t527 - t477 * t514 + t479 * t515) * t527 + (t474 * t527 - t476 * t514 + t478 * t515) * t526) * qJD(3)) * t523 / 0.2e1 + t487 * ((t423 * t511 + t425 * t489 + t427 * t490) * t459 + (t422 * t511 + t424 * t489 + t426 * t490) * t458 + (t511 * t434 + t489 * t435 + t490 * t436) * t487) / 0.2e1 + t459 * ((t496 * t423 + t462 * t425 + t463 * t427) * t459 + (t422 * t496 + t424 * t462 + t426 * t463) * t458 + (t434 * t496 + t435 * t462 + t436 * t463) * t487) / 0.2e1 + t458 * ((t423 * t494 + t425 * t460 + t427 * t461) * t459 + (t494 * t422 + t460 * t424 + t461 * t426) * t458 + (t434 * t494 + t435 * t460 + t436 * t461) * t487) / 0.2e1 + ((t500 * t526 - t501 * t512 + t502 * t513) * t530 + ((t475 * t526 - t477 * t512 + t479 * t513) * t527 + (t474 * t526 - t476 * t512 + t478 * t513) * t526) * qJD(3)) * t522 / 0.2e1 + t530 * ((t533 * t500 - t524 * t501 + t525 * t502) * t530 + ((t475 * t533 - t477 * t524 + t479 * t525) * t527 + (t474 * t533 - t476 * t524 + t478 * t525) * t526) * qJD(3)) / 0.2e1 + m(7) * (t413 ^ 2 + t414 ^ 2 + t415 ^ 2) / 0.2e1 + m(6) * (t416 ^ 2 + t417 ^ 2 + t418 ^ 2) / 0.2e1 + m(5) * (t419 ^ 2 + t420 ^ 2 + t421 ^ 2) / 0.2e1 + m(4) * (t430 ^ 2 + t431 ^ 2 + t432 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t547 ^ 2 + t504 ^ 2 + t505 ^ 2) / 0.2e1 + ((t592 * t493 + t590 * t494 + t591 * t512) * t516 + (t597 * t493 + t593 * t494 + t595 * t512) * t499 + (t598 * t493 + t594 * t494 + t596 * t512) * t498) * t498 / 0.2e1 + ((t592 * t495 + t590 * t496 + t591 * t514) * t516 + (t597 * t495 + t593 * t496 + t595 * t514) * t499 + (t598 * t495 + t594 * t496 + t596 * t514) * t498) * t499 / 0.2e1 + ((t592 * t510 + t590 * t511 + t591 * t524) * t516 + (t597 * t510 + t593 * t511 + t595 * t524) * t499 + (t598 * t510 + t594 * t511 + t596 * t524) * t498) * t516 / 0.2e1 + (Icges(2,3) + (Icges(3,5) * t547 + (Icges(3,1) * t545 + Icges(3,4) * t580) * t546) * t578 + t546 * t580 * (Icges(3,6) * t547 + (Icges(3,4) * t545 + Icges(3,2) * t580) * t546) + t547 * (Icges(3,3) * t547 + (Icges(3,5) * t545 + Icges(3,6) * t580) * t546) + m(2) * (t540 ^ 2 + t541 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
