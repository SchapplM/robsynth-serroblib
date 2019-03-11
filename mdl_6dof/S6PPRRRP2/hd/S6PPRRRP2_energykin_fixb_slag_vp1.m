% Calculate kinetic energy for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:56:00
% EndTime: 2019-03-08 18:56:02
% DurationCPUTime: 2.14s
% Computational Cost: add. (4930->265), mult. (13631->409), div. (0->0), fcn. (17804->14), ass. (0->133)
t615 = Icges(6,1) + Icges(7,1);
t614 = -Icges(6,4) + Icges(7,5);
t613 = Icges(7,4) + Icges(6,5);
t612 = Icges(6,2) + Icges(7,3);
t611 = Icges(7,2) + Icges(6,3);
t610 = -Icges(6,6) + Icges(7,6);
t609 = rSges(7,1) + pkin(5);
t608 = rSges(7,3) + qJ(6);
t558 = sin(pkin(11));
t560 = cos(pkin(11));
t588 = sin(pkin(12));
t592 = cos(pkin(6));
t577 = t592 * t588;
t590 = cos(pkin(12));
t549 = t558 * t590 + t560 * t577;
t563 = sin(qJ(3));
t579 = t592 * t590;
t570 = t558 * t588 - t560 * t579;
t591 = cos(pkin(7));
t566 = t570 * t591;
t559 = sin(pkin(6));
t589 = sin(pkin(7));
t580 = t559 * t589;
t595 = cos(qJ(3));
t528 = t549 * t595 + (-t560 * t580 - t566) * t563;
t581 = t559 * t591;
t542 = -t560 * t581 + t570 * t589;
t562 = sin(qJ(4));
t594 = cos(qJ(4));
t511 = t528 * t594 + t542 * t562;
t575 = t595 * t580;
t527 = t549 * t563 + t560 * t575 + t566 * t595;
t561 = sin(qJ(5));
t593 = cos(qJ(5));
t487 = t511 * t561 - t527 * t593;
t488 = t511 * t593 + t527 * t561;
t510 = t528 * t562 - t542 * t594;
t607 = t487 * t612 + t488 * t614 + t510 * t610;
t550 = -t558 * t577 + t560 * t590;
t569 = t558 * t579 + t560 * t588;
t565 = t569 * t591;
t530 = t550 * t595 + (t558 * t580 - t565) * t563;
t543 = t558 * t581 + t569 * t589;
t513 = t530 * t594 + t543 * t562;
t529 = t550 * t563 - t558 * t575 + t565 * t595;
t489 = t513 * t561 - t529 * t593;
t490 = t513 * t593 + t529 * t561;
t512 = t530 * t562 - t543 * t594;
t606 = t489 * t612 + t490 * t614 + t512 * t610;
t605 = t487 * t610 + t488 * t613 + t510 * t611;
t604 = t489 * t610 + t490 * t613 + t512 * t611;
t603 = t614 * t487 + t488 * t615 + t613 * t510;
t602 = t614 * t489 + t490 * t615 + t613 * t512;
t576 = t591 * t590;
t578 = t592 * t589;
t541 = t563 * t578 + (t563 * t576 + t588 * t595) * t559;
t548 = -t580 * t590 + t591 * t592;
t532 = t541 * t594 + t548 * t562;
t540 = -t578 * t595 + (t563 * t588 - t576 * t595) * t559;
t514 = t532 * t561 - t540 * t593;
t515 = t532 * t593 + t540 * t561;
t531 = t541 * t562 - t548 * t594;
t601 = t514 * t612 + t515 * t614 + t531 * t610;
t600 = t514 * t610 + t515 * t613 + t531 * t611;
t599 = t614 * t514 + t515 * t615 + t613 * t531;
t587 = rSges(7,2) * t510 + t608 * t487 + t609 * t488;
t586 = rSges(7,2) * t512 + t608 * t489 + t609 * t490;
t585 = rSges(7,2) * t531 + t608 * t514 + t609 * t515;
t538 = qJD(3) * t542;
t516 = qJD(4) * t527 + t538;
t539 = qJD(3) * t543;
t517 = qJD(4) * t529 + t539;
t547 = qJD(3) * t548;
t533 = qJD(4) * t540 + t547;
t584 = qJD(2) * t559;
t555 = qJD(2) * t592 + qJD(1);
t582 = t560 * t584;
t505 = pkin(3) * t528 + pkin(9) * t527;
t522 = pkin(3) * t541 + pkin(9) * t540;
t554 = t558 * t584;
t574 = -t505 * t547 + t522 * t538 + t554;
t506 = pkin(3) * t530 + pkin(9) * t529;
t573 = t505 * t539 - t506 * t538 + t555;
t572 = t506 * t547 - t522 * t539 - t582;
t482 = pkin(4) * t511 + pkin(10) * t510;
t507 = pkin(4) * t532 + pkin(10) * t531;
t571 = -t482 * t533 + t516 * t507 + t574;
t483 = pkin(4) * t513 + pkin(10) * t512;
t568 = t517 * t482 - t483 * t516 + t573;
t567 = t533 * t483 - t507 * t517 + t572;
t521 = rSges(4,1) * t541 - rSges(4,2) * t540 + rSges(4,3) * t548;
t520 = Icges(4,1) * t541 - Icges(4,4) * t540 + Icges(4,5) * t548;
t519 = Icges(4,4) * t541 - Icges(4,2) * t540 + Icges(4,6) * t548;
t518 = Icges(4,5) * t541 - Icges(4,6) * t540 + Icges(4,3) * t548;
t508 = qJD(5) * t531 + t533;
t503 = rSges(5,1) * t532 - rSges(5,2) * t531 + rSges(5,3) * t540;
t502 = Icges(5,1) * t532 - Icges(5,4) * t531 + Icges(5,5) * t540;
t501 = Icges(5,4) * t532 - Icges(5,2) * t531 + Icges(5,6) * t540;
t500 = Icges(5,5) * t532 - Icges(5,6) * t531 + Icges(5,3) * t540;
t498 = rSges(4,1) * t530 - rSges(4,2) * t529 + rSges(4,3) * t543;
t497 = rSges(4,1) * t528 - rSges(4,2) * t527 + rSges(4,3) * t542;
t496 = Icges(4,1) * t530 - Icges(4,4) * t529 + Icges(4,5) * t543;
t495 = Icges(4,1) * t528 - Icges(4,4) * t527 + Icges(4,5) * t542;
t494 = Icges(4,4) * t530 - Icges(4,2) * t529 + Icges(4,6) * t543;
t493 = Icges(4,4) * t528 - Icges(4,2) * t527 + Icges(4,6) * t542;
t492 = Icges(4,5) * t530 - Icges(4,6) * t529 + Icges(4,3) * t543;
t491 = Icges(4,5) * t528 - Icges(4,6) * t527 + Icges(4,3) * t542;
t486 = qJD(5) * t512 + t517;
t485 = qJD(5) * t510 + t516;
t479 = rSges(6,1) * t515 - rSges(6,2) * t514 + rSges(6,3) * t531;
t471 = rSges(5,1) * t513 - rSges(5,2) * t512 + rSges(5,3) * t529;
t470 = rSges(5,1) * t511 - rSges(5,2) * t510 + rSges(5,3) * t527;
t469 = Icges(5,1) * t513 - Icges(5,4) * t512 + Icges(5,5) * t529;
t468 = Icges(5,1) * t511 - Icges(5,4) * t510 + Icges(5,5) * t527;
t467 = Icges(5,4) * t513 - Icges(5,2) * t512 + Icges(5,6) * t529;
t466 = Icges(5,4) * t511 - Icges(5,2) * t510 + Icges(5,6) * t527;
t465 = Icges(5,5) * t513 - Icges(5,6) * t512 + Icges(5,3) * t529;
t464 = Icges(5,5) * t511 - Icges(5,6) * t510 + Icges(5,3) * t527;
t462 = -t582 + (t498 * t548 - t521 * t543) * qJD(3);
t461 = t554 + (-t497 * t548 + t521 * t542) * qJD(3);
t458 = (t497 * t543 - t498 * t542) * qJD(3) + t555;
t457 = rSges(6,1) * t490 - rSges(6,2) * t489 + rSges(6,3) * t512;
t455 = rSges(6,1) * t488 - rSges(6,2) * t487 + rSges(6,3) * t510;
t441 = t471 * t533 - t503 * t517 + t572;
t440 = -t470 * t533 + t503 * t516 + t574;
t439 = t470 * t517 - t471 * t516 + t573;
t438 = t457 * t508 - t479 * t486 + t567;
t437 = -t455 * t508 + t479 * t485 + t571;
t436 = t455 * t486 - t457 * t485 + t568;
t435 = qJD(6) * t487 - t486 * t585 + t508 * t586 + t567;
t434 = qJD(6) * t489 + t485 * t585 - t508 * t587 + t571;
t433 = qJD(6) * t514 - t485 * t586 + t486 * t587 + t568;
t1 = m(7) * (t433 ^ 2 + t434 ^ 2 + t435 ^ 2) / 0.2e1 + m(6) * (t436 ^ 2 + t437 ^ 2 + t438 ^ 2) / 0.2e1 + m(5) * (t439 ^ 2 + t440 ^ 2 + t441 ^ 2) / 0.2e1 + m(4) * (t458 ^ 2 + t461 ^ 2 + t462 ^ 2) / 0.2e1 + m(3) * (t555 ^ 2 + (t558 ^ 2 + t560 ^ 2) * qJD(2) ^ 2 * t559 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t516 * ((t465 * t527 - t467 * t510 + t469 * t511) * t517 + (t464 * t527 - t466 * t510 + t468 * t511) * t516 + (t500 * t527 - t501 * t510 + t502 * t511) * t533) / 0.2e1 + t533 * ((t465 * t540 - t467 * t531 + t469 * t532) * t517 + (t464 * t540 - t466 * t531 + t468 * t532) * t516 + (t500 * t540 - t501 * t531 + t502 * t532) * t533) / 0.2e1 + t517 * ((t465 * t529 - t467 * t512 + t469 * t513) * t517 + (t464 * t529 - t466 * t512 + t468 * t513) * t516 + (t500 * t529 - t501 * t512 + t502 * t513) * t533) / 0.2e1 + ((t487 * t601 + t488 * t599 + t510 * t600) * t508 + (t487 * t606 + t488 * t602 + t510 * t604) * t486 + (t607 * t487 + t603 * t488 + t605 * t510) * t485) * t485 / 0.2e1 + ((t489 * t601 + t490 * t599 + t512 * t600) * t508 + (t606 * t489 + t602 * t490 + t604 * t512) * t486 + (t489 * t607 + t603 * t490 + t605 * t512) * t485) * t486 / 0.2e1 + ((t601 * t514 + t599 * t515 + t600 * t531) * t508 + (t514 * t606 + t515 * t602 + t531 * t604) * t486 + (t514 * t607 + t603 * t515 + t605 * t531) * t485) * t508 / 0.2e1 + (t543 * ((t492 * t543 - t494 * t529 + t496 * t530) * t543 + (t491 * t543 - t493 * t529 + t495 * t530) * t542 + (t518 * t543 - t519 * t529 + t530 * t520) * t548) + t542 * ((t492 * t542 - t494 * t527 + t496 * t528) * t543 + (t491 * t542 - t493 * t527 + t495 * t528) * t542 + (t518 * t542 - t519 * t527 + t520 * t528) * t548) + t548 * ((t492 * t548 - t494 * t540 + t496 * t541) * t543 + (t491 * t548 - t493 * t540 + t495 * t541) * t542 + (t518 * t548 - t519 * t540 + t520 * t541) * t548)) * qJD(3) ^ 2 / 0.2e1;
T  = t1;
