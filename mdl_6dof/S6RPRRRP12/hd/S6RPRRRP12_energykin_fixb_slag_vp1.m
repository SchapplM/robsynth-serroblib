% Calculate kinetic energy for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP12_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP12_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:43:58
% EndTime: 2019-03-09 06:44:01
% DurationCPUTime: 2.58s
% Computational Cost: add. (5044->292), mult. (13829->449), div. (0->0), fcn. (17994->14), ass. (0->143)
t630 = Icges(6,1) + Icges(7,1);
t629 = -Icges(6,4) + Icges(7,5);
t628 = Icges(7,4) + Icges(6,5);
t627 = Icges(6,2) + Icges(7,3);
t626 = Icges(7,2) + Icges(6,3);
t625 = -Icges(6,6) + Icges(7,6);
t624 = rSges(7,1) + pkin(5);
t623 = rSges(7,3) + qJ(6);
t568 = cos(pkin(6));
t603 = cos(pkin(12));
t605 = sin(qJ(1));
t589 = t605 * t603;
t566 = sin(pkin(12));
t572 = cos(qJ(1));
t599 = t572 * t566;
t555 = t568 * t599 + t589;
t571 = sin(qJ(3));
t591 = t572 * t603;
t594 = t605 * t566;
t581 = -t568 * t591 + t594;
t604 = cos(pkin(7));
t576 = t581 * t604;
t567 = sin(pkin(6));
t602 = sin(pkin(7));
t592 = t567 * t602;
t608 = cos(qJ(3));
t533 = t555 * t608 + (-t572 * t592 - t576) * t571;
t593 = t567 * t604;
t547 = -t572 * t593 + t581 * t602;
t570 = sin(qJ(4));
t607 = cos(qJ(4));
t514 = t533 * t607 + t547 * t570;
t590 = t608 * t602;
t600 = t567 * t572;
t532 = t555 * t571 + t576 * t608 + t590 * t600;
t569 = sin(qJ(5));
t606 = cos(qJ(5));
t488 = t514 * t569 - t532 * t606;
t489 = t514 * t606 + t532 * t569;
t513 = t533 * t570 - t547 * t607;
t622 = t627 * t488 + t629 * t489 + t625 * t513;
t556 = -t568 * t594 + t591;
t582 = t568 * t589 + t599;
t613 = t582 * t604 - t605 * t592;
t535 = t556 * t608 - t571 * t613;
t548 = t582 * t602 + t593 * t605;
t516 = t535 * t607 + t548 * t570;
t534 = t556 * t571 + t608 * t613;
t490 = t516 * t569 - t534 * t606;
t491 = t516 * t606 + t534 * t569;
t515 = t535 * t570 - t548 * t607;
t621 = t627 * t490 + t629 * t491 + t625 * t515;
t620 = t625 * t488 + t628 * t489 + t626 * t513;
t619 = t625 * t490 + t628 * t491 + t626 * t515;
t618 = t629 * t488 + t630 * t489 + t628 * t513;
t617 = t629 * t490 + t630 * t491 + t628 * t515;
t587 = t604 * t603;
t546 = t568 * t602 * t571 + (t566 * t608 + t571 * t587) * t567;
t554 = t568 * t604 - t592 * t603;
t531 = t546 * t607 + t554 * t570;
t601 = t566 * t567;
t545 = -t567 * t587 * t608 - t568 * t590 + t571 * t601;
t511 = t531 * t569 - t545 * t606;
t512 = t531 * t606 + t545 * t569;
t530 = t546 * t570 - t554 * t607;
t616 = t627 * t511 + t629 * t512 + t625 * t530;
t615 = t625 * t511 + t628 * t512 + t626 * t530;
t614 = t629 * t511 + t630 * t512 + t628 * t530;
t598 = rSges(7,2) * t513 + t623 * t488 + t624 * t489;
t597 = rSges(7,2) * t515 + t623 * t490 + t624 * t491;
t596 = rSges(7,2) * t530 + t623 * t511 + t624 * t512;
t543 = qJD(3) * t547;
t517 = qJD(4) * t532 + t543;
t544 = qJD(3) * t548;
t518 = qJD(4) * t534 + t544;
t551 = qJD(3) * t554 + qJD(1);
t595 = t567 * t605;
t536 = qJD(4) * t545 + t551;
t588 = -qJD(2) * t600 + qJD(1) * (t572 * pkin(1) + qJ(2) * t595);
t558 = pkin(1) * t605 - qJ(2) * t600;
t564 = qJD(2) * t595;
t586 = t564 + (-t555 * pkin(2) - pkin(9) * t547 - t558) * qJD(1);
t584 = qJD(1) * (t556 * pkin(2) + pkin(9) * t548) + t588;
t507 = pkin(3) * t533 + pkin(10) * t532;
t508 = pkin(3) * t535 + pkin(10) * t534;
t565 = qJD(2) * t568;
t583 = t507 * t544 - t508 * t543 + t565;
t525 = pkin(3) * t546 + pkin(10) * t545;
t580 = -t551 * t507 + t525 * t543 + t586;
t484 = pkin(4) * t514 + pkin(11) * t513;
t485 = pkin(4) * t516 + pkin(11) * t515;
t579 = t518 * t484 - t517 * t485 + t583;
t578 = t551 * t508 - t525 * t544 + t584;
t506 = pkin(4) * t531 + pkin(11) * t530;
t575 = -t536 * t484 + t517 * t506 + t580;
t574 = t536 * t485 - t506 * t518 + t578;
t562 = t572 * rSges(2,1) - rSges(2,2) * t605;
t561 = rSges(2,1) * t605 + t572 * rSges(2,2);
t524 = qJD(1) * (t556 * rSges(3,1) - rSges(3,2) * t582 + rSges(3,3) * t595) + t588;
t523 = t564 + (-t555 * rSges(3,1) + rSges(3,2) * t581 + rSges(3,3) * t600 - t558) * qJD(1);
t522 = rSges(4,1) * t546 - rSges(4,2) * t545 + rSges(4,3) * t554;
t521 = Icges(4,1) * t546 - Icges(4,4) * t545 + Icges(4,5) * t554;
t520 = Icges(4,4) * t546 - Icges(4,2) * t545 + Icges(4,6) * t554;
t519 = Icges(4,5) * t546 - Icges(4,6) * t545 + Icges(4,3) * t554;
t509 = qJD(5) * t530 + t536;
t503 = rSges(4,1) * t535 - rSges(4,2) * t534 + rSges(4,3) * t548;
t502 = rSges(4,1) * t533 - rSges(4,2) * t532 + rSges(4,3) * t547;
t501 = Icges(4,1) * t535 - Icges(4,4) * t534 + Icges(4,5) * t548;
t500 = Icges(4,1) * t533 - Icges(4,4) * t532 + Icges(4,5) * t547;
t499 = Icges(4,4) * t535 - Icges(4,2) * t534 + Icges(4,6) * t548;
t498 = Icges(4,4) * t533 - Icges(4,2) * t532 + Icges(4,6) * t547;
t497 = Icges(4,5) * t535 - Icges(4,6) * t534 + Icges(4,3) * t548;
t496 = Icges(4,5) * t533 - Icges(4,6) * t532 + Icges(4,3) * t547;
t495 = rSges(5,1) * t531 - rSges(5,2) * t530 + rSges(5,3) * t545;
t494 = Icges(5,1) * t531 - Icges(5,4) * t530 + Icges(5,5) * t545;
t493 = Icges(5,4) * t531 - Icges(5,2) * t530 + Icges(5,6) * t545;
t492 = Icges(5,5) * t531 - Icges(5,6) * t530 + Icges(5,3) * t545;
t487 = qJD(5) * t515 + t518;
t486 = qJD(5) * t513 + t517;
t480 = rSges(5,1) * t516 - rSges(5,2) * t515 + rSges(5,3) * t534;
t479 = rSges(5,1) * t514 - rSges(5,2) * t513 + rSges(5,3) * t532;
t478 = Icges(5,1) * t516 - Icges(5,4) * t515 + Icges(5,5) * t534;
t477 = Icges(5,1) * t514 - Icges(5,4) * t513 + Icges(5,5) * t532;
t476 = Icges(5,4) * t516 - Icges(5,2) * t515 + Icges(5,6) * t534;
t475 = Icges(5,4) * t514 - Icges(5,2) * t513 + Icges(5,6) * t532;
t474 = Icges(5,5) * t516 - Icges(5,6) * t515 + Icges(5,3) * t534;
t473 = Icges(5,5) * t514 - Icges(5,6) * t513 + Icges(5,3) * t532;
t472 = rSges(6,1) * t512 - rSges(6,2) * t511 + rSges(6,3) * t530;
t461 = t503 * t551 - t522 * t544 + t584;
t460 = -t551 * t502 + t522 * t543 + t586;
t459 = t565 + (t502 * t548 - t503 * t547) * qJD(3);
t458 = rSges(6,1) * t491 - rSges(6,2) * t490 + rSges(6,3) * t515;
t456 = rSges(6,1) * t489 - rSges(6,2) * t488 + rSges(6,3) * t513;
t442 = t480 * t536 - t495 * t518 + t578;
t441 = -t536 * t479 + t517 * t495 + t580;
t440 = t479 * t518 - t480 * t517 + t583;
t439 = t458 * t509 - t472 * t487 + t574;
t438 = -t509 * t456 + t486 * t472 + t575;
t437 = t456 * t487 - t458 * t486 + t579;
t436 = qJD(6) * t488 - t487 * t596 + t509 * t597 + t574;
t435 = qJD(6) * t490 + t486 * t596 - t509 * t598 + t575;
t434 = qJD(6) * t511 - t486 * t597 + t487 * t598 + t579;
t1 = ((t519 * t548 - t520 * t534 + t521 * t535) * t551 + ((t497 * t548 - t499 * t534 + t501 * t535) * t548 + (t496 * t548 - t498 * t534 + t500 * t535) * t547) * qJD(3)) * t544 / 0.2e1 + t536 * ((t474 * t545 - t476 * t530 + t478 * t531) * t518 + (t473 * t545 - t475 * t530 + t477 * t531) * t517 + (t545 * t492 - t530 * t493 + t531 * t494) * t536) / 0.2e1 + t518 * ((t534 * t474 - t515 * t476 + t516 * t478) * t518 + (t473 * t534 - t475 * t515 + t477 * t516) * t517 + (t492 * t534 - t493 * t515 + t494 * t516) * t536) / 0.2e1 + t517 * ((t474 * t532 - t476 * t513 + t478 * t514) * t518 + (t532 * t473 - t513 * t475 + t514 * t477) * t517 + (t492 * t532 - t493 * t513 + t494 * t514) * t536) / 0.2e1 + t551 * ((t554 * t519 - t545 * t520 + t546 * t521) * t551 + ((t497 * t554 - t499 * t545 + t501 * t546) * t548 + (t496 * t554 - t498 * t545 + t500 * t546) * t547) * qJD(3)) / 0.2e1 + ((t519 * t547 - t520 * t532 + t521 * t533) * t551 + ((t497 * t547 - t499 * t532 + t501 * t533) * t548 + (t496 * t547 - t498 * t532 + t500 * t533) * t547) * qJD(3)) * t543 / 0.2e1 + m(7) * (t434 ^ 2 + t435 ^ 2 + t436 ^ 2) / 0.2e1 + m(6) * (t437 ^ 2 + t438 ^ 2 + t439 ^ 2) / 0.2e1 + m(5) * (t440 ^ 2 + t441 ^ 2 + t442 ^ 2) / 0.2e1 + m(4) * (t459 ^ 2 + t460 ^ 2 + t461 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t568 ^ 2 + t523 ^ 2 + t524 ^ 2) / 0.2e1 + ((t488 * t616 + t489 * t614 + t513 * t615) * t509 + (t488 * t621 + t489 * t617 + t513 * t619) * t487 + (t622 * t488 + t618 * t489 + t620 * t513) * t486) * t486 / 0.2e1 + ((t490 * t616 + t491 * t614 + t515 * t615) * t509 + (t621 * t490 + t617 * t491 + t619 * t515) * t487 + (t490 * t622 + t618 * t491 + t620 * t515) * t486) * t487 / 0.2e1 + ((t616 * t511 + t614 * t512 + t615 * t530) * t509 + (t511 * t621 + t512 * t617 + t530 * t619) * t487 + (t511 * t622 + t618 * t512 + t620 * t530) * t486) * t509 / 0.2e1 + (Icges(2,3) + (Icges(3,5) * t568 + (Icges(3,1) * t566 + Icges(3,4) * t603) * t567) * t601 + t567 * t603 * (Icges(3,6) * t568 + (Icges(3,4) * t566 + Icges(3,2) * t603) * t567) + t568 * (Icges(3,3) * t568 + (Icges(3,5) * t566 + Icges(3,6) * t603) * t567) + m(2) * (t561 ^ 2 + t562 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
