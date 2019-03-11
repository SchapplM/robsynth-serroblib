% Calculate kinetic energy for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:42:00
% EndTime: 2019-03-08 18:42:03
% DurationCPUTime: 2.62s
% Computational Cost: add. (5374->324), mult. (14573->500), div. (0->0), fcn. (19076->16), ass. (0->153)
t611 = Icges(4,3) + Icges(5,3);
t576 = sin(qJ(3));
t578 = cos(qJ(3));
t605 = sin(pkin(13));
t606 = cos(pkin(13));
t560 = -t576 * t605 + t578 * t606;
t572 = cos(pkin(7));
t548 = t560 * t572;
t566 = sin(pkin(12));
t567 = sin(pkin(11));
t570 = cos(pkin(12));
t571 = cos(pkin(11));
t573 = cos(pkin(6));
t598 = t571 * t573;
t555 = -t566 * t567 + t570 * t598;
t556 = t566 * t598 + t567 * t570;
t559 = -t576 * t606 - t578 * t605;
t569 = sin(pkin(6));
t568 = sin(pkin(7));
t583 = t568 * t560;
t581 = t569 * t583;
t516 = t555 * t548 + t556 * t559 - t571 * t581;
t547 = t559 * t568;
t549 = t559 * t572;
t600 = t569 * t571;
t517 = t547 * t600 - t549 * t555 + t556 * t560;
t603 = t567 * t573;
t557 = -t566 * t571 - t570 * t603;
t558 = -t566 * t603 + t570 * t571;
t518 = t557 * t548 + t558 * t559 + t567 * t581;
t604 = t567 * t569;
t519 = -t547 * t604 - t549 * t557 + t558 * t560;
t601 = t569 * t570;
t529 = t569 * t566 * t559 + t548 * t601 + t573 * t583;
t530 = -t573 * t547 + (-t549 * t570 + t560 * t566) * t569;
t586 = t555 * t572 - t568 * t600;
t531 = -t556 * t576 + t578 * t586;
t532 = t556 * t578 + t576 * t586;
t585 = t557 * t572 + t568 * t604;
t533 = -t558 * t576 + t578 * t585;
t534 = t558 * t578 + t576 * t585;
t541 = t573 * t568 * t578 + (t570 * t572 * t578 - t566 * t576) * t569;
t597 = t572 * t576;
t602 = t568 * t576;
t542 = t573 * t602 + (t566 * t578 + t570 * t597) * t569;
t599 = t569 * t572;
t544 = -t557 * t568 + t567 * t599;
t554 = -t568 * t601 + t572 * t573;
t587 = t555 * t568 + t571 * t599;
t610 = -(Icges(4,5) * t532 + Icges(5,5) * t517 + Icges(4,6) * t531 + Icges(5,6) * t516 - t611 * t587) * t587 + (Icges(4,5) * t542 + Icges(5,5) * t530 + Icges(4,6) * t541 + Icges(5,6) * t529 + t611 * t554) * t554 + (Icges(4,5) * t534 + Icges(5,5) * t519 + Icges(4,6) * t533 + Icges(5,6) * t518 + t611 * t544) * t544;
t609 = cos(qJ(5));
t608 = pkin(3) * t578;
t607 = pkin(8) + qJ(4);
t539 = qJD(3) * t587;
t511 = -qJD(5) * t516 - t539;
t540 = qJD(3) * t544;
t512 = -qJD(5) * t518 + t540;
t551 = qJD(3) * t554;
t522 = -qJD(5) * t529 + t551;
t596 = qJD(2) * t569;
t562 = qJD(2) * t573 + qJD(1);
t552 = pkin(3) * t602 + t572 * t607;
t553 = pkin(3) * t597 - t568 * t607;
t527 = (-pkin(8) * t572 + t552) * t573 + ((pkin(8) * t568 + t553) * t570 + t608 * t566) * t569;
t561 = t567 * t596;
t594 = qJD(4) * t544 - t527 * t539 + t561;
t593 = t571 * t596;
t505 = pkin(8) * t587 - t552 * t600 + t555 * t553 + t556 * t608;
t589 = qJD(4) * t554 + t505 * t540 + t562;
t506 = -pkin(8) * t544 + t552 * t604 + t557 * t553 + t558 * t608;
t588 = -qJD(4) * t587 + t506 * t551 - t593;
t483 = pkin(4) * t517 - pkin(9) * t516;
t504 = pkin(4) * t530 - pkin(9) * t529;
t584 = -t504 * t539 + (-t483 - t505) * t551 + t594;
t484 = pkin(4) * t519 - pkin(9) * t518;
t582 = t483 * t540 - (-t484 - t506) * t539 + t589;
t580 = t484 * t551 + (-t504 - t527) * t540 + t588;
t577 = cos(qJ(6));
t575 = sin(qJ(5));
t574 = sin(qJ(6));
t526 = rSges(4,1) * t542 + rSges(4,2) * t541 + rSges(4,3) * t554;
t525 = Icges(4,1) * t542 + Icges(4,4) * t541 + Icges(4,5) * t554;
t524 = Icges(4,4) * t542 + Icges(4,2) * t541 + Icges(4,6) * t554;
t521 = t530 * t609 + t554 * t575;
t520 = t530 * t575 - t554 * t609;
t510 = t519 * t609 + t544 * t575;
t509 = t519 * t575 - t544 * t609;
t508 = t517 * t609 - t575 * t587;
t507 = t517 * t575 + t587 * t609;
t503 = rSges(4,1) * t534 + rSges(4,2) * t533 + rSges(4,3) * t544;
t502 = rSges(4,1) * t532 + rSges(4,2) * t531 - rSges(4,3) * t587;
t501 = Icges(4,1) * t534 + Icges(4,4) * t533 + Icges(4,5) * t544;
t500 = Icges(4,1) * t532 + Icges(4,4) * t531 - Icges(4,5) * t587;
t499 = Icges(4,4) * t534 + Icges(4,2) * t533 + Icges(4,6) * t544;
t498 = Icges(4,4) * t532 + Icges(4,2) * t531 - Icges(4,6) * t587;
t494 = rSges(5,1) * t530 + rSges(5,2) * t529 + rSges(5,3) * t554;
t493 = Icges(5,1) * t530 + Icges(5,4) * t529 + Icges(5,5) * t554;
t492 = Icges(5,4) * t530 + Icges(5,2) * t529 + Icges(5,6) * t554;
t489 = t521 * t577 - t529 * t574;
t488 = -t521 * t574 - t529 * t577;
t486 = qJD(6) * t520 + t522;
t485 = pkin(5) * t521 + pkin(10) * t520;
t480 = rSges(5,1) * t519 + rSges(5,2) * t518 + rSges(5,3) * t544;
t479 = rSges(5,1) * t517 + rSges(5,2) * t516 - rSges(5,3) * t587;
t478 = Icges(5,1) * t519 + Icges(5,4) * t518 + Icges(5,5) * t544;
t477 = Icges(5,1) * t517 + Icges(5,4) * t516 - Icges(5,5) * t587;
t476 = Icges(5,4) * t519 + Icges(5,2) * t518 + Icges(5,6) * t544;
t475 = Icges(5,4) * t517 + Icges(5,2) * t516 - Icges(5,6) * t587;
t472 = t510 * t577 - t518 * t574;
t471 = -t510 * t574 - t518 * t577;
t470 = t508 * t577 - t516 * t574;
t469 = -t508 * t574 - t516 * t577;
t468 = rSges(6,1) * t521 - rSges(6,2) * t520 - rSges(6,3) * t529;
t467 = Icges(6,1) * t521 - Icges(6,4) * t520 - Icges(6,5) * t529;
t466 = Icges(6,4) * t521 - Icges(6,2) * t520 - Icges(6,6) * t529;
t465 = Icges(6,5) * t521 - Icges(6,6) * t520 - Icges(6,3) * t529;
t464 = qJD(6) * t509 + t512;
t463 = qJD(6) * t507 + t511;
t462 = pkin(5) * t510 + pkin(10) * t509;
t461 = pkin(5) * t508 + pkin(10) * t507;
t460 = -t593 + (t503 * t554 - t526 * t544) * qJD(3);
t459 = t561 + (-t502 * t554 - t526 * t587) * qJD(3);
t458 = rSges(6,1) * t510 - rSges(6,2) * t509 - rSges(6,3) * t518;
t457 = rSges(6,1) * t508 - rSges(6,2) * t507 - rSges(6,3) * t516;
t456 = Icges(6,1) * t510 - Icges(6,4) * t509 - Icges(6,5) * t518;
t455 = Icges(6,1) * t508 - Icges(6,4) * t507 - Icges(6,5) * t516;
t454 = Icges(6,4) * t510 - Icges(6,2) * t509 - Icges(6,6) * t518;
t453 = Icges(6,4) * t508 - Icges(6,2) * t507 - Icges(6,6) * t516;
t452 = Icges(6,5) * t510 - Icges(6,6) * t509 - Icges(6,3) * t518;
t451 = Icges(6,5) * t508 - Icges(6,6) * t507 - Icges(6,3) * t516;
t450 = (t502 * t544 + t503 * t587) * qJD(3) + t562;
t449 = rSges(7,1) * t489 + rSges(7,2) * t488 + rSges(7,3) * t520;
t448 = Icges(7,1) * t489 + Icges(7,4) * t488 + Icges(7,5) * t520;
t447 = Icges(7,4) * t489 + Icges(7,2) * t488 + Icges(7,6) * t520;
t446 = Icges(7,5) * t489 + Icges(7,6) * t488 + Icges(7,3) * t520;
t445 = rSges(7,1) * t472 + rSges(7,2) * t471 + rSges(7,3) * t509;
t444 = rSges(7,1) * t470 + rSges(7,2) * t469 + rSges(7,3) * t507;
t443 = Icges(7,1) * t472 + Icges(7,4) * t471 + Icges(7,5) * t509;
t442 = Icges(7,1) * t470 + Icges(7,4) * t469 + Icges(7,5) * t507;
t441 = Icges(7,4) * t472 + Icges(7,2) * t471 + Icges(7,6) * t509;
t440 = Icges(7,4) * t470 + Icges(7,2) * t469 + Icges(7,6) * t507;
t439 = Icges(7,5) * t472 + Icges(7,6) * t471 + Icges(7,3) * t509;
t438 = Icges(7,5) * t470 + Icges(7,6) * t469 + Icges(7,3) * t507;
t437 = (t480 * t554 + (-t494 - t527) * t544) * qJD(3) + t588;
t436 = (-t494 * t587 + (-t479 - t505) * t554) * qJD(3) + t594;
t435 = (t479 * t544 - (-t480 - t506) * t587) * qJD(3) + t589;
t434 = t522 * t458 - t512 * t468 + t580;
t433 = -t522 * t457 + t511 * t468 + t584;
t432 = t512 * t457 - t511 * t458 + t582;
t431 = t486 * t445 - t464 * t449 + t522 * t462 - t512 * t485 + t580;
t430 = -t486 * t444 + t463 * t449 - t522 * t461 + t511 * t485 + t584;
t429 = t464 * t444 - t463 * t445 + t512 * t461 - t511 * t462 + t582;
t1 = t512 * ((-t518 * t452 - t509 * t454 + t510 * t456) * t512 + (-t451 * t518 - t453 * t509 + t455 * t510) * t511 + (-t465 * t518 - t466 * t509 + t467 * t510) * t522) / 0.2e1 + t511 * ((-t452 * t516 - t454 * t507 + t456 * t508) * t512 + (-t516 * t451 - t507 * t453 + t508 * t455) * t511 + (-t465 * t516 - t466 * t507 + t467 * t508) * t522) / 0.2e1 + t522 * ((-t452 * t529 - t454 * t520 + t456 * t521) * t512 + (-t451 * t529 - t453 * t520 + t455 * t521) * t511 + (-t529 * t465 - t520 * t466 + t521 * t467) * t522) / 0.2e1 + t464 * ((t509 * t439 + t471 * t441 + t472 * t443) * t464 + (t438 * t509 + t440 * t471 + t442 * t472) * t463 + (t446 * t509 + t447 * t471 + t448 * t472) * t486) / 0.2e1 + t463 * ((t439 * t507 + t441 * t469 + t443 * t470) * t464 + (t507 * t438 + t469 * t440 + t470 * t442) * t463 + (t446 * t507 + t447 * t469 + t448 * t470) * t486) / 0.2e1 + t486 * ((t439 * t520 + t441 * t488 + t443 * t489) * t464 + (t438 * t520 + t440 * t488 + t442 * t489) * t463 + (t520 * t446 + t488 * t447 + t489 * t448) * t486) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t562 ^ 2 + (t567 ^ 2 + t571 ^ 2) * qJD(2) ^ 2 * t569 ^ 2) / 0.2e1 + m(4) * (t450 ^ 2 + t459 ^ 2 + t460 ^ 2) / 0.2e1 + m(5) * (t435 ^ 2 + t436 ^ 2 + t437 ^ 2) / 0.2e1 + m(6) * (t432 ^ 2 + t433 ^ 2 + t434 ^ 2) / 0.2e1 + m(7) * (t429 ^ 2 + t430 ^ 2 + t431 ^ 2) / 0.2e1 + (-((t492 * t516 + t493 * t517 + t524 * t531 + t525 * t532) * t554 + (t476 * t516 + t478 * t517 + t499 * t531 + t501 * t532) * t544 - (t516 * t475 + t517 * t477 + t531 * t498 + t532 * t500 + t610) * t587) * t587 + ((t492 * t518 + t493 * t519 + t524 * t533 + t525 * t534) * t554 - (t475 * t518 + t477 * t519 + t498 * t533 + t500 * t534) * t587 + (t518 * t476 + t519 * t478 + t533 * t499 + t534 * t501 + t610) * t544) * t544 + ((t476 * t529 + t478 * t530 + t499 * t541 + t501 * t542) * t544 - (t475 * t529 + t477 * t530 + t498 * t541 + t500 * t542) * t587 + (t529 * t492 + t530 * t493 + t541 * t524 + t542 * t525 + t610) * t554) * t554) * qJD(3) ^ 2 / 0.2e1;
T  = t1;
