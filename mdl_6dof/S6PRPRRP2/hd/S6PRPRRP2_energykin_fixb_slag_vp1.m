% Calculate kinetic energy for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:59:55
% EndTime: 2019-03-08 19:59:58
% DurationCPUTime: 3.28s
% Computational Cost: add. (3431->317), mult. (8756->481), div. (0->0), fcn. (11165->12), ass. (0->151)
t611 = Icges(6,1) + Icges(7,1);
t610 = -Icges(6,4) + Icges(7,5);
t609 = Icges(7,4) + Icges(6,5);
t608 = Icges(6,2) + Icges(7,3);
t607 = Icges(7,2) + Icges(6,3);
t606 = -Icges(6,6) + Icges(7,6);
t543 = sin(pkin(11));
t550 = sin(qJ(2));
t551 = cos(qJ(2));
t581 = cos(pkin(11));
t532 = -t551 * t543 - t550 * t581;
t544 = sin(pkin(10));
t546 = cos(pkin(10));
t533 = -t550 * t543 + t551 * t581;
t547 = cos(pkin(6));
t557 = t547 * t533;
t510 = t544 * t532 + t546 * t557;
t526 = t532 * t547;
t511 = -t526 * t546 + t533 * t544;
t545 = sin(pkin(6));
t579 = t545 * t546;
t459 = Icges(4,5) * t511 + Icges(4,6) * t510 - Icges(4,3) * t579;
t512 = t532 * t546 - t544 * t557;
t513 = t526 * t544 + t533 * t546;
t580 = t544 * t545;
t460 = Icges(4,5) * t513 + Icges(4,6) * t512 + Icges(4,3) * t580;
t576 = t547 * t551;
t528 = -t544 * t550 + t546 * t576;
t577 = t547 * t550;
t529 = t544 * t551 + t546 * t577;
t495 = Icges(3,5) * t529 + Icges(3,6) * t528 - Icges(3,3) * t579;
t530 = -t544 * t576 - t546 * t550;
t531 = -t544 * t577 + t546 * t551;
t496 = Icges(3,5) * t531 + Icges(3,6) * t530 + Icges(3,3) * t580;
t524 = t533 * t545;
t525 = t532 * t545;
t589 = (-Icges(4,5) * t525 + Icges(4,6) * t524 + (Icges(3,5) * t550 + Icges(3,6) * t551) * t545 + (Icges(4,3) + Icges(3,3)) * t547) * t547;
t605 = -t589 + (t459 + t495) * t579 - (t496 + t460) * t580;
t604 = rSges(7,1) + pkin(5);
t603 = rSges(7,3) + qJ(6);
t549 = sin(qJ(4));
t578 = t545 * t549;
t584 = cos(qJ(4));
t483 = t511 * t584 - t546 * t578;
t548 = sin(qJ(5));
t583 = cos(qJ(5));
t455 = t483 * t548 + t510 * t583;
t456 = t483 * t583 - t510 * t548;
t567 = t545 * t584;
t482 = t511 * t549 + t546 * t567;
t601 = t455 * t608 + t456 * t610 + t482 * t606;
t485 = t513 * t584 + t544 * t578;
t457 = t485 * t548 + t512 * t583;
t458 = t485 * t583 - t512 * t548;
t484 = t513 * t549 - t544 * t567;
t600 = t457 * t608 + t458 * t610 + t484 * t606;
t599 = t455 * t606 + t456 * t609 + t482 * t607;
t598 = t457 * t606 + t458 * t609 + t484 * t607;
t597 = t610 * t455 + t456 * t611 + t609 * t482;
t596 = t610 * t457 + t458 * t611 + t609 * t484;
t515 = -t525 * t584 + t547 * t549;
t480 = t515 * t548 + t524 * t583;
t481 = t515 * t583 - t524 * t548;
t514 = -t525 * t549 - t547 * t584;
t595 = t480 * t608 + t481 * t610 + t514 * t606;
t594 = t480 * t606 + t481 * t609 + t514 * t607;
t593 = t610 * t480 + t481 * t611 + t609 * t514;
t588 = qJD(2) ^ 2;
t582 = pkin(2) * t551;
t574 = rSges(7,2) * t482 + t603 * t455 + t456 * t604;
t573 = rSges(7,2) * t484 + t603 * t457 + t458 * t604;
t572 = rSges(7,2) * t514 + t603 * t480 + t481 * t604;
t534 = pkin(2) * t545 * t550 + qJ(3) * t547;
t571 = pkin(3) * t525 + pkin(8) * t524 - t534;
t570 = qJD(2) * t545;
t538 = t544 * t570;
t486 = -qJD(4) * t512 + t538;
t542 = qJD(2) * t547;
t516 = -qJD(4) * t524 + t542;
t569 = qJD(3) * t546;
t566 = t546 * t570;
t564 = pkin(2) * t577 - qJ(3) * t545;
t562 = (rSges(4,1) * t525 - rSges(4,2) * t524 - rSges(4,3) * t547 - t534) * t545;
t505 = t544 * t582 + t546 * t564;
t506 = -t544 * t564 + t546 * t582;
t561 = qJD(3) * t547 + t505 * t538 + t506 * t566 + qJD(1);
t487 = -qJD(4) * t510 - t566;
t470 = pkin(3) * t511 - pkin(8) * t510;
t471 = pkin(3) * t513 - pkin(8) * t512;
t558 = t470 * t538 + t471 * t566 + t561;
t450 = pkin(4) * t483 + pkin(9) * t482;
t451 = pkin(4) * t485 + pkin(9) * t484;
t556 = t486 * t450 - t451 * t487 + t558;
t490 = t506 * t542;
t555 = t471 * t542 + t490 + (qJD(2) * t544 * t571 - t569) * t545;
t537 = qJD(3) * t580;
t554 = t537 + ((-t470 - t505) * t547 + t571 * t579) * qJD(2);
t478 = pkin(4) * t515 + pkin(9) * t514;
t553 = t516 * t451 - t478 * t486 + t555;
t552 = -t450 * t516 + t487 * t478 + t554;
t523 = t547 * rSges(3,3) + (rSges(3,1) * t550 + rSges(3,2) * t551) * t545;
t522 = Icges(3,5) * t547 + (Icges(3,1) * t550 + Icges(3,4) * t551) * t545;
t521 = Icges(3,6) * t547 + (Icges(3,4) * t550 + Icges(3,2) * t551) * t545;
t502 = rSges(3,1) * t531 + rSges(3,2) * t530 + rSges(3,3) * t580;
t501 = rSges(3,1) * t529 + rSges(3,2) * t528 - rSges(3,3) * t579;
t500 = Icges(3,1) * t531 + Icges(3,4) * t530 + Icges(3,5) * t580;
t499 = Icges(3,1) * t529 + Icges(3,4) * t528 - Icges(3,5) * t579;
t498 = Icges(3,4) * t531 + Icges(3,2) * t530 + Icges(3,6) * t580;
t497 = Icges(3,4) * t529 + Icges(3,2) * t528 - Icges(3,6) * t579;
t493 = -Icges(4,1) * t525 + Icges(4,4) * t524 + Icges(4,5) * t547;
t492 = -Icges(4,4) * t525 + Icges(4,2) * t524 + Icges(4,6) * t547;
t479 = qJD(5) * t514 + t516;
t477 = (-t501 * t547 - t523 * t579) * qJD(2);
t476 = (t502 * t547 - t523 * t580) * qJD(2);
t475 = rSges(5,1) * t515 - rSges(5,2) * t514 - rSges(5,3) * t524;
t474 = Icges(5,1) * t515 - Icges(5,4) * t514 - Icges(5,5) * t524;
t473 = Icges(5,4) * t515 - Icges(5,2) * t514 - Icges(5,6) * t524;
t472 = Icges(5,5) * t515 - Icges(5,6) * t514 - Icges(5,3) * t524;
t466 = rSges(4,1) * t513 + rSges(4,2) * t512 + rSges(4,3) * t580;
t465 = rSges(4,1) * t511 + rSges(4,2) * t510 - rSges(4,3) * t579;
t464 = Icges(4,1) * t513 + Icges(4,4) * t512 + Icges(4,5) * t580;
t463 = Icges(4,1) * t511 + Icges(4,4) * t510 - Icges(4,5) * t579;
t462 = Icges(4,4) * t513 + Icges(4,2) * t512 + Icges(4,6) * t580;
t461 = Icges(4,4) * t511 + Icges(4,2) * t510 - Icges(4,6) * t579;
t454 = qJD(1) + (t501 * t544 + t502 * t546) * t570;
t453 = qJD(5) * t482 + t487;
t452 = qJD(5) * t484 + t486;
t446 = rSges(6,1) * t481 - rSges(6,2) * t480 + rSges(6,3) * t514;
t438 = rSges(5,1) * t485 - rSges(5,2) * t484 - rSges(5,3) * t512;
t437 = rSges(5,1) * t483 - rSges(5,2) * t482 - rSges(5,3) * t510;
t436 = Icges(5,1) * t485 - Icges(5,4) * t484 - Icges(5,5) * t512;
t435 = Icges(5,1) * t483 - Icges(5,4) * t482 - Icges(5,5) * t510;
t434 = Icges(5,4) * t485 - Icges(5,2) * t484 - Icges(5,6) * t512;
t433 = Icges(5,4) * t483 - Icges(5,2) * t482 - Icges(5,6) * t510;
t432 = Icges(5,5) * t485 - Icges(5,6) * t484 - Icges(5,3) * t512;
t431 = Icges(5,5) * t483 - Icges(5,6) * t482 - Icges(5,3) * t510;
t427 = t537 + ((-t465 - t505) * t547 + t546 * t562) * qJD(2);
t426 = -t545 * t569 + t490 + (t466 * t547 + t544 * t562) * qJD(2);
t425 = rSges(6,1) * t458 - rSges(6,2) * t457 + rSges(6,3) * t484;
t423 = rSges(6,1) * t456 - rSges(6,2) * t455 + rSges(6,3) * t482;
t409 = (t465 * t544 + t466 * t546) * t570 + t561;
t408 = -t437 * t516 + t475 * t487 + t554;
t407 = t438 * t516 - t475 * t486 + t555;
t406 = t437 * t486 - t438 * t487 + t558;
t405 = -t423 * t479 + t446 * t453 + t552;
t404 = t425 * t479 - t446 * t452 + t553;
t403 = t423 * t452 - t425 * t453 + t556;
t402 = qJD(6) * t457 + t453 * t572 - t479 * t574 + t552;
t401 = qJD(6) * t455 - t452 * t572 + t479 * t573 + t553;
t400 = qJD(6) * t480 + t452 * t574 - t453 * t573 + t556;
t1 = t486 * ((-t432 * t512 - t434 * t484 + t436 * t485) * t486 + (-t431 * t512 - t433 * t484 + t435 * t485) * t487 + (-t472 * t512 - t473 * t484 + t474 * t485) * t516) / 0.2e1 + t487 * ((-t432 * t510 - t434 * t482 + t436 * t483) * t486 + (-t431 * t510 - t433 * t482 + t435 * t483) * t487 + (-t472 * t510 - t473 * t482 + t474 * t483) * t516) / 0.2e1 + t516 * ((-t432 * t524 - t434 * t514 + t436 * t515) * t486 + (-t431 * t524 - t433 * t514 + t435 * t515) * t487 + (-t472 * t524 - t473 * t514 + t474 * t515) * t516) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t454 ^ 2 + t476 ^ 2 + t477 ^ 2) / 0.2e1 + m(4) * (t409 ^ 2 + t426 ^ 2 + t427 ^ 2) / 0.2e1 + m(5) * (t406 ^ 2 + t407 ^ 2 + t408 ^ 2) / 0.2e1 + m(6) * (t403 ^ 2 + t404 ^ 2 + t405 ^ 2) / 0.2e1 + m(7) * (t400 ^ 2 + t401 ^ 2 + t402 ^ 2) / 0.2e1 + ((t457 * t595 + t458 * t593 + t484 * t594) * t479 + (t457 * t601 + t597 * t458 + t599 * t484) * t453 + (t600 * t457 + t596 * t458 + t598 * t484) * t452) * t452 / 0.2e1 + ((t455 * t595 + t456 * t593 + t482 * t594) * t479 + (t601 * t455 + t597 * t456 + t599 * t482) * t453 + (t455 * t600 + t456 * t596 + t482 * t598) * t452) * t453 / 0.2e1 + ((t595 * t480 + t593 * t481 + t594 * t514) * t479 + (t480 * t601 + t597 * t481 + t599 * t514) * t453 + (t480 * t600 + t481 * t596 + t514 * t598) * t452) * t479 / 0.2e1 - ((t462 * t510 + t464 * t511 + t498 * t528 + t500 * t529) * t580 + (t492 * t510 + t493 * t511 + t521 * t528 + t522 * t529) * t547 + (-t461 * t510 - t463 * t511 - t497 * t528 - t499 * t529 + t605) * t579) * t588 * t579 / 0.2e1 + (((t492 * t524 - t493 * t525 + t589) * t547 + ((t460 * t547 + t462 * t524 - t464 * t525) * t544 - (t459 * t547 + t461 * t524 - t463 * t525) * t546 + (-t495 * t546 + t496 * t544 + t521 * t551 + t522 * t550) * t547 + ((t498 * t551 + t500 * t550) * t544 - (t497 * t551 + t499 * t550) * t546) * t545) * t545) * t547 + ((-t461 * t512 - t463 * t513 - t497 * t530 - t499 * t531) * t579 + (t492 * t512 + t493 * t513 + t521 * t530 + t522 * t531) * t547 + (t462 * t512 + t464 * t513 + t498 * t530 + t500 * t531 - t605) * t580) * t580) * t588 / 0.2e1;
T  = t1;
