% Calculate kinetic energy for
% S6PPRRRP1
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
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:18
% EndTime: 2019-03-08 18:52:21
% DurationCPUTime: 2.10s
% Computational Cost: add. (4987->270), mult. (13746->416), div. (0->0), fcn. (17946->14), ass. (0->136)
t605 = Icges(6,1) + Icges(7,1);
t604 = Icges(6,4) + Icges(7,4);
t603 = Icges(6,5) + Icges(7,5);
t602 = Icges(6,2) + Icges(7,2);
t601 = Icges(6,6) + Icges(7,6);
t600 = Icges(6,3) + Icges(7,3);
t599 = rSges(7,3) + qJ(6);
t543 = sin(pkin(11));
t545 = cos(pkin(11));
t578 = sin(pkin(12));
t582 = cos(pkin(6));
t564 = t582 * t578;
t580 = cos(pkin(12));
t533 = t543 * t580 + t545 * t564;
t549 = sin(qJ(3));
t566 = t582 * t580;
t557 = t543 * t578 - t545 * t566;
t581 = cos(pkin(7));
t553 = t557 * t581;
t544 = sin(pkin(6));
t579 = sin(pkin(7));
t567 = t544 * t579;
t586 = cos(qJ(3));
t513 = t533 * t586 + (-t545 * t567 - t553) * t549;
t568 = t544 * t581;
t526 = -t545 * t568 + t557 * t579;
t548 = sin(qJ(4));
t585 = cos(qJ(4));
t498 = t513 * t585 + t526 * t548;
t562 = t586 * t567;
t512 = t533 * t549 + t545 * t562 + t553 * t586;
t547 = sin(qJ(5));
t550 = cos(qJ(5));
t474 = -t498 * t547 + t512 * t550;
t577 = t512 * t547;
t475 = t498 * t550 + t577;
t497 = t513 * t548 - t526 * t585;
t598 = t601 * t474 + t603 * t475 + t600 * t497;
t534 = -t543 * t564 + t545 * t580;
t556 = t543 * t566 + t545 * t578;
t552 = t556 * t581;
t515 = t534 * t586 + (t543 * t567 - t552) * t549;
t527 = t543 * t568 + t556 * t579;
t500 = t515 * t585 + t527 * t548;
t514 = t534 * t549 - t543 * t562 + t552 * t586;
t476 = -t500 * t547 + t514 * t550;
t576 = t514 * t547;
t477 = t500 * t550 + t576;
t499 = t515 * t548 - t527 * t585;
t597 = t601 * t476 + t603 * t477 + t600 * t499;
t596 = t602 * t474 + t604 * t475 + t601 * t497;
t595 = t602 * t476 + t604 * t477 + t601 * t499;
t594 = t604 * t474 + t605 * t475 + t603 * t497;
t593 = t604 * t476 + t605 * t477 + t603 * t499;
t563 = t581 * t580;
t565 = t582 * t579;
t525 = t549 * t565 + (t549 * t563 + t578 * t586) * t544;
t532 = -t567 * t580 + t581 * t582;
t517 = t525 * t585 + t532 * t548;
t524 = -t565 * t586 + (t549 * t578 - t563 * t586) * t544;
t501 = -t517 * t547 + t524 * t550;
t575 = t524 * t547;
t502 = t517 * t550 + t575;
t516 = t525 * t548 - t532 * t585;
t592 = t601 * t501 + t603 * t502 + t600 * t516;
t591 = t602 * t501 + t604 * t502 + t601 * t516;
t590 = t604 * t501 + t605 * t502 + t603 * t516;
t584 = pkin(5) * t550;
t574 = rSges(7,1) * t475 + rSges(7,2) * t474 + pkin(5) * t577 + t599 * t497 + t498 * t584;
t573 = rSges(7,1) * t477 + rSges(7,2) * t476 + pkin(5) * t576 + t599 * t499 + t500 * t584;
t572 = rSges(7,1) * t502 + rSges(7,2) * t501 + pkin(5) * t575 + t599 * t516 + t517 * t584;
t522 = qJD(3) * t526;
t503 = qJD(4) * t512 + t522;
t523 = qJD(3) * t527;
t504 = qJD(4) * t514 + t523;
t531 = qJD(3) * t532;
t518 = qJD(4) * t524 + t531;
t571 = qJD(2) * t544;
t539 = qJD(2) * t582 + qJD(1);
t569 = t545 * t571;
t492 = pkin(3) * t513 + pkin(9) * t512;
t509 = pkin(3) * t525 + pkin(9) * t524;
t538 = t543 * t571;
t561 = -t492 * t531 + t509 * t522 + t538;
t493 = pkin(3) * t515 + pkin(9) * t514;
t560 = t492 * t523 - t493 * t522 + t539;
t559 = t493 * t531 - t509 * t523 - t569;
t470 = t498 * pkin(4) + t497 * pkin(10);
t494 = t517 * pkin(4) + t516 * pkin(10);
t558 = -t518 * t470 + t503 * t494 + t561;
t471 = t500 * pkin(4) + t499 * pkin(10);
t555 = t504 * t470 - t503 * t471 + t560;
t554 = t518 * t471 - t504 * t494 + t559;
t508 = rSges(4,1) * t525 - rSges(4,2) * t524 + rSges(4,3) * t532;
t507 = Icges(4,1) * t525 - Icges(4,4) * t524 + Icges(4,5) * t532;
t506 = Icges(4,4) * t525 - Icges(4,2) * t524 + Icges(4,6) * t532;
t505 = Icges(4,5) * t525 - Icges(4,6) * t524 + Icges(4,3) * t532;
t495 = qJD(5) * t516 + t518;
t490 = rSges(5,1) * t517 - rSges(5,2) * t516 + rSges(5,3) * t524;
t489 = Icges(5,1) * t517 - Icges(5,4) * t516 + Icges(5,5) * t524;
t488 = Icges(5,4) * t517 - Icges(5,2) * t516 + Icges(5,6) * t524;
t487 = Icges(5,5) * t517 - Icges(5,6) * t516 + Icges(5,3) * t524;
t485 = rSges(4,1) * t515 - rSges(4,2) * t514 + rSges(4,3) * t527;
t484 = rSges(4,1) * t513 - rSges(4,2) * t512 + rSges(4,3) * t526;
t483 = Icges(4,1) * t515 - Icges(4,4) * t514 + Icges(4,5) * t527;
t482 = Icges(4,1) * t513 - Icges(4,4) * t512 + Icges(4,5) * t526;
t481 = Icges(4,4) * t515 - Icges(4,2) * t514 + Icges(4,6) * t527;
t480 = Icges(4,4) * t513 - Icges(4,2) * t512 + Icges(4,6) * t526;
t479 = Icges(4,5) * t515 - Icges(4,6) * t514 + Icges(4,3) * t527;
t478 = Icges(4,5) * t513 - Icges(4,6) * t512 + Icges(4,3) * t526;
t473 = qJD(5) * t499 + t504;
t472 = qJD(5) * t497 + t503;
t467 = rSges(6,1) * t502 + rSges(6,2) * t501 + rSges(6,3) * t516;
t459 = rSges(5,1) * t500 - rSges(5,2) * t499 + rSges(5,3) * t514;
t458 = rSges(5,1) * t498 - rSges(5,2) * t497 + rSges(5,3) * t512;
t457 = Icges(5,1) * t500 - Icges(5,4) * t499 + Icges(5,5) * t514;
t456 = Icges(5,1) * t498 - Icges(5,4) * t497 + Icges(5,5) * t512;
t455 = Icges(5,4) * t500 - Icges(5,2) * t499 + Icges(5,6) * t514;
t454 = Icges(5,4) * t498 - Icges(5,2) * t497 + Icges(5,6) * t512;
t453 = Icges(5,5) * t500 - Icges(5,6) * t499 + Icges(5,3) * t514;
t452 = Icges(5,5) * t498 - Icges(5,6) * t497 + Icges(5,3) * t512;
t449 = -t569 + (t485 * t532 - t508 * t527) * qJD(3);
t448 = t538 + (-t484 * t532 + t508 * t526) * qJD(3);
t447 = (t484 * t527 - t485 * t526) * qJD(3) + t539;
t446 = rSges(6,1) * t477 + rSges(6,2) * t476 + rSges(6,3) * t499;
t444 = rSges(6,1) * t475 + rSges(6,2) * t474 + rSges(6,3) * t497;
t428 = t459 * t518 - t490 * t504 + t559;
t427 = -t458 * t518 + t490 * t503 + t561;
t426 = t458 * t504 - t459 * t503 + t560;
t425 = t446 * t495 - t467 * t473 + t554;
t424 = -t444 * t495 + t467 * t472 + t558;
t423 = t444 * t473 - t446 * t472 + t555;
t422 = qJD(6) * t497 - t473 * t572 + t495 * t573 + t554;
t421 = qJD(6) * t499 + t472 * t572 - t495 * t574 + t558;
t420 = qJD(6) * t516 - t472 * t573 + t473 * t574 + t555;
t1 = m(6) * (t423 ^ 2 + t424 ^ 2 + t425 ^ 2) / 0.2e1 + m(7) * (t420 ^ 2 + t421 ^ 2 + t422 ^ 2) / 0.2e1 + m(5) * (t426 ^ 2 + t427 ^ 2 + t428 ^ 2) / 0.2e1 + m(4) * (t447 ^ 2 + t448 ^ 2 + t449 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t539 ^ 2 + (t543 ^ 2 + t545 ^ 2) * qJD(2) ^ 2 * t544 ^ 2) / 0.2e1 + t518 * ((t453 * t524 - t455 * t516 + t457 * t517) * t504 + (t452 * t524 - t454 * t516 + t456 * t517) * t503 + (t524 * t487 - t516 * t488 + t517 * t489) * t518) / 0.2e1 + t504 * ((t514 * t453 - t499 * t455 + t500 * t457) * t504 + (t452 * t514 - t454 * t499 + t456 * t500) * t503 + (t487 * t514 - t488 * t499 + t489 * t500) * t518) / 0.2e1 + t503 * ((t453 * t512 - t455 * t497 + t457 * t498) * t504 + (t512 * t452 - t497 * t454 + t498 * t456) * t503 + (t487 * t512 - t488 * t497 + t489 * t498) * t518) / 0.2e1 + ((t474 * t591 + t475 * t590 + t497 * t592) * t495 + (t474 * t595 + t475 * t593 + t497 * t597) * t473 + (t596 * t474 + t594 * t475 + t598 * t497) * t472) * t472 / 0.2e1 + ((t476 * t591 + t477 * t590 + t499 * t592) * t495 + (t595 * t476 + t593 * t477 + t597 * t499) * t473 + (t596 * t476 + t594 * t477 + t499 * t598) * t472) * t473 / 0.2e1 + ((t591 * t501 + t590 * t502 + t592 * t516) * t495 + (t501 * t595 + t502 * t593 + t516 * t597) * t473 + (t596 * t501 + t594 * t502 + t598 * t516) * t472) * t495 / 0.2e1 + (t526 * ((t479 * t526 - t481 * t512 + t483 * t513) * t527 + (t526 * t478 - t512 * t480 + t513 * t482) * t526 + (t505 * t526 - t506 * t512 + t507 * t513) * t532) + t532 * ((t479 * t532 - t481 * t524 + t483 * t525) * t527 + (t478 * t532 - t480 * t524 + t482 * t525) * t526 + (t532 * t505 - t524 * t506 + t525 * t507) * t532) + t527 * ((t527 * t479 - t514 * t481 + t515 * t483) * t527 + (t478 * t527 - t480 * t514 + t482 * t515) * t526 + (t505 * t527 - t506 * t514 + t507 * t515) * t532)) * qJD(3) ^ 2 / 0.2e1;
T  = t1;
