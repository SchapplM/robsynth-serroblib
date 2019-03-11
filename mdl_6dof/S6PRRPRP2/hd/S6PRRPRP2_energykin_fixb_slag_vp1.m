% Calculate kinetic energy for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:22
% EndTime: 2019-03-08 21:28:24
% DurationCPUTime: 2.67s
% Computational Cost: add. (3116->319), mult. (5728->473), div. (0->0), fcn. (6945->12), ass. (0->150)
t597 = Icges(6,1) + Icges(7,1);
t596 = -Icges(6,4) + Icges(7,5);
t595 = Icges(7,4) + Icges(6,5);
t594 = Icges(6,2) + Icges(7,3);
t593 = Icges(7,2) + Icges(6,3);
t592 = -Icges(6,6) + Icges(7,6);
t591 = Icges(4,3) + Icges(5,3);
t590 = rSges(7,1) + pkin(5);
t589 = rSges(7,3) + qJ(6);
t528 = sin(pkin(10));
t530 = cos(pkin(10));
t537 = cos(qJ(2));
t531 = cos(pkin(6));
t535 = sin(qJ(2));
t559 = t531 * t535;
t514 = t528 * t537 + t530 * t559;
t552 = qJ(3) + pkin(11);
t527 = sin(t552);
t547 = cos(t552);
t529 = sin(pkin(6));
t565 = t529 * t530;
t490 = t514 * t547 - t527 * t565;
t558 = t531 * t537;
t513 = t528 * t535 - t530 * t558;
t533 = sin(qJ(5));
t569 = cos(qJ(5));
t461 = t490 * t533 - t513 * t569;
t462 = t490 * t569 + t513 * t533;
t546 = t529 * t547;
t489 = t514 * t527 + t530 * t546;
t588 = t594 * t461 + t596 * t462 + t592 * t489;
t516 = -t528 * t559 + t530 * t537;
t566 = t528 * t529;
t492 = t516 * t547 + t527 * t566;
t515 = t528 * t558 + t530 * t535;
t463 = t492 * t533 - t515 * t569;
t464 = t492 * t569 + t515 * t533;
t491 = t516 * t527 - t528 * t546;
t587 = t594 * t463 + t596 * t464 + t592 * t491;
t586 = t592 * t461 + t595 * t462 + t593 * t489;
t585 = t592 * t463 + t595 * t464 + t593 * t491;
t584 = t596 * t461 + t597 * t462 + t595 * t489;
t583 = t596 * t463 + t597 * t464 + t595 * t491;
t505 = t531 * t527 + t535 * t546;
t561 = t529 * t537;
t493 = t505 * t533 + t561 * t569;
t494 = t505 * t569 - t533 * t561;
t563 = t529 * t535;
t504 = t527 * t563 - t531 * t547;
t582 = t594 * t493 + t596 * t494 + t592 * t504;
t581 = t592 * t493 + t595 * t494 + t593 * t504;
t580 = t596 * t493 + t597 * t494 + t595 * t504;
t534 = sin(qJ(3));
t536 = cos(qJ(3));
t562 = t529 * t536;
t495 = -t514 * t534 - t530 * t562;
t564 = t529 * t534;
t550 = t530 * t564;
t496 = t514 * t536 - t550;
t579 = Icges(4,5) * t496 + Icges(5,5) * t490 + Icges(4,6) * t495 - Icges(5,6) * t489 + t591 * t513;
t497 = -t516 * t534 + t528 * t562;
t551 = t528 * t564;
t498 = t516 * t536 + t551;
t578 = Icges(4,5) * t498 + Icges(5,5) * t492 + Icges(4,6) * t497 - Icges(5,6) * t491 + t591 * t515;
t517 = t531 * t536 - t534 * t563;
t560 = t531 * t534;
t518 = t535 * t562 + t560;
t577 = Icges(4,5) * t518 + Icges(5,5) * t505 + Icges(4,6) * t517 - Icges(5,6) * t504 - t591 * t561;
t576 = qJD(2) ^ 2;
t568 = pkin(3) * t536;
t557 = rSges(7,2) * t489 + t589 * t461 + t462 * t590;
t556 = rSges(7,2) * t491 + t589 * t463 + t464 * t590;
t555 = rSges(7,2) * t504 + t589 * t493 + t494 * t590;
t554 = qJD(2) * t529;
t524 = t528 * t554;
t499 = qJD(3) * t515 + t524;
t526 = qJD(2) * t531;
t549 = t530 * t554;
t486 = pkin(2) * t514 + pkin(8) * t513;
t487 = pkin(2) * t516 + pkin(8) * t515;
t548 = t486 * t524 + t487 * t549 + qJD(1);
t500 = qJD(3) * t513 - t549;
t520 = -qJD(3) * t561 + t526;
t519 = (pkin(2) * t535 - pkin(8) * t537) * t529;
t545 = t487 * t526 - t519 * t524;
t442 = pkin(3) * t551 + qJ(4) * t515 + t516 * t568;
t544 = qJD(4) * t513 + t520 * t442 + t545;
t543 = (-t486 * t531 - t519 * t565) * qJD(2);
t441 = -pkin(3) * t550 + qJ(4) * t513 + t514 * t568;
t542 = -qJD(4) * t561 + t499 * t441 + t548;
t483 = pkin(3) * t560 + (-qJ(4) * t537 + t535 * t568) * t529;
t541 = qJD(4) * t515 + t500 * t483 + t543;
t456 = pkin(4) * t492 + pkin(9) * t491;
t475 = pkin(4) * t505 + pkin(9) * t504;
t540 = t520 * t456 + (-t475 - t483) * t499 + t544;
t455 = pkin(4) * t490 + pkin(9) * t489;
t539 = t499 * t455 + (-t442 - t456) * t500 + t542;
t538 = t500 * t475 + (-t441 - t455) * t520 + t541;
t506 = t531 * rSges(3,3) + (rSges(3,1) * t535 + rSges(3,2) * t537) * t529;
t503 = Icges(3,5) * t531 + (Icges(3,1) * t535 + Icges(3,4) * t537) * t529;
t502 = Icges(3,6) * t531 + (Icges(3,4) * t535 + Icges(3,2) * t537) * t529;
t501 = Icges(3,3) * t531 + (Icges(3,5) * t535 + Icges(3,6) * t537) * t529;
t488 = qJD(5) * t504 + t520;
t484 = t518 * rSges(4,1) + t517 * rSges(4,2) - rSges(4,3) * t561;
t482 = Icges(4,1) * t518 + Icges(4,4) * t517 - Icges(4,5) * t561;
t481 = Icges(4,4) * t518 + Icges(4,2) * t517 - Icges(4,6) * t561;
t477 = rSges(3,1) * t516 - rSges(3,2) * t515 + rSges(3,3) * t566;
t476 = rSges(3,1) * t514 - rSges(3,2) * t513 - rSges(3,3) * t565;
t474 = Icges(3,1) * t516 - Icges(3,4) * t515 + Icges(3,5) * t566;
t473 = Icges(3,1) * t514 - Icges(3,4) * t513 - Icges(3,5) * t565;
t472 = Icges(3,4) * t516 - Icges(3,2) * t515 + Icges(3,6) * t566;
t471 = Icges(3,4) * t514 - Icges(3,2) * t513 - Icges(3,6) * t565;
t470 = Icges(3,5) * t516 - Icges(3,6) * t515 + Icges(3,3) * t566;
t469 = Icges(3,5) * t514 - Icges(3,6) * t513 - Icges(3,3) * t565;
t468 = t505 * rSges(5,1) - t504 * rSges(5,2) - rSges(5,3) * t561;
t467 = Icges(5,1) * t505 - Icges(5,4) * t504 - Icges(5,5) * t561;
t466 = Icges(5,4) * t505 - Icges(5,2) * t504 - Icges(5,6) * t561;
t460 = qJD(5) * t489 + t500;
t459 = qJD(5) * t491 + t499;
t453 = (-t476 * t531 - t506 * t565) * qJD(2);
t452 = (t477 * t531 - t506 * t566) * qJD(2);
t451 = rSges(4,1) * t498 + rSges(4,2) * t497 + rSges(4,3) * t515;
t450 = rSges(4,1) * t496 + rSges(4,2) * t495 + rSges(4,3) * t513;
t449 = Icges(4,1) * t498 + Icges(4,4) * t497 + Icges(4,5) * t515;
t448 = Icges(4,1) * t496 + Icges(4,4) * t495 + Icges(4,5) * t513;
t447 = Icges(4,4) * t498 + Icges(4,2) * t497 + Icges(4,6) * t515;
t446 = Icges(4,4) * t496 + Icges(4,2) * t495 + Icges(4,6) * t513;
t440 = rSges(5,1) * t492 - rSges(5,2) * t491 + rSges(5,3) * t515;
t439 = rSges(5,1) * t490 - rSges(5,2) * t489 + rSges(5,3) * t513;
t438 = Icges(5,1) * t492 - Icges(5,4) * t491 + Icges(5,5) * t515;
t437 = Icges(5,1) * t490 - Icges(5,4) * t489 + Icges(5,5) * t513;
t436 = Icges(5,4) * t492 - Icges(5,2) * t491 + Icges(5,6) * t515;
t435 = Icges(5,4) * t490 - Icges(5,2) * t489 + Icges(5,6) * t513;
t432 = rSges(6,1) * t494 - rSges(6,2) * t493 + rSges(6,3) * t504;
t422 = qJD(1) + (t476 * t528 + t477 * t530) * t554;
t418 = rSges(6,1) * t464 - rSges(6,2) * t463 + rSges(6,3) * t491;
t416 = rSges(6,1) * t462 - rSges(6,2) * t461 + rSges(6,3) * t489;
t402 = -t450 * t520 + t484 * t500 + t543;
t401 = t451 * t520 - t484 * t499 + t545;
t400 = t450 * t499 - t451 * t500 + t548;
t399 = t468 * t500 + (-t439 - t441) * t520 + t541;
t398 = t440 * t520 + (-t468 - t483) * t499 + t544;
t397 = t499 * t439 + (-t440 - t442) * t500 + t542;
t396 = -t416 * t488 + t432 * t460 + t538;
t395 = t418 * t488 - t432 * t459 + t540;
t394 = t459 * t416 - t460 * t418 + t539;
t393 = qJD(6) * t463 + t460 * t555 - t488 * t557 + t538;
t392 = qJD(6) * t461 - t459 * t555 + t488 * t556 + t540;
t391 = qJD(6) * t493 + t459 * t557 - t460 * t556 + t539;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t422 ^ 2 + t452 ^ 2 + t453 ^ 2) / 0.2e1 + m(4) * (t400 ^ 2 + t401 ^ 2 + t402 ^ 2) / 0.2e1 + m(5) * (t397 ^ 2 + t398 ^ 2 + t399 ^ 2) / 0.2e1 + m(6) * (t394 ^ 2 + t395 ^ 2 + t396 ^ 2) / 0.2e1 + m(7) * (t391 ^ 2 + t392 ^ 2 + t393 ^ 2) / 0.2e1 - t576 * ((-t470 * t565 - t472 * t513 + t474 * t514) * t566 - (-t469 * t565 - t471 * t513 + t473 * t514) * t565 + (-t501 * t565 - t502 * t513 + t503 * t514) * t531) * t565 / 0.2e1 + ((t463 * t582 + t464 * t580 + t491 * t581) * t488 + (t463 * t588 + t584 * t464 + t586 * t491) * t460 + (t587 * t463 + t583 * t464 + t585 * t491) * t459) * t459 / 0.2e1 + ((t461 * t582 + t462 * t580 + t489 * t581) * t488 + (t588 * t461 + t584 * t462 + t586 * t489) * t460 + (t461 * t587 + t462 * t583 + t489 * t585) * t459) * t460 / 0.2e1 + ((t582 * t493 + t580 * t494 + t581 * t504) * t488 + (t493 * t588 + t584 * t494 + t586 * t504) * t460 + (t493 * t587 + t494 * t583 + t504 * t585) * t459) * t488 / 0.2e1 + ((-t466 * t491 + t467 * t492 + t481 * t497 + t482 * t498 + t515 * t577) * t520 + (-t435 * t491 + t437 * t492 + t446 * t497 + t448 * t498 + t515 * t579) * t500 + (-t436 * t491 + t438 * t492 + t447 * t497 + t449 * t498 + t578 * t515) * t499) * t499 / 0.2e1 + ((-t466 * t489 + t467 * t490 + t481 * t495 + t482 * t496 + t513 * t577) * t520 + (-t435 * t489 + t437 * t490 + t446 * t495 + t448 * t496 + t579 * t513) * t500 + (-t436 * t489 + t438 * t490 + t447 * t495 + t449 * t496 + t513 * t578) * t499) * t500 / 0.2e1 + ((-t504 * t466 + t505 * t467 + t517 * t481 + t518 * t482 - t577 * t561) * t520 + (-t504 * t435 + t505 * t437 + t517 * t446 + t518 * t448 - t579 * t561) * t500 + (-t504 * t436 + t505 * t438 + t517 * t447 + t518 * t449 - t578 * t561) * t499) * t520 / 0.2e1 + (t531 * (t531 ^ 2 * t501 + (((t472 * t537 + t474 * t535) * t528 - (t471 * t537 + t473 * t535) * t530) * t529 + (-t469 * t530 + t470 * t528 + t502 * t537 + t503 * t535) * t531) * t529) + ((t470 * t566 - t472 * t515 + t474 * t516) * t566 - (t469 * t566 - t471 * t515 + t473 * t516) * t565 + (t501 * t566 - t502 * t515 + t503 * t516) * t531) * t566) * t576 / 0.2e1;
T  = t1;
