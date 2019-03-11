% Calculate kinetic energy for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:03:49
% EndTime: 2019-03-09 17:03:52
% DurationCPUTime: 3.01s
% Computational Cost: add. (3185->327), mult. (5780->481), div. (0->0), fcn. (6979->12), ass. (0->153)
t602 = Icges(6,1) + Icges(7,1);
t601 = -Icges(6,4) + Icges(7,5);
t600 = Icges(7,4) + Icges(6,5);
t599 = Icges(6,2) + Icges(7,3);
t598 = Icges(7,2) + Icges(6,3);
t597 = -Icges(6,6) + Icges(7,6);
t596 = Icges(4,3) + Icges(5,3);
t595 = rSges(7,1) + pkin(5);
t594 = rSges(7,3) + qJ(6);
t535 = sin(qJ(2));
t536 = sin(qJ(1));
t538 = cos(qJ(2));
t539 = cos(qJ(1));
t570 = cos(pkin(6));
t553 = t539 * t570;
t513 = t535 * t536 - t538 * t553;
t514 = t535 * t553 + t536 * t538;
t531 = sin(pkin(6));
t565 = t531 * t539;
t473 = Icges(3,5) * t514 - Icges(3,6) * t513 - Icges(3,3) * t565;
t554 = t536 * t570;
t515 = t539 * t535 + t538 * t554;
t516 = -t535 * t554 + t539 * t538;
t568 = t531 * t536;
t474 = Icges(3,5) * t516 - Icges(3,6) * t515 + Icges(3,3) * t568;
t593 = t531 * (t473 * t539 - t474 * t536);
t558 = qJ(3) + pkin(11);
t530 = sin(t558);
t551 = cos(t558);
t490 = t514 * t551 - t530 * t565;
t533 = sin(qJ(5));
t573 = cos(qJ(5));
t459 = t490 * t533 - t513 * t573;
t460 = t490 * t573 + t513 * t533;
t550 = t531 * t551;
t489 = t514 * t530 + t539 * t550;
t592 = t599 * t459 + t601 * t460 + t597 * t489;
t492 = t516 * t551 + t530 * t568;
t461 = t492 * t533 - t515 * t573;
t462 = t492 * t573 + t515 * t533;
t491 = t516 * t530 - t536 * t550;
t591 = t599 * t461 + t601 * t462 + t597 * t491;
t590 = t597 * t459 + t600 * t460 + t598 * t489;
t589 = t597 * t461 + t600 * t462 + t598 * t491;
t588 = t601 * t459 + t602 * t460 + t600 * t489;
t587 = t601 * t461 + t602 * t462 + t600 * t491;
t503 = t530 * t570 + t535 * t550;
t566 = t531 * t538;
t487 = t503 * t533 + t566 * t573;
t488 = t503 * t573 - t533 * t566;
t569 = t531 * t535;
t502 = t530 * t569 - t551 * t570;
t586 = t599 * t487 + t601 * t488 + t597 * t502;
t585 = t597 * t487 + t600 * t488 + t598 * t502;
t584 = t601 * t487 + t602 * t488 + t600 * t502;
t534 = sin(qJ(3));
t537 = cos(qJ(3));
t493 = -t514 * t534 - t537 * t565;
t556 = t534 * t565;
t494 = t514 * t537 - t556;
t583 = Icges(4,5) * t494 + Icges(5,5) * t490 + Icges(4,6) * t493 - Icges(5,6) * t489 + t596 * t513;
t567 = t531 * t537;
t495 = -t516 * t534 + t536 * t567;
t557 = t534 * t568;
t496 = t516 * t537 + t557;
t582 = Icges(4,5) * t496 + Icges(5,5) * t492 + Icges(4,6) * t495 - Icges(5,6) * t491 + t596 * t515;
t511 = -t534 * t569 + t537 * t570;
t552 = t570 * t534;
t512 = t535 * t567 + t552;
t581 = Icges(4,5) * t512 + Icges(5,5) * t503 + Icges(4,6) * t511 - Icges(5,6) * t502 - t596 * t566;
t572 = pkin(3) * t537;
t564 = rSges(7,2) * t489 + t594 * t459 + t460 * t595;
t563 = rSges(7,2) * t491 + t594 * t461 + t462 * t595;
t562 = rSges(7,2) * t502 + t594 * t487 + t488 * t595;
t485 = pkin(2) * t514 + pkin(9) * t513;
t486 = pkin(2) * t516 + pkin(9) * t515;
t559 = qJD(2) * t531;
t526 = t536 * t559;
t555 = t539 * t559;
t561 = t485 * t526 + t486 * t555;
t497 = qJD(3) * t515 + t526;
t560 = qJD(1) * (pkin(1) * t536 - pkin(8) * t565);
t527 = qJD(2) * t570 + qJD(1);
t498 = qJD(3) * t513 - t555;
t518 = -qJD(3) * t566 + t527;
t517 = (pkin(2) * t535 - pkin(9) * t538) * t531;
t519 = qJD(1) * (pkin(1) * t539 + pkin(8) * t568);
t548 = t527 * t486 - t517 * t526 + t519;
t441 = -pkin(3) * t556 + qJ(4) * t513 + t514 * t572;
t547 = -qJD(4) * t566 + t497 * t441 + t561;
t442 = pkin(3) * t557 + qJ(4) * t515 + t516 * t572;
t546 = qJD(4) * t513 + t518 * t442 + t548;
t545 = -t485 * t527 - t517 * t555 - t560;
t455 = pkin(4) * t490 + pkin(10) * t489;
t456 = pkin(4) * t492 + pkin(10) * t491;
t544 = t497 * t455 + (-t442 - t456) * t498 + t547;
t471 = pkin(3) * t552 + (-qJ(4) * t538 + t535 * t572) * t531;
t543 = qJD(4) * t515 + t498 * t471 + t545;
t467 = pkin(4) * t503 + pkin(10) * t502;
t542 = t518 * t456 + (-t467 - t471) * t497 + t546;
t541 = t498 * t467 + (-t441 - t455) * t518 + t543;
t523 = rSges(2,1) * t539 - rSges(2,2) * t536;
t522 = rSges(2,1) * t536 + rSges(2,2) * t539;
t504 = t570 * rSges(3,3) + (rSges(3,1) * t535 + rSges(3,2) * t538) * t531;
t501 = Icges(3,5) * t570 + (Icges(3,1) * t535 + Icges(3,4) * t538) * t531;
t500 = Icges(3,6) * t570 + (Icges(3,4) * t535 + Icges(3,2) * t538) * t531;
t499 = Icges(3,3) * t570 + (Icges(3,5) * t535 + Icges(3,6) * t538) * t531;
t484 = qJD(5) * t502 + t518;
t481 = rSges(3,1) * t516 - rSges(3,2) * t515 + rSges(3,3) * t568;
t480 = rSges(3,1) * t514 - rSges(3,2) * t513 - rSges(3,3) * t565;
t478 = Icges(3,1) * t516 - Icges(3,4) * t515 + Icges(3,5) * t568;
t477 = Icges(3,1) * t514 - Icges(3,4) * t513 - Icges(3,5) * t565;
t476 = Icges(3,4) * t516 - Icges(3,2) * t515 + Icges(3,6) * t568;
t475 = Icges(3,4) * t514 - Icges(3,2) * t513 - Icges(3,6) * t565;
t472 = rSges(4,1) * t512 + rSges(4,2) * t511 - rSges(4,3) * t566;
t470 = Icges(4,1) * t512 + Icges(4,4) * t511 - Icges(4,5) * t566;
t469 = Icges(4,4) * t512 + Icges(4,2) * t511 - Icges(4,6) * t566;
t466 = rSges(5,1) * t503 - rSges(5,2) * t502 - rSges(5,3) * t566;
t465 = Icges(5,1) * t503 - Icges(5,4) * t502 - Icges(5,5) * t566;
t464 = Icges(5,4) * t503 - Icges(5,2) * t502 - Icges(5,6) * t566;
t458 = qJD(5) * t489 + t498;
t457 = qJD(5) * t491 + t497;
t451 = rSges(4,1) * t496 + rSges(4,2) * t495 + rSges(4,3) * t515;
t450 = rSges(4,1) * t494 + rSges(4,2) * t493 + rSges(4,3) * t513;
t449 = Icges(4,1) * t496 + Icges(4,4) * t495 + Icges(4,5) * t515;
t448 = Icges(4,1) * t494 + Icges(4,4) * t493 + Icges(4,5) * t513;
t447 = Icges(4,4) * t496 + Icges(4,2) * t495 + Icges(4,6) * t515;
t446 = Icges(4,4) * t494 + Icges(4,2) * t493 + Icges(4,6) * t513;
t440 = rSges(5,1) * t492 - rSges(5,2) * t491 + rSges(5,3) * t515;
t439 = rSges(5,1) * t490 - rSges(5,2) * t489 + rSges(5,3) * t513;
t438 = Icges(5,1) * t492 - Icges(5,4) * t491 + Icges(5,5) * t515;
t437 = Icges(5,1) * t490 - Icges(5,4) * t489 + Icges(5,5) * t513;
t436 = Icges(5,4) * t492 - Icges(5,2) * t491 + Icges(5,6) * t515;
t435 = Icges(5,4) * t490 - Icges(5,2) * t489 + Icges(5,6) * t513;
t432 = t481 * t527 - t504 * t526 + t519;
t431 = -t480 * t527 - t504 * t555 - t560;
t430 = rSges(6,1) * t488 - rSges(6,2) * t487 + rSges(6,3) * t502;
t420 = (t480 * t536 + t481 * t539) * t559;
t416 = rSges(6,1) * t462 - rSges(6,2) * t461 + rSges(6,3) * t491;
t414 = rSges(6,1) * t460 - rSges(6,2) * t459 + rSges(6,3) * t489;
t400 = t451 * t518 - t472 * t497 + t548;
t399 = -t450 * t518 + t472 * t498 + t545;
t398 = t450 * t497 - t451 * t498 + t561;
t397 = t440 * t518 + (-t466 - t471) * t497 + t546;
t396 = t466 * t498 + (-t439 - t441) * t518 + t543;
t395 = t439 * t497 + (-t440 - t442) * t498 + t547;
t394 = t416 * t484 - t430 * t457 + t542;
t393 = -t414 * t484 + t430 * t458 + t541;
t392 = t414 * t457 - t416 * t458 + t544;
t391 = qJD(6) * t459 - t457 * t562 + t484 * t563 + t542;
t390 = qJD(6) * t461 + t458 * t562 - t484 * t564 + t541;
t389 = qJD(6) * t487 + t457 * t564 - t458 * t563 + t544;
t1 = m(3) * (t420 ^ 2 + t431 ^ 2 + t432 ^ 2) / 0.2e1 + t527 * ((t570 * t474 + (t476 * t538 + t478 * t535) * t531) * t526 - (t570 * t473 + (t475 * t538 + t477 * t535) * t531) * t555 + (t570 * t499 + (t500 * t538 + t501 * t535) * t531) * t527) / 0.2e1 + m(7) * (t389 ^ 2 + t390 ^ 2 + t391 ^ 2) / 0.2e1 + m(6) * (t392 ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + m(5) * (t395 ^ 2 + t396 ^ 2 + t397 ^ 2) / 0.2e1 + m(4) * (t398 ^ 2 + t399 ^ 2 + t400 ^ 2) / 0.2e1 - ((-t499 * t565 - t500 * t513 + t501 * t514) * t527 + ((-t476 * t513 + t478 * t514) * t536 + (t513 * t475 - t514 * t477 + t593) * t539) * t559) * t555 / 0.2e1 + ((t499 * t568 - t500 * t515 + t501 * t516) * t527 + (-(-t475 * t515 + t477 * t516) * t539 + (-t515 * t476 + t516 * t478 - t593) * t536) * t559) * t526 / 0.2e1 + ((t461 * t586 + t462 * t584 + t491 * t585) * t484 + (t461 * t592 + t462 * t588 + t491 * t590) * t458 + (t591 * t461 + t587 * t462 + t589 * t491) * t457) * t457 / 0.2e1 + ((t586 * t459 + t460 * t584 + t585 * t489) * t484 + (t592 * t459 + t588 * t460 + t590 * t489) * t458 + (t459 * t591 + t460 * t587 + t489 * t589) * t457) * t458 / 0.2e1 + ((t586 * t487 + t584 * t488 + t585 * t502) * t484 + (t487 * t592 + t488 * t588 + t502 * t590) * t458 + (t487 * t591 + t488 * t587 + t502 * t589) * t457) * t484 / 0.2e1 + ((-t464 * t491 + t465 * t492 + t469 * t495 + t470 * t496 + t515 * t581) * t518 + (-t435 * t491 + t437 * t492 + t446 * t495 + t448 * t496 + t515 * t583) * t498 + (-t436 * t491 + t438 * t492 + t447 * t495 + t449 * t496 + t582 * t515) * t497) * t497 / 0.2e1 + ((-t464 * t489 + t465 * t490 + t469 * t493 + t470 * t494 + t513 * t581) * t518 + (-t435 * t489 + t437 * t490 + t446 * t493 + t448 * t494 + t583 * t513) * t498 + (-t436 * t489 + t438 * t490 + t447 * t493 + t449 * t494 + t513 * t582) * t497) * t498 / 0.2e1 + ((-t464 * t502 + t465 * t503 + t469 * t511 + t470 * t512 - t581 * t566) * t518 + (-t435 * t502 + t437 * t503 + t446 * t511 + t448 * t512 - t566 * t583) * t498 + (-t436 * t502 + t438 * t503 + t447 * t511 + t449 * t512 - t566 * t582) * t497) * t518 / 0.2e1 + (Icges(2,3) + m(2) * (t522 ^ 2 + t523 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
