% Calculate kinetic energy for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:08
% EndTime: 2019-03-09 00:01:11
% DurationCPUTime: 2.82s
% Computational Cost: add. (3221->325), mult. (5938->503), div. (0->0), fcn. (7197->12), ass. (0->155)
t589 = Icges(6,1) + Icges(7,1);
t588 = -Icges(6,4) + Icges(7,5);
t587 = Icges(7,4) + Icges(6,5);
t586 = Icges(6,2) + Icges(7,3);
t585 = Icges(7,2) + Icges(6,3);
t584 = -Icges(6,6) + Icges(7,6);
t583 = rSges(7,1) + pkin(5);
t582 = rSges(7,3) + qJ(6);
t527 = sin(pkin(11));
t529 = cos(pkin(11));
t535 = cos(qJ(2));
t530 = cos(pkin(6));
t533 = sin(qJ(2));
t558 = t530 * t533;
t513 = t527 * t535 + t529 * t558;
t556 = qJ(3) + qJ(4);
t526 = sin(t556);
t546 = cos(t556);
t528 = sin(pkin(6));
t564 = t528 * t529;
t488 = t513 * t546 - t526 * t564;
t557 = t530 * t535;
t512 = t527 * t533 - t529 * t557;
t531 = sin(qJ(5));
t568 = cos(qJ(5));
t457 = t488 * t531 - t512 * t568;
t458 = t488 * t568 + t512 * t531;
t545 = t528 * t546;
t487 = t513 * t526 + t529 * t545;
t581 = t457 * t586 + t458 * t588 + t584 * t487;
t515 = -t527 * t558 + t529 * t535;
t565 = t527 * t528;
t490 = t515 * t546 + t526 * t565;
t514 = t527 * t557 + t529 * t533;
t459 = t490 * t531 - t514 * t568;
t460 = t490 * t568 + t514 * t531;
t489 = t515 * t526 - t527 * t545;
t580 = t459 * t586 + t460 * t588 + t584 * t489;
t579 = t584 * t457 + t458 * t587 + t585 * t487;
t578 = t584 * t459 + t460 * t587 + t585 * t489;
t577 = t588 * t457 + t458 * t589 + t587 * t487;
t576 = t588 * t459 + t460 * t589 + t587 * t489;
t505 = t530 * t526 + t533 * t545;
t560 = t528 * t535;
t491 = t505 * t531 + t560 * t568;
t492 = t505 * t568 - t531 * t560;
t562 = t528 * t533;
t504 = t526 * t562 - t530 * t546;
t575 = t491 * t586 + t492 * t588 + t584 * t504;
t574 = t584 * t491 + t492 * t587 + t585 * t504;
t573 = t588 * t491 + t492 * t589 + t587 * t504;
t572 = qJD(2) ^ 2;
t534 = cos(qJ(3));
t567 = pkin(3) * t534;
t532 = sin(qJ(3));
t563 = t528 * t532;
t561 = t528 * t534;
t559 = t530 * t532;
t555 = rSges(7,2) * t487 + t582 * t457 + t583 * t458;
t554 = rSges(7,2) * t489 + t582 * t459 + t583 * t460;
t553 = rSges(7,2) * t504 + t582 * t491 + t583 * t492;
t552 = qJD(2) * t528;
t523 = t527 * t552;
t497 = qJD(3) * t514 + t523;
t525 = qJD(2) * t530;
t550 = t527 * t563;
t549 = t529 * t563;
t465 = qJD(4) * t514 + t497;
t548 = t529 * t552;
t485 = t513 * pkin(2) + t512 * pkin(8);
t486 = t515 * pkin(2) + t514 * pkin(8);
t547 = t485 * t523 + t486 * t548 + qJD(1);
t498 = qJD(3) * t512 - t548;
t518 = (pkin(2) * t533 - pkin(8) * t535) * t528;
t544 = t486 * t525 - t518 * t523;
t466 = qJD(4) * t512 + t498;
t499 = t525 + (-qJD(3) - qJD(4)) * t560;
t439 = -pkin(3) * t549 + pkin(9) * t512 + t513 * t567;
t440 = pkin(3) * t550 + pkin(9) * t514 + t515 * t567;
t543 = t497 * t439 - t440 * t498 + t547;
t542 = (-t485 * t530 - t518 * t564) * qJD(2);
t482 = pkin(3) * t559 + (-pkin(9) * t535 + t533 * t567) * t528;
t519 = -qJD(3) * t560 + t525;
t541 = t519 * t440 - t482 * t497 + t544;
t454 = pkin(4) * t488 + pkin(10) * t487;
t455 = pkin(4) * t490 + pkin(10) * t489;
t540 = t465 * t454 - t455 * t466 + t543;
t539 = -t439 * t519 + t498 * t482 + t542;
t481 = pkin(4) * t505 + pkin(10) * t504;
t538 = t499 * t455 - t465 * t481 + t541;
t537 = -t454 * t499 + t466 * t481 + t539;
t517 = t533 * t561 + t559;
t516 = t530 * t534 - t532 * t562;
t503 = rSges(3,3) * t530 + (rSges(3,1) * t533 + rSges(3,2) * t535) * t528;
t502 = Icges(3,5) * t530 + (Icges(3,1) * t533 + Icges(3,4) * t535) * t528;
t501 = Icges(3,6) * t530 + (Icges(3,4) * t533 + Icges(3,2) * t535) * t528;
t500 = Icges(3,3) * t530 + (Icges(3,5) * t533 + Icges(3,6) * t535) * t528;
t496 = t515 * t534 + t550;
t495 = -t515 * t532 + t527 * t561;
t494 = t513 * t534 - t549;
t493 = -t513 * t532 - t529 * t561;
t483 = qJD(5) * t504 + t499;
t480 = rSges(4,1) * t517 + rSges(4,2) * t516 - rSges(4,3) * t560;
t479 = Icges(4,1) * t517 + Icges(4,4) * t516 - Icges(4,5) * t560;
t478 = Icges(4,4) * t517 + Icges(4,2) * t516 - Icges(4,6) * t560;
t477 = Icges(4,5) * t517 + Icges(4,6) * t516 - Icges(4,3) * t560;
t474 = rSges(3,1) * t515 - rSges(3,2) * t514 + rSges(3,3) * t565;
t473 = rSges(3,1) * t513 - rSges(3,2) * t512 - rSges(3,3) * t564;
t472 = Icges(3,1) * t515 - Icges(3,4) * t514 + Icges(3,5) * t565;
t471 = Icges(3,1) * t513 - Icges(3,4) * t512 - Icges(3,5) * t564;
t470 = Icges(3,4) * t515 - Icges(3,2) * t514 + Icges(3,6) * t565;
t469 = Icges(3,4) * t513 - Icges(3,2) * t512 - Icges(3,6) * t564;
t468 = Icges(3,5) * t515 - Icges(3,6) * t514 + Icges(3,3) * t565;
t467 = Icges(3,5) * t513 - Icges(3,6) * t512 - Icges(3,3) * t564;
t464 = rSges(5,1) * t505 - rSges(5,2) * t504 - rSges(5,3) * t560;
t463 = Icges(5,1) * t505 - Icges(5,4) * t504 - Icges(5,5) * t560;
t462 = Icges(5,4) * t505 - Icges(5,2) * t504 - Icges(5,6) * t560;
t461 = Icges(5,5) * t505 - Icges(5,6) * t504 - Icges(5,3) * t560;
t452 = (-t473 * t530 - t503 * t564) * qJD(2);
t451 = (t474 * t530 - t503 * t565) * qJD(2);
t450 = qJD(5) * t487 + t466;
t449 = qJD(5) * t489 + t465;
t448 = rSges(4,1) * t496 + rSges(4,2) * t495 + rSges(4,3) * t514;
t447 = rSges(4,1) * t494 + rSges(4,2) * t493 + rSges(4,3) * t512;
t446 = Icges(4,1) * t496 + Icges(4,4) * t495 + Icges(4,5) * t514;
t445 = Icges(4,1) * t494 + Icges(4,4) * t493 + Icges(4,5) * t512;
t444 = Icges(4,4) * t496 + Icges(4,2) * t495 + Icges(4,6) * t514;
t443 = Icges(4,4) * t494 + Icges(4,2) * t493 + Icges(4,6) * t512;
t442 = Icges(4,5) * t496 + Icges(4,6) * t495 + Icges(4,3) * t514;
t441 = Icges(4,5) * t494 + Icges(4,6) * t493 + Icges(4,3) * t512;
t438 = rSges(5,1) * t490 - rSges(5,2) * t489 + rSges(5,3) * t514;
t437 = rSges(5,1) * t488 - rSges(5,2) * t487 + rSges(5,3) * t512;
t435 = Icges(5,1) * t490 - Icges(5,4) * t489 + Icges(5,5) * t514;
t434 = Icges(5,1) * t488 - Icges(5,4) * t487 + Icges(5,5) * t512;
t433 = Icges(5,4) * t490 - Icges(5,2) * t489 + Icges(5,6) * t514;
t432 = Icges(5,4) * t488 - Icges(5,2) * t487 + Icges(5,6) * t512;
t431 = Icges(5,5) * t490 - Icges(5,6) * t489 + Icges(5,3) * t514;
t430 = Icges(5,5) * t488 - Icges(5,6) * t487 + Icges(5,3) * t512;
t429 = rSges(6,1) * t492 - rSges(6,2) * t491 + rSges(6,3) * t504;
t419 = qJD(1) + (t473 * t527 + t474 * t529) * t552;
t414 = rSges(6,1) * t460 - rSges(6,2) * t459 + rSges(6,3) * t489;
t412 = rSges(6,1) * t458 - rSges(6,2) * t457 + rSges(6,3) * t487;
t398 = -t447 * t519 + t480 * t498 + t542;
t397 = t448 * t519 - t480 * t497 + t544;
t396 = t447 * t497 - t448 * t498 + t547;
t395 = -t437 * t499 + t464 * t466 + t539;
t394 = t438 * t499 - t464 * t465 + t541;
t393 = t437 * t465 - t438 * t466 + t543;
t392 = -t412 * t483 + t429 * t450 + t537;
t391 = t414 * t483 - t429 * t449 + t538;
t390 = t412 * t449 - t414 * t450 + t540;
t389 = qJD(6) * t459 + t450 * t553 - t483 * t555 + t537;
t388 = qJD(6) * t457 - t449 * t553 + t483 * t554 + t538;
t387 = qJD(6) * t491 + t449 * t555 - t450 * t554 + t540;
t1 = -t572 * ((-t468 * t564 - t470 * t512 + t472 * t513) * t565 - (-t467 * t564 - t469 * t512 + t471 * t513) * t564 + (-t500 * t564 - t501 * t512 + t502 * t513) * t530) * t564 / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t497 * ((t442 * t514 + t444 * t495 + t446 * t496) * t497 + (t441 * t514 + t443 * t495 + t445 * t496) * t498 + (t477 * t514 + t478 * t495 + t479 * t496) * t519) / 0.2e1 + t498 * ((t442 * t512 + t444 * t493 + t446 * t494) * t497 + (t441 * t512 + t443 * t493 + t445 * t494) * t498 + (t477 * t512 + t478 * t493 + t479 * t494) * t519) / 0.2e1 + t519 * ((-t442 * t560 + t444 * t516 + t446 * t517) * t497 + (-t441 * t560 + t443 * t516 + t445 * t517) * t498 + (-t477 * t560 + t478 * t516 + t479 * t517) * t519) / 0.2e1 + t465 * ((t431 * t514 - t433 * t489 + t435 * t490) * t465 + (t430 * t514 - t432 * t489 + t434 * t490) * t466 + (t461 * t514 - t462 * t489 + t463 * t490) * t499) / 0.2e1 + t466 * ((t431 * t512 - t433 * t487 + t435 * t488) * t465 + (t430 * t512 - t432 * t487 + t434 * t488) * t466 + (t461 * t512 - t462 * t487 + t463 * t488) * t499) / 0.2e1 + t499 * ((-t431 * t560 - t433 * t504 + t435 * t505) * t465 + (-t430 * t560 - t432 * t504 + t434 * t505) * t466 + (-t461 * t560 - t462 * t504 + t463 * t505) * t499) / 0.2e1 + m(4) * (t396 ^ 2 + t397 ^ 2 + t398 ^ 2) / 0.2e1 + m(5) * (t393 ^ 2 + t394 ^ 2 + t395 ^ 2) / 0.2e1 + m(6) * (t390 ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(7) * (t387 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 + m(3) * (t419 ^ 2 + t451 ^ 2 + t452 ^ 2) / 0.2e1 + ((t459 * t575 + t460 * t573 + t489 * t574) * t483 + (t459 * t581 + t577 * t460 + t579 * t489) * t450 + (t459 * t580 + t460 * t576 + t489 * t578) * t449) * t449 / 0.2e1 + ((t457 * t575 + t458 * t573 + t487 * t574) * t483 + (t457 * t581 + t577 * t458 + t579 * t487) * t450 + (t457 * t580 + t458 * t576 + t487 * t578) * t449) * t450 / 0.2e1 + ((t491 * t575 + t492 * t573 + t504 * t574) * t483 + (t491 * t581 + t577 * t492 + t579 * t504) * t450 + (t491 * t580 + t492 * t576 + t504 * t578) * t449) * t483 / 0.2e1 + (((t468 * t565 - t470 * t514 + t472 * t515) * t565 - (t467 * t565 - t469 * t514 + t471 * t515) * t564 + (t500 * t565 - t501 * t514 + t502 * t515) * t530) * t565 + t530 * (t530 ^ 2 * t500 + (((t470 * t535 + t472 * t533) * t527 - (t469 * t535 + t471 * t533) * t529) * t528 + (-t467 * t529 + t468 * t527 + t501 * t535 + t502 * t533) * t530) * t528)) * t572 / 0.2e1;
T  = t1;
