% Calculate kinetic energy for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:47:14
% EndTime: 2019-03-10 01:47:17
% DurationCPUTime: 3.09s
% Computational Cost: add. (3290->333), mult. (5990->511), div. (0->0), fcn. (7231->12), ass. (0->158)
t594 = Icges(6,1) + Icges(7,1);
t593 = -Icges(6,4) + Icges(7,5);
t592 = Icges(7,4) + Icges(6,5);
t591 = Icges(6,2) + Icges(7,3);
t590 = Icges(7,2) + Icges(6,3);
t589 = -Icges(6,6) + Icges(7,6);
t588 = rSges(7,1) + pkin(5);
t587 = rSges(7,3) + qJ(6);
t533 = sin(qJ(2));
t534 = sin(qJ(1));
t536 = cos(qJ(2));
t537 = cos(qJ(1));
t569 = cos(pkin(6));
t551 = t537 * t569;
t512 = t533 * t534 - t536 * t551;
t513 = t533 * t551 + t534 * t536;
t530 = sin(pkin(6));
t564 = t530 * t537;
t471 = Icges(3,5) * t513 - Icges(3,6) * t512 - Icges(3,3) * t564;
t552 = t534 * t569;
t514 = t537 * t533 + t536 * t552;
t515 = -t533 * t552 + t537 * t536;
t567 = t530 * t534;
t472 = Icges(3,5) * t515 - Icges(3,6) * t514 + Icges(3,3) * t567;
t586 = t530 * (t471 * t537 - t472 * t534);
t563 = qJ(3) + qJ(4);
t529 = sin(t563);
t553 = cos(t563);
t488 = t513 * t553 - t529 * t564;
t531 = sin(qJ(5));
t572 = cos(qJ(5));
t455 = t488 * t531 - t512 * t572;
t456 = t488 * t572 + t512 * t531;
t549 = t530 * t553;
t487 = t513 * t529 + t537 * t549;
t585 = t591 * t455 + t593 * t456 + t589 * t487;
t490 = t515 * t553 + t529 * t567;
t457 = t490 * t531 - t514 * t572;
t458 = t490 * t572 + t514 * t531;
t489 = t515 * t529 - t534 * t549;
t584 = t591 * t457 + t593 * t458 + t589 * t489;
t583 = t589 * t455 + t592 * t456 + t590 * t487;
t582 = t589 * t457 + t592 * t458 + t590 * t489;
t581 = t593 * t455 + t594 * t456 + t592 * t487;
t580 = t593 * t457 + t594 * t458 + t592 * t489;
t503 = t529 * t569 + t533 * t549;
t565 = t530 * t536;
t485 = t503 * t531 + t565 * t572;
t486 = t503 * t572 - t531 * t565;
t568 = t530 * t533;
t502 = t529 * t568 - t553 * t569;
t579 = t591 * t485 + t593 * t486 + t589 * t502;
t578 = t589 * t485 + t592 * t486 + t590 * t502;
t577 = t593 * t485 + t594 * t486 + t592 * t502;
t535 = cos(qJ(3));
t571 = pkin(3) * t535;
t566 = t530 * t535;
t562 = rSges(7,2) * t487 + t587 * t455 + t456 * t588;
t561 = rSges(7,2) * t489 + t587 * t457 + t458 * t588;
t560 = rSges(7,2) * t502 + t587 * t485 + t486 * t588;
t483 = pkin(2) * t513 + pkin(9) * t512;
t484 = pkin(2) * t515 + pkin(9) * t514;
t557 = qJD(2) * t530;
t525 = t534 * t557;
t554 = t537 * t557;
t559 = t483 * t525 + t484 * t554;
t495 = qJD(3) * t514 + t525;
t558 = qJD(1) * (pkin(1) * t534 - pkin(8) * t564);
t526 = qJD(2) * t569 + qJD(1);
t532 = sin(qJ(3));
t556 = t532 * t567;
t555 = t532 * t564;
t467 = qJD(4) * t514 + t495;
t550 = t569 * t532;
t496 = qJD(3) * t512 - t554;
t439 = -pkin(3) * t555 + pkin(10) * t512 + t513 * t571;
t440 = pkin(3) * t556 + pkin(10) * t514 + t515 * t571;
t547 = t495 * t439 - t440 * t496 + t559;
t468 = qJD(4) * t512 + t496;
t516 = (pkin(2) * t533 - pkin(9) * t536) * t530;
t518 = qJD(1) * (pkin(1) * t537 + pkin(8) * t567);
t546 = t526 * t484 - t516 * t525 + t518;
t497 = (-qJD(3) - qJD(4)) * t565 + t526;
t453 = pkin(4) * t488 + pkin(11) * t487;
t454 = pkin(4) * t490 + pkin(11) * t489;
t545 = t467 * t453 - t454 * t468 + t547;
t544 = -t483 * t526 - t516 * t554 - t558;
t477 = pkin(3) * t550 + (-pkin(10) * t536 + t533 * t571) * t530;
t517 = -qJD(3) * t565 + t526;
t543 = t517 * t440 - t477 * t495 + t546;
t542 = -t439 * t517 + t496 * t477 + t544;
t470 = pkin(4) * t503 + pkin(11) * t502;
t541 = t497 * t454 - t467 * t470 + t543;
t540 = -t453 * t497 + t468 * t470 + t542;
t522 = rSges(2,1) * t537 - rSges(2,2) * t534;
t521 = rSges(2,1) * t534 + rSges(2,2) * t537;
t511 = t533 * t566 + t550;
t510 = -t532 * t568 + t535 * t569;
t501 = t569 * rSges(3,3) + (rSges(3,1) * t533 + rSges(3,2) * t536) * t530;
t500 = Icges(3,5) * t569 + (Icges(3,1) * t533 + Icges(3,4) * t536) * t530;
t499 = Icges(3,6) * t569 + (Icges(3,4) * t533 + Icges(3,2) * t536) * t530;
t498 = Icges(3,3) * t569 + (Icges(3,5) * t533 + Icges(3,6) * t536) * t530;
t494 = t515 * t535 + t556;
t493 = -t515 * t532 + t534 * t566;
t492 = t513 * t535 - t555;
t491 = -t513 * t532 - t535 * t564;
t480 = rSges(3,1) * t515 - rSges(3,2) * t514 + rSges(3,3) * t567;
t479 = rSges(3,1) * t513 - rSges(3,2) * t512 - rSges(3,3) * t564;
t476 = Icges(3,1) * t515 - Icges(3,4) * t514 + Icges(3,5) * t567;
t475 = Icges(3,1) * t513 - Icges(3,4) * t512 - Icges(3,5) * t564;
t474 = Icges(3,4) * t515 - Icges(3,2) * t514 + Icges(3,6) * t567;
t473 = Icges(3,4) * t513 - Icges(3,2) * t512 - Icges(3,6) * t564;
t469 = rSges(4,1) * t511 + rSges(4,2) * t510 - rSges(4,3) * t565;
t466 = Icges(4,1) * t511 + Icges(4,4) * t510 - Icges(4,5) * t565;
t465 = Icges(4,4) * t511 + Icges(4,2) * t510 - Icges(4,6) * t565;
t464 = Icges(4,5) * t511 + Icges(4,6) * t510 - Icges(4,3) * t565;
t463 = qJD(5) * t502 + t497;
t462 = rSges(5,1) * t503 - rSges(5,2) * t502 - rSges(5,3) * t565;
t461 = Icges(5,1) * t503 - Icges(5,4) * t502 - Icges(5,5) * t565;
t460 = Icges(5,4) * t503 - Icges(5,2) * t502 - Icges(5,6) * t565;
t459 = Icges(5,5) * t503 - Icges(5,6) * t502 - Icges(5,3) * t565;
t450 = qJD(5) * t487 + t468;
t449 = qJD(5) * t489 + t467;
t448 = rSges(4,1) * t494 + rSges(4,2) * t493 + rSges(4,3) * t514;
t447 = rSges(4,1) * t492 + rSges(4,2) * t491 + rSges(4,3) * t512;
t446 = Icges(4,1) * t494 + Icges(4,4) * t493 + Icges(4,5) * t514;
t445 = Icges(4,1) * t492 + Icges(4,4) * t491 + Icges(4,5) * t512;
t444 = Icges(4,4) * t494 + Icges(4,2) * t493 + Icges(4,6) * t514;
t443 = Icges(4,4) * t492 + Icges(4,2) * t491 + Icges(4,6) * t512;
t442 = Icges(4,5) * t494 + Icges(4,6) * t493 + Icges(4,3) * t514;
t441 = Icges(4,5) * t492 + Icges(4,6) * t491 + Icges(4,3) * t512;
t438 = rSges(5,1) * t490 - rSges(5,2) * t489 + rSges(5,3) * t514;
t437 = rSges(5,1) * t488 - rSges(5,2) * t487 + rSges(5,3) * t512;
t436 = Icges(5,1) * t490 - Icges(5,4) * t489 + Icges(5,5) * t514;
t435 = Icges(5,1) * t488 - Icges(5,4) * t487 + Icges(5,5) * t512;
t434 = Icges(5,4) * t490 - Icges(5,2) * t489 + Icges(5,6) * t514;
t433 = Icges(5,4) * t488 - Icges(5,2) * t487 + Icges(5,6) * t512;
t432 = Icges(5,5) * t490 - Icges(5,6) * t489 + Icges(5,3) * t514;
t431 = Icges(5,5) * t488 - Icges(5,6) * t487 + Icges(5,3) * t512;
t429 = rSges(6,1) * t486 - rSges(6,2) * t485 + rSges(6,3) * t502;
t420 = t480 * t526 - t501 * t525 + t518;
t419 = -t479 * t526 - t501 * t554 - t558;
t417 = (t479 * t534 + t480 * t537) * t557;
t412 = rSges(6,1) * t458 - rSges(6,2) * t457 + rSges(6,3) * t489;
t410 = rSges(6,1) * t456 - rSges(6,2) * t455 + rSges(6,3) * t487;
t396 = t448 * t517 - t469 * t495 + t546;
t395 = -t447 * t517 + t469 * t496 + t544;
t394 = t447 * t495 - t448 * t496 + t559;
t393 = t438 * t497 - t462 * t467 + t543;
t392 = -t437 * t497 + t462 * t468 + t542;
t391 = t437 * t467 - t438 * t468 + t547;
t390 = t412 * t463 - t429 * t449 + t541;
t389 = -t410 * t463 + t429 * t450 + t540;
t388 = t410 * t449 - t412 * t450 + t545;
t387 = qJD(6) * t455 - t449 * t560 + t463 * t561 + t541;
t386 = qJD(6) * t457 + t450 * t560 - t463 * t562 + t540;
t385 = qJD(6) * t485 + t449 * t562 - t450 * t561 + t545;
t1 = -((-t498 * t564 - t499 * t512 + t500 * t513) * t526 + ((-t474 * t512 + t476 * t513) * t534 + (t512 * t473 - t513 * t475 + t586) * t537) * t557) * t554 / 0.2e1 + ((t498 * t567 - t499 * t514 + t500 * t515) * t526 + (-(-t473 * t514 + t475 * t515) * t537 + (-t514 * t474 + t515 * t476 - t586) * t534) * t557) * t525 / 0.2e1 + t496 * ((t442 * t512 + t444 * t491 + t446 * t492) * t495 + (t441 * t512 + t443 * t491 + t445 * t492) * t496 + (t464 * t512 + t465 * t491 + t466 * t492) * t517) / 0.2e1 + t517 * ((-t442 * t565 + t444 * t510 + t446 * t511) * t495 + (-t441 * t565 + t443 * t510 + t445 * t511) * t496 + (-t464 * t565 + t465 * t510 + t466 * t511) * t517) / 0.2e1 + t495 * ((t442 * t514 + t444 * t493 + t446 * t494) * t495 + (t441 * t514 + t443 * t493 + t445 * t494) * t496 + (t464 * t514 + t465 * t493 + t466 * t494) * t517) / 0.2e1 + m(4) * (t394 ^ 2 + t395 ^ 2 + t396 ^ 2) / 0.2e1 + m(3) * (t417 ^ 2 + t419 ^ 2 + t420 ^ 2) / 0.2e1 + m(7) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(6) * (t388 ^ 2 + t389 ^ 2 + t390 ^ 2) / 0.2e1 + m(5) * (t391 ^ 2 + t392 ^ 2 + t393 ^ 2) / 0.2e1 + t467 * ((t514 * t432 - t489 * t434 + t490 * t436) * t467 + (t431 * t514 - t433 * t489 + t435 * t490) * t468 + (t459 * t514 - t460 * t489 + t461 * t490) * t497) / 0.2e1 + t468 * ((t432 * t512 - t434 * t487 + t436 * t488) * t467 + (t512 * t431 - t487 * t433 + t488 * t435) * t468 + (t459 * t512 - t460 * t487 + t461 * t488) * t497) / 0.2e1 + t497 * ((-t432 * t565 - t434 * t502 + t436 * t503) * t467 + (-t431 * t565 - t433 * t502 + t435 * t503) * t468 + (-t459 * t565 - t460 * t502 + t461 * t503) * t497) / 0.2e1 + t526 * ((t569 * t472 + (t474 * t536 + t476 * t533) * t530) * t525 - (t569 * t471 + (t473 * t536 + t475 * t533) * t530) * t554 + (t569 * t498 + (t499 * t536 + t500 * t533) * t530) * t526) / 0.2e1 + ((t457 * t579 + t458 * t577 + t489 * t578) * t463 + (t457 * t585 + t458 * t581 + t489 * t583) * t450 + (t584 * t457 + t580 * t458 + t582 * t489) * t449) * t449 / 0.2e1 + ((t455 * t579 + t456 * t577 + t487 * t578) * t463 + (t585 * t455 + t581 * t456 + t583 * t487) * t450 + (t455 * t584 + t456 * t580 + t487 * t582) * t449) * t450 / 0.2e1 + ((t579 * t485 + t577 * t486 + t578 * t502) * t463 + (t485 * t585 + t486 * t581 + t502 * t583) * t450 + (t485 * t584 + t486 * t580 + t502 * t582) * t449) * t463 / 0.2e1 + (m(2) * (t521 ^ 2 + t522 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
