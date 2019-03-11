% Calculate kinetic energy for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:40
% EndTime: 2019-03-08 19:17:43
% DurationCPUTime: 3.04s
% Computational Cost: add. (2638->321), mult. (6734->489), div. (0->0), fcn. (8398->12), ass. (0->152)
t577 = Icges(4,1) + Icges(5,2);
t576 = Icges(4,4) + Icges(5,6);
t575 = Icges(5,4) - Icges(4,5);
t574 = Icges(5,5) - Icges(4,6);
t573 = Icges(4,2) + Icges(5,3);
t515 = sin(pkin(6));
t516 = cos(pkin(10));
t552 = t515 * t516;
t514 = sin(pkin(10));
t553 = t514 * t515;
t520 = sin(qJ(2));
t522 = cos(qJ(2));
t554 = sin(pkin(11));
t555 = cos(pkin(11));
t506 = -t520 * t554 + t522 * t555;
t498 = t506 * t515;
t528 = t520 * t555 + t522 * t554;
t499 = t528 * t515;
t517 = cos(pkin(6));
t571 = Icges(5,1) + Icges(4,3) + Icges(3,3);
t559 = ((Icges(3,5) * t520 + Icges(3,6) * t522) * t515 - t575 * t499 - t574 * t498 + t571 * t517) * t517;
t526 = t517 * t506;
t483 = -t514 * t526 - t516 * t528;
t525 = t517 * t528;
t484 = t506 * t516 - t514 * t525;
t549 = t517 * t522;
t503 = -t514 * t549 - t516 * t520;
t550 = t517 * t520;
t504 = -t514 * t550 + t516 * t522;
t561 = Icges(3,5) * t504 + Icges(3,6) * t503 - t483 * t574 - t484 * t575 + t553 * t571;
t481 = -t514 * t528 + t516 * t526;
t482 = t514 * t506 + t516 * t525;
t501 = -t514 * t520 + t516 * t549;
t502 = t514 * t522 + t516 * t550;
t562 = -Icges(3,5) * t502 - Icges(3,6) * t501 + t481 * t574 + t482 * t575 + t552 * t571;
t572 = -t552 * t562 - t553 * t561 - t559;
t569 = t481 * t573 + t482 * t576 + t552 * t574;
t568 = -t483 * t573 - t484 * t576 + t553 * t574;
t567 = -t481 * t576 - t482 * t577 - t552 * t575;
t566 = t483 * t576 + t484 * t577 - t553 * t575;
t432 = pkin(3) * t484 - qJ(4) * t483;
t513 = qJD(2) * t517;
t565 = -qJD(4) * t481 + t432 * t513;
t564 = -t498 * t573 - t499 * t576 + t517 * t574;
t563 = t498 * t576 + t499 * t577 - t517 * t575;
t558 = qJD(2) ^ 2;
t557 = cos(qJ(5));
t556 = pkin(2) * t522;
t519 = sin(qJ(5));
t551 = t515 * t519;
t431 = pkin(3) * t482 - qJ(4) * t481;
t538 = pkin(2) * t550 - qJ(3) * t515;
t474 = t514 * t556 + t516 * t538;
t548 = -t431 - t474;
t507 = pkin(2) * t515 * t520 + qJ(3) * t517;
t547 = -pkin(3) * t499 + qJ(4) * t498 - t507;
t509 = qJD(3) * t553;
t546 = -qJD(4) * t483 + t509;
t545 = qJD(2) * t515;
t510 = t514 * t545;
t450 = qJD(5) * t484 + t510;
t487 = qJD(5) * t499 + t513;
t544 = qJD(3) * t516;
t542 = -pkin(4) * t517 - pkin(8) * t499 + t547;
t541 = t515 * t557;
t540 = t516 * t545;
t535 = (-rSges(4,1) * t499 - rSges(4,2) * t498 - rSges(4,3) * t517 - t507) * t515;
t475 = -t514 * t538 + t516 * t556;
t534 = qJD(3) * t517 + t474 * t510 + t475 * t540 + qJD(1);
t531 = (-rSges(5,1) * t517 + rSges(5,2) * t499 + rSges(5,3) * t498 + t547) * t515;
t451 = qJD(5) * t482 - t540;
t456 = t475 * t513;
t530 = -t515 * t544 + t456;
t529 = -qJD(4) * t498 + t431 * t510 + t432 * t540 + t534;
t454 = pkin(4) * t553 + pkin(8) * t484;
t455 = -pkin(4) * t552 + pkin(8) * t482;
t527 = t454 * t540 + t455 * t510 + t529;
t524 = t454 * t513 + t456 + (qJD(2) * t514 * t542 - t544) * t515 + t565;
t523 = ((-t455 + t548) * t517 + t542 * t552) * qJD(2) + t546;
t521 = cos(qJ(6));
t518 = sin(qJ(6));
t497 = t517 * rSges(3,3) + (rSges(3,1) * t520 + rSges(3,2) * t522) * t515;
t496 = Icges(3,5) * t517 + (Icges(3,1) * t520 + Icges(3,4) * t522) * t515;
t495 = Icges(3,6) * t517 + (Icges(3,4) * t520 + Icges(3,2) * t522) * t515;
t486 = -t498 * t519 + t517 * t557;
t485 = t498 * t557 + t517 * t519;
t472 = rSges(3,1) * t504 + rSges(3,2) * t503 + rSges(3,3) * t553;
t471 = rSges(3,1) * t502 + rSges(3,2) * t501 - rSges(3,3) * t552;
t470 = Icges(3,1) * t504 + Icges(3,4) * t503 + Icges(3,5) * t553;
t469 = Icges(3,1) * t502 + Icges(3,4) * t501 - Icges(3,5) * t552;
t468 = Icges(3,4) * t504 + Icges(3,2) * t503 + Icges(3,6) * t553;
t467 = Icges(3,4) * t502 + Icges(3,2) * t501 - Icges(3,6) * t552;
t449 = -t481 * t519 - t516 * t541;
t448 = -t481 * t557 + t516 * t551;
t447 = -t483 * t519 + t514 * t541;
t446 = t483 * t557 + t514 * t551;
t444 = t486 * t521 + t499 * t518;
t443 = -t486 * t518 + t499 * t521;
t440 = qJD(6) * t485 + t487;
t439 = pkin(5) * t486 + pkin(9) * t485;
t438 = (-t471 * t517 - t497 * t552) * qJD(2);
t437 = (t472 * t517 - t497 * t553) * qJD(2);
t436 = rSges(6,1) * t486 - rSges(6,2) * t485 + rSges(6,3) * t499;
t435 = Icges(6,1) * t486 - Icges(6,4) * t485 + Icges(6,5) * t499;
t434 = Icges(6,4) * t486 - Icges(6,2) * t485 + Icges(6,6) * t499;
t433 = Icges(6,5) * t486 - Icges(6,6) * t485 + Icges(6,3) * t499;
t429 = rSges(4,1) * t484 + rSges(4,2) * t483 + rSges(4,3) * t553;
t428 = rSges(4,1) * t482 + rSges(4,2) * t481 - rSges(4,3) * t552;
t427 = -rSges(5,1) * t552 - rSges(5,2) * t482 - rSges(5,3) * t481;
t426 = rSges(5,1) * t553 - rSges(5,2) * t484 - rSges(5,3) * t483;
t411 = t449 * t521 + t482 * t518;
t410 = -t449 * t518 + t482 * t521;
t409 = t447 * t521 + t484 * t518;
t408 = -t447 * t518 + t484 * t521;
t407 = qJD(1) + (t471 * t514 + t472 * t516) * t545;
t406 = -qJD(6) * t448 + t451;
t405 = qJD(6) * t446 + t450;
t404 = pkin(5) * t449 - pkin(9) * t448;
t403 = pkin(5) * t447 + pkin(9) * t446;
t402 = rSges(7,1) * t444 + rSges(7,2) * t443 + rSges(7,3) * t485;
t401 = Icges(7,1) * t444 + Icges(7,4) * t443 + Icges(7,5) * t485;
t400 = Icges(7,4) * t444 + Icges(7,2) * t443 + Icges(7,6) * t485;
t399 = Icges(7,5) * t444 + Icges(7,6) * t443 + Icges(7,3) * t485;
t398 = rSges(6,1) * t449 + rSges(6,2) * t448 + rSges(6,3) * t482;
t397 = rSges(6,1) * t447 - rSges(6,2) * t446 + rSges(6,3) * t484;
t396 = Icges(6,1) * t449 + Icges(6,4) * t448 + Icges(6,5) * t482;
t395 = Icges(6,1) * t447 - Icges(6,4) * t446 + Icges(6,5) * t484;
t394 = Icges(6,4) * t449 + Icges(6,2) * t448 + Icges(6,6) * t482;
t393 = Icges(6,4) * t447 - Icges(6,2) * t446 + Icges(6,6) * t484;
t392 = Icges(6,5) * t449 + Icges(6,6) * t448 + Icges(6,3) * t482;
t391 = Icges(6,5) * t447 - Icges(6,6) * t446 + Icges(6,3) * t484;
t390 = t509 + ((-t428 - t474) * t517 + t516 * t535) * qJD(2);
t389 = (t429 * t517 + t514 * t535) * qJD(2) + t530;
t388 = rSges(7,1) * t411 + rSges(7,2) * t410 - rSges(7,3) * t448;
t387 = rSges(7,1) * t409 + rSges(7,2) * t408 + rSges(7,3) * t446;
t386 = Icges(7,1) * t411 + Icges(7,4) * t410 - Icges(7,5) * t448;
t385 = Icges(7,1) * t409 + Icges(7,4) * t408 + Icges(7,5) * t446;
t384 = Icges(7,4) * t411 + Icges(7,2) * t410 - Icges(7,6) * t448;
t383 = Icges(7,4) * t409 + Icges(7,2) * t408 + Icges(7,6) * t446;
t382 = Icges(7,5) * t411 + Icges(7,6) * t410 - Icges(7,3) * t448;
t381 = Icges(7,5) * t409 + Icges(7,6) * t408 + Icges(7,3) * t446;
t380 = (t428 * t514 + t429 * t516) * t545 + t534;
t379 = ((-t427 + t548) * t517 + t516 * t531) * qJD(2) + t546;
t378 = (t426 * t517 + t514 * t531) * qJD(2) + t530 + t565;
t377 = (t426 * t516 + t427 * t514) * t545 + t529;
t376 = -t398 * t487 + t436 * t451 + t523;
t375 = t397 * t487 - t436 * t450 + t524;
t374 = -t397 * t451 + t398 * t450 + t527;
t373 = -t388 * t440 + t402 * t406 - t404 * t487 + t439 * t451 + t523;
t372 = t387 * t440 - t402 * t405 + t403 * t487 - t439 * t450 + t524;
t371 = -t387 * t406 + t388 * t405 - t403 * t451 + t404 * t450 + t527;
t1 = t405 * ((t446 * t381 + t408 * t383 + t409 * t385) * t405 + (t382 * t446 + t384 * t408 + t386 * t409) * t406 + (t399 * t446 + t400 * t408 + t401 * t409) * t440) / 0.2e1 + t406 * ((-t381 * t448 + t383 * t410 + t385 * t411) * t405 + (-t448 * t382 + t410 * t384 + t411 * t386) * t406 + (-t399 * t448 + t400 * t410 + t401 * t411) * t440) / 0.2e1 + t440 * ((t381 * t485 + t383 * t443 + t385 * t444) * t405 + (t382 * t485 + t384 * t443 + t386 * t444) * t406 + (t485 * t399 + t443 * t400 + t444 * t401) * t440) / 0.2e1 + t451 * ((t391 * t482 + t393 * t448 + t395 * t449) * t450 + (t482 * t392 + t448 * t394 + t449 * t396) * t451 + (t433 * t482 + t434 * t448 + t435 * t449) * t487) / 0.2e1 + t487 * ((t391 * t499 - t393 * t485 + t395 * t486) * t450 + (t392 * t499 - t394 * t485 + t396 * t486) * t451 + (t433 * t499 - t434 * t485 + t435 * t486) * t487) / 0.2e1 + t450 * ((t484 * t391 - t446 * t393 + t447 * t395) * t450 + (t392 * t484 - t394 * t446 + t396 * t447) * t451 + (t433 * t484 - t434 * t446 + t435 * t447) * t487) / 0.2e1 + m(4) * (t380 ^ 2 + t389 ^ 2 + t390 ^ 2) / 0.2e1 + m(5) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(6) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(3) * (t407 ^ 2 + t437 ^ 2 + t438 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(7) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 - ((t468 * t501 + t470 * t502 - t481 * t568 + t482 * t566) * t553 + (-t481 * t564 + t482 * t563 + t495 * t501 + t496 * t502) * t517 + (-t467 * t501 - t469 * t502 - t481 * t569 + t482 * t567 + t572) * t552) * t558 * t552 / 0.2e1 + (((-t564 * t498 + t563 * t499 + t559) * t517 + ((t495 * t522 + t496 * t520) * t517 + (-(t467 * t522 + t469 * t520) * t515 + t567 * t499 - t569 * t498 + t562 * t517) * t516 + ((t468 * t522 + t470 * t520) * t515 + t566 * t499 - t568 * t498 + t561 * t517) * t514) * t515) * t517 + ((-t467 * t503 - t469 * t504 - t483 * t569 + t484 * t567) * t552 + (-t483 * t564 + t484 * t563 + t495 * t503 + t496 * t504) * t517 + (t468 * t503 + t470 * t504 - t483 * t568 + t484 * t566 - t572) * t553) * t553) * t558 / 0.2e1;
T  = t1;
