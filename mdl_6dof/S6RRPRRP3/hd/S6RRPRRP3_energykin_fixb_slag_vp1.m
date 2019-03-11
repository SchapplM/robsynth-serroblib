% Calculate kinetic energy for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:18
% EndTime: 2019-03-09 11:47:20
% DurationCPUTime: 2.74s
% Computational Cost: add. (2037->288), mult. (2234->438), div. (0->0), fcn. (2191->10), ass. (0->158)
t581 = Icges(6,1) + Icges(7,1);
t580 = Icges(6,4) + Icges(7,4);
t579 = -Icges(7,5) - Icges(6,5);
t578 = Icges(6,2) + Icges(7,2);
t577 = -Icges(7,6) - Icges(6,6);
t576 = Icges(3,3) + Icges(4,3);
t575 = -Icges(7,3) - Icges(6,3);
t477 = qJ(2) + pkin(10);
t469 = sin(t477);
t470 = cos(t477);
t481 = sin(qJ(2));
t484 = cos(qJ(2));
t574 = Icges(3,5) * t484 + Icges(4,5) * t470 - Icges(3,6) * t481 - Icges(4,6) * t469;
t485 = cos(qJ(1));
t483 = cos(qJ(4));
t546 = t483 * pkin(4);
t492 = pkin(9) * t469 + t470 * t546;
t480 = sin(qJ(4));
t482 = sin(qJ(1));
t531 = t480 * t482;
t383 = pkin(4) * t531 + t485 * t492;
t397 = -pkin(9) * t470 + t469 * t546;
t472 = qJD(2) * t482;
t519 = qJD(4) * t469;
t447 = t485 * t519 + t472;
t466 = -qJD(4) * t470 + qJD(1);
t573 = t466 * t383 - t397 * t447;
t478 = qJ(4) + qJ(5);
t474 = cos(t478);
t532 = t474 * t485;
t473 = sin(t478);
t535 = t473 * t482;
t427 = -t470 * t535 - t532;
t533 = t474 * t482;
t534 = t473 * t485;
t428 = t470 * t533 - t534;
t537 = t469 * t482;
t572 = -t577 * t427 - t579 * t428 - t575 * t537;
t429 = -t470 * t534 + t533;
t430 = t470 * t532 + t535;
t536 = t469 * t485;
t571 = -t577 * t429 - t579 * t430 - t575 * t536;
t570 = t578 * t427 + t580 * t428 - t577 * t537;
t569 = t578 * t429 + t580 * t430 - t577 * t536;
t568 = t580 * t427 + t581 * t428 - t579 * t537;
t567 = t580 * t429 + t581 * t430 - t579 * t536;
t530 = t480 * t485;
t382 = -pkin(4) * t530 + t482 * t492;
t520 = qJD(2) * t485;
t448 = t482 * t519 - t520;
t566 = -t382 * t466 + t448 * t397;
t565 = t575 * t470 + (t577 * t473 - t579 * t474) * t469;
t564 = t577 * t470 + (-t578 * t473 + t580 * t474) * t469;
t563 = t579 * t470 + (-t580 * t473 + t581 * t474) * t469;
t562 = t574 * t482 - t576 * t485;
t561 = t576 * t482 + t574 * t485;
t560 = Icges(3,5) * t481 + Icges(4,5) * t469 + Icges(3,6) * t484 + Icges(4,6) * t470;
t541 = Icges(4,4) * t469;
t450 = Icges(4,2) * t470 + t541;
t540 = Icges(4,4) * t470;
t451 = Icges(4,1) * t469 + t540;
t543 = Icges(3,4) * t481;
t458 = Icges(3,2) * t484 + t543;
t542 = Icges(3,4) * t484;
t459 = Icges(3,1) * t481 + t542;
t559 = -t450 * t469 + t451 * t470 - t458 * t481 + t459 * t484;
t504 = -Icges(4,2) * t469 + t540;
t420 = Icges(4,6) * t482 + t485 * t504;
t506 = Icges(4,1) * t470 - t541;
t422 = Icges(4,5) * t482 + t485 * t506;
t505 = -Icges(3,2) * t481 + t542;
t435 = Icges(3,6) * t482 + t485 * t505;
t507 = Icges(3,1) * t484 - t543;
t437 = Icges(3,5) * t482 + t485 * t507;
t558 = -t420 * t469 + t422 * t470 - t435 * t481 + t437 * t484;
t419 = -Icges(4,6) * t485 + t482 * t504;
t421 = -Icges(4,5) * t485 + t482 * t506;
t434 = -Icges(3,6) * t485 + t482 * t505;
t436 = -Icges(3,5) * t485 + t482 * t507;
t557 = t419 * t469 - t421 * t470 + t434 * t481 - t436 * t484;
t550 = pkin(2) * t481;
t547 = pkin(2) * t484;
t529 = t482 * t483;
t528 = t483 * t485;
t522 = pkin(5) * t474;
t490 = qJ(6) * t469 + t470 * t522;
t515 = pkin(5) * t473;
t527 = rSges(7,1) * t428 + rSges(7,2) * t427 + rSges(7,3) * t537 + t482 * t490 - t485 * t515;
t526 = rSges(7,1) * t430 + rSges(7,2) * t429 + rSges(7,3) * t536 + t482 * t515 + t485 * t490;
t525 = (-qJ(6) - rSges(7,3)) * t470 + (rSges(7,1) * t474 - rSges(7,2) * t473 + t522) * t469;
t415 = -qJ(3) * t485 + t482 * t547;
t416 = qJ(3) * t482 + t485 * t547;
t524 = t415 * t472 + t416 * t520;
t465 = pkin(1) * t482 - pkin(7) * t485;
t523 = -t415 - t465;
t518 = qJD(5) * t469;
t513 = pkin(3) * t470 + pkin(8) * t469;
t440 = t513 * t482;
t441 = t513 * t485;
t514 = t440 * t472 + t441 * t520 + t524;
t455 = qJD(1) * (pkin(1) * t485 + pkin(7) * t482);
t512 = qJD(1) * t416 - qJD(3) * t485 + t455;
t511 = rSges(3,1) * t484 - rSges(3,2) * t481;
t510 = rSges(4,1) * t470 - rSges(4,2) * t469;
t509 = qJD(2) * (-rSges(4,1) * t469 - rSges(4,2) * t470 - t550);
t508 = (-pkin(3) * t469 + pkin(8) * t470 - t550) * qJD(2);
t495 = qJD(1) * t441 + t512;
t471 = qJD(3) * t482;
t494 = t471 + (-t440 + t523) * qJD(1);
t493 = t447 * t382 - t383 * t448 + t514;
t491 = qJD(6) * t469 + t508;
t489 = t482 * t508 + t495;
t488 = t485 * t508 + t494;
t464 = rSges(2,1) * t485 - rSges(2,2) * t482;
t463 = rSges(2,1) * t482 + rSges(2,2) * t485;
t462 = rSges(3,1) * t481 + rSges(3,2) * t484;
t446 = qJD(1) + (-qJD(4) - qJD(5)) * t470;
t445 = t470 * t528 + t531;
t444 = -t470 * t530 + t529;
t443 = t470 * t529 - t530;
t442 = -t470 * t531 - t528;
t439 = rSges(3,3) * t482 + t485 * t511;
t438 = -rSges(3,3) * t485 + t482 * t511;
t424 = rSges(4,3) * t482 + t485 * t510;
t423 = -rSges(4,3) * t485 + t482 * t510;
t414 = t482 * t518 + t448;
t413 = t485 * t518 + t447;
t412 = -rSges(5,3) * t470 + (rSges(5,1) * t483 - rSges(5,2) * t480) * t469;
t410 = -Icges(5,5) * t470 + (Icges(5,1) * t483 - Icges(5,4) * t480) * t469;
t409 = -Icges(5,6) * t470 + (Icges(5,4) * t483 - Icges(5,2) * t480) * t469;
t408 = -Icges(5,3) * t470 + (Icges(5,5) * t483 - Icges(5,6) * t480) * t469;
t405 = -rSges(6,3) * t470 + (rSges(6,1) * t474 - rSges(6,2) * t473) * t469;
t395 = qJD(1) * t439 - t462 * t472 + t455;
t394 = -t462 * t520 + (-t438 - t465) * qJD(1);
t393 = (t438 * t482 + t439 * t485) * qJD(2);
t392 = rSges(5,1) * t445 + rSges(5,2) * t444 + rSges(5,3) * t536;
t391 = rSges(5,1) * t443 + rSges(5,2) * t442 + rSges(5,3) * t537;
t389 = Icges(5,1) * t445 + Icges(5,4) * t444 + Icges(5,5) * t536;
t388 = Icges(5,1) * t443 + Icges(5,4) * t442 + Icges(5,5) * t537;
t387 = Icges(5,4) * t445 + Icges(5,2) * t444 + Icges(5,6) * t536;
t386 = Icges(5,4) * t443 + Icges(5,2) * t442 + Icges(5,6) * t537;
t385 = Icges(5,5) * t445 + Icges(5,6) * t444 + Icges(5,3) * t536;
t384 = Icges(5,5) * t443 + Icges(5,6) * t442 + Icges(5,3) * t537;
t381 = rSges(6,1) * t430 + rSges(6,2) * t429 + rSges(6,3) * t536;
t379 = rSges(6,1) * t428 + rSges(6,2) * t427 + rSges(6,3) * t537;
t361 = qJD(1) * t424 + t482 * t509 + t512;
t360 = t471 + t485 * t509 + (-t423 + t523) * qJD(1);
t359 = (t423 * t482 + t424 * t485) * qJD(2) + t524;
t358 = t392 * t466 - t412 * t447 + t489;
t357 = -t391 * t466 + t412 * t448 + t488;
t356 = t391 * t447 - t392 * t448 + t514;
t355 = t381 * t446 - t405 * t413 + t489 + t573;
t354 = -t379 * t446 + t405 * t414 + t488 + t566;
t353 = t379 * t413 - t381 * t414 + t493;
t352 = -t413 * t525 + t446 * t526 + t482 * t491 + t495 + t573;
t351 = t414 * t525 - t446 * t527 + t485 * t491 + t494 + t566;
t350 = -qJD(6) * t470 + t413 * t527 - t414 * t526 + t493;
t1 = m(6) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + m(7) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + m(5) * (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + m(4) * (t359 ^ 2 + t360 ^ 2 + t361 ^ 2) / 0.2e1 + m(3) * (t393 ^ 2 + t394 ^ 2 + t395 ^ 2) / 0.2e1 + t448 * ((t385 * t537 + t387 * t442 + t389 * t443) * t447 + (t384 * t537 + t442 * t386 + t443 * t388) * t448 + (t408 * t537 + t409 * t442 + t410 * t443) * t466) / 0.2e1 + t466 * ((-t384 * t448 - t385 * t447 - t408 * t466) * t470 + ((-t387 * t480 + t389 * t483) * t447 + (-t386 * t480 + t388 * t483) * t448 + (-t409 * t480 + t410 * t483) * t466) * t469) / 0.2e1 + t447 * ((t385 * t536 + t444 * t387 + t445 * t389) * t447 + (t384 * t536 + t386 * t444 + t388 * t445) * t448 + (t408 * t536 + t409 * t444 + t410 * t445) * t466) / 0.2e1 + ((t429 * t564 + t430 * t563 + t536 * t565) * t446 + (t429 * t570 + t430 * t568 + t536 * t572) * t414 + (t569 * t429 + t567 * t430 + t571 * t536) * t413) * t413 / 0.2e1 + ((t427 * t564 + t428 * t563 + t537 * t565) * t446 + (t570 * t427 + t568 * t428 + t572 * t537) * t414 + (t427 * t569 + t428 * t567 + t537 * t571) * t413) * t414 / 0.2e1 + ((-t413 * t571 - t414 * t572 - t446 * t565) * t470 + ((-t473 * t564 + t474 * t563) * t446 + (-t473 * t570 + t474 * t568) * t414 + (-t473 * t569 + t474 * t567) * t413) * t469) * t446 / 0.2e1 + (m(2) * (t463 ^ 2 + t464 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t419 * t470 - t421 * t469 - t434 * t484 - t436 * t481) * t485 + (t420 * t470 + t422 * t469 + t435 * t484 + t437 * t481) * t482) * qJD(2) + (t470 * t450 + t469 * t451 + t484 * t458 + t481 * t459) * qJD(1)) * qJD(1) / 0.2e1 + ((t561 * t482 ^ 2 + (t557 * t485 + (t558 - t562) * t482) * t485) * qJD(2) + (t482 * t560 + t485 * t559) * qJD(1)) * t472 / 0.2e1 - ((t562 * t485 ^ 2 + (t558 * t482 + (t557 - t561) * t485) * t482) * qJD(2) + (t482 * t559 - t485 * t560) * qJD(1)) * t520 / 0.2e1;
T  = t1;
