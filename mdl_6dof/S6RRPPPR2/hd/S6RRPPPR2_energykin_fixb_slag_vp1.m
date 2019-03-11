% Calculate kinetic energy for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:09:54
% EndTime: 2019-03-09 08:09:57
% DurationCPUTime: 3.53s
% Computational Cost: add. (1458->284), mult. (1790->419), div. (0->0), fcn. (1691->10), ass. (0->156)
t576 = Icges(4,4) + Icges(5,6);
t575 = Icges(4,1) + Icges(5,2);
t574 = Icges(4,2) + Icges(5,3);
t471 = qJ(2) + pkin(9);
t468 = cos(t471);
t573 = t576 * t468;
t466 = sin(t471);
t572 = t576 * t466;
t571 = -Icges(5,4) + Icges(4,5);
t570 = Icges(5,5) - Icges(4,6);
t569 = t574 * t466 - t573;
t568 = -t575 * t468 + t572;
t477 = sin(qJ(1));
t479 = cos(qJ(1));
t567 = -t569 * t477 + t570 * t479;
t566 = t570 * t477 + t569 * t479;
t565 = t568 * t477 + t571 * t479;
t564 = t571 * t477 - t568 * t479;
t563 = -t574 * t468 - t572;
t562 = t575 * t466 + t573;
t561 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t476 = sin(qJ(2));
t478 = cos(qJ(2));
t560 = Icges(3,5) * t478 - Icges(3,6) * t476 + t570 * t466 + t571 * t468;
t559 = t561 * t477 + t560 * t479;
t558 = t560 * t477 - t561 * t479;
t557 = Icges(3,5) * t476 + Icges(3,6) * t478 + t571 * t466 - t570 * t468;
t544 = Icges(3,4) * t476;
t451 = Icges(3,2) * t478 + t544;
t543 = Icges(3,4) * t478;
t452 = Icges(3,1) * t476 + t543;
t556 = -t451 * t476 + t452 * t478 + t563 * t466 + t562 * t468;
t505 = -Icges(3,2) * t476 + t543;
t422 = Icges(3,6) * t477 + t479 * t505;
t507 = Icges(3,1) * t478 - t544;
t424 = Icges(3,5) * t477 + t479 * t507;
t555 = -t422 * t476 + t424 * t478 + t566 * t466 + t564 * t468;
t421 = -Icges(3,6) * t479 + t477 * t505;
t423 = -Icges(3,5) * t479 + t477 * t507;
t554 = t421 * t476 - t423 * t478 + t567 * t466 + t565 * t468;
t550 = pkin(2) * t476;
t472 = sin(pkin(10));
t549 = pkin(5) * t472;
t547 = pkin(2) * t478;
t473 = cos(pkin(10));
t546 = pkin(5) * t473;
t470 = pkin(10) + qJ(6);
t465 = sin(t470);
t538 = t465 * t477;
t537 = t465 * t479;
t467 = cos(t470);
t536 = t467 * t477;
t535 = t467 * t479;
t534 = t468 * t477;
t533 = t468 * t479;
t532 = t472 * t477;
t531 = t472 * t479;
t530 = t473 * t477;
t529 = t473 * t479;
t394 = -qJ(3) * t479 + t477 * t547;
t395 = qJ(3) * t477 + t479 * t547;
t523 = qJD(2) * t479;
t524 = qJD(2) * t477;
t527 = t394 * t524 + t395 * t523;
t460 = pkin(1) * t477 - pkin(7) * t479;
t526 = -t394 - t460;
t469 = qJD(3) * t477;
t522 = qJD(4) * t466;
t525 = t479 * t522 + t469;
t521 = qJD(5) * t468;
t520 = qJD(6) * t468;
t509 = pkin(3) * t468 + qJ(4) * t466;
t425 = t509 * t477;
t519 = -t425 + t526;
t518 = t479 * t521 + t525;
t515 = -pkin(3) * t466 + qJ(4) * t468 - t550;
t439 = -pkin(4) * t479 + qJ(5) * t534;
t514 = -t439 + t519;
t449 = qJD(1) * (pkin(1) * t479 + pkin(7) * t477);
t513 = qJD(1) * t395 - qJD(3) * t479 + t449;
t512 = rSges(3,1) * t478 - rSges(3,2) * t476;
t511 = rSges(4,1) * t468 - rSges(4,2) * t466;
t510 = -rSges(5,2) * t468 + rSges(5,3) * t466;
t508 = qJD(2) * (-rSges(4,1) * t466 - rSges(4,2) * t468 - t550);
t489 = -qJ(5) * t466 + t515;
t488 = qJD(2) * (rSges(5,2) * t466 + rSges(5,3) * t468 + t515);
t426 = t509 * t479;
t487 = qJD(1) * t426 + t477 * t522 + t513;
t486 = -qJD(4) * t468 + t425 * t524 + t426 * t523 + t527;
t485 = pkin(8) * t468 + t466 * t549;
t484 = qJD(2) * (-rSges(6,3) * t466 - (-rSges(6,1) * t472 - rSges(6,2) * t473) * t468 + t489);
t483 = qJD(2) * (-pkin(8) * t466 + t468 * t549 + t489);
t438 = pkin(4) * t477 + qJ(5) * t533;
t482 = qJD(1) * t438 + t477 * t521 + t487;
t481 = qJD(5) * t466 + t438 * t523 + t439 * t524 + t486;
t461 = qJD(6) * t466 + qJD(1);
t459 = rSges(2,1) * t479 - rSges(2,2) * t477;
t458 = rSges(2,1) * t477 + rSges(2,2) * t479;
t457 = rSges(3,1) * t476 + rSges(3,2) * t478;
t437 = t477 * t520 - t523;
t436 = t479 * t520 + t524;
t434 = t466 * t532 - t529;
t433 = t466 * t530 + t531;
t432 = t466 * t531 + t530;
t431 = t466 * t529 - t532;
t428 = rSges(3,3) * t477 + t479 * t512;
t427 = -rSges(3,3) * t479 + t477 * t512;
t417 = t466 * t538 - t535;
t416 = t466 * t536 + t537;
t415 = t466 * t537 + t536;
t414 = t466 * t535 - t538;
t413 = -rSges(5,1) * t479 + t477 * t510;
t412 = rSges(5,1) * t477 + t479 * t510;
t411 = rSges(4,3) * t477 + t479 * t511;
t410 = -rSges(4,3) * t479 + t477 * t511;
t390 = Icges(6,5) * t466 + (-Icges(6,1) * t472 - Icges(6,4) * t473) * t468;
t389 = Icges(6,6) * t466 + (-Icges(6,4) * t472 - Icges(6,2) * t473) * t468;
t388 = Icges(6,3) * t466 + (-Icges(6,5) * t472 - Icges(6,6) * t473) * t468;
t385 = rSges(7,3) * t466 + (-rSges(7,1) * t465 - rSges(7,2) * t467) * t468;
t384 = Icges(7,5) * t466 + (-Icges(7,1) * t465 - Icges(7,4) * t467) * t468;
t383 = Icges(7,6) * t466 + (-Icges(7,4) * t465 - Icges(7,2) * t467) * t468;
t382 = Icges(7,3) * t466 + (-Icges(7,5) * t465 - Icges(7,6) * t467) * t468;
t381 = t477 * t485 - t479 * t546;
t380 = t477 * t546 + t479 * t485;
t379 = qJD(1) * t428 - t457 * t524 + t449;
t378 = -t457 * t523 + (-t427 - t460) * qJD(1);
t377 = (t427 * t477 + t428 * t479) * qJD(2);
t376 = rSges(6,1) * t434 + rSges(6,2) * t433 + rSges(6,3) * t534;
t375 = rSges(6,1) * t432 + rSges(6,2) * t431 + rSges(6,3) * t533;
t374 = Icges(6,1) * t434 + Icges(6,4) * t433 + Icges(6,5) * t534;
t373 = Icges(6,1) * t432 + Icges(6,4) * t431 + Icges(6,5) * t533;
t372 = Icges(6,4) * t434 + Icges(6,2) * t433 + Icges(6,6) * t534;
t371 = Icges(6,4) * t432 + Icges(6,2) * t431 + Icges(6,6) * t533;
t370 = Icges(6,5) * t434 + Icges(6,6) * t433 + Icges(6,3) * t534;
t369 = Icges(6,5) * t432 + Icges(6,6) * t431 + Icges(6,3) * t533;
t368 = rSges(7,1) * t417 + rSges(7,2) * t416 + rSges(7,3) * t534;
t367 = rSges(7,1) * t415 + rSges(7,2) * t414 + rSges(7,3) * t533;
t366 = Icges(7,1) * t417 + Icges(7,4) * t416 + Icges(7,5) * t534;
t365 = Icges(7,1) * t415 + Icges(7,4) * t414 + Icges(7,5) * t533;
t364 = Icges(7,4) * t417 + Icges(7,2) * t416 + Icges(7,6) * t534;
t363 = Icges(7,4) * t415 + Icges(7,2) * t414 + Icges(7,6) * t533;
t362 = Icges(7,5) * t417 + Icges(7,6) * t416 + Icges(7,3) * t534;
t361 = Icges(7,5) * t415 + Icges(7,6) * t414 + Icges(7,3) * t533;
t360 = qJD(1) * t411 + t477 * t508 + t513;
t359 = t469 + t479 * t508 + (-t410 + t526) * qJD(1);
t358 = (t410 * t477 + t411 * t479) * qJD(2) + t527;
t357 = qJD(1) * t412 + t477 * t488 + t487;
t356 = t479 * t488 + (-t413 + t519) * qJD(1) + t525;
t355 = (t412 * t479 + t413 * t477) * qJD(2) + t486;
t354 = qJD(1) * t375 + t477 * t484 + t482;
t353 = t479 * t484 + (-t376 + t514) * qJD(1) + t518;
t352 = (t375 * t479 + t376 * t477) * qJD(2) + t481;
t351 = qJD(1) * t380 + t367 * t461 - t385 * t436 + t477 * t483 + t482;
t350 = -t368 * t461 + t385 * t437 + t479 * t483 + (-t381 + t514) * qJD(1) + t518;
t349 = -t367 * t437 + t368 * t436 + (t380 * t479 + t381 * t477) * qJD(2) + t481;
t1 = t437 * ((t361 * t534 + t363 * t416 + t365 * t417) * t436 + (t362 * t534 + t416 * t364 + t417 * t366) * t437 + (t382 * t534 + t383 * t416 + t384 * t417) * t461) / 0.2e1 + t461 * ((t361 * t436 + t362 * t437 + t382 * t461) * t466 + ((-t363 * t467 - t365 * t465) * t436 + (-t364 * t467 - t366 * t465) * t437 + (-t383 * t467 - t384 * t465) * t461) * t468) / 0.2e1 + t436 * ((t361 * t533 + t414 * t363 + t415 * t365) * t436 + (t362 * t533 + t364 * t414 + t366 * t415) * t437 + (t382 * t533 + t383 * t414 + t384 * t415) * t461) / 0.2e1 + m(4) * (t358 ^ 2 + t359 ^ 2 + t360 ^ 2) / 0.2e1 + m(3) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(7) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + m(5) * (t355 ^ 2 + t356 ^ 2 + t357 ^ 2) / 0.2e1 + m(6) * (t352 ^ 2 + t353 ^ 2 + t354 ^ 2) / 0.2e1 + (m(2) * (t458 ^ 2 + t459 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t421 * t478 - t423 * t476 + (t372 * t473 + t374 * t472 - t567) * t468 + (-t370 + t565) * t466) * t479 + (t422 * t478 + t424 * t476 + (-t371 * t473 - t373 * t472 - t566) * t468 + (t369 + t564) * t466) * t477) * qJD(2) + (t478 * t451 + t476 * t452 + (-t473 * t389 - t472 * t390 - t563) * t468 + (t388 + t562) * t466) * qJD(1)) * qJD(1) / 0.2e1 + (((-t370 * t533 - t372 * t431 - t374 * t432 + t554 * t479) * t479 + (t369 * t533 + t371 * t431 + t373 * t432 + (t555 - t558) * t479 + t559 * t477) * t477) * qJD(2) + (t388 * t533 + t389 * t431 + t390 * t432 + t477 * t557 + t479 * t556) * qJD(1)) * t524 / 0.2e1 - (((t369 * t534 + t371 * t433 + t373 * t434 + t555 * t477) * t477 + (-t370 * t534 - t372 * t433 - t374 * t434 + (t554 - t559) * t477 + t558 * t479) * t479) * qJD(2) + (t388 * t534 + t389 * t433 + t390 * t434 + t477 * t556 - t479 * t557) * qJD(1)) * t523 / 0.2e1;
T  = t1;
