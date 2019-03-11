% Calculate kinetic energy for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:10
% EndTime: 2019-03-08 19:51:12
% DurationCPUTime: 2.73s
% Computational Cost: add. (1720->285), mult. (4425->427), div. (0->0), fcn. (5224->10), ass. (0->137)
t566 = Icges(3,1) + Icges(4,2);
t565 = Icges(5,1) + Icges(6,2);
t564 = Icges(6,1) + Icges(5,3);
t563 = Icges(3,4) + Icges(4,6);
t562 = Icges(5,4) + Icges(6,6);
t561 = Icges(6,4) - Icges(5,5);
t560 = Icges(3,5) - Icges(4,4);
t559 = Icges(6,5) - Icges(5,6);
t558 = Icges(3,2) + Icges(4,3);
t557 = Icges(5,2) + Icges(6,3);
t556 = Icges(3,6) - Icges(4,5);
t555 = Icges(3,3) + Icges(4,1);
t500 = cos(pkin(6));
t498 = sin(pkin(6));
t499 = cos(pkin(10));
t527 = t498 * t499;
t497 = sin(pkin(10));
t528 = t497 * t498;
t502 = sin(qJ(2));
t504 = cos(qJ(2));
t536 = t555 * t500 + (t502 * t560 + t504 * t556) * t498;
t524 = t500 * t504;
t482 = t497 * t502 - t499 * t524;
t525 = t500 * t502;
t483 = t497 * t504 + t499 * t525;
t542 = -t482 * t556 + t483 * t560 - t527 * t555;
t484 = t497 * t524 + t499 * t502;
t485 = -t497 * t525 + t499 * t504;
t543 = -t484 * t556 + t485 * t560 + t528 * t555;
t554 = -t536 * t500 + t527 * t542 - t528 * t543;
t529 = sin(qJ(4));
t518 = t498 * t529;
t530 = cos(qJ(4));
t458 = t482 * t530 + t499 * t518;
t519 = t498 * t530;
t459 = -t482 * t529 + t499 * t519;
t553 = -t458 * t557 + t459 * t562 + t483 * t559;
t456 = -t484 * t530 + t497 * t518;
t457 = t484 * t529 + t497 * t519;
t552 = t456 * t557 - t457 * t562 + t485 * t559;
t551 = t456 * t559 - t457 * t561 + t485 * t564;
t550 = -t458 * t559 + t459 * t561 + t483 * t564;
t549 = -t458 * t562 + t459 * t565 + t483 * t561;
t548 = t456 * t562 - t457 * t565 + t485 * t561;
t547 = t484 * t558 - t485 * t563 - t528 * t556;
t546 = t482 * t558 - t483 * t563 + t527 * t556;
t545 = -t484 * t563 + t485 * t566 + t528 * t560;
t544 = t482 * t563 - t483 * t566 + t527 * t560;
t486 = t500 * t529 + t504 * t519;
t487 = t500 * t530 - t504 * t518;
t526 = t498 * t502;
t541 = t486 * t557 - t487 * t562 + t526 * t559;
t540 = t486 * t562 - t487 * t565 + t526 * t561;
t539 = t486 * t559 - t487 * t561 + t526 * t564;
t538 = t556 * t500 + (t502 * t563 + t504 * t558) * t498;
t537 = t560 * t500 + (t502 * t566 + t504 * t563) * t498;
t534 = qJD(2) ^ 2;
t448 = pkin(2) * t485 + qJ(3) * t484;
t496 = qJD(2) * t500;
t523 = qJD(3) * t482 + t448 * t496;
t522 = qJD(2) * t498;
t494 = t497 * t522;
t462 = qJD(4) * t485 + t494;
t489 = qJD(4) * t526 + t496;
t521 = qJD(3) * t504;
t517 = t499 * t522;
t447 = pkin(2) * t483 + qJ(3) * t482;
t516 = t447 * t494 + t448 * t517 + qJD(1);
t488 = (pkin(2) * t502 - qJ(3) * t504) * t498;
t514 = (-t500 * rSges(4,1) - (-rSges(4,2) * t502 - rSges(4,3) * t504) * t498 - t488) * t498;
t513 = (-pkin(3) * t500 - pkin(8) * t526 - t488) * t498;
t463 = qJD(4) * t483 - t517;
t464 = pkin(3) * t528 + pkin(8) * t485;
t465 = -pkin(3) * t527 + pkin(8) * t483;
t510 = t464 * t517 + t465 * t494 - t498 * t521 + t516;
t509 = qJD(2) * t497 * t513 + t464 * t496 + t523;
t411 = -pkin(4) * t459 - qJ(5) * t458;
t508 = qJD(5) * t486 + t411 * t462 + t510;
t410 = pkin(4) * t457 + qJ(5) * t456;
t507 = -qJD(5) * t458 + t410 * t489 + t509;
t480 = qJD(3) * t484;
t506 = t480 + ((-t447 - t465) * t500 + t499 * t513) * qJD(2);
t449 = pkin(4) * t487 + qJ(5) * t486;
t505 = qJD(5) * t456 + t449 * t463 + t506;
t503 = cos(qJ(6));
t501 = sin(qJ(6));
t473 = t500 * rSges(3,3) + (rSges(3,1) * t502 + rSges(3,2) * t504) * t498;
t466 = pkin(5) * t526 + pkin(9) * t487;
t461 = t486 * t501 + t503 * t526;
t460 = t486 * t503 - t501 * t526;
t454 = qJD(6) * t487 + t489;
t445 = rSges(5,1) * t487 - rSges(5,2) * t486 + rSges(5,3) * t526;
t444 = rSges(6,1) * t526 - rSges(6,2) * t487 + rSges(6,3) * t486;
t437 = rSges(3,1) * t485 - rSges(3,2) * t484 + rSges(3,3) * t528;
t436 = rSges(3,1) * t483 - rSges(3,2) * t482 - rSges(3,3) * t527;
t435 = -rSges(4,1) * t527 - rSges(4,2) * t483 + rSges(4,3) * t482;
t434 = rSges(4,1) * t528 - rSges(4,2) * t485 + rSges(4,3) * t484;
t419 = pkin(5) * t485 + pkin(9) * t457;
t418 = pkin(5) * t483 - pkin(9) * t459;
t417 = t456 * t501 + t485 * t503;
t416 = t456 * t503 - t485 * t501;
t415 = -t458 * t501 + t483 * t503;
t414 = -t458 * t503 - t483 * t501;
t413 = -qJD(6) * t459 + t463;
t412 = qJD(6) * t457 + t462;
t407 = (-t436 * t500 - t473 * t527) * qJD(2);
t406 = (t437 * t500 - t473 * t528) * qJD(2);
t405 = rSges(7,1) * t461 + rSges(7,2) * t460 + rSges(7,3) * t487;
t404 = Icges(7,1) * t461 + Icges(7,4) * t460 + Icges(7,5) * t487;
t403 = Icges(7,4) * t461 + Icges(7,2) * t460 + Icges(7,6) * t487;
t402 = Icges(7,5) * t461 + Icges(7,6) * t460 + Icges(7,3) * t487;
t401 = rSges(6,1) * t485 - rSges(6,2) * t457 + rSges(6,3) * t456;
t400 = rSges(6,1) * t483 + rSges(6,2) * t459 - rSges(6,3) * t458;
t399 = -rSges(5,1) * t459 + rSges(5,2) * t458 + rSges(5,3) * t483;
t398 = rSges(5,1) * t457 - rSges(5,2) * t456 + rSges(5,3) * t485;
t384 = qJD(1) + (t436 * t497 + t437 * t499) * t522;
t383 = rSges(7,1) * t417 + rSges(7,2) * t416 + rSges(7,3) * t457;
t382 = rSges(7,1) * t415 + rSges(7,2) * t414 - rSges(7,3) * t459;
t381 = Icges(7,1) * t417 + Icges(7,4) * t416 + Icges(7,5) * t457;
t380 = Icges(7,1) * t415 + Icges(7,4) * t414 - Icges(7,5) * t459;
t379 = Icges(7,4) * t417 + Icges(7,2) * t416 + Icges(7,6) * t457;
t378 = Icges(7,4) * t415 + Icges(7,2) * t414 - Icges(7,6) * t459;
t377 = Icges(7,5) * t417 + Icges(7,6) * t416 + Icges(7,3) * t457;
t376 = Icges(7,5) * t415 + Icges(7,6) * t414 - Icges(7,3) * t459;
t375 = t480 + ((-t435 - t447) * t500 + t499 * t514) * qJD(2);
t374 = (t434 * t500 + t497 * t514) * qJD(2) + t523;
t373 = (-t521 + (t434 * t499 + t435 * t497) * qJD(2)) * t498 + t516;
t372 = -t399 * t489 + t445 * t463 + t506;
t371 = t398 * t489 - t445 * t462 + t509;
t370 = -t398 * t463 + t399 * t462 + t510;
t369 = t444 * t463 + (-t400 - t411) * t489 + t505;
t368 = t401 * t489 + (-t444 - t449) * t462 + t507;
t367 = t462 * t400 + (-t401 - t410) * t463 + t508;
t366 = -t382 * t454 + t405 * t413 + t463 * t466 + (-t411 - t418) * t489 + t505;
t365 = t383 * t454 - t405 * t412 + t419 * t489 + (-t449 - t466) * t462 + t507;
t364 = t412 * t382 - t413 * t383 + t462 * t418 + (-t410 - t419) * t463 + t508;
t1 = m(7) * (t364 ^ 2 + t365 ^ 2 + t366 ^ 2) / 0.2e1 + m(6) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(5) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(3) * (t384 ^ 2 + t406 ^ 2 + t407 ^ 2) / 0.2e1 + m(4) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t412 * ((t377 * t457 + t379 * t416 + t381 * t417) * t412 + (t376 * t457 + t378 * t416 + t380 * t417) * t413 + (t402 * t457 + t403 * t416 + t404 * t417) * t454) / 0.2e1 + t413 * ((-t377 * t459 + t379 * t414 + t381 * t415) * t412 + (-t376 * t459 + t378 * t414 + t380 * t415) * t413 + (-t402 * t459 + t403 * t414 + t404 * t415) * t454) / 0.2e1 + t454 * ((t377 * t487 + t379 * t460 + t381 * t461) * t412 + (t376 * t487 + t378 * t460 + t380 * t461) * t413 + (t402 * t487 + t403 * t460 + t404 * t461) * t454) / 0.2e1 + ((t456 * t541 - t457 * t540 + t485 * t539) * t489 + (t456 * t553 - t457 * t549 + t485 * t550) * t463 + (t552 * t456 - t548 * t457 + t551 * t485) * t462) * t462 / 0.2e1 + ((-t458 * t541 + t459 * t540 + t483 * t539) * t489 + (-t553 * t458 + t549 * t459 + t550 * t483) * t463 + (-t458 * t552 + t459 * t548 + t483 * t551) * t462) * t463 / 0.2e1 + ((t541 * t486 - t540 * t487 + t539 * t526) * t489 + (t486 * t553 - t487 * t549 + t526 * t550) * t463 + (t486 * t552 - t487 * t548 + t526 * t551) * t462) * t489 / 0.2e1 - ((t482 * t547 + t483 * t545) * t528 + (-t482 * t538 + t483 * t537) * t500 + (-t482 * t546 + t483 * t544 + t554) * t527) * t534 * t527 / 0.2e1 + ((t536 * t500 ^ 2 + (((t502 * t544 + t504 * t546) * t499 + (t502 * t545 - t504 * t547) * t497) * t498 + (t497 * t543 - t499 * t542 + t502 * t537 + t504 * t538) * t500) * t498) * t500 + ((-t484 * t546 + t485 * t544) * t527 + (-t484 * t538 + t485 * t537) * t500 + (t484 * t547 + t485 * t545 - t554) * t528) * t528) * t534 / 0.2e1;
T  = t1;
