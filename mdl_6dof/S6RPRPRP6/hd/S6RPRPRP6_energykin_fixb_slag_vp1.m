% Calculate kinetic energy for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:17:57
% EndTime: 2019-03-09 03:17:59
% DurationCPUTime: 2.51s
% Computational Cost: add. (1236->208), mult. (1578->314), div. (0->0), fcn. (1512->8), ass. (0->120)
t550 = Icges(4,4) + Icges(5,6);
t549 = Icges(4,1) + Icges(5,2);
t548 = -Icges(4,2) - Icges(5,3);
t443 = pkin(9) + qJ(3);
t441 = cos(t443);
t547 = t550 * t441;
t440 = sin(t443);
t546 = t550 * t440;
t545 = -Icges(5,4) + Icges(4,5);
t544 = Icges(5,5) - Icges(4,6);
t543 = t548 * t440 + t547;
t542 = -t549 * t441 + t546;
t541 = Icges(5,1) + Icges(4,3);
t540 = Icges(6,1) + Icges(7,1);
t539 = Icges(6,4) + Icges(7,4);
t538 = Icges(7,5) + Icges(6,5);
t537 = Icges(6,2) + Icges(7,2);
t536 = Icges(7,6) + Icges(6,6);
t535 = Icges(7,3) + Icges(6,3);
t449 = sin(qJ(1));
t451 = cos(qJ(1));
t534 = t449 * t543 + t451 * t544;
t533 = -t449 * t544 + t451 * t543;
t532 = t542 * t449 + t451 * t545;
t531 = t449 * t545 - t542 * t451;
t530 = t548 * t441 - t546;
t529 = t549 * t440 + t547;
t528 = t544 * t440 + t441 * t545;
t450 = cos(qJ(5));
t492 = t450 * t451;
t448 = sin(qJ(5));
t495 = t448 * t449;
t413 = t440 * t492 - t495;
t493 = t449 * t450;
t494 = t448 * t451;
t414 = t440 * t494 + t493;
t496 = t441 * t451;
t527 = t413 * t536 + t414 * t538 + t496 * t535;
t415 = t440 * t493 + t494;
t416 = t440 * t495 - t492;
t497 = t441 * t449;
t526 = t415 * t536 + t416 * t538 + t497 * t535;
t525 = t537 * t413 + t414 * t539 + t536 * t496;
t524 = t537 * t415 + t416 * t539 + t536 * t497;
t523 = t539 * t413 + t540 * t414 + t538 * t496;
t522 = t539 * t415 + t540 * t416 + t538 * t497;
t521 = (-t448 * t538 - t450 * t536) * t441 + t535 * t440;
t520 = (-t448 * t539 - t537 * t450) * t441 + t536 * t440;
t519 = (-t540 * t448 - t539 * t450) * t441 + t538 * t440;
t518 = t528 * t449 - t541 * t451;
t517 = t541 * t449 + t528 * t451;
t516 = t440 * t545 - t544 * t441;
t515 = t440 * t530 + t441 * t529;
t514 = -t440 * t533 + t441 * t531;
t513 = t440 * t534 + t441 * t532;
t506 = pkin(5) * t448;
t445 = cos(pkin(9));
t504 = pkin(2) * t445;
t503 = pkin(5) * t450;
t456 = qJ(6) * t441 + t440 * t506;
t490 = rSges(7,1) * t414 + rSges(7,2) * t413 + rSges(7,3) * t496 + t449 * t503 + t451 * t456;
t489 = rSges(7,1) * t416 + rSges(7,2) * t415 + rSges(7,3) * t497 + t449 * t456 - t451 * t503;
t488 = (-rSges(7,1) * t448 - rSges(7,2) * t450 - t506) * t441 + (rSges(7,3) + qJ(6)) * t440;
t434 = pkin(1) * t449 - qJ(2) * t451;
t487 = pkin(7) * t451 - t449 * t504 - t434;
t442 = qJD(2) * t449;
t483 = qJD(4) * t440;
t486 = t451 * t483 + t442;
t485 = qJD(3) * t449;
t484 = qJD(3) * t451;
t482 = qJD(5) * t441;
t472 = pkin(3) * t441 + qJ(4) * t440;
t409 = t472 * t449;
t481 = -t409 + t487;
t426 = pkin(3) * t440 - qJ(4) * t441;
t478 = qJD(3) * (rSges(5,2) * t440 + rSges(5,3) * t441 - t426);
t431 = qJD(1) * (pkin(1) * t451 + qJ(2) * t449);
t477 = -qJD(2) * t451 + qJD(1) * (pkin(7) * t449 + t451 * t504) + t431;
t410 = t472 * t451;
t476 = -qJD(4) * t441 + t409 * t485 + t410 * t484;
t444 = sin(pkin(9));
t475 = rSges(3,1) * t445 - rSges(3,2) * t444;
t474 = rSges(4,1) * t441 - rSges(4,2) * t440;
t473 = -rSges(5,2) * t441 + rSges(5,3) * t440;
t471 = (-pkin(8) * t440 - t426) * qJD(3);
t458 = qJD(1) * t410 + t449 * t483 + t477;
t429 = pkin(4) * t449 + pkin(8) * t496;
t430 = -pkin(4) * t451 + pkin(8) * t497;
t457 = t429 * t484 + t430 * t485 + t476;
t455 = qJD(1) * t429 + t458;
t454 = qJD(6) * t441 + t471;
t453 = (-t430 + t481) * qJD(1) + t486;
t437 = qJD(5) * t440 + qJD(1);
t436 = rSges(2,1) * t451 - rSges(2,2) * t449;
t435 = rSges(2,1) * t449 + rSges(2,2) * t451;
t428 = rSges(4,1) * t440 + rSges(4,2) * t441;
t419 = t449 * t482 - t484;
t418 = t451 * t482 + t485;
t407 = -rSges(5,1) * t451 + t449 * t473;
t406 = rSges(5,1) * t449 + t451 * t473;
t405 = rSges(4,3) * t449 + t451 * t474;
t404 = -rSges(4,3) * t451 + t449 * t474;
t387 = rSges(6,3) * t440 + (-rSges(6,1) * t448 - rSges(6,2) * t450) * t441;
t378 = qJD(1) * t449 * rSges(3,3) + t431 + (qJD(1) * t475 - qJD(2)) * t451;
t377 = t442 + (t451 * rSges(3,3) - t449 * t475 - t434) * qJD(1);
t374 = rSges(6,1) * t416 + rSges(6,2) * t415 + rSges(6,3) * t497;
t372 = rSges(6,1) * t414 + rSges(6,2) * t413 + rSges(6,3) * t496;
t358 = (t404 * t449 + t405 * t451) * qJD(3);
t357 = qJD(1) * t405 - t428 * t485 + t477;
t356 = -t428 * t484 + t442 + (-t404 + t487) * qJD(1);
t355 = (t406 * t451 + t407 * t449) * qJD(3) + t476;
t354 = qJD(1) * t406 + t449 * t478 + t458;
t353 = t451 * t478 + (-t407 + t481) * qJD(1) + t486;
t352 = t372 * t437 - t387 * t418 + t449 * t471 + t455;
t351 = -t374 * t437 + t387 * t419 + t451 * t471 + t453;
t350 = -t372 * t419 + t374 * t418 + t457;
t349 = -t418 * t488 + t437 * t490 + t449 * t454 + t455;
t348 = t419 * t488 - t437 * t489 + t451 * t454 + t453;
t347 = qJD(6) * t440 + t418 * t489 - t419 * t490 + t457;
t1 = m(5) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + m(3) * (t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(4) * (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + m(7) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + m(6) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + ((t520 * t413 + t519 * t414 + t521 * t496) * t437 + (t524 * t413 + t522 * t414 + t526 * t496) * t419 + (t525 * t413 + t523 * t414 + t527 * t496) * t418) * t418 / 0.2e1 + ((t520 * t415 + t519 * t416 + t521 * t497) * t437 + (t524 * t415 + t522 * t416 + t526 * t497) * t419 + (t525 * t415 + t523 * t416 + t527 * t497) * t418) * t419 / 0.2e1 + (((-t519 * t448 - t520 * t450) * t437 + (-t522 * t448 - t524 * t450) * t419 + (-t523 * t448 - t525 * t450) * t418) * t441 + (t527 * t418 + t526 * t419 + t521 * t437) * t440) * t437 / 0.2e1 + (((t440 * t532 - t441 * t534) * t451 + (t440 * t531 + t441 * t533) * t449) * qJD(3) + (t440 * t529 - t441 * t530) * qJD(1)) * qJD(1) / 0.2e1 + ((t517 * t449 ^ 2 + (t513 * t451 + (t514 - t518) * t449) * t451) * qJD(3) + (t516 * t449 + t515 * t451) * qJD(1)) * t485 / 0.2e1 - ((t518 * t451 ^ 2 + (t514 * t449 + (t513 - t517) * t451) * t449) * qJD(3) + (t515 * t449 - t516 * t451) * qJD(1)) * t484 / 0.2e1 + (Icges(3,2) * t445 ^ 2 + (Icges(3,1) * t444 + 0.2e1 * Icges(3,4) * t445) * t444 + Icges(2,3) + m(2) * (t435 ^ 2 + t436 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
