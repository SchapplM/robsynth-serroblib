% Calculate kinetic energy for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:33
% EndTime: 2019-03-09 16:43:36
% DurationCPUTime: 3.15s
% Computational Cost: add. (1404->240), mult. (1886->358), div. (0->0), fcn. (1797->8), ass. (0->142)
t565 = Icges(4,4) + Icges(5,6);
t564 = Icges(4,1) + Icges(5,2);
t563 = -Icges(4,2) - Icges(5,3);
t460 = qJ(2) + qJ(3);
t459 = cos(t460);
t562 = t565 * t459;
t458 = sin(t460);
t561 = t565 * t458;
t560 = Icges(5,4) - Icges(4,5);
t559 = Icges(5,5) - Icges(4,6);
t558 = t458 * t563 + t562;
t557 = t459 * t564 - t561;
t556 = Icges(5,1) + Icges(4,3);
t555 = Icges(6,1) + Icges(7,1);
t554 = Icges(6,4) - Icges(7,5);
t553 = Icges(7,4) + Icges(6,5);
t552 = Icges(6,2) + Icges(7,3);
t551 = Icges(7,6) - Icges(6,6);
t550 = Icges(6,3) + Icges(7,2);
t463 = sin(qJ(1));
t466 = cos(qJ(1));
t549 = t463 * t558 + t466 * t559;
t548 = -t463 * t559 + t466 * t558;
t547 = t557 * t463 + t466 * t560;
t546 = -t463 * t560 + t557 * t466;
t545 = t459 * t563 - t561;
t544 = t458 * t564 + t562;
t543 = t559 * t458 - t459 * t560;
t542 = rSges(7,1) + pkin(5);
t541 = rSges(7,3) + qJ(6);
t464 = cos(qJ(5));
t508 = t464 * t466;
t461 = sin(qJ(5));
t511 = t461 * t463;
t424 = -t458 * t508 + t511;
t509 = t463 * t464;
t510 = t461 * t466;
t425 = t458 * t510 + t509;
t512 = t459 * t466;
t540 = t552 * t424 - t554 * t425 + t551 * t512;
t426 = t458 * t509 + t510;
t427 = t458 * t511 - t508;
t513 = t459 * t463;
t539 = -t552 * t426 - t554 * t427 + t551 * t513;
t538 = t551 * t424 + t553 * t425 + t550 * t512;
t537 = -t551 * t426 + t553 * t427 + t550 * t513;
t536 = -t554 * t424 + t555 * t425 + t553 * t512;
t535 = t554 * t426 + t555 * t427 + t553 * t513;
t534 = (t554 * t461 + t552 * t464) * t459 + t551 * t458;
t533 = (-t553 * t461 + t551 * t464) * t459 + t550 * t458;
t532 = (-t555 * t461 - t554 * t464) * t459 + t553 * t458;
t457 = qJD(2) * t463;
t443 = qJD(3) * t463 + t457;
t444 = (-qJD(2) - qJD(3)) * t466;
t531 = (-t549 * t458 + t547 * t459) * t444 + (-t548 * t458 + t546 * t459) * t443 + (t545 * t458 + t544 * t459) * qJD(1);
t530 = (t543 * t463 - t556 * t466) * t444 + (t556 * t463 + t543 * t466) * t443 + (-t458 * t560 - t559 * t459) * qJD(1);
t523 = pkin(9) * t458;
t465 = cos(qJ(2));
t521 = pkin(2) * t465;
t462 = sin(qJ(2));
t519 = Icges(3,4) * t462;
t518 = Icges(3,4) * t465;
t507 = rSges(7,2) * t512 + t541 * t424 + t542 * t425;
t506 = rSges(7,2) * t513 - t541 * t426 + t542 * t427;
t392 = -pkin(8) * t466 + t463 * t521;
t393 = pkin(8) * t463 + t466 * t521;
t502 = qJD(2) * t466;
t505 = t392 * t457 + t393 * t502;
t504 = rSges(7,2) * t458 + (-t542 * t461 + t541 * t464) * t459;
t453 = pkin(1) * t463 - pkin(7) * t466;
t503 = -t392 - t453;
t501 = qJD(4) * t458;
t500 = qJD(4) * t459;
t499 = qJD(5) * t459;
t498 = pkin(2) * qJD(2) * t462;
t491 = pkin(3) * t459 + qJ(4) * t458;
t421 = t491 * t463;
t497 = t443 * t421 + t505;
t496 = -t421 + t503;
t495 = t466 * t498;
t494 = rSges(3,1) * t465 - rSges(3,2) * t462;
t493 = rSges(4,1) * t459 - rSges(4,2) * t458;
t492 = -rSges(5,2) * t459 + rSges(5,3) * t458;
t490 = Icges(3,1) * t465 - t519;
t488 = -Icges(3,2) * t462 + t518;
t485 = Icges(3,5) * t465 - Icges(3,6) * t462;
t413 = -Icges(3,6) * t466 + t463 * t488;
t415 = -Icges(3,5) * t466 + t463 * t490;
t481 = t413 * t462 - t415 * t465;
t414 = Icges(3,6) * t463 + t466 * t488;
t416 = Icges(3,5) * t463 + t466 * t490;
t480 = -t414 * t462 + t416 * t465;
t446 = Icges(3,2) * t465 + t519;
t447 = Icges(3,1) * t462 + t518;
t479 = -t446 * t462 + t447 * t465;
t441 = qJD(1) * (pkin(1) * t466 + pkin(7) * t463);
t478 = qJD(1) * t393 - t463 * t498 + t441;
t438 = pkin(3) * t458 - qJ(4) * t459;
t477 = t444 * t438 + t466 * t501 - t495;
t423 = t491 * t466;
t430 = pkin(4) * t463 + pkin(9) * t512;
t431 = -pkin(4) * t466 + pkin(9) * t513;
t474 = t443 * t431 + (-t423 - t430) * t444 + t497;
t473 = qJD(1) * t423 + t463 * t501 + t478;
t472 = t444 * t523 + (-t431 + t496) * qJD(1) + t477;
t471 = qJD(1) * t430 + (-t438 - t523) * t443 + t473;
t454 = qJD(5) * t458 + qJD(1);
t450 = rSges(2,1) * t466 - rSges(2,2) * t463;
t449 = rSges(2,1) * t463 + rSges(2,2) * t466;
t448 = rSges(3,1) * t462 + rSges(3,2) * t465;
t445 = Icges(3,5) * t462 + Icges(3,6) * t465;
t440 = rSges(4,1) * t458 + rSges(4,2) * t459;
t439 = -rSges(5,2) * t458 - rSges(5,3) * t459;
t420 = rSges(3,3) * t463 + t466 * t494;
t419 = -rSges(3,3) * t466 + t463 * t494;
t418 = t463 * t499 + t444;
t417 = t466 * t499 + t443;
t412 = Icges(3,3) * t463 + t466 * t485;
t411 = -Icges(3,3) * t466 + t463 * t485;
t409 = -rSges(5,1) * t466 + t463 * t492;
t408 = rSges(5,1) * t463 + t466 * t492;
t407 = rSges(4,3) * t463 + t466 * t493;
t406 = -rSges(4,3) * t466 + t463 * t493;
t389 = rSges(6,3) * t458 + (-rSges(6,1) * t461 - rSges(6,2) * t464) * t459;
t375 = qJD(1) * t420 - t448 * t457 + t441;
t374 = -t448 * t502 + (-t419 - t453) * qJD(1);
t373 = rSges(6,1) * t427 + rSges(6,2) * t426 + rSges(6,3) * t513;
t371 = rSges(6,1) * t425 - rSges(6,2) * t424 + rSges(6,3) * t512;
t357 = (t419 * t463 + t420 * t466) * qJD(2);
t356 = qJD(1) * t407 - t440 * t443 + t478;
t355 = -t495 + t440 * t444 + (-t406 + t503) * qJD(1);
t354 = t406 * t443 - t407 * t444 + t505;
t353 = qJD(1) * t408 + (-t438 - t439) * t443 + t473;
t352 = t439 * t444 + (-t409 + t496) * qJD(1) + t477;
t351 = -t500 + t409 * t443 + (-t408 - t423) * t444 + t497;
t350 = t371 * t454 - t389 * t417 + t471;
t349 = -t373 * t454 + t389 * t418 + t472;
t348 = -t371 * t418 + t373 * t417 + t474 - t500;
t347 = -qJD(6) * t426 - t417 * t504 + t454 * t507 + t471;
t346 = qJD(6) * t424 + t418 * t504 - t454 * t506 + t472;
t345 = (qJD(6) * t464 - qJD(4)) * t459 - t507 * t418 + t506 * t417 + t474;
t1 = ((t463 * t445 + t466 * t479) * qJD(1) + (t463 ^ 2 * t412 + (t481 * t466 + (-t411 + t480) * t463) * t466) * qJD(2)) * t457 / 0.2e1 - ((-t466 * t445 + t463 * t479) * qJD(1) + (t466 ^ 2 * t411 + (t480 * t463 + (-t412 + t481) * t466) * t463) * qJD(2)) * t502 / 0.2e1 + m(7) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + m(6) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + m(5) * (t351 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + m(4) * (t354 ^ 2 + t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(3) * (t357 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + ((t424 * t534 + t425 * t532 + t512 * t533) * t454 + (t424 * t539 + t425 * t535 + t512 * t537) * t418 + (t424 * t540 + t536 * t425 + t538 * t512) * t417) * t417 / 0.2e1 + ((-t426 * t534 + t427 * t532 + t513 * t533) * t454 + (-t426 * t539 + t427 * t535 + t513 * t537) * t418 + (-t426 * t540 + t536 * t427 + t538 * t513) * t417) * t418 / 0.2e1 + (t463 * t530 + t466 * t531) * t443 / 0.2e1 + (t463 * t531 - t466 * t530) * t444 / 0.2e1 + (((-t461 * t532 + t464 * t534) * t454 + (-t461 * t535 + t464 * t539) * t418 + (-t536 * t461 + t464 * t540) * t417) * t459 + (t417 * t538 + t418 * t537 + t454 * t533) * t458) * t454 / 0.2e1 + (Icges(2,3) + m(2) * (t449 ^ 2 + t450 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t465 * t414 + t462 * t416) * t463 - (t413 * t465 + t415 * t462) * t466) * qJD(2) + (t547 * t458 + t549 * t459) * t444 + (t546 * t458 + t548 * t459) * t443 + (t465 * t446 + t462 * t447 + t544 * t458 - t545 * t459) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
