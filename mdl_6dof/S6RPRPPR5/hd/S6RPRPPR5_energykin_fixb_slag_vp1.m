% Calculate kinetic energy for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:02
% EndTime: 2019-03-09 02:50:05
% DurationCPUTime: 2.96s
% Computational Cost: add. (1334->255), mult. (1543->384), div. (0->0), fcn. (1477->10), ass. (0->138)
t547 = Icges(4,4) + Icges(5,6);
t546 = Icges(4,1) + Icges(5,2);
t545 = Icges(4,2) + Icges(5,3);
t454 = pkin(9) + qJ(3);
t451 = cos(t454);
t544 = t547 * t451;
t449 = sin(t454);
t543 = t547 * t449;
t542 = -Icges(5,4) + Icges(4,5);
t541 = Icges(5,5) - Icges(4,6);
t540 = t545 * t449 - t544;
t539 = -t546 * t451 + t543;
t538 = Icges(5,1) + Icges(4,3);
t461 = sin(qJ(1));
t462 = cos(qJ(1));
t537 = -t540 * t461 + t541 * t462;
t536 = t541 * t461 + t540 * t462;
t535 = t542 * t461 - t539 * t462;
t534 = t539 * t461 + t542 * t462;
t533 = -t545 * t451 - t543;
t532 = t546 * t449 + t544;
t531 = t541 * t449 + t542 * t451;
t530 = t531 * t461 - t462 * t538;
t529 = t461 * t538 + t531 * t462;
t528 = t542 * t449 - t541 * t451;
t527 = t449 * t533 + t451 * t532;
t526 = t449 * t537 + t451 * t534;
t525 = t449 * t536 + t451 * t535;
t455 = sin(pkin(10));
t521 = pkin(5) * t455;
t458 = cos(pkin(9));
t519 = pkin(2) * t458;
t457 = cos(pkin(10));
t518 = pkin(5) * t457;
t453 = pkin(10) + qJ(6);
t448 = sin(t453);
t513 = t448 * t461;
t512 = t448 * t462;
t450 = cos(t453);
t511 = t450 * t461;
t510 = t450 * t462;
t509 = t451 * t461;
t508 = t451 * t462;
t507 = t455 * t461;
t506 = t455 * t462;
t505 = t457 * t461;
t504 = t457 * t462;
t441 = pkin(1) * t461 - qJ(2) * t462;
t501 = pkin(7) * t462 - t519 * t461 - t441;
t452 = qJD(2) * t461;
t497 = qJD(4) * t449;
t500 = t462 * t497 + t452;
t499 = qJD(3) * t461;
t498 = qJD(3) * t462;
t496 = qJD(5) * t451;
t495 = qJD(6) * t451;
t482 = pkin(3) * t451 + qJ(4) * t449;
t414 = t482 * t461;
t494 = -t414 + t501;
t493 = t462 * t496 + t500;
t433 = pkin(3) * t449 - qJ(4) * t451;
t490 = -qJ(5) * t449 - t433;
t489 = qJD(3) * (rSges(5,2) * t449 + rSges(5,3) * t451 - t433);
t426 = -pkin(4) * t462 + qJ(5) * t509;
t488 = -t426 + t494;
t436 = qJD(1) * (pkin(1) * t462 + qJ(2) * t461);
t487 = -qJD(2) * t462 + qJD(1) * (pkin(7) * t461 + t519 * t462) + t436;
t415 = t482 * t462;
t486 = -qJD(4) * t451 + t414 * t499 + t415 * t498;
t456 = sin(pkin(9));
t485 = rSges(3,1) * t458 - rSges(3,2) * t456;
t484 = rSges(4,1) * t451 - rSges(4,2) * t449;
t483 = -rSges(5,2) * t451 + rSges(5,3) * t449;
t469 = qJD(3) * (-rSges(6,3) * t449 - (-rSges(6,1) * t455 - rSges(6,2) * t457) * t451 + t490);
t468 = qJD(3) * (-pkin(8) * t449 + t451 * t521 + t490);
t467 = qJD(1) * t415 + t461 * t497 + t487;
t466 = pkin(8) * t451 + t449 * t521;
t425 = pkin(4) * t461 + qJ(5) * t508;
t465 = qJD(5) * t449 + t425 * t498 + t426 * t499 + t486;
t464 = qJD(1) * t425 + t461 * t496 + t467;
t444 = qJD(6) * t449 + qJD(1);
t443 = rSges(2,1) * t462 - rSges(2,2) * t461;
t442 = rSges(2,1) * t461 + rSges(2,2) * t462;
t435 = rSges(4,1) * t449 + rSges(4,2) * t451;
t424 = t461 * t495 - t498;
t423 = t462 * t495 + t499;
t421 = t449 * t507 - t504;
t420 = t449 * t505 + t506;
t419 = t449 * t506 + t505;
t418 = t449 * t504 - t507;
t412 = t449 * t513 - t510;
t411 = t449 * t511 + t512;
t410 = t449 * t512 + t511;
t409 = t449 * t510 - t513;
t408 = -rSges(5,1) * t462 + t483 * t461;
t407 = rSges(5,1) * t461 + t483 * t462;
t406 = rSges(4,3) * t461 + t484 * t462;
t405 = -rSges(4,3) * t462 + t484 * t461;
t386 = Icges(6,5) * t449 + (-Icges(6,1) * t455 - Icges(6,4) * t457) * t451;
t385 = Icges(6,6) * t449 + (-Icges(6,4) * t455 - Icges(6,2) * t457) * t451;
t384 = Icges(6,3) * t449 + (-Icges(6,5) * t455 - Icges(6,6) * t457) * t451;
t383 = rSges(7,3) * t449 + (-rSges(7,1) * t448 - rSges(7,2) * t450) * t451;
t382 = Icges(7,5) * t449 + (-Icges(7,1) * t448 - Icges(7,4) * t450) * t451;
t381 = Icges(7,6) * t449 + (-Icges(7,4) * t448 - Icges(7,2) * t450) * t451;
t380 = Icges(7,3) * t449 + (-Icges(7,5) * t448 - Icges(7,6) * t450) * t451;
t379 = qJD(1) * t461 * rSges(3,3) + t436 + (qJD(1) * t485 - qJD(2)) * t462;
t378 = t452 + (t462 * rSges(3,3) - t485 * t461 - t441) * qJD(1);
t377 = t466 * t461 - t518 * t462;
t376 = t518 * t461 + t466 * t462;
t375 = rSges(6,1) * t421 + rSges(6,2) * t420 + rSges(6,3) * t509;
t374 = rSges(6,1) * t419 + rSges(6,2) * t418 + rSges(6,3) * t508;
t373 = Icges(6,1) * t421 + Icges(6,4) * t420 + Icges(6,5) * t509;
t372 = Icges(6,1) * t419 + Icges(6,4) * t418 + Icges(6,5) * t508;
t371 = Icges(6,4) * t421 + Icges(6,2) * t420 + Icges(6,6) * t509;
t370 = Icges(6,4) * t419 + Icges(6,2) * t418 + Icges(6,6) * t508;
t369 = Icges(6,5) * t421 + Icges(6,6) * t420 + Icges(6,3) * t509;
t368 = Icges(6,5) * t419 + Icges(6,6) * t418 + Icges(6,3) * t508;
t367 = (t405 * t461 + t406 * t462) * qJD(3);
t366 = rSges(7,1) * t412 + rSges(7,2) * t411 + rSges(7,3) * t509;
t365 = rSges(7,1) * t410 + rSges(7,2) * t409 + rSges(7,3) * t508;
t364 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t509;
t363 = Icges(7,1) * t410 + Icges(7,4) * t409 + Icges(7,5) * t508;
t362 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t509;
t361 = Icges(7,4) * t410 + Icges(7,2) * t409 + Icges(7,6) * t508;
t360 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t509;
t359 = Icges(7,5) * t410 + Icges(7,6) * t409 + Icges(7,3) * t508;
t358 = qJD(1) * t406 - t435 * t499 + t487;
t357 = -t435 * t498 + t452 + (-t405 + t501) * qJD(1);
t356 = (t407 * t462 + t408 * t461) * qJD(3) + t486;
t355 = qJD(1) * t407 + t461 * t489 + t467;
t354 = t462 * t489 + (-t408 + t494) * qJD(1) + t500;
t353 = (t374 * t462 + t375 * t461) * qJD(3) + t465;
t352 = qJD(1) * t374 + t461 * t469 + t464;
t351 = t462 * t469 + (-t375 + t488) * qJD(1) + t493;
t350 = qJD(1) * t376 + t365 * t444 - t383 * t423 + t461 * t468 + t464;
t349 = -t366 * t444 + t383 * t424 + t462 * t468 + (-t377 + t488) * qJD(1) + t493;
t348 = -t365 * t424 + t366 * t423 + (t376 * t462 + t377 * t461) * qJD(3) + t465;
t1 = t424 * ((t359 * t509 + t361 * t411 + t363 * t412) * t423 + (t360 * t509 + t362 * t411 + t364 * t412) * t424 + (t380 * t509 + t381 * t411 + t382 * t412) * t444) / 0.2e1 + t444 * ((t359 * t423 + t360 * t424 + t380 * t444) * t449 + ((-t361 * t450 - t363 * t448) * t423 + (-t362 * t450 - t364 * t448) * t424 + (-t381 * t450 - t382 * t448) * t444) * t451) / 0.2e1 + t423 * ((t359 * t508 + t361 * t409 + t363 * t410) * t423 + (t360 * t508 + t362 * t409 + t364 * t410) * t424 + (t380 * t508 + t381 * t409 + t382 * t410) * t444) / 0.2e1 + m(5) * (t354 ^ 2 + t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(6) * (t351 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + m(7) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + m(4) * (t357 ^ 2 + t358 ^ 2 + t367 ^ 2) / 0.2e1 + m(3) * (t378 ^ 2 + t379 ^ 2) / 0.2e1 + (Icges(3,2) * t458 ^ 2 + (Icges(3,1) * t456 + 0.2e1 * Icges(3,4) * t458) * t456 + Icges(2,3) + m(2) * (t442 ^ 2 + t443 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((((t371 * t457 + t373 * t455 - t537) * t451 + (-t369 + t534) * t449) * t462 + ((-t370 * t457 - t372 * t455 - t536) * t451 + (t368 + t535) * t449) * t461) * qJD(3) + ((-t385 * t457 - t386 * t455 - t533) * t451 + (t384 + t532) * t449) * qJD(1)) * qJD(1) / 0.2e1 + (((-t369 * t508 - t371 * t418 - t373 * t419 + t526 * t462) * t462 + (t368 * t508 + t370 * t418 + t372 * t419 + (t525 - t530) * t462 + t529 * t461) * t461) * qJD(3) + (t384 * t508 + t385 * t418 + t386 * t419 + t461 * t528 + t462 * t527) * qJD(1)) * t499 / 0.2e1 - (((t368 * t509 + t370 * t420 + t372 * t421 + t525 * t461) * t461 + (-t369 * t509 - t371 * t420 - t373 * t421 + (t526 - t529) * t461 + t530 * t462) * t462) * qJD(3) + (t384 * t509 + t385 * t420 + t386 * t421 + t461 * t527 - t528 * t462) * qJD(1)) * t498 / 0.2e1;
T  = t1;
