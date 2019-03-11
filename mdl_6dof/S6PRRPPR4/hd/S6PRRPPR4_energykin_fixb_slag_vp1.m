% Calculate kinetic energy for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:30
% EndTime: 2019-03-08 21:12:33
% DurationCPUTime: 2.55s
% Computational Cost: add. (2738->324), mult. (7055->475), div. (0->0), fcn. (8839->12), ass. (0->146)
t557 = Icges(5,1) + Icges(6,1);
t556 = -Icges(5,4) + Icges(6,5);
t555 = Icges(6,4) + Icges(5,5);
t554 = Icges(5,2) + Icges(6,3);
t553 = -Icges(5,6) + Icges(6,6);
t552 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t506 = sin(pkin(10));
t508 = cos(pkin(10));
t514 = cos(qJ(2));
t509 = cos(pkin(6));
t512 = sin(qJ(2));
t532 = t509 * t512;
t493 = t506 * t514 + t508 * t532;
t511 = sin(qJ(3));
t507 = sin(pkin(6));
t533 = t508 * t507;
t538 = cos(qJ(3));
t478 = t493 * t538 - t511 * t533;
t531 = t509 * t514;
t492 = t506 * t512 - t508 * t531;
t505 = sin(pkin(11));
t537 = cos(pkin(11));
t448 = t478 * t505 - t492 * t537;
t449 = t478 * t537 + t492 * t505;
t525 = t507 * t538;
t477 = t493 * t511 + t508 * t525;
t551 = t554 * t448 + t556 * t449 + t553 * t477;
t495 = -t506 * t532 + t508 * t514;
t535 = t507 * t511;
t480 = t495 * t538 + t506 * t535;
t494 = t506 * t531 + t508 * t512;
t450 = t480 * t505 - t494 * t537;
t451 = t480 * t537 + t494 * t505;
t479 = t495 * t511 - t506 * t525;
t550 = t554 * t450 + t556 * t451 + t553 * t479;
t549 = t556 * t448 + t557 * t449 + t555 * t477;
t548 = t556 * t450 + t557 * t451 + t555 * t479;
t497 = t509 * t511 + t512 * t525;
t534 = t507 * t514;
t475 = t497 * t505 + t534 * t537;
t476 = t497 * t537 - t505 * t534;
t496 = -t509 * t538 + t512 * t535;
t547 = t554 * t475 + t556 * t476 + t553 * t496;
t546 = t556 * t475 + t557 * t476 + t555 * t496;
t545 = -Icges(4,4) * t478 - Icges(4,6) * t492 + t553 * t448 + t555 * t449 + t552 * t477;
t544 = -Icges(4,4) * t480 - Icges(4,6) * t494 + t553 * t450 + t555 * t451 + t552 * t479;
t543 = -Icges(4,4) * t497 + Icges(4,6) * t534 + t553 * t475 + t555 * t476 + t552 * t496;
t542 = qJD(2) ^ 2;
t536 = t506 * t507;
t412 = pkin(4) * t449 + qJ(5) * t448;
t442 = pkin(3) * t478 + qJ(4) * t477;
t530 = -t412 - t442;
t413 = pkin(4) * t451 + qJ(5) * t450;
t443 = pkin(3) * t480 + qJ(4) * t479;
t529 = -t413 - t443;
t441 = pkin(4) * t476 + qJ(5) * t475;
t470 = pkin(3) * t497 + qJ(4) * t496;
t528 = -t441 - t470;
t527 = qJD(2) * t507;
t502 = t506 * t527;
t481 = qJD(3) * t494 + t502;
t504 = qJD(2) * t509;
t524 = t508 * t527;
t468 = pkin(2) * t493 + pkin(8) * t492;
t469 = pkin(2) * t495 + pkin(8) * t494;
t523 = t468 * t502 + t469 * t524 + qJD(1);
t482 = qJD(3) * t492 - t524;
t499 = -qJD(3) * t534 + t504;
t522 = qJD(4) * t496 + t481 * t442 + t523;
t498 = (pkin(2) * t512 - pkin(8) * t514) * t507;
t521 = t469 * t504 - t498 * t502;
t520 = qJD(5) * t475 + t481 * t412 + t522;
t519 = qJD(4) * t477 + t499 * t443 + t521;
t518 = (-t468 * t509 - t498 * t533) * qJD(2);
t517 = qJD(5) * t448 + t499 * t413 + t519;
t516 = qJD(4) * t479 + t482 * t470 + t518;
t515 = qJD(5) * t450 + t482 * t441 + t516;
t513 = cos(qJ(6));
t510 = sin(qJ(6));
t486 = t509 * rSges(3,3) + (rSges(3,1) * t512 + rSges(3,2) * t514) * t507;
t485 = Icges(3,5) * t509 + (Icges(3,1) * t512 + Icges(3,4) * t514) * t507;
t484 = Icges(3,6) * t509 + (Icges(3,4) * t512 + Icges(3,2) * t514) * t507;
t483 = Icges(3,3) * t509 + (Icges(3,5) * t512 + Icges(3,6) * t514) * t507;
t474 = -qJD(6) * t496 + t499;
t466 = t497 * rSges(4,1) - t496 * rSges(4,2) - rSges(4,3) * t534;
t465 = Icges(4,1) * t497 - Icges(4,4) * t496 - Icges(4,5) * t534;
t463 = Icges(4,5) * t497 - Icges(4,6) * t496 - Icges(4,3) * t534;
t460 = rSges(3,1) * t495 - rSges(3,2) * t494 + rSges(3,3) * t536;
t459 = rSges(3,1) * t493 - rSges(3,2) * t492 - rSges(3,3) * t533;
t458 = Icges(3,1) * t495 - Icges(3,4) * t494 + Icges(3,5) * t536;
t457 = Icges(3,1) * t493 - Icges(3,4) * t492 - Icges(3,5) * t533;
t456 = Icges(3,4) * t495 - Icges(3,2) * t494 + Icges(3,6) * t536;
t455 = Icges(3,4) * t493 - Icges(3,2) * t492 - Icges(3,6) * t533;
t454 = Icges(3,5) * t495 - Icges(3,6) * t494 + Icges(3,3) * t536;
t453 = Icges(3,5) * t493 - Icges(3,6) * t492 - Icges(3,3) * t533;
t452 = pkin(5) * t476 - pkin(9) * t496;
t445 = -qJD(6) * t477 + t482;
t444 = -qJD(6) * t479 + t481;
t439 = t475 * t510 + t476 * t513;
t438 = t475 * t513 - t476 * t510;
t436 = (-t459 * t509 - t486 * t533) * qJD(2);
t435 = (t460 * t509 - t486 * t536) * qJD(2);
t434 = rSges(5,1) * t476 - rSges(5,2) * t475 + rSges(5,3) * t496;
t433 = rSges(6,1) * t476 + rSges(6,2) * t496 + rSges(6,3) * t475;
t432 = rSges(4,1) * t480 - rSges(4,2) * t479 + rSges(4,3) * t494;
t431 = rSges(4,1) * t478 - rSges(4,2) * t477 + rSges(4,3) * t492;
t424 = Icges(4,1) * t480 - Icges(4,4) * t479 + Icges(4,5) * t494;
t423 = Icges(4,1) * t478 - Icges(4,4) * t477 + Icges(4,5) * t492;
t420 = Icges(4,5) * t480 - Icges(4,6) * t479 + Icges(4,3) * t494;
t419 = Icges(4,5) * t478 - Icges(4,6) * t477 + Icges(4,3) * t492;
t418 = pkin(5) * t451 - pkin(9) * t479;
t417 = pkin(5) * t449 - pkin(9) * t477;
t414 = qJD(1) + (t459 * t506 + t460 * t508) * t527;
t411 = t450 * t510 + t451 * t513;
t410 = t450 * t513 - t451 * t510;
t409 = t448 * t510 + t449 * t513;
t408 = t448 * t513 - t449 * t510;
t405 = rSges(5,1) * t451 - rSges(5,2) * t450 + rSges(5,3) * t479;
t404 = rSges(6,1) * t451 + rSges(6,2) * t479 + rSges(6,3) * t450;
t403 = rSges(5,1) * t449 - rSges(5,2) * t448 + rSges(5,3) * t477;
t402 = rSges(6,1) * t449 + rSges(6,2) * t477 + rSges(6,3) * t448;
t389 = rSges(7,1) * t439 + rSges(7,2) * t438 - rSges(7,3) * t496;
t388 = Icges(7,1) * t439 + Icges(7,4) * t438 - Icges(7,5) * t496;
t387 = Icges(7,4) * t439 + Icges(7,2) * t438 - Icges(7,6) * t496;
t386 = Icges(7,5) * t439 + Icges(7,6) * t438 - Icges(7,3) * t496;
t385 = -t431 * t499 + t466 * t482 + t518;
t384 = t432 * t499 - t466 * t481 + t521;
t383 = rSges(7,1) * t411 + rSges(7,2) * t410 - rSges(7,3) * t479;
t382 = rSges(7,1) * t409 + rSges(7,2) * t408 - rSges(7,3) * t477;
t381 = Icges(7,1) * t411 + Icges(7,4) * t410 - Icges(7,5) * t479;
t380 = Icges(7,1) * t409 + Icges(7,4) * t408 - Icges(7,5) * t477;
t379 = Icges(7,4) * t411 + Icges(7,2) * t410 - Icges(7,6) * t479;
t378 = Icges(7,4) * t409 + Icges(7,2) * t408 - Icges(7,6) * t477;
t377 = Icges(7,5) * t411 + Icges(7,6) * t410 - Icges(7,3) * t479;
t376 = Icges(7,5) * t409 + Icges(7,6) * t408 - Icges(7,3) * t477;
t375 = t431 * t481 - t432 * t482 + t523;
t374 = t434 * t482 + (-t403 - t442) * t499 + t516;
t373 = t405 * t499 + (-t434 - t470) * t481 + t519;
t372 = t403 * t481 + (-t405 - t443) * t482 + t522;
t371 = t433 * t482 + (-t402 + t530) * t499 + t515;
t370 = t404 * t499 + (-t433 + t528) * t481 + t517;
t369 = t402 * t481 + (-t404 + t529) * t482 + t520;
t368 = -t382 * t474 + t389 * t445 + t452 * t482 + (-t417 + t530) * t499 + t515;
t367 = t383 * t474 - t389 * t444 + t418 * t499 + (-t452 + t528) * t481 + t517;
t366 = t382 * t444 - t383 * t445 + t417 * t481 + (-t418 + t529) * t482 + t520;
t1 = -t542 * ((-t454 * t533 - t456 * t492 + t458 * t493) * t536 - (-t453 * t533 - t455 * t492 + t457 * t493) * t533 + (-t483 * t533 - t484 * t492 + t485 * t493) * t509) * t533 / 0.2e1 + t445 * ((-t377 * t477 + t379 * t408 + t381 * t409) * t444 + (-t477 * t376 + t408 * t378 + t409 * t380) * t445 + (-t386 * t477 + t387 * t408 + t388 * t409) * t474) / 0.2e1 + t474 * ((-t377 * t496 + t379 * t438 + t381 * t439) * t444 + (-t376 * t496 + t378 * t438 + t380 * t439) * t445 + (-t496 * t386 + t438 * t387 + t439 * t388) * t474) / 0.2e1 + t444 * ((-t479 * t377 + t410 * t379 + t411 * t381) * t444 + (-t376 * t479 + t378 * t410 + t380 * t411) * t445 + (-t386 * t479 + t387 * t410 + t388 * t411) * t474) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t414 ^ 2 + t435 ^ 2 + t436 ^ 2) / 0.2e1 + m(4) * (t375 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(5) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(7) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + (((t454 * t536 - t456 * t494 + t458 * t495) * t536 - (t453 * t536 - t455 * t494 + t457 * t495) * t533 + (t483 * t536 - t484 * t494 + t485 * t495) * t509) * t536 + t509 * (t509 ^ 2 * t483 + (((t456 * t514 + t458 * t512) * t506 - (t455 * t514 + t457 * t512) * t508) * t507 + (-t453 * t508 + t454 * t506 + t484 * t514 + t485 * t512) * t509) * t507)) * t542 / 0.2e1 + ((t547 * t450 + t546 * t451 + t463 * t494 + t465 * t480 + t543 * t479) * t499 + (t419 * t494 + t423 * t480 + t551 * t450 + t549 * t451 + t545 * t479) * t482 + (t494 * t420 + t480 * t424 + t550 * t450 + t548 * t451 + t544 * t479) * t481) * t481 / 0.2e1 + ((t547 * t448 + t546 * t449 + t463 * t492 + t465 * t478 + t543 * t477) * t499 + (t492 * t419 + t478 * t423 + t551 * t448 + t549 * t449 + t545 * t477) * t482 + (t420 * t492 + t424 * t478 + t550 * t448 + t548 * t449 + t544 * t477) * t481) * t482 / 0.2e1 + ((-t463 * t534 + t497 * t465 + t547 * t475 + t546 * t476 + t543 * t496) * t499 + (-t419 * t534 + t497 * t423 + t551 * t475 + t549 * t476 + t545 * t496) * t482 + (-t420 * t534 + t497 * t424 + t550 * t475 + t548 * t476 + t544 * t496) * t481) * t499 / 0.2e1;
T  = t1;
