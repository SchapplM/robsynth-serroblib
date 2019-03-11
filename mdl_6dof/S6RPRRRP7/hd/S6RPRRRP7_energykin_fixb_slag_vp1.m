% Calculate kinetic energy for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:37
% EndTime: 2019-03-09 06:18:39
% DurationCPUTime: 2.21s
% Computational Cost: add. (1866->250), mult. (1954->394), div. (0->0), fcn. (1954->10), ass. (0->135)
t527 = Icges(6,1) + Icges(7,1);
t526 = -Icges(6,4) + Icges(7,5);
t525 = Icges(7,4) + Icges(6,5);
t524 = Icges(6,2) + Icges(7,3);
t523 = -Icges(7,6) + Icges(6,6);
t522 = -Icges(6,3) - Icges(7,2);
t521 = rSges(7,1) + pkin(5);
t520 = rSges(7,3) + qJ(6);
t456 = pkin(10) + qJ(3);
t451 = cos(t456);
t457 = qJ(4) + qJ(5);
t455 = cos(t457);
t464 = cos(qJ(1));
t497 = t455 * t464;
t454 = sin(t457);
t462 = sin(qJ(1));
t499 = t454 * t462;
t421 = t451 * t499 + t497;
t494 = t462 * t455;
t498 = t454 * t464;
t422 = t451 * t494 - t498;
t450 = sin(t456);
t501 = t450 * t462;
t519 = t524 * t421 + t526 * t422 - t523 * t501;
t423 = t451 * t498 - t494;
t424 = t451 * t497 + t499;
t500 = t450 * t464;
t518 = t524 * t423 + t526 * t424 - t523 * t500;
t517 = -t523 * t421 + t525 * t422 - t522 * t501;
t516 = -t523 * t423 + t525 * t424 - t522 * t500;
t515 = t526 * t421 + t527 * t422 + t525 * t501;
t514 = t526 * t423 + t527 * t424 + t525 * t500;
t513 = t523 * t451 + (t524 * t454 + t526 * t455) * t450;
t512 = t522 * t451 + (-t523 * t454 + t525 * t455) * t450;
t511 = -t525 * t451 + (t526 * t454 + t527 * t455) * t450;
t459 = cos(pkin(10));
t506 = pkin(2) * t459;
t463 = cos(qJ(4));
t505 = pkin(4) * t463;
t503 = Icges(4,4) * t450;
t502 = Icges(4,4) * t451;
t461 = sin(qJ(4));
t496 = t461 * t462;
t495 = t461 * t464;
t493 = t462 * t463;
t492 = t463 * t464;
t490 = rSges(7,2) * t501 + t520 * t421 + t521 * t422;
t489 = rSges(7,2) * t500 + t520 * t423 + t521 * t424;
t488 = -rSges(7,2) * t451 + (t520 * t454 + t521 * t455) * t450;
t443 = pkin(1) * t462 - qJ(2) * t464;
t487 = pkin(7) * t464 - t462 * t506 - t443;
t482 = pkin(3) * t451 + pkin(8) * t450;
t426 = t482 * t462;
t427 = t482 * t464;
t452 = qJD(3) * t462;
t485 = qJD(3) * t464;
t486 = t426 * t452 + t427 * t485;
t484 = qJD(4) * t450;
t433 = t464 * t484 + t452;
t483 = qJD(5) * t450;
t434 = t462 * t484 - t485;
t440 = qJD(1) * (pkin(1) * t464 + qJ(2) * t462);
t481 = -qJD(2) * t464 + qJD(1) * (pkin(7) * t462 + t464 * t506) + t440;
t458 = sin(pkin(10));
t480 = rSges(3,1) * t459 - rSges(3,2) * t458;
t479 = rSges(4,1) * t451 - rSges(4,2) * t450;
t478 = Icges(4,1) * t451 - t503;
t477 = -Icges(4,2) * t450 + t502;
t476 = Icges(4,5) * t451 - Icges(4,6) * t450;
t412 = -Icges(4,6) * t464 + t462 * t477;
t414 = -Icges(4,5) * t464 + t462 * t478;
t475 = t412 * t450 - t414 * t451;
t413 = Icges(4,6) * t462 + t464 * t477;
t415 = Icges(4,5) * t462 + t464 * t478;
t474 = -t413 * t450 + t415 * t451;
t436 = Icges(4,2) * t451 + t503;
t437 = Icges(4,1) * t450 + t502;
t473 = -t436 * t450 + t437 * t451;
t471 = pkin(9) * t450 + t451 * t505;
t378 = -pkin(4) * t495 + t462 * t471;
t379 = pkin(4) * t496 + t464 * t471;
t472 = t433 * t378 - t379 * t434 + t486;
t439 = pkin(3) * t450 - pkin(8) * t451;
t470 = qJD(1) * t427 - t439 * t452 + t481;
t453 = qJD(2) * t462;
t469 = t453 + (-t426 + t487) * qJD(1) - t439 * t485;
t393 = -pkin(9) * t451 + t450 * t505;
t446 = -qJD(4) * t451 + qJD(1);
t468 = t446 * t379 - t393 * t433 + t470;
t467 = -t378 * t446 + t434 * t393 + t469;
t445 = rSges(2,1) * t464 - rSges(2,2) * t462;
t444 = rSges(2,1) * t462 + rSges(2,2) * t464;
t438 = rSges(4,1) * t450 + rSges(4,2) * t451;
t435 = Icges(4,5) * t450 + Icges(4,6) * t451;
t432 = qJD(1) + (-qJD(4) - qJD(5)) * t451;
t431 = t451 * t492 + t496;
t430 = -t451 * t495 + t493;
t429 = t451 * t493 - t495;
t428 = -t451 * t496 - t492;
t417 = rSges(4,3) * t462 + t464 * t479;
t416 = -rSges(4,3) * t464 + t462 * t479;
t411 = Icges(4,3) * t462 + t464 * t476;
t410 = -Icges(4,3) * t464 + t462 * t476;
t409 = t462 * t483 + t434;
t408 = t464 * t483 + t433;
t406 = -rSges(5,3) * t451 + (rSges(5,1) * t463 - rSges(5,2) * t461) * t450;
t405 = -Icges(5,5) * t451 + (Icges(5,1) * t463 - Icges(5,4) * t461) * t450;
t404 = -Icges(5,6) * t451 + (Icges(5,4) * t463 - Icges(5,2) * t461) * t450;
t403 = -Icges(5,3) * t451 + (Icges(5,5) * t463 - Icges(5,6) * t461) * t450;
t401 = -rSges(6,3) * t451 + (rSges(6,1) * t455 - rSges(6,2) * t454) * t450;
t392 = qJD(1) * t462 * rSges(3,3) + t440 + (qJD(1) * t480 - qJD(2)) * t464;
t391 = t453 + (t464 * rSges(3,3) - t462 * t480 - t443) * qJD(1);
t388 = rSges(5,1) * t431 + rSges(5,2) * t430 + rSges(5,3) * t500;
t387 = rSges(5,1) * t429 + rSges(5,2) * t428 + rSges(5,3) * t501;
t385 = Icges(5,1) * t431 + Icges(5,4) * t430 + Icges(5,5) * t500;
t384 = Icges(5,1) * t429 + Icges(5,4) * t428 + Icges(5,5) * t501;
t383 = Icges(5,4) * t431 + Icges(5,2) * t430 + Icges(5,6) * t500;
t382 = Icges(5,4) * t429 + Icges(5,2) * t428 + Icges(5,6) * t501;
t381 = Icges(5,5) * t431 + Icges(5,6) * t430 + Icges(5,3) * t500;
t380 = Icges(5,5) * t429 + Icges(5,6) * t428 + Icges(5,3) * t501;
t377 = rSges(6,1) * t424 - rSges(6,2) * t423 + rSges(6,3) * t500;
t375 = rSges(6,1) * t422 - rSges(6,2) * t421 + rSges(6,3) * t501;
t361 = (t416 * t462 + t417 * t464) * qJD(3);
t358 = qJD(1) * t417 - t438 * t452 + t481;
t357 = -t438 * t485 + t453 + (-t416 + t487) * qJD(1);
t356 = t387 * t433 - t388 * t434 + t486;
t355 = t388 * t446 - t406 * t433 + t470;
t354 = -t387 * t446 + t406 * t434 + t469;
t353 = t377 * t432 - t401 * t408 + t468;
t352 = -t375 * t432 + t401 * t409 + t467;
t351 = t375 * t408 - t377 * t409 + t472;
t350 = qJD(6) * t421 - t408 * t488 + t432 * t489 + t468;
t349 = qJD(6) * t423 + t409 * t488 - t432 * t490 + t467;
t348 = qJD(6) * t450 * t454 + t408 * t490 - t409 * t489 + t472;
t1 = m(3) * (t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(4) * (t357 ^ 2 + t358 ^ 2 + t361 ^ 2) / 0.2e1 + m(5) * (t354 ^ 2 + t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(6) * (t351 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + m(7) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 - ((-t464 * t435 + t462 * t473) * qJD(1) + (t464 ^ 2 * t410 + (t474 * t462 + (-t411 + t475) * t464) * t462) * qJD(3)) * t485 / 0.2e1 + qJD(1) * ((t451 * t436 + t450 * t437) * qJD(1) + ((t413 * t451 + t415 * t450) * t462 - (t412 * t451 + t414 * t450) * t464) * qJD(3)) / 0.2e1 + t433 * ((t381 * t500 + t383 * t430 + t431 * t385) * t433 + (t380 * t500 + t382 * t430 + t384 * t431) * t434 + (t403 * t500 + t404 * t430 + t405 * t431) * t446) / 0.2e1 + t434 * ((t381 * t501 + t383 * t428 + t385 * t429) * t433 + (t380 * t501 + t382 * t428 + t384 * t429) * t434 + (t403 * t501 + t404 * t428 + t405 * t429) * t446) / 0.2e1 + t446 * ((-t380 * t434 - t381 * t433 - t403 * t446) * t451 + ((-t383 * t461 + t385 * t463) * t433 + (-t382 * t461 + t384 * t463) * t434 + (-t404 * t461 + t405 * t463) * t446) * t450) / 0.2e1 + ((t462 * t435 + t464 * t473) * qJD(1) + (t462 ^ 2 * t411 + (t475 * t464 + (-t410 + t474) * t462) * t464) * qJD(3)) * t452 / 0.2e1 + ((t513 * t423 + t511 * t424 + t512 * t500) * t432 + (t519 * t423 + t515 * t424 + t517 * t500) * t409 + (t518 * t423 + t514 * t424 + t516 * t500) * t408) * t408 / 0.2e1 + ((t513 * t421 + t511 * t422 + t512 * t501) * t432 + (t519 * t421 + t515 * t422 + t517 * t501) * t409 + (t518 * t421 + t514 * t422 + t516 * t501) * t408) * t409 / 0.2e1 + ((-t516 * t408 - t517 * t409 - t512 * t432) * t451 + ((t513 * t454 + t511 * t455) * t432 + (t519 * t454 + t515 * t455) * t409 + (t518 * t454 + t514 * t455) * t408) * t450) * t432 / 0.2e1 + (m(2) * (t444 ^ 2 + t445 ^ 2) + Icges(2,3) + Icges(3,2) * t459 ^ 2 + (Icges(3,1) * t458 + 0.2e1 * Icges(3,4) * t459) * t458) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
