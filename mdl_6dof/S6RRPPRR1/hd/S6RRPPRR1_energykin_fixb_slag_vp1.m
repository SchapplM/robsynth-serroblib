% Calculate kinetic energy for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:46
% EndTime: 2019-03-09 08:45:50
% DurationCPUTime: 3.20s
% Computational Cost: add. (1649->274), mult. (2210->417), div. (0->0), fcn. (2289->10), ass. (0->147)
t547 = Icges(4,4) - Icges(5,5);
t546 = Icges(4,1) + Icges(5,1);
t545 = Icges(4,2) + Icges(5,3);
t457 = qJ(2) + pkin(10);
t452 = sin(t457);
t544 = t547 * t452;
t453 = cos(t457);
t543 = t547 * t453;
t542 = Icges(5,4) + Icges(4,5);
t541 = Icges(4,6) - Icges(5,6);
t540 = t545 * t452 - t543;
t539 = t546 * t453 - t544;
t462 = sin(qJ(1));
t465 = cos(qJ(1));
t538 = t462 * t540 + t465 * t541;
t537 = -t462 * t541 + t465 * t540;
t536 = -t539 * t462 + t465 * t542;
t535 = t462 * t542 + t539 * t465;
t534 = -t545 * t453 - t544;
t533 = t546 * t452 + t543;
t532 = Icges(5,2) + Icges(3,3) + Icges(4,3);
t461 = sin(qJ(2));
t464 = cos(qJ(2));
t531 = Icges(3,5) * t464 - Icges(3,6) * t461 - t541 * t452 + t453 * t542;
t530 = t462 * t532 + t465 * t531;
t529 = t462 * t531 - t465 * t532;
t528 = Icges(3,5) * t461 + Icges(3,6) * t464 + t452 * t542 + t541 * t453;
t516 = Icges(3,4) * t461;
t443 = Icges(3,2) * t464 + t516;
t515 = Icges(3,4) * t464;
t444 = Icges(3,1) * t461 + t515;
t527 = -t443 * t461 + t444 * t464 + t452 * t534 + t453 * t533;
t488 = -Icges(3,2) * t461 + t515;
t414 = Icges(3,6) * t462 + t465 * t488;
t491 = Icges(3,1) * t464 - t516;
t416 = Icges(3,5) * t462 + t465 * t491;
t526 = -t414 * t461 + t416 * t464 + t452 * t537 + t453 * t535;
t413 = -Icges(3,6) * t465 + t462 * t488;
t415 = -Icges(3,5) * t465 + t462 * t491;
t525 = t413 * t461 - t415 * t464 - t452 * t538 + t536 * t453;
t521 = cos(qJ(5));
t520 = pkin(2) * t461;
t518 = pkin(2) * t464;
t510 = t453 * t462;
t509 = t453 * t465;
t385 = -qJ(3) * t465 + t462 * t518;
t386 = qJ(3) * t462 + t465 * t518;
t456 = qJD(2) * t462;
t505 = qJD(2) * t465;
t508 = t385 * t456 + t386 * t505;
t450 = pkin(1) * t462 - pkin(7) * t465;
t507 = -t385 - t450;
t455 = qJD(3) * t462;
t504 = qJD(4) * t452;
t506 = t465 * t504 + t455;
t493 = pkin(3) * t453 + qJ(4) * t452;
t417 = t493 * t462;
t503 = -t417 + t507;
t502 = t452 * t521;
t499 = -pkin(3) * t452 + qJ(4) * t453 - t520;
t441 = qJD(5) * t465 - t505;
t440 = -qJD(5) * t462 + t456;
t432 = pkin(4) * t510 + pkin(8) * t465;
t498 = -t432 + t503;
t437 = qJD(1) * (pkin(1) * t465 + pkin(7) * t462);
t497 = qJD(1) * t386 - qJD(3) * t465 + t437;
t496 = rSges(3,1) * t464 - rSges(3,2) * t461;
t495 = rSges(4,1) * t453 - rSges(4,2) * t452;
t494 = rSges(5,1) * t453 + rSges(5,3) * t452;
t492 = qJD(2) * (-rSges(4,1) * t452 - rSges(4,2) * t453 - t520);
t473 = qJD(2) * (-rSges(5,1) * t452 + rSges(5,3) * t453 + t499);
t460 = sin(qJ(5));
t423 = t452 * t460 + t453 * t521;
t418 = t493 * t465;
t472 = qJD(1) * t418 + t462 * t504 + t497;
t471 = -qJD(4) * t453 + t417 * t456 + t418 * t505 + t508;
t470 = qJD(2) * (-pkin(4) * t452 + t499);
t433 = pkin(4) * t509 - pkin(8) * t462;
t469 = t432 * t456 + t433 * t505 + t471;
t468 = t465 * t470 + t506;
t467 = qJD(1) * t433 + t462 * t470 + t472;
t463 = cos(qJ(6));
t459 = sin(qJ(6));
t449 = rSges(2,1) * t465 - rSges(2,2) * t462;
t448 = rSges(2,1) * t462 + rSges(2,2) * t465;
t447 = rSges(3,1) * t461 + rSges(3,2) * t464;
t424 = -t453 * t460 + t502;
t420 = rSges(3,3) * t462 + t465 * t496;
t419 = -rSges(3,3) * t465 + t462 * t496;
t409 = qJD(6) * t423 + qJD(1);
t408 = t423 * t465;
t407 = t460 * t509 - t465 * t502;
t406 = t423 * t462;
t405 = t460 * t510 - t462 * t502;
t404 = rSges(4,3) * t462 + t465 * t495;
t403 = rSges(5,2) * t462 + t465 * t494;
t402 = -rSges(4,3) * t465 + t462 * t495;
t401 = -rSges(5,2) * t465 + t462 * t494;
t381 = t408 * t463 - t459 * t462;
t380 = -t408 * t459 - t462 * t463;
t379 = t406 * t463 + t459 * t465;
t378 = -t406 * t459 + t463 * t465;
t377 = qJD(6) * t405 + t441;
t376 = qJD(6) * t407 + t440;
t375 = pkin(5) * t424 + pkin(9) * t423;
t374 = rSges(6,1) * t424 - rSges(6,2) * t423;
t373 = Icges(6,1) * t424 - Icges(6,4) * t423;
t372 = Icges(6,4) * t424 - Icges(6,2) * t423;
t371 = Icges(6,5) * t424 - Icges(6,6) * t423;
t370 = qJD(1) * t420 - t447 * t456 + t437;
t369 = -t447 * t505 + (-t419 - t450) * qJD(1);
t368 = (t419 * t462 + t420 * t465) * qJD(2);
t367 = pkin(5) * t408 + pkin(9) * t407;
t366 = pkin(5) * t406 + pkin(9) * t405;
t365 = rSges(6,1) * t408 - rSges(6,2) * t407 - rSges(6,3) * t462;
t364 = rSges(6,1) * t406 - rSges(6,2) * t405 + rSges(6,3) * t465;
t363 = Icges(6,1) * t408 - Icges(6,4) * t407 - Icges(6,5) * t462;
t362 = Icges(6,1) * t406 - Icges(6,4) * t405 + Icges(6,5) * t465;
t361 = Icges(6,4) * t408 - Icges(6,2) * t407 - Icges(6,6) * t462;
t360 = Icges(6,4) * t406 - Icges(6,2) * t405 + Icges(6,6) * t465;
t359 = Icges(6,5) * t408 - Icges(6,6) * t407 - Icges(6,3) * t462;
t358 = Icges(6,5) * t406 - Icges(6,6) * t405 + Icges(6,3) * t465;
t357 = rSges(7,3) * t423 + (rSges(7,1) * t463 - rSges(7,2) * t459) * t424;
t356 = Icges(7,5) * t423 + (Icges(7,1) * t463 - Icges(7,4) * t459) * t424;
t355 = Icges(7,6) * t423 + (Icges(7,4) * t463 - Icges(7,2) * t459) * t424;
t354 = Icges(7,3) * t423 + (Icges(7,5) * t463 - Icges(7,6) * t459) * t424;
t353 = rSges(7,1) * t381 + rSges(7,2) * t380 + rSges(7,3) * t407;
t352 = rSges(7,1) * t379 + rSges(7,2) * t378 + rSges(7,3) * t405;
t351 = Icges(7,1) * t381 + Icges(7,4) * t380 + Icges(7,5) * t407;
t350 = Icges(7,1) * t379 + Icges(7,4) * t378 + Icges(7,5) * t405;
t349 = Icges(7,4) * t381 + Icges(7,2) * t380 + Icges(7,6) * t407;
t348 = Icges(7,4) * t379 + Icges(7,2) * t378 + Icges(7,6) * t405;
t347 = Icges(7,5) * t381 + Icges(7,6) * t380 + Icges(7,3) * t407;
t346 = Icges(7,5) * t379 + Icges(7,6) * t378 + Icges(7,3) * t405;
t345 = qJD(1) * t404 + t462 * t492 + t497;
t344 = t455 + t465 * t492 + (-t402 + t507) * qJD(1);
t343 = (t402 * t462 + t404 * t465) * qJD(2) + t508;
t342 = qJD(1) * t403 + t462 * t473 + t472;
t341 = t465 * t473 + (-t401 + t503) * qJD(1) + t506;
t340 = (t401 * t462 + t403 * t465) * qJD(2) + t471;
t339 = qJD(1) * t365 - t374 * t440 + t467;
t338 = t374 * t441 + (-t364 + t498) * qJD(1) + t468;
t337 = t364 * t440 - t365 * t441 + t469;
t336 = qJD(1) * t367 + t353 * t409 - t357 * t376 - t375 * t440 + t467;
t335 = -t352 * t409 + t357 * t377 + t375 * t441 + (-t366 + t498) * qJD(1) + t468;
t334 = t352 * t376 - t353 * t377 + t366 * t440 - t367 * t441 + t469;
t1 = m(6) * (t337 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + m(4) * (t343 ^ 2 + t344 ^ 2 + t345 ^ 2) / 0.2e1 + m(3) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(7) * (t334 ^ 2 + t335 ^ 2 + t336 ^ 2) / 0.2e1 + m(5) * (t340 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 + t377 * ((t347 * t405 + t349 * t378 + t351 * t379) * t376 + (t405 * t346 + t378 * t348 + t379 * t350) * t377 + (t354 * t405 + t355 * t378 + t356 * t379) * t409) / 0.2e1 + t409 * ((t346 * t377 + t347 * t376 + t354 * t409) * t423 + ((-t349 * t459 + t351 * t463) * t376 + (-t348 * t459 + t350 * t463) * t377 + (-t355 * t459 + t356 * t463) * t409) * t424) / 0.2e1 + t376 * ((t407 * t347 + t380 * t349 + t381 * t351) * t376 + (t346 * t407 + t348 * t380 + t350 * t381) * t377 + (t354 * t407 + t355 * t380 + t356 * t381) * t409) / 0.2e1 + t440 * ((-t462 * t359 - t407 * t361 + t408 * t363) * t440 + (-t358 * t462 - t360 * t407 + t362 * t408) * t441 + (-t371 * t462 - t372 * t407 + t373 * t408) * qJD(1)) / 0.2e1 + t441 * ((t359 * t465 - t361 * t405 + t363 * t406) * t440 + (t465 * t358 - t405 * t360 + t406 * t362) * t441 + (t371 * t465 - t372 * t405 + t373 * t406) * qJD(1)) / 0.2e1 + (Icges(2,3) + m(2) * (t448 ^ 2 + t449 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t530 * t462 ^ 2 + (t525 * t465 + (t526 - t529) * t462) * t465) * qJD(2) + (t528 * t462 + t527 * t465) * qJD(1)) * t456 / 0.2e1 - ((t529 * t465 ^ 2 + (t526 * t462 + (t525 - t530) * t465) * t462) * qJD(2) + (t462 * t527 - t465 * t528) * qJD(1)) * t505 / 0.2e1 + ((-t361 * t423 + t363 * t424) * t440 + (-t360 * t423 + t362 * t424) * t441 + ((-t464 * t413 - t461 * t415 + t536 * t452 + t453 * t538) * t465 + (t464 * t414 + t461 * t416 + t452 * t535 - t537 * t453) * t462) * qJD(2) + (-t423 * t372 + t424 * t373 + t464 * t443 + t461 * t444 + t533 * t452 - t534 * t453) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
