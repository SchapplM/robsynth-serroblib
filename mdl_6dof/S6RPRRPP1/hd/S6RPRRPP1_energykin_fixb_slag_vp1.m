% Calculate kinetic energy for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:27:52
% EndTime: 2019-03-09 04:27:54
% DurationCPUTime: 1.93s
% Computational Cost: add. (1835->236), mult. (1839->351), div. (0->0), fcn. (1838->10), ass. (0->124)
t486 = Icges(6,1) + Icges(7,1);
t485 = -Icges(6,4) + Icges(7,5);
t484 = Icges(7,4) + Icges(6,5);
t483 = Icges(6,2) + Icges(7,3);
t482 = -Icges(7,6) + Icges(6,6);
t481 = -Icges(6,3) - Icges(7,2) - Icges(5,3);
t480 = rSges(7,1) + pkin(5);
t479 = rSges(7,3) + qJ(6);
t419 = qJ(4) + pkin(10);
t415 = sin(t419);
t417 = cos(t419);
t420 = qJ(1) + pkin(9);
t418 = cos(t420);
t416 = sin(t420);
t426 = cos(qJ(3));
t458 = t416 * t426;
t372 = t415 * t458 + t417 * t418;
t373 = -t415 * t418 + t417 * t458;
t423 = sin(qJ(3));
t459 = t416 * t423;
t478 = t483 * t372 + t485 * t373 - t482 * t459;
t455 = t418 * t426;
t374 = t415 * t455 - t416 * t417;
t375 = t415 * t416 + t417 * t455;
t456 = t418 * t423;
t477 = t483 * t374 + t485 * t375 - t482 * t456;
t476 = t485 * t372 + t486 * t373 + t484 * t459;
t475 = t485 * t374 + t486 * t375 + t484 * t456;
t474 = t482 * t426 + (t483 * t415 + t485 * t417) * t423;
t473 = -t484 * t426 + (t485 * t415 + t486 * t417) * t423;
t425 = cos(qJ(4));
t422 = sin(qJ(4));
t454 = t422 * t426;
t387 = -t416 * t454 - t418 * t425;
t453 = t425 * t426;
t457 = t418 * t422;
t388 = t416 * t453 - t457;
t472 = Icges(5,5) * t388 + Icges(5,6) * t387 - t482 * t372 + t484 * t373 - t481 * t459;
t389 = t416 * t425 - t418 * t454;
t460 = t416 * t422;
t390 = t418 * t453 + t460;
t471 = Icges(5,5) * t390 + Icges(5,6) * t389 - t482 * t374 + t484 * t375 - t481 * t456;
t470 = t481 * t426 + (Icges(5,5) * t425 - Icges(5,6) * t422 - t482 * t415 + t484 * t417) * t423;
t424 = sin(qJ(1));
t465 = pkin(1) * t424;
t464 = pkin(4) * t425;
t462 = Icges(4,4) * t423;
t461 = Icges(4,4) * t426;
t452 = rSges(7,2) * t459 + t479 * t372 + t480 * t373;
t451 = rSges(7,2) * t456 + t479 * t374 + t480 * t375;
t450 = -rSges(7,2) * t426 + (t479 * t415 + t480 * t417) * t423;
t427 = cos(qJ(1));
t414 = qJD(1) * t427 * pkin(1);
t449 = qJD(1) * (pkin(2) * t418 + pkin(7) * t416) + t414;
t448 = qJD(3) * t416;
t447 = qJD(3) * t418;
t446 = qJD(4) * t423;
t445 = qJD(5) * t423;
t442 = pkin(3) * t426 + pkin(8) * t423;
t396 = t442 * t416;
t397 = t442 * t418;
t444 = t396 * t448 + t397 * t447 + qJD(2);
t443 = -pkin(2) * t416 + pkin(7) * t418 - t465;
t441 = rSges(4,1) * t426 - rSges(4,2) * t423;
t440 = Icges(4,1) * t426 - t462;
t439 = -Icges(4,2) * t423 + t461;
t438 = Icges(4,5) * t426 - Icges(4,6) * t423;
t363 = -Icges(4,6) * t418 + t439 * t416;
t365 = -Icges(4,5) * t418 + t440 * t416;
t437 = t363 * t423 - t365 * t426;
t364 = Icges(4,6) * t416 + t439 * t418;
t366 = Icges(4,5) * t416 + t440 * t418;
t436 = -t364 * t423 + t366 * t426;
t404 = Icges(4,2) * t426 + t462;
t405 = Icges(4,1) * t423 + t461;
t435 = -t404 * t423 + t405 * t426;
t411 = pkin(3) * t423 - pkin(8) * t426;
t434 = qJD(1) * t397 - t411 * t448 + t449;
t432 = qJ(5) * t423 + t464 * t426;
t354 = -pkin(4) * t457 + t432 * t416;
t398 = t418 * t446 + t448;
t433 = -qJD(5) * t426 + t398 * t354 + t444;
t355 = pkin(4) * t460 + t432 * t418;
t412 = -qJD(4) * t426 + qJD(1);
t431 = t412 * t355 + t416 * t445 + t434;
t430 = (-t396 + t443) * qJD(1) - t411 * t447;
t371 = -qJ(5) * t426 + t464 * t423;
t399 = t416 * t446 - t447;
t429 = t399 * t371 + t418 * t445 + t430;
t410 = rSges(2,1) * t427 - rSges(2,2) * t424;
t409 = rSges(2,1) * t424 + rSges(2,2) * t427;
t408 = rSges(4,1) * t423 + rSges(4,2) * t426;
t403 = Icges(4,5) * t423 + Icges(4,6) * t426;
t395 = -rSges(5,3) * t426 + (rSges(5,1) * t425 - rSges(5,2) * t422) * t423;
t393 = -Icges(5,5) * t426 + (Icges(5,1) * t425 - Icges(5,4) * t422) * t423;
t392 = -Icges(5,6) * t426 + (Icges(5,4) * t425 - Icges(5,2) * t422) * t423;
t385 = t414 + qJD(1) * (rSges(3,1) * t418 - rSges(3,2) * t416);
t384 = (-rSges(3,1) * t416 - rSges(3,2) * t418 - t465) * qJD(1);
t383 = -rSges(6,3) * t426 + (rSges(6,1) * t417 - rSges(6,2) * t415) * t423;
t368 = rSges(4,3) * t416 + t441 * t418;
t367 = -rSges(4,3) * t418 + t441 * t416;
t362 = Icges(4,3) * t416 + t438 * t418;
t361 = -Icges(4,3) * t418 + t438 * t416;
t357 = rSges(5,1) * t390 + rSges(5,2) * t389 + rSges(5,3) * t456;
t356 = rSges(5,1) * t388 + rSges(5,2) * t387 + rSges(5,3) * t459;
t353 = Icges(5,1) * t390 + Icges(5,4) * t389 + Icges(5,5) * t456;
t352 = Icges(5,1) * t388 + Icges(5,4) * t387 + Icges(5,5) * t459;
t351 = Icges(5,4) * t390 + Icges(5,2) * t389 + Icges(5,6) * t456;
t350 = Icges(5,4) * t388 + Icges(5,2) * t387 + Icges(5,6) * t459;
t346 = rSges(6,1) * t375 - rSges(6,2) * t374 + rSges(6,3) * t456;
t344 = rSges(6,1) * t373 - rSges(6,2) * t372 + rSges(6,3) * t459;
t330 = qJD(1) * t368 - t408 * t448 + t449;
t329 = -t408 * t447 + (-t367 + t443) * qJD(1);
t328 = qJD(2) + (t367 * t416 + t368 * t418) * qJD(3);
t326 = t357 * t412 - t395 * t398 + t434;
t325 = -t356 * t412 + t395 * t399 + t430;
t324 = t356 * t398 - t357 * t399 + t444;
t323 = t346 * t412 + (-t371 - t383) * t398 + t431;
t322 = t383 * t399 + (-t344 - t354) * t412 + t429;
t321 = t344 * t398 + (-t346 - t355) * t399 + t433;
t320 = qJD(6) * t372 + t451 * t412 + (-t371 - t450) * t398 + t431;
t319 = qJD(6) * t374 + t450 * t399 + (-t354 - t452) * t412 + t429;
t318 = qJD(6) * t415 * t423 + t452 * t398 + (-t355 - t451) * t399 + t433;
t1 = ((t416 * t403 + t435 * t418) * qJD(1) + (t416 ^ 2 * t362 + (t437 * t418 + (-t361 + t436) * t416) * t418) * qJD(3)) * t448 / 0.2e1 + m(3) * (qJD(2) ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(4) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(5) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(6) * (t321 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(7) * (t318 ^ 2 + t319 ^ 2 + t320 ^ 2) / 0.2e1 - ((-t418 * t403 + t435 * t416) * qJD(1) + (t418 ^ 2 * t361 + (t436 * t416 + (-t362 + t437) * t418) * t416) * qJD(3)) * t447 / 0.2e1 + qJD(1) * ((t426 * t404 + t423 * t405) * qJD(1) + ((t364 * t426 + t366 * t423) * t416 - (t363 * t426 + t365 * t423) * t418) * qJD(3)) / 0.2e1 + ((t374 * t474 + t375 * t473 + t389 * t392 + t390 * t393 + t456 * t470) * t412 + (t350 * t389 + t352 * t390 + t374 * t478 + t476 * t375 + t472 * t456) * t399 + (t389 * t351 + t390 * t353 + t477 * t374 + t475 * t375 + t471 * t456) * t398) * t398 / 0.2e1 + ((t372 * t474 + t373 * t473 + t387 * t392 + t388 * t393 + t459 * t470) * t412 + (t387 * t350 + t388 * t352 + t478 * t372 + t476 * t373 + t472 * t459) * t399 + (t351 * t387 + t353 * t388 + t372 * t477 + t373 * t475 + t459 * t471) * t398) * t399 / 0.2e1 + ((-t398 * t471 - t399 * t472 - t412 * t470) * t426 + ((-t392 * t422 + t393 * t425 + t415 * t474 + t417 * t473) * t412 + (-t350 * t422 + t352 * t425 + t415 * t478 + t476 * t417) * t399 + (-t351 * t422 + t353 * t425 + t415 * t477 + t417 * t475) * t398) * t423) * t412 / 0.2e1 + (m(2) * (t409 ^ 2 + t410 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
