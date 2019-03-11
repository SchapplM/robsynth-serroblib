% Calculate kinetic energy for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:14
% EndTime: 2019-03-09 21:05:17
% DurationCPUTime: 2.61s
% Computational Cost: add. (1621->265), mult. (2568->405), div. (0->0), fcn. (2605->8), ass. (0->138)
t510 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t509 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t508 = Icges(7,5) - Icges(6,4) - Icges(5,5);
t507 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t506 = -Icges(6,6) + Icges(7,6) + Icges(5,6);
t505 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t504 = rSges(7,1) + pkin(5);
t503 = rSges(7,3) + qJ(6);
t440 = qJ(3) + qJ(4);
t438 = sin(t440);
t439 = cos(t440);
t443 = sin(qJ(1));
t445 = cos(qJ(2));
t446 = cos(qJ(1));
t474 = t445 * t446;
t406 = t438 * t474 - t443 * t439;
t407 = t438 * t443 + t439 * t474;
t360 = pkin(4) * t407 + qJ(5) * t406;
t475 = t443 * t445;
t404 = t438 * t475 + t439 * t446;
t421 = qJD(1) + (-qJD(3) - qJD(4)) * t445;
t502 = qJD(5) * t404 + t421 * t360;
t442 = sin(qJ(2));
t444 = cos(qJ(3));
t485 = pkin(3) * t444;
t452 = pkin(9) * t442 + t445 * t485;
t441 = sin(qJ(3));
t479 = t441 * t443;
t370 = pkin(3) * t479 + t446 * t452;
t375 = -pkin(9) * t445 + t442 * t485;
t437 = qJD(2) * t443;
t467 = qJD(3) * t442;
t419 = t446 * t467 + t437;
t435 = -qJD(3) * t445 + qJD(1);
t501 = t435 * t370 - t375 * t419;
t468 = qJD(2) * t446;
t420 = t443 * t467 - t468;
t466 = qJD(4) * t442;
t391 = t443 * t466 + t420;
t410 = (pkin(4) * t439 + qJ(5) * t438) * t442;
t500 = qJD(5) * t406 + t391 * t410;
t478 = t441 * t446;
t369 = -pkin(3) * t478 + t443 * t452;
t499 = -t369 * t435 + t420 * t375;
t405 = -t438 * t446 + t439 * t475;
t477 = t442 * t443;
t498 = t506 * t404 + t508 * t405 - t505 * t477;
t476 = t442 * t446;
t497 = t506 * t406 + t508 * t407 - t505 * t476;
t496 = t507 * t404 + t509 * t405 - t506 * t477;
t495 = t507 * t406 + t509 * t407 - t506 * t476;
t494 = t509 * t404 + t510 * t405 - t508 * t477;
t493 = t509 * t406 + t510 * t407 - t508 * t476;
t492 = t505 * t445 + (t506 * t438 + t508 * t439) * t442;
t491 = t506 * t445 + (t507 * t438 + t509 * t439) * t442;
t490 = t508 * t445 + (t509 * t438 + t510 * t439) * t442;
t483 = Icges(3,4) * t442;
t482 = Icges(3,4) * t445;
t473 = rSges(7,2) * t404 + t405 * t504 - t503 * t477;
t472 = rSges(7,2) * t406 + t407 * t504 - t503 * t476;
t471 = t503 * t445 + (rSges(7,2) * t438 + t439 * t504) * t442;
t464 = pkin(2) * t445 + pkin(8) * t442;
t416 = t464 * t443;
t417 = t464 * t446;
t470 = t416 * t437 + t417 * t468;
t422 = qJD(1) * (pkin(1) * t446 + pkin(7) * t443);
t469 = qJD(1) * t417 + t422;
t431 = pkin(1) * t443 - pkin(7) * t446;
t465 = (-t416 - t431) * qJD(1);
t463 = rSges(3,1) * t445 - rSges(3,2) * t442;
t462 = Icges(3,1) * t445 - t483;
t461 = -Icges(3,2) * t442 + t482;
t460 = Icges(3,5) * t445 - Icges(3,6) * t442;
t396 = -Icges(3,6) * t446 + t443 * t461;
t399 = -Icges(3,5) * t446 + t443 * t462;
t459 = t396 * t442 - t399 * t445;
t397 = Icges(3,6) * t443 + t446 * t461;
t400 = Icges(3,5) * t443 + t446 * t462;
t458 = -t397 * t442 + t400 * t445;
t424 = Icges(3,2) * t445 + t483;
t425 = Icges(3,1) * t442 + t482;
t457 = -t424 * t442 + t425 * t445;
t430 = pkin(2) * t442 - pkin(8) * t445;
t456 = -qJD(2) * t430 - qJD(6) * t442;
t455 = t419 * t369 - t370 * t420 + t470;
t454 = -t430 * t437 + t469;
t359 = pkin(4) * t405 + qJ(5) * t404;
t390 = t446 * t466 + t419;
t453 = qJD(5) * t442 * t438 + t390 * t359 + t455;
t451 = -t430 * t468 + t465;
t450 = t454 + t501;
t449 = t451 + t499;
t428 = rSges(2,1) * t446 - rSges(2,2) * t443;
t427 = rSges(2,1) * t443 + rSges(2,2) * t446;
t426 = rSges(3,1) * t442 + rSges(3,2) * t445;
t423 = Icges(3,5) * t442 + Icges(3,6) * t445;
t415 = t444 * t474 + t479;
t414 = -t441 * t474 + t443 * t444;
t413 = t444 * t475 - t478;
t412 = -t441 * t475 - t444 * t446;
t403 = rSges(3,3) * t443 + t446 * t463;
t402 = -rSges(3,3) * t446 + t443 * t463;
t401 = -rSges(4,3) * t445 + (rSges(4,1) * t444 - rSges(4,2) * t441) * t442;
t398 = -Icges(4,5) * t445 + (Icges(4,1) * t444 - Icges(4,4) * t441) * t442;
t395 = -Icges(4,6) * t445 + (Icges(4,4) * t444 - Icges(4,2) * t441) * t442;
t394 = Icges(3,3) * t443 + t446 * t460;
t393 = -Icges(3,3) * t446 + t443 * t460;
t392 = -Icges(4,3) * t445 + (Icges(4,5) * t444 - Icges(4,6) * t441) * t442;
t387 = -rSges(5,3) * t445 + (rSges(5,1) * t439 - rSges(5,2) * t438) * t442;
t386 = -rSges(6,2) * t445 + (rSges(6,1) * t439 + rSges(6,3) * t438) * t442;
t368 = rSges(4,1) * t415 + rSges(4,2) * t414 + rSges(4,3) * t476;
t367 = rSges(4,1) * t413 + rSges(4,2) * t412 + rSges(4,3) * t477;
t366 = Icges(4,1) * t415 + Icges(4,4) * t414 + Icges(4,5) * t476;
t365 = Icges(4,1) * t413 + Icges(4,4) * t412 + Icges(4,5) * t477;
t364 = Icges(4,4) * t415 + Icges(4,2) * t414 + Icges(4,6) * t476;
t363 = Icges(4,4) * t413 + Icges(4,2) * t412 + Icges(4,6) * t477;
t362 = Icges(4,5) * t415 + Icges(4,6) * t414 + Icges(4,3) * t476;
t361 = Icges(4,5) * t413 + Icges(4,6) * t412 + Icges(4,3) * t477;
t358 = qJD(1) * t403 - t426 * t437 + t422;
t357 = -t426 * t468 + (-t402 - t431) * qJD(1);
t356 = (t402 * t443 + t403 * t446) * qJD(2);
t354 = rSges(5,1) * t407 - rSges(5,2) * t406 + rSges(5,3) * t476;
t353 = rSges(6,1) * t407 + rSges(6,2) * t476 + rSges(6,3) * t406;
t351 = rSges(5,1) * t405 - rSges(5,2) * t404 + rSges(5,3) * t477;
t350 = rSges(6,1) * t405 + rSges(6,2) * t477 + rSges(6,3) * t404;
t327 = t368 * t435 - t401 * t419 + t454;
t326 = -t367 * t435 + t401 * t420 + t451;
t325 = t367 * t419 - t368 * t420 + t470;
t324 = t354 * t421 - t387 * t390 + t450;
t323 = -t351 * t421 + t387 * t391 + t449;
t322 = t351 * t390 - t354 * t391 + t455;
t321 = t353 * t421 + (-t386 - t410) * t390 + t450 + t502;
t320 = t386 * t391 + (-t350 - t359) * t421 + t449 + t500;
t319 = t350 * t390 + (-t353 - t360) * t391 + t453;
t318 = t456 * t443 + t472 * t421 + (-t410 - t471) * t390 + t469 + t501 + t502;
t317 = t456 * t446 + t471 * t391 + t465 + (-t359 - t473) * t421 + t499 + t500;
t316 = qJD(6) * t445 + t473 * t390 + (-t360 - t472) * t391 + t453;
t1 = -((-t423 * t446 + t443 * t457) * qJD(1) + (t393 * t446 ^ 2 + (t458 * t443 + (-t394 + t459) * t446) * t443) * qJD(2)) * t468 / 0.2e1 + ((t423 * t443 + t446 * t457) * qJD(1) + (t394 * t443 ^ 2 + (t459 * t446 + (-t393 + t458) * t443) * t446) * qJD(2)) * t437 / 0.2e1 + m(3) * (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + m(4) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(5) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + m(6) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(7) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + qJD(1) * ((t424 * t445 + t425 * t442) * qJD(1) + ((t397 * t445 + t400 * t442) * t443 - (t396 * t445 + t399 * t442) * t446) * qJD(2)) / 0.2e1 + t419 * ((t362 * t476 + t364 * t414 + t366 * t415) * t419 + (t361 * t476 + t363 * t414 + t365 * t415) * t420 + (t392 * t476 + t395 * t414 + t398 * t415) * t435) / 0.2e1 + t420 * ((t362 * t477 + t364 * t412 + t366 * t413) * t419 + (t361 * t477 + t363 * t412 + t365 * t413) * t420 + (t392 * t477 + t395 * t412 + t398 * t413) * t435) / 0.2e1 + t435 * ((-t361 * t420 - t362 * t419 - t392 * t435) * t445 + ((-t364 * t441 + t366 * t444) * t419 + (-t363 * t441 + t365 * t444) * t420 + (-t395 * t441 + t398 * t444) * t435) * t442) / 0.2e1 + (Icges(2,3) + m(2) * (t427 ^ 2 + t428 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t491 * t406 + t490 * t407 - t492 * t476) * t421 + (t496 * t406 + t494 * t407 - t498 * t476) * t391 + (t495 * t406 + t493 * t407 - t497 * t476) * t390) * t390 / 0.2e1 + ((t491 * t404 + t490 * t405 - t492 * t477) * t421 + (t496 * t404 + t494 * t405 - t498 * t477) * t391 + (t495 * t404 + t493 * t405 - t497 * t477) * t390) * t391 / 0.2e1 + ((t497 * t390 + t498 * t391 + t492 * t421) * t445 + ((t491 * t438 + t490 * t439) * t421 + (t496 * t438 + t494 * t439) * t391 + (t495 * t438 + t493 * t439) * t390) * t442) * t421 / 0.2e1;
T  = t1;
