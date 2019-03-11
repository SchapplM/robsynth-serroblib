% Calculate kinetic energy for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:12
% EndTime: 2019-03-09 09:54:15
% DurationCPUTime: 2.42s
% Computational Cost: add. (1532->261), mult. (2413->391), div. (0->0), fcn. (2452->8), ass. (0->133)
t512 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t511 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t510 = -Icges(5,5) - Icges(7,5) + Icges(6,4);
t509 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t508 = Icges(5,6) - Icges(6,5) - Icges(7,4);
t507 = -Icges(5,3) - Icges(7,1) - Icges(6,1);
t506 = rSges(7,1) + pkin(5);
t505 = rSges(7,3) + qJ(6);
t441 = pkin(9) + qJ(4);
t439 = sin(t441);
t440 = cos(t441);
t448 = cos(qJ(1));
t446 = sin(qJ(1));
t447 = cos(qJ(2));
t482 = t446 * t447;
t398 = t439 * t482 + t440 * t448;
t480 = t448 * t439;
t399 = t440 * t482 - t480;
t445 = sin(qJ(2));
t484 = t445 * t446;
t504 = t511 * t398 + t399 * t512 - t510 * t484;
t400 = -t446 * t440 + t447 * t480;
t481 = t447 * t448;
t401 = t439 * t446 + t440 * t481;
t483 = t445 * t448;
t503 = t511 * t400 + t401 * t512 - t510 * t483;
t502 = t398 * t509 + t399 * t511 - t484 * t508;
t501 = t400 * t509 + t401 * t511 - t483 * t508;
t500 = -t398 * t508 - t399 * t510 - t484 * t507;
t499 = -t400 * t508 - t401 * t510 - t483 * t507;
t498 = t507 * t447 + (-t439 * t508 - t440 * t510) * t445;
t497 = t508 * t447 + (t439 * t509 + t440 * t511) * t445;
t496 = t510 * t447 + (t511 * t439 + t440 * t512) * t445;
t443 = cos(pkin(9));
t490 = pkin(3) * t443;
t489 = Icges(3,4) * t445;
t488 = Icges(3,4) * t447;
t487 = t440 * t445;
t442 = sin(pkin(9));
t486 = t442 * t446;
t485 = t442 * t448;
t478 = rSges(7,2) * t398 + t505 * t399 + t506 * t484;
t477 = rSges(7,2) * t400 + t505 * t401 + t506 * t483;
t476 = (rSges(7,2) * t439 + rSges(7,3) * t440) * t445 + qJ(6) * t487 - t506 * t447;
t463 = pkin(2) * t447 + qJ(3) * t445;
t418 = t463 * t446;
t432 = pkin(1) * t446 - pkin(7) * t448;
t475 = -t418 - t432;
t474 = qJD(2) * t446;
t473 = qJD(2) * t448;
t472 = qJD(3) * t445;
t471 = qJD(4) * t445;
t419 = t463 * t448;
t423 = qJD(1) * (pkin(1) * t448 + pkin(7) * t446);
t470 = qJD(1) * t419 + t446 * t472 + t423;
t428 = pkin(2) * t445 - qJ(3) * t447;
t467 = qJD(2) * (pkin(8) * t447 - t445 * t490 - t428);
t466 = qJD(2) * (rSges(4,3) * t447 - (rSges(4,1) * t443 - rSges(4,2) * t442) * t445 - t428);
t465 = -qJD(3) * t447 + t418 * t474 + t419 * t473;
t464 = rSges(3,1) * t447 - rSges(3,2) * t445;
t462 = Icges(3,1) * t447 - t489;
t461 = -Icges(3,2) * t445 + t488;
t460 = Icges(3,5) * t447 - Icges(3,6) * t445;
t404 = -Icges(3,6) * t448 + t446 * t461;
t406 = -Icges(3,5) * t448 + t446 * t462;
t459 = t404 * t445 - t406 * t447;
t405 = Icges(3,6) * t446 + t448 * t461;
t407 = Icges(3,5) * t446 + t448 * t462;
t458 = -t405 * t445 + t407 * t447;
t425 = Icges(3,2) * t447 + t489;
t426 = Icges(3,1) * t445 + t488;
t457 = -t425 * t445 + t426 * t447;
t455 = pkin(8) * t445 + t447 * t490;
t374 = -pkin(3) * t485 + t446 * t455;
t375 = pkin(3) * t486 + t448 * t455;
t456 = t374 * t474 + t375 * t473 + t465;
t364 = pkin(4) * t399 + qJ(5) * t398;
t421 = t448 * t471 + t474;
t454 = qJD(5) * t445 * t439 + t421 * t364 + t456;
t453 = qJD(1) * t375 + t446 * t467 + t470;
t365 = pkin(4) * t401 + qJ(5) * t400;
t437 = -qJD(4) * t447 + qJD(1);
t452 = qJD(5) * t398 + t437 * t365 + t453;
t436 = t448 * t472;
t451 = t436 + (-t374 + t475) * qJD(1) + t448 * t467;
t410 = (pkin(4) * t440 + qJ(5) * t439) * t445;
t422 = t446 * t471 - t473;
t450 = qJD(5) * t400 + t422 * t410 + t451;
t431 = rSges(2,1) * t448 - rSges(2,2) * t446;
t430 = rSges(2,1) * t446 + rSges(2,2) * t448;
t429 = rSges(3,1) * t445 + rSges(3,2) * t447;
t424 = Icges(3,5) * t445 + Icges(3,6) * t447;
t417 = t443 * t481 + t486;
t416 = -t442 * t481 + t443 * t446;
t415 = t443 * t482 - t485;
t414 = -t442 * t482 - t443 * t448;
t412 = rSges(3,3) * t446 + t448 * t464;
t411 = -rSges(3,3) * t448 + t446 * t464;
t403 = Icges(3,3) * t446 + t448 * t460;
t402 = -Icges(3,3) * t448 + t446 * t460;
t396 = -Icges(4,5) * t447 + (Icges(4,1) * t443 - Icges(4,4) * t442) * t445;
t395 = -Icges(4,6) * t447 + (Icges(4,4) * t443 - Icges(4,2) * t442) * t445;
t394 = -Icges(4,3) * t447 + (Icges(4,5) * t443 - Icges(4,6) * t442) * t445;
t391 = -rSges(6,1) * t447 + (-rSges(6,2) * t440 + rSges(6,3) * t439) * t445;
t389 = -rSges(5,3) * t447 + (rSges(5,1) * t440 - rSges(5,2) * t439) * t445;
t373 = rSges(4,1) * t417 + rSges(4,2) * t416 + rSges(4,3) * t483;
t372 = rSges(4,1) * t415 + rSges(4,2) * t414 + rSges(4,3) * t484;
t371 = Icges(4,1) * t417 + Icges(4,4) * t416 + Icges(4,5) * t483;
t370 = Icges(4,1) * t415 + Icges(4,4) * t414 + Icges(4,5) * t484;
t369 = Icges(4,4) * t417 + Icges(4,2) * t416 + Icges(4,6) * t483;
t368 = Icges(4,4) * t415 + Icges(4,2) * t414 + Icges(4,6) * t484;
t367 = Icges(4,5) * t417 + Icges(4,6) * t416 + Icges(4,3) * t483;
t366 = Icges(4,5) * t415 + Icges(4,6) * t414 + Icges(4,3) * t484;
t362 = qJD(1) * t412 - t429 * t474 + t423;
t361 = -t429 * t473 + (-t411 - t432) * qJD(1);
t358 = (t411 * t446 + t412 * t448) * qJD(2);
t357 = rSges(5,1) * t401 - rSges(5,2) * t400 + rSges(5,3) * t483;
t356 = rSges(5,1) * t399 - rSges(5,2) * t398 + rSges(5,3) * t484;
t355 = rSges(6,1) * t483 - rSges(6,2) * t401 + rSges(6,3) * t400;
t353 = rSges(6,1) * t484 - rSges(6,2) * t399 + rSges(6,3) * t398;
t331 = qJD(1) * t373 + t446 * t466 + t470;
t330 = t436 + t448 * t466 + (-t372 + t475) * qJD(1);
t329 = (t372 * t446 + t373 * t448) * qJD(2) + t465;
t328 = t357 * t437 - t389 * t421 + t453;
t327 = -t356 * t437 + t389 * t422 + t451;
t326 = t356 * t421 - t357 * t422 + t456;
t325 = t355 * t437 + (-t391 - t410) * t421 + t452;
t324 = t391 * t422 + (-t353 - t364) * t437 + t450;
t323 = t353 * t421 + (-t355 - t365) * t422 + t454;
t322 = qJD(6) * t399 + t477 * t437 + (-t410 - t476) * t421 + t452;
t321 = qJD(6) * t401 + t476 * t422 + (-t364 - t478) * t437 + t450;
t320 = qJD(6) * t487 + t478 * t421 + (-t365 - t477) * t422 + t454;
t1 = m(4) * (t329 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + m(3) * (t358 ^ 2 + t361 ^ 2 + t362 ^ 2) / 0.2e1 + m(6) * (t323 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + m(7) * (t320 ^ 2 + t321 ^ 2 + t322 ^ 2) / 0.2e1 + m(5) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + (m(2) * (t430 ^ 2 + t431 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t366 * t448 - t367 * t446) * t447 + ((-t369 * t442 + t371 * t443) * t446 - (-t368 * t442 + t370 * t443) * t448) * t445 + (t447 * t405 + t407 * t445) * t446 - (t404 * t447 + t406 * t445) * t448) * qJD(2) + ((-t394 + t425) * t447 + (-t395 * t442 + t396 * t443 + t426) * t445) * qJD(1)) * qJD(1) / 0.2e1 + (((-t366 * t483 - t368 * t416 - t370 * t417 + t459 * t448) * t448 + ((-t402 + t458) * t448 + t367 * t483 + t416 * t369 + t417 * t371 + t446 * t403) * t446) * qJD(2) + (t394 * t483 + t395 * t416 + t396 * t417 + t446 * t424 + t448 * t457) * qJD(1)) * t474 / 0.2e1 - (((-t366 * t484 - t414 * t368 - t415 * t370 + t448 * t402) * t448 + (t367 * t484 + t369 * t414 + t371 * t415 + (-t403 + t459) * t448 + t458 * t446) * t446) * qJD(2) + (t394 * t484 + t395 * t414 + t396 * t415 - t448 * t424 + t446 * t457) * qJD(1)) * t473 / 0.2e1 + ((t497 * t400 + t496 * t401 + t498 * t483) * t437 + (t502 * t400 + t504 * t401 + t500 * t483) * t422 + (t501 * t400 + t503 * t401 + t499 * t483) * t421) * t421 / 0.2e1 + ((t497 * t398 + t496 * t399 + t498 * t484) * t437 + (t502 * t398 + t504 * t399 + t500 * t484) * t422 + (t501 * t398 + t503 * t399 + t499 * t484) * t421) * t422 / 0.2e1 + ((-t499 * t421 - t500 * t422 - t498 * t437) * t447 + ((t497 * t439 + t496 * t440) * t437 + (t502 * t439 + t504 * t440) * t422 + (t501 * t439 + t503 * t440) * t421) * t445) * t437 / 0.2e1;
T  = t1;
