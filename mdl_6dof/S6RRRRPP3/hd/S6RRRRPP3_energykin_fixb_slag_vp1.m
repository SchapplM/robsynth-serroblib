% Calculate kinetic energy for
% S6RRRRPP3
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:41
% EndTime: 2019-03-09 20:53:43
% DurationCPUTime: 2.22s
% Computational Cost: add. (1613->236), mult. (2177->363), div. (0->0), fcn. (2160->8), ass. (0->133)
t509 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t508 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t507 = -Icges(5,5) - Icges(7,5) + Icges(6,4);
t506 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t505 = Icges(5,6) - Icges(6,5) - Icges(7,4);
t504 = -Icges(5,3) - Icges(7,1) - Icges(6,1);
t503 = rSges(7,1) + pkin(5);
t502 = rSges(7,3) + qJ(6);
t436 = qJ(2) + qJ(3);
t435 = cos(t436);
t440 = cos(qJ(4));
t442 = cos(qJ(1));
t476 = t440 * t442;
t437 = sin(qJ(4));
t439 = sin(qJ(1));
t478 = t437 * t439;
t407 = t435 * t478 + t476;
t475 = t442 * t437;
t477 = t439 * t440;
t408 = t435 * t477 - t475;
t434 = sin(t436);
t481 = t434 * t439;
t501 = t508 * t407 + t509 * t408 - t507 * t481;
t409 = t435 * t475 - t477;
t410 = t435 * t476 + t478;
t479 = t434 * t442;
t500 = t508 * t409 + t509 * t410 - t507 * t479;
t499 = t506 * t407 + t508 * t408 - t505 * t481;
t498 = t506 * t409 + t508 * t410 - t505 * t479;
t497 = -t505 * t407 - t507 * t408 - t504 * t481;
t496 = -t505 * t409 - t507 * t410 - t504 * t479;
t495 = t504 * t435 + (-t505 * t437 - t507 * t440) * t434;
t494 = t505 * t435 + (t506 * t437 + t508 * t440) * t434;
t493 = t507 * t435 + (t508 * t437 + t509 * t440) * t434;
t441 = cos(qJ(2));
t487 = pkin(2) * t441;
t438 = sin(qJ(2));
t485 = Icges(3,4) * t438;
t484 = Icges(3,4) * t441;
t483 = Icges(4,4) * t434;
t482 = Icges(4,4) * t435;
t480 = t434 * t440;
t474 = rSges(7,2) * t407 + t502 * t408 + t503 * t481;
t473 = rSges(7,2) * t409 + t502 * t410 + t503 * t479;
t380 = -pkin(8) * t442 + t439 * t487;
t381 = pkin(8) * t439 + t442 * t487;
t433 = qJD(2) * t439;
t469 = qJD(2) * t442;
t472 = t380 * t433 + t381 * t469;
t471 = (rSges(7,2) * t437 + rSges(7,3) * t440) * t434 + qJ(6) * t480 - t503 * t435;
t427 = pkin(1) * t439 - pkin(7) * t442;
t470 = -t380 - t427;
t418 = qJD(3) * t439 + t433;
t468 = qJD(4) * t434;
t467 = pkin(2) * qJD(2) * t438;
t419 = (-qJD(2) - qJD(3)) * t442;
t466 = t442 * t467;
t465 = pkin(3) * t435 + pkin(9) * t434;
t464 = rSges(3,1) * t441 - rSges(3,2) * t438;
t463 = rSges(4,1) * t435 - rSges(4,2) * t434;
t462 = Icges(3,1) * t441 - t485;
t461 = Icges(4,1) * t435 - t483;
t460 = -Icges(3,2) * t438 + t484;
t459 = -Icges(4,2) * t434 + t482;
t458 = Icges(3,5) * t441 - Icges(3,6) * t438;
t457 = Icges(4,5) * t435 - Icges(4,6) * t434;
t394 = -Icges(3,6) * t442 + t439 * t460;
t396 = -Icges(3,5) * t442 + t439 * t462;
t456 = t394 * t438 - t396 * t441;
t395 = Icges(3,6) * t439 + t442 * t460;
t397 = Icges(3,5) * t439 + t442 * t462;
t455 = -t395 * t438 + t397 * t441;
t421 = Icges(3,2) * t441 + t485;
t422 = Icges(3,1) * t438 + t484;
t454 = -t421 * t438 + t422 * t441;
t405 = t465 * t439;
t406 = t465 * t442;
t453 = t418 * t405 - t406 * t419 + t472;
t417 = qJD(1) * (pkin(1) * t442 + pkin(7) * t439);
t452 = qJD(1) * t381 - t439 * t467 + t417;
t360 = pkin(4) * t408 + qJ(5) * t407;
t398 = t442 * t468 + t418;
t451 = qJD(5) * t434 * t437 + t398 * t360 + t453;
t450 = qJD(1) * (Icges(4,5) * t434 + Icges(4,6) * t435) + (-Icges(4,3) * t442 + t439 * t457) * t419 + (Icges(4,3) * t439 + t442 * t457) * t418;
t416 = pkin(3) * t434 - pkin(9) * t435;
t449 = qJD(1) * t406 - t416 * t418 + t452;
t448 = t419 * t416 + (-t405 + t470) * qJD(1) - t466;
t361 = pkin(4) * t410 + qJ(5) * t409;
t428 = -qJD(4) * t435 + qJD(1);
t447 = qJD(5) * t407 + t428 * t361 + t449;
t399 = t439 * t468 + t419;
t404 = (pkin(4) * t440 + qJ(5) * t437) * t434;
t446 = qJD(5) * t409 + t399 * t404 + t448;
t385 = -Icges(4,6) * t442 + t439 * t459;
t386 = Icges(4,6) * t439 + t442 * t459;
t387 = -Icges(4,5) * t442 + t439 * t461;
t388 = Icges(4,5) * t439 + t442 * t461;
t413 = Icges(4,2) * t435 + t483;
t414 = Icges(4,1) * t434 + t482;
t445 = (-t386 * t434 + t388 * t435) * t418 + (-t385 * t434 + t387 * t435) * t419 + (-t413 * t434 + t414 * t435) * qJD(1);
t425 = rSges(2,1) * t442 - rSges(2,2) * t439;
t424 = rSges(2,1) * t439 + rSges(2,2) * t442;
t423 = rSges(3,1) * t438 + rSges(3,2) * t441;
t420 = Icges(3,5) * t438 + Icges(3,6) * t441;
t415 = rSges(4,1) * t434 + rSges(4,2) * t435;
t403 = rSges(3,3) * t439 + t442 * t464;
t402 = -rSges(3,3) * t442 + t439 * t464;
t393 = Icges(3,3) * t439 + t442 * t458;
t392 = -Icges(3,3) * t442 + t439 * t458;
t390 = rSges(4,3) * t439 + t442 * t463;
t389 = -rSges(4,3) * t442 + t439 * t463;
t379 = -rSges(6,1) * t435 + (-rSges(6,2) * t440 + rSges(6,3) * t437) * t434;
t377 = -rSges(5,3) * t435 + (rSges(5,1) * t440 - rSges(5,2) * t437) * t434;
t358 = qJD(1) * t403 - t423 * t433 + t417;
t357 = -t423 * t469 + (-t402 - t427) * qJD(1);
t356 = rSges(5,1) * t410 - rSges(5,2) * t409 + rSges(5,3) * t479;
t355 = rSges(5,1) * t408 - rSges(5,2) * t407 + rSges(5,3) * t481;
t354 = rSges(6,1) * t479 - rSges(6,2) * t410 + rSges(6,3) * t409;
t352 = rSges(6,1) * t481 - rSges(6,2) * t408 + rSges(6,3) * t407;
t331 = (t402 * t439 + t403 * t442) * qJD(2);
t329 = qJD(1) * t390 - t415 * t418 + t452;
t328 = -t466 + t415 * t419 + (-t389 + t470) * qJD(1);
t327 = t389 * t418 - t390 * t419 + t472;
t326 = t356 * t428 - t377 * t398 + t449;
t325 = -t355 * t428 + t377 * t399 + t448;
t324 = t355 * t398 - t356 * t399 + t453;
t323 = t354 * t428 + (-t379 - t404) * t398 + t447;
t322 = t379 * t399 + (-t352 - t360) * t428 + t446;
t321 = t352 * t398 + (-t354 - t361) * t399 + t451;
t320 = qJD(6) * t408 + t473 * t428 + (-t404 - t471) * t398 + t447;
t319 = qJD(6) * t410 + t471 * t399 + (-t360 - t474) * t428 + t446;
t318 = qJD(6) * t480 + t474 * t398 + (-t361 - t473) * t399 + t451;
t1 = t418 * (t450 * t439 + t445 * t442) / 0.2e1 + t419 * (t445 * t439 - t450 * t442) / 0.2e1 + m(3) * (t331 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + m(6) * (t321 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(7) * (t318 ^ 2 + t319 ^ 2 + t320 ^ 2) / 0.2e1 + m(5) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(4) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 - ((-t442 * t420 + t439 * t454) * qJD(1) + (t442 ^ 2 * t392 + (t455 * t439 + (-t393 + t456) * t442) * t439) * qJD(2)) * t469 / 0.2e1 + ((t439 * t420 + t442 * t454) * qJD(1) + (t439 ^ 2 * t393 + (t456 * t442 + (-t392 + t455) * t439) * t442) * qJD(2)) * t433 / 0.2e1 + (m(2) * (t424 ^ 2 + t425 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t386 * t435 + t388 * t434) * t418 + (t385 * t435 + t387 * t434) * t419 + ((t395 * t441 + t397 * t438) * t439 - (t394 * t441 + t396 * t438) * t442) * qJD(2) + (t435 * t413 + t434 * t414 + t421 * t441 + t422 * t438) * qJD(1)) * qJD(1) / 0.2e1 + ((t409 * t494 + t410 * t493 + t479 * t495) * t428 + (t499 * t409 + t410 * t501 + t497 * t479) * t399 + (t498 * t409 + t500 * t410 + t496 * t479) * t398) * t398 / 0.2e1 + ((t407 * t494 + t408 * t493 + t481 * t495) * t428 + (t499 * t407 + t501 * t408 + t497 * t481) * t399 + (t407 * t498 + t408 * t500 + t481 * t496) * t398) * t399 / 0.2e1 + ((-t398 * t496 - t399 * t497 - t428 * t495) * t435 + ((t437 * t494 + t440 * t493) * t428 + (t499 * t437 + t440 * t501) * t399 + (t437 * t498 + t440 * t500) * t398) * t434) * t428 / 0.2e1;
T  = t1;
