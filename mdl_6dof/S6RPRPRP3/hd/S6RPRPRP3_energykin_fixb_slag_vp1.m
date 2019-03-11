% Calculate kinetic energy for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:07:46
% EndTime: 2019-03-09 03:07:48
% DurationCPUTime: 2.12s
% Computational Cost: add. (1790->243), mult. (1794->370), div. (0->0), fcn. (1793->10), ass. (0->128)
t494 = Icges(6,1) + Icges(7,1);
t493 = -Icges(6,4) + Icges(7,5);
t492 = Icges(7,4) + Icges(6,5);
t491 = Icges(6,2) + Icges(7,3);
t490 = -Icges(7,6) + Icges(6,6);
t489 = -Icges(6,3) - Icges(7,2);
t488 = rSges(7,1) + pkin(5);
t487 = rSges(7,3) + qJ(6);
t422 = pkin(10) + qJ(5);
t418 = sin(t422);
t420 = cos(t422);
t423 = qJ(1) + pkin(9);
t421 = cos(t423);
t419 = sin(t423);
t429 = cos(qJ(3));
t466 = t419 * t429;
t375 = t418 * t466 + t420 * t421;
t376 = -t418 * t421 + t420 * t466;
t427 = sin(qJ(3));
t467 = t419 * t427;
t486 = t491 * t375 + t493 * t376 - t490 * t467;
t463 = t421 * t429;
t377 = t418 * t463 - t419 * t420;
t378 = t418 * t419 + t420 * t463;
t464 = t421 * t427;
t485 = t491 * t377 + t493 * t378 - t490 * t464;
t484 = -t490 * t375 + t492 * t376 - t489 * t467;
t483 = -t490 * t377 + t492 * t378 - t489 * t464;
t482 = t493 * t375 + t494 * t376 + t492 * t467;
t481 = t493 * t377 + t494 * t378 + t492 * t464;
t480 = t490 * t429 + (t491 * t418 + t493 * t420) * t427;
t479 = t489 * t429 + (-t490 * t418 + t492 * t420) * t427;
t478 = -t492 * t429 + (t493 * t418 + t494 * t420) * t427;
t428 = sin(qJ(1));
t473 = pkin(1) * t428;
t425 = cos(pkin(10));
t471 = pkin(4) * t425;
t470 = Icges(4,4) * t427;
t469 = Icges(4,4) * t429;
t424 = sin(pkin(10));
t468 = t419 * t424;
t465 = t421 * t424;
t462 = t424 * t429;
t461 = t425 * t429;
t459 = rSges(7,2) * t467 + t487 * t375 + t488 * t376;
t458 = rSges(7,2) * t464 + t487 * t377 + t488 * t378;
t457 = -rSges(7,2) * t429 + (t487 * t418 + t488 * t420) * t427;
t430 = cos(qJ(1));
t417 = qJD(1) * t430 * pkin(1);
t456 = qJD(1) * (pkin(2) * t421 + pkin(7) * t419) + t417;
t455 = qJD(3) * t419;
t454 = qJD(3) * t421;
t453 = qJD(4) * t427;
t452 = qJD(5) * t427;
t449 = -pkin(2) * t419 + pkin(7) * t421 - t473;
t411 = pkin(3) * t427 - qJ(4) * t429;
t448 = qJD(3) * (pkin(8) * t429 - t427 * t471 - t411);
t447 = qJD(3) * (rSges(5,3) * t429 - (rSges(5,1) * t425 - rSges(5,2) * t424) * t427 - t411);
t443 = pkin(3) * t429 + qJ(4) * t427;
t400 = t443 * t421;
t446 = qJD(1) * t400 + t419 * t453 + t456;
t398 = t443 * t419;
t445 = -t398 + t449;
t444 = rSges(4,1) * t429 - rSges(4,2) * t427;
t442 = Icges(4,1) * t429 - t470;
t441 = -Icges(4,2) * t427 + t469;
t440 = Icges(4,5) * t429 - Icges(4,6) * t427;
t366 = -Icges(4,6) * t421 + t419 * t441;
t368 = -Icges(4,5) * t421 + t419 * t442;
t439 = t366 * t427 - t368 * t429;
t367 = Icges(4,6) * t419 + t421 * t441;
t369 = Icges(4,5) * t419 + t421 * t442;
t438 = -t367 * t427 + t369 * t429;
t407 = Icges(4,2) * t429 + t470;
t408 = Icges(4,1) * t427 + t469;
t437 = -t407 * t427 + t408 * t429;
t436 = -qJD(4) * t429 + t398 * t455 + t400 * t454 + qJD(2);
t434 = pkin(8) * t427 + t429 * t471;
t360 = -pkin(4) * t465 + t419 * t434;
t361 = pkin(4) * t468 + t421 * t434;
t435 = t360 * t455 + t361 * t454 + t436;
t433 = qJD(1) * t361 + t419 * t448 + t446;
t410 = t421 * t453;
t432 = t410 + (-t360 + t445) * qJD(1) + t421 * t448;
t415 = -qJD(5) * t429 + qJD(1);
t414 = rSges(2,1) * t430 - rSges(2,2) * t428;
t413 = rSges(2,1) * t428 + rSges(2,2) * t430;
t412 = rSges(4,1) * t427 + rSges(4,2) * t429;
t406 = Icges(4,5) * t427 + Icges(4,6) * t429;
t402 = t419 * t452 - t454;
t401 = t421 * t452 + t455;
t396 = -Icges(5,5) * t429 + (Icges(5,1) * t425 - Icges(5,4) * t424) * t427;
t395 = -Icges(5,6) * t429 + (Icges(5,4) * t425 - Icges(5,2) * t424) * t427;
t394 = -Icges(5,3) * t429 + (Icges(5,5) * t425 - Icges(5,6) * t424) * t427;
t393 = t421 * t461 + t468;
t392 = t419 * t425 - t421 * t462;
t391 = t419 * t461 - t465;
t390 = -t419 * t462 - t421 * t425;
t388 = t417 + qJD(1) * (rSges(3,1) * t421 - rSges(3,2) * t419);
t387 = (-rSges(3,1) * t419 - rSges(3,2) * t421 - t473) * qJD(1);
t386 = -rSges(6,3) * t429 + (rSges(6,1) * t420 - rSges(6,2) * t418) * t427;
t373 = rSges(4,3) * t419 + t421 * t444;
t372 = -rSges(4,3) * t421 + t419 * t444;
t365 = Icges(4,3) * t419 + t421 * t440;
t364 = -Icges(4,3) * t421 + t419 * t440;
t359 = rSges(5,1) * t393 + rSges(5,2) * t392 + rSges(5,3) * t464;
t358 = rSges(5,1) * t391 + rSges(5,2) * t390 + rSges(5,3) * t467;
t357 = Icges(5,1) * t393 + Icges(5,4) * t392 + Icges(5,5) * t464;
t356 = Icges(5,1) * t391 + Icges(5,4) * t390 + Icges(5,5) * t467;
t355 = Icges(5,4) * t393 + Icges(5,2) * t392 + Icges(5,6) * t464;
t354 = Icges(5,4) * t391 + Icges(5,2) * t390 + Icges(5,6) * t467;
t353 = Icges(5,5) * t393 + Icges(5,6) * t392 + Icges(5,3) * t464;
t352 = Icges(5,5) * t391 + Icges(5,6) * t390 + Icges(5,3) * t467;
t348 = rSges(6,1) * t378 - rSges(6,2) * t377 + rSges(6,3) * t464;
t346 = rSges(6,1) * t376 - rSges(6,2) * t375 + rSges(6,3) * t467;
t332 = qJD(1) * t373 - t412 * t455 + t456;
t331 = -t412 * t454 + (-t372 + t449) * qJD(1);
t330 = qJD(2) + (t372 * t419 + t373 * t421) * qJD(3);
t329 = qJD(1) * t359 + t419 * t447 + t446;
t328 = t410 + t421 * t447 + (-t358 + t445) * qJD(1);
t327 = (t358 * t419 + t359 * t421) * qJD(3) + t436;
t326 = t348 * t415 - t386 * t401 + t433;
t325 = -t346 * t415 + t386 * t402 + t432;
t324 = t346 * t401 - t348 * t402 + t435;
t323 = qJD(6) * t375 - t401 * t457 + t415 * t458 + t433;
t322 = qJD(6) * t377 + t402 * t457 - t415 * t459 + t432;
t321 = qJD(6) * t418 * t427 + t401 * t459 - t402 * t458 + t435;
t1 = m(7) * (t321 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t387 ^ 2 + t388 ^ 2) / 0.2e1 + m(4) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(5) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + m(6) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + ((t480 * t377 + t478 * t378 + t479 * t464) * t415 + (t486 * t377 + t482 * t378 + t484 * t464) * t402 + (t485 * t377 + t481 * t378 + t483 * t464) * t401) * t401 / 0.2e1 + ((t480 * t375 + t478 * t376 + t479 * t467) * t415 + (t486 * t375 + t482 * t376 + t484 * t467) * t402 + (t485 * t375 + t481 * t376 + t483 * t467) * t401) * t402 / 0.2e1 + ((-t483 * t401 - t484 * t402 - t479 * t415) * t429 + ((t480 * t418 + t478 * t420) * t415 + (t486 * t418 + t482 * t420) * t402 + (t485 * t418 + t481 * t420) * t401) * t427) * t415 / 0.2e1 + (((t367 * t429 + t369 * t427) * t419 - (t366 * t429 + t368 * t427) * t421 + (t352 * t421 - t353 * t419) * t429 + ((-t355 * t424 + t357 * t425) * t419 - (-t354 * t424 + t356 * t425) * t421) * t427) * qJD(3) + ((t407 - t394) * t429 + (-t395 * t424 + t396 * t425 + t408) * t427) * qJD(1)) * qJD(1) / 0.2e1 + (((-t352 * t464 - t354 * t392 - t356 * t393 + t439 * t421) * t421 + ((-t364 + t438) * t421 + t353 * t464 + t355 * t392 + t357 * t393 + t365 * t419) * t419) * qJD(3) + (t392 * t395 + t393 * t396 + t394 * t464 + t419 * t406 + t421 * t437) * qJD(1)) * t455 / 0.2e1 - (((-t352 * t467 - t354 * t390 - t356 * t391 + t364 * t421) * t421 + ((-t365 + t439) * t421 + t353 * t467 + t355 * t390 + t357 * t391 + t438 * t419) * t419) * qJD(3) + (t390 * t395 + t391 * t396 + t394 * t467 - t421 * t406 + t419 * t437) * qJD(1)) * t454 / 0.2e1 + (m(2) * (t413 ^ 2 + t414 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
