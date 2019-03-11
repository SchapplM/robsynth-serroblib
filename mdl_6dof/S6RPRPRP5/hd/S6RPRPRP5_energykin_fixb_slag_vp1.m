% Calculate kinetic energy for
% S6RPRPRP5
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
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:14:25
% EndTime: 2019-03-09 03:14:27
% DurationCPUTime: 2.17s
% Computational Cost: add. (1740->250), mult. (1849->381), div. (0->0), fcn. (1849->10), ass. (0->132)
t520 = Icges(6,1) + Icges(7,1);
t519 = -Icges(6,4) + Icges(7,5);
t518 = Icges(7,4) + Icges(6,5);
t517 = Icges(6,2) + Icges(7,3);
t516 = -Icges(7,6) + Icges(6,6);
t515 = -Icges(6,3) - Icges(7,2);
t514 = rSges(7,1) + pkin(5);
t513 = rSges(7,3) + qJ(6);
t443 = pkin(10) + qJ(5);
t440 = cos(t443);
t444 = pkin(9) + qJ(3);
t441 = cos(t444);
t452 = cos(qJ(1));
t438 = sin(t443);
t451 = sin(qJ(1));
t494 = t438 * t451;
t410 = t440 * t452 + t441 * t494;
t486 = t451 * t440;
t411 = -t438 * t452 + t441 * t486;
t439 = sin(t444);
t493 = t439 * t451;
t512 = t517 * t410 + t519 * t411 - t516 * t493;
t491 = t441 * t452;
t412 = t438 * t491 - t486;
t413 = t440 * t491 + t494;
t492 = t439 * t452;
t511 = t517 * t412 + t519 * t413 - t516 * t492;
t510 = -t516 * t410 + t518 * t411 - t515 * t493;
t509 = -t516 * t412 + t518 * t413 - t515 * t492;
t508 = t519 * t410 + t520 * t411 + t518 * t493;
t507 = t519 * t412 + t520 * t413 + t518 * t492;
t506 = t516 * t441 + (t517 * t438 + t519 * t440) * t439;
t505 = t515 * t441 + (-t516 * t438 + t518 * t440) * t439;
t504 = -t518 * t441 + (t519 * t438 + t520 * t440) * t439;
t448 = cos(pkin(9));
t498 = pkin(2) * t448;
t447 = cos(pkin(10));
t497 = pkin(4) * t447;
t496 = Icges(4,4) * t439;
t495 = Icges(4,4) * t441;
t445 = sin(pkin(10));
t490 = t445 * t451;
t489 = t445 * t452;
t488 = t447 * t451;
t487 = t447 * t452;
t483 = rSges(7,2) * t493 + t513 * t410 + t514 * t411;
t482 = rSges(7,2) * t492 + t513 * t412 + t514 * t413;
t481 = -rSges(7,2) * t441 + (t513 * t438 + t514 * t440) * t439;
t431 = pkin(1) * t451 - qJ(2) * t452;
t480 = pkin(7) * t452 - t451 * t498 - t431;
t442 = qJD(2) * t451;
t476 = qJD(4) * t439;
t479 = t452 * t476 + t442;
t478 = qJD(3) * t451;
t477 = qJD(3) * t452;
t475 = qJD(5) * t439;
t465 = pkin(3) * t441 + qJ(4) * t439;
t415 = t465 * t451;
t474 = -t415 + t480;
t426 = pkin(3) * t439 - qJ(4) * t441;
t471 = qJD(3) * (pkin(8) * t441 - t439 * t497 - t426);
t470 = qJD(3) * (rSges(5,3) * t441 - (rSges(5,1) * t447 - rSges(5,2) * t445) * t439 - t426);
t428 = qJD(1) * (pkin(1) * t452 + qJ(2) * t451);
t469 = -qJD(2) * t452 + qJD(1) * (pkin(7) * t451 + t452 * t498) + t428;
t416 = t465 * t452;
t468 = -qJD(4) * t441 + t415 * t478 + t416 * t477;
t446 = sin(pkin(9));
t467 = rSges(3,1) * t448 - rSges(3,2) * t446;
t466 = rSges(4,1) * t441 - rSges(4,2) * t439;
t464 = Icges(4,1) * t441 - t496;
t463 = -Icges(4,2) * t439 + t495;
t462 = Icges(4,5) * t441 - Icges(4,6) * t439;
t401 = -Icges(4,6) * t452 + t451 * t463;
t403 = -Icges(4,5) * t452 + t451 * t464;
t461 = t401 * t439 - t403 * t441;
t402 = Icges(4,6) * t451 + t452 * t463;
t404 = Icges(4,5) * t451 + t452 * t464;
t460 = -t402 * t439 + t404 * t441;
t424 = Icges(4,2) * t441 + t496;
t425 = Icges(4,1) * t439 + t495;
t459 = -t424 * t439 + t425 * t441;
t458 = qJD(1) * t416 + t451 * t476 + t469;
t456 = pkin(8) * t439 + t441 * t497;
t370 = -pkin(4) * t489 + t451 * t456;
t371 = pkin(4) * t490 + t452 * t456;
t457 = t370 * t478 + t371 * t477 + t468;
t455 = qJD(1) * t371 + t451 * t471 + t458;
t454 = (-t370 + t474) * qJD(1) + t452 * t471 + t479;
t434 = -qJD(5) * t441 + qJD(1);
t433 = rSges(2,1) * t452 - rSges(2,2) * t451;
t432 = rSges(2,1) * t451 + rSges(2,2) * t452;
t427 = rSges(4,1) * t439 + rSges(4,2) * t441;
t423 = Icges(4,5) * t439 + Icges(4,6) * t441;
t422 = t451 * t475 - t477;
t421 = t452 * t475 + t478;
t420 = t441 * t487 + t490;
t419 = -t441 * t489 + t488;
t418 = t441 * t488 - t489;
t417 = -t441 * t490 - t487;
t409 = rSges(4,3) * t451 + t452 * t466;
t408 = -rSges(4,3) * t452 + t451 * t466;
t400 = Icges(4,3) * t451 + t452 * t462;
t399 = -Icges(4,3) * t452 + t451 * t462;
t395 = -Icges(5,5) * t441 + (Icges(5,1) * t447 - Icges(5,4) * t445) * t439;
t394 = -Icges(5,6) * t441 + (Icges(5,4) * t447 - Icges(5,2) * t445) * t439;
t393 = -Icges(5,3) * t441 + (Icges(5,5) * t447 - Icges(5,6) * t445) * t439;
t392 = -rSges(6,3) * t441 + (rSges(6,1) * t440 - rSges(6,2) * t438) * t439;
t383 = qJD(1) * t451 * rSges(3,3) + t428 + (qJD(1) * t467 - qJD(2)) * t452;
t382 = t442 + (t452 * rSges(3,3) - t451 * t467 - t431) * qJD(1);
t379 = rSges(5,1) * t420 + rSges(5,2) * t419 + rSges(5,3) * t492;
t378 = rSges(5,1) * t418 + rSges(5,2) * t417 + rSges(5,3) * t493;
t377 = Icges(5,1) * t420 + Icges(5,4) * t419 + Icges(5,5) * t492;
t376 = Icges(5,1) * t418 + Icges(5,4) * t417 + Icges(5,5) * t493;
t375 = Icges(5,4) * t420 + Icges(5,2) * t419 + Icges(5,6) * t492;
t374 = Icges(5,4) * t418 + Icges(5,2) * t417 + Icges(5,6) * t493;
t373 = Icges(5,5) * t420 + Icges(5,6) * t419 + Icges(5,3) * t492;
t372 = Icges(5,5) * t418 + Icges(5,6) * t417 + Icges(5,3) * t493;
t366 = (t408 * t451 + t409 * t452) * qJD(3);
t365 = rSges(6,1) * t413 - rSges(6,2) * t412 + rSges(6,3) * t492;
t363 = rSges(6,1) * t411 - rSges(6,2) * t410 + rSges(6,3) * t493;
t349 = qJD(1) * t409 - t427 * t478 + t469;
t348 = -t427 * t477 + t442 + (-t408 + t480) * qJD(1);
t347 = (t378 * t451 + t379 * t452) * qJD(3) + t468;
t346 = qJD(1) * t379 + t451 * t470 + t458;
t345 = t452 * t470 + (-t378 + t474) * qJD(1) + t479;
t344 = t365 * t434 - t392 * t421 + t455;
t343 = -t363 * t434 + t392 * t422 + t454;
t342 = t363 * t421 - t365 * t422 + t457;
t341 = qJD(6) * t410 - t421 * t481 + t434 * t482 + t455;
t340 = qJD(6) * t412 + t422 * t481 - t434 * t483 + t454;
t339 = qJD(6) * t438 * t439 + t421 * t483 - t422 * t482 + t457;
t1 = m(3) * (t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(4) * (t348 ^ 2 + t349 ^ 2 + t366 ^ 2) / 0.2e1 + m(5) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + m(6) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(7) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + ((t506 * t412 + t504 * t413 + t505 * t492) * t434 + (t512 * t412 + t508 * t413 + t510 * t492) * t422 + (t511 * t412 + t507 * t413 + t509 * t492) * t421) * t421 / 0.2e1 + ((t506 * t410 + t504 * t411 + t505 * t493) * t434 + (t512 * t410 + t508 * t411 + t510 * t493) * t422 + (t511 * t410 + t507 * t411 + t509 * t493) * t421) * t422 / 0.2e1 + ((-t509 * t421 - t510 * t422 - t505 * t434) * t441 + ((t506 * t438 + t504 * t440) * t434 + (t512 * t438 + t508 * t440) * t422 + (t511 * t438 + t507 * t440) * t421) * t439) * t434 / 0.2e1 + (((t402 * t441 + t404 * t439) * t451 - (t401 * t441 + t439 * t403) * t452 + (t372 * t452 - t373 * t451) * t441 + ((-t375 * t445 + t377 * t447) * t451 - (-t374 * t445 + t376 * t447) * t452) * t439) * qJD(3) + ((t424 - t393) * t441 + (-t394 * t445 + t395 * t447 + t425) * t439) * qJD(1)) * qJD(1) / 0.2e1 + (((-t372 * t492 - t374 * t419 - t376 * t420 + t461 * t452) * t452 + (t373 * t492 + t375 * t419 + t377 * t420 + (-t399 + t460) * t452 + t400 * t451) * t451) * qJD(3) + (t393 * t492 + t394 * t419 + t395 * t420 + t451 * t423 + t452 * t459) * qJD(1)) * t478 / 0.2e1 - (((-t372 * t493 - t374 * t417 - t376 * t418 + t399 * t452) * t452 + ((-t400 + t461) * t452 + t373 * t493 + t375 * t417 + t377 * t418 + t460 * t451) * t451) * qJD(3) + (t393 * t493 + t394 * t417 + t395 * t418 - t452 * t423 + t451 * t459) * qJD(1)) * t477 / 0.2e1 + (m(2) * (t432 ^ 2 + t433 ^ 2) + Icges(2,3) + Icges(3,2) * t448 ^ 2 + (Icges(3,1) * t446 + 0.2e1 * Icges(3,4) * t448) * t446) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
