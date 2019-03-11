% Calculate kinetic energy for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:19
% EndTime: 2019-03-09 08:21:21
% DurationCPUTime: 2.31s
% Computational Cost: add. (1075->266), mult. (2645->394), div. (0->0), fcn. (2844->8), ass. (0->134)
t510 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t509 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t508 = Icges(4,4) - Icges(6,5) + Icges(5,6);
t507 = Icges(4,6) + Icges(6,4) - Icges(5,5);
t506 = Icges(6,6) - Icges(4,5) + Icges(5,4);
t505 = -Icges(4,3) - Icges(6,2) - Icges(5,1);
t445 = sin(pkin(9));
t449 = sin(qJ(1));
t451 = cos(qJ(2));
t486 = t449 * t451;
t446 = cos(pkin(9));
t452 = cos(qJ(1));
t489 = t446 * t452;
t418 = t445 * t486 + t489;
t485 = t452 * t445;
t419 = t446 * t486 - t485;
t448 = sin(qJ(2));
t488 = t448 * t449;
t504 = -t418 * t509 + t419 * t508 + t488 * t507;
t420 = -t449 * t446 + t451 * t485;
t421 = t445 * t449 + t451 * t489;
t487 = t448 * t452;
t503 = t420 * t509 - t421 * t508 - t487 * t507;
t502 = t508 * t418 - t419 * t510 + t506 * t488;
t501 = -t508 * t420 + t421 * t510 - t506 * t487;
t500 = t418 * t507 + t419 * t506 + t488 * t505;
t499 = -t420 * t507 - t421 * t506 - t487 * t505;
t498 = t506 * t451 + (-t508 * t445 + t446 * t510) * t448;
t497 = t505 * t451 + (-t445 * t507 - t446 * t506) * t448;
t496 = t507 * t451 + (t445 * t509 - t446 * t508) * t448;
t493 = Icges(3,4) * t448;
t492 = Icges(3,4) * t451;
t491 = t445 * t448;
t490 = t446 * t448;
t479 = qJD(3) * t448;
t441 = t452 * t479;
t484 = qJD(4) * t420 + t441;
t433 = pkin(2) * t448 - qJ(3) * t451;
t483 = -(pkin(3) * t446 + qJ(4) * t445) * t448 - t433;
t466 = pkin(2) * t451 + qJ(3) * t448;
t423 = t466 * t449;
t437 = pkin(1) * t449 - pkin(7) * t452;
t482 = -t423 - t437;
t481 = qJD(2) * t449;
t480 = qJD(2) * t452;
t478 = qJD(6) * t448;
t380 = pkin(3) * t419 + qJ(4) * t418;
t477 = -t380 + t482;
t476 = qJD(5) * t421 + t484;
t424 = t466 * t452;
t429 = qJD(1) * (pkin(1) * t452 + pkin(7) * t449);
t475 = qJD(1) * t424 + t449 * t479 + t429;
t474 = pkin(4) * t451 - qJ(5) * t490 + t483;
t471 = qJD(2) * (rSges(4,3) * t451 - (rSges(4,1) * t446 - rSges(4,2) * t445) * t448 - t433);
t385 = pkin(4) * t488 + qJ(5) * t419;
t470 = -t385 + t477;
t469 = qJD(2) * (rSges(5,1) * t451 - (-rSges(5,2) * t446 + rSges(5,3) * t445) * t448 + t483);
t468 = -qJD(3) * t451 + t423 * t481 + t424 * t480;
t467 = rSges(3,1) * t451 - rSges(3,2) * t448;
t381 = pkin(3) * t421 + qJ(4) * t420;
t465 = qJD(1) * t381 + qJD(4) * t418 + t475;
t464 = Icges(3,1) * t451 - t493;
t463 = -Icges(3,2) * t448 + t492;
t462 = Icges(3,5) * t451 - Icges(3,6) * t448;
t403 = -Icges(3,6) * t452 + t449 * t463;
t405 = -Icges(3,5) * t452 + t449 * t464;
t461 = t403 * t448 - t405 * t451;
t404 = Icges(3,6) * t449 + t452 * t463;
t406 = Icges(3,5) * t449 + t452 * t464;
t460 = -t404 * t448 + t406 * t451;
t431 = Icges(3,2) * t451 + t493;
t432 = Icges(3,1) * t448 + t492;
t459 = -t431 * t448 + t432 * t451;
t458 = qJD(2) * (-rSges(6,2) * t451 - (rSges(6,1) * t445 + rSges(6,3) * t446) * t448 + t474);
t457 = qJD(2) * (-pkin(5) * t491 + pkin(8) * t451 + t474);
t386 = pkin(4) * t487 + qJ(5) * t421;
t456 = qJD(1) * t386 + qJD(5) * t419 + t465;
t455 = qJD(4) * t491 + t380 * t481 + t381 * t480 + t468;
t454 = qJD(5) * t490 + t385 * t481 + t386 * t480 + t455;
t450 = cos(qJ(6));
t447 = sin(qJ(6));
t442 = -qJD(6) * t451 + qJD(1);
t436 = rSges(2,1) * t452 - rSges(2,2) * t449;
t435 = rSges(2,1) * t449 + rSges(2,2) * t452;
t434 = rSges(3,1) * t448 + rSges(3,2) * t451;
t430 = Icges(3,5) * t448 + Icges(3,6) * t451;
t426 = t449 * t478 - t480;
t425 = t452 * t478 + t481;
t412 = (t445 * t450 + t446 * t447) * t448;
t411 = (-t445 * t447 + t446 * t450) * t448;
t410 = rSges(3,3) * t449 + t452 * t467;
t409 = -rSges(3,3) * t452 + t449 * t467;
t402 = Icges(3,3) * t449 + t452 * t462;
t401 = -Icges(3,3) * t452 + t449 * t462;
t388 = pkin(5) * t420 + pkin(8) * t487;
t387 = pkin(5) * t418 + pkin(8) * t488;
t378 = t420 * t450 + t421 * t447;
t377 = -t420 * t447 + t421 * t450;
t376 = t418 * t450 + t419 * t447;
t375 = -t418 * t447 + t419 * t450;
t372 = rSges(4,1) * t421 - rSges(4,2) * t420 + rSges(4,3) * t487;
t371 = rSges(6,1) * t420 - rSges(6,2) * t487 + rSges(6,3) * t421;
t370 = rSges(4,1) * t419 - rSges(4,2) * t418 + rSges(4,3) * t488;
t369 = rSges(6,1) * t418 - rSges(6,2) * t488 + rSges(6,3) * t419;
t368 = rSges(5,1) * t487 - rSges(5,2) * t421 + rSges(5,3) * t420;
t367 = rSges(5,1) * t488 - rSges(5,2) * t419 + rSges(5,3) * t418;
t348 = rSges(7,1) * t412 + rSges(7,2) * t411 - rSges(7,3) * t451;
t347 = Icges(7,1) * t412 + Icges(7,4) * t411 - Icges(7,5) * t451;
t346 = Icges(7,4) * t412 + Icges(7,2) * t411 - Icges(7,6) * t451;
t345 = Icges(7,5) * t412 + Icges(7,6) * t411 - Icges(7,3) * t451;
t344 = qJD(1) * t410 - t434 * t481 + t429;
t343 = -t434 * t480 + (-t409 - t437) * qJD(1);
t342 = (t409 * t449 + t410 * t452) * qJD(2);
t341 = rSges(7,1) * t378 + rSges(7,2) * t377 + rSges(7,3) * t487;
t340 = rSges(7,1) * t376 + rSges(7,2) * t375 + rSges(7,3) * t488;
t339 = Icges(7,1) * t378 + Icges(7,4) * t377 + Icges(7,5) * t487;
t338 = Icges(7,1) * t376 + Icges(7,4) * t375 + Icges(7,5) * t488;
t337 = Icges(7,4) * t378 + Icges(7,2) * t377 + Icges(7,6) * t487;
t336 = Icges(7,4) * t376 + Icges(7,2) * t375 + Icges(7,6) * t488;
t335 = Icges(7,5) * t378 + Icges(7,6) * t377 + Icges(7,3) * t487;
t334 = Icges(7,5) * t376 + Icges(7,6) * t375 + Icges(7,3) * t488;
t333 = qJD(1) * t372 + t449 * t471 + t475;
t332 = t441 + t452 * t471 + (-t370 + t482) * qJD(1);
t331 = (t370 * t449 + t372 * t452) * qJD(2) + t468;
t330 = qJD(1) * t368 + t449 * t469 + t465;
t329 = t452 * t469 + (-t367 + t477) * qJD(1) + t484;
t328 = (t367 * t449 + t368 * t452) * qJD(2) + t455;
t327 = qJD(1) * t371 + t449 * t458 + t456;
t326 = t452 * t458 + (-t369 + t470) * qJD(1) + t476;
t325 = (t369 * t449 + t371 * t452) * qJD(2) + t454;
t324 = qJD(1) * t388 + t341 * t442 - t348 * t425 + t449 * t457 + t456;
t323 = -t340 * t442 + t348 * t426 + t452 * t457 + (-t387 + t470) * qJD(1) + t476;
t322 = t340 * t425 - t341 * t426 + (t387 * t449 + t388 * t452) * qJD(2) + t454;
t1 = m(6) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(7) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + m(3) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(4) * (t331 ^ 2 + t332 ^ 2 + t333 ^ 2) / 0.2e1 + m(5) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + t442 * ((-t335 * t451 + t337 * t411 + t339 * t412) * t425 + (-t334 * t451 + t336 * t411 + t338 * t412) * t426 + (-t451 * t345 + t411 * t346 + t412 * t347) * t442) / 0.2e1 + t425 * ((t335 * t487 + t377 * t337 + t378 * t339) * t425 + (t334 * t487 + t336 * t377 + t338 * t378) * t426 + (t345 * t487 + t346 * t377 + t347 * t378) * t442) / 0.2e1 + t426 * ((t335 * t488 + t337 * t375 + t339 * t376) * t425 + (t334 * t488 + t375 * t336 + t376 * t338) * t426 + (t345 * t488 + t346 * t375 + t347 * t376) * t442) / 0.2e1 + (m(2) * (t435 ^ 2 + t436 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((((-t403 - t500) * t452 + (t404 - t499) * t449) * t451 + ((t445 * t504 + t502 * t446 - t405) * t452 + (t445 * t503 + t446 * t501 + t406) * t449) * t448) * qJD(2) + ((t431 - t497) * t451 + (t445 * t496 + t446 * t498 + t432) * t448) * qJD(1)) * qJD(1) / 0.2e1 + (((t504 * t420 + t502 * t421 + t461 * t452 + t500 * t487) * t452 + ((-t401 + t460) * t452 + t449 * t402 + t499 * t487 + t501 * t421 + t503 * t420) * t449) * qJD(2) + (t420 * t496 + t421 * t498 + t449 * t430 + t452 * t459 + t487 * t497) * qJD(1)) * t481 / 0.2e1 - (((t452 * t401 + t504 * t418 + t502 * t419 + t500 * t488) * t452 + (t460 * t449 + (-t402 + t461) * t452 + t499 * t488 + t501 * t419 + t503 * t418) * t449) * qJD(2) + (t418 * t496 + t419 * t498 - t452 * t430 + t449 * t459 + t488 * t497) * qJD(1)) * t480 / 0.2e1;
T  = t1;
