% Calculate kinetic energy for
% S6RPRRRP5
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:44
% EndTime: 2019-03-09 06:10:45
% DurationCPUTime: 1.89s
% Computational Cost: add. (1745->225), mult. (1650->354), div. (0->0), fcn. (1594->10), ass. (0->130)
t511 = Icges(6,1) + Icges(7,1);
t510 = -Icges(6,4) + Icges(7,5);
t509 = Icges(7,4) + Icges(6,5);
t508 = Icges(6,2) + Icges(7,3);
t507 = -Icges(7,6) + Icges(6,6);
t506 = -Icges(6,3) - Icges(7,2);
t505 = rSges(7,1) + pkin(5);
t504 = rSges(7,3) + qJ(6);
t437 = pkin(10) + qJ(3);
t431 = qJ(4) + t437;
t426 = cos(t431);
t443 = cos(qJ(5));
t444 = cos(qJ(1));
t479 = t443 * t444;
t441 = sin(qJ(5));
t442 = sin(qJ(1));
t482 = t441 * t442;
t404 = t426 * t482 + t479;
t480 = t442 * t443;
t481 = t441 * t444;
t405 = t426 * t480 - t481;
t425 = sin(t431);
t484 = t425 * t442;
t503 = t508 * t404 + t510 * t405 - t507 * t484;
t406 = t426 * t481 - t480;
t407 = t426 * t479 + t482;
t483 = t425 * t444;
t502 = t508 * t406 + t510 * t407 - t507 * t483;
t501 = -t507 * t404 + t509 * t405 - t506 * t484;
t500 = -t507 * t406 + t509 * t407 - t506 * t483;
t499 = t510 * t404 + t511 * t405 + t509 * t484;
t498 = t510 * t406 + t511 * t407 + t509 * t483;
t497 = t507 * t426 + (t508 * t441 + t510 * t443) * t425;
t496 = t506 * t426 + (-t507 * t441 + t509 * t443) * t425;
t495 = -t509 * t426 + (t510 * t441 + t511 * t443) * t425;
t439 = cos(pkin(10));
t489 = t439 * pkin(2);
t429 = sin(t437);
t488 = Icges(4,4) * t429;
t430 = cos(t437);
t487 = Icges(4,4) * t430;
t486 = Icges(5,4) * t425;
t485 = Icges(5,4) * t426;
t477 = rSges(7,2) * t484 + t504 * t404 + t505 * t405;
t476 = rSges(7,2) * t483 + t504 * t406 + t505 * t407;
t472 = pkin(3) * t430;
t366 = -pkin(8) * t444 + t442 * t472;
t367 = pkin(8) * t442 + t444 * t472;
t433 = qJD(3) * t442;
t470 = qJD(3) * t444;
t475 = t366 * t433 + t367 * t470;
t474 = -rSges(7,2) * t426 + (t504 * t441 + t505 * t443) * t425;
t422 = pkin(1) * t442 - qJ(2) * t444;
t473 = pkin(7) * t444 - t442 * t489 - t422;
t420 = qJD(4) * t442 + t433;
t469 = qJD(5) * t425;
t468 = pkin(3) * qJD(3) * t429;
t467 = -t366 + t473;
t421 = (-qJD(3) - qJD(4)) * t444;
t466 = pkin(4) * t426 + pkin(9) * t425;
t418 = qJD(1) * (pkin(1) * t444 + qJ(2) * t442);
t465 = -qJD(2) * t444 + qJD(1) * (pkin(7) * t442 + t444 * t489) + t418;
t438 = sin(pkin(10));
t464 = rSges(3,1) * t439 - rSges(3,2) * t438;
t463 = rSges(4,1) * t430 - rSges(4,2) * t429;
t462 = rSges(5,1) * t426 - rSges(5,2) * t425;
t461 = Icges(4,1) * t430 - t488;
t460 = Icges(5,1) * t426 - t486;
t459 = -Icges(4,2) * t429 + t487;
t458 = -Icges(5,2) * t425 + t485;
t457 = Icges(4,5) * t430 - Icges(4,6) * t429;
t456 = Icges(5,5) * t426 - Icges(5,6) * t425;
t392 = -Icges(4,6) * t444 + t442 * t459;
t394 = -Icges(4,5) * t444 + t442 * t461;
t455 = t392 * t429 - t394 * t430;
t393 = Icges(4,6) * t442 + t444 * t459;
t395 = Icges(4,5) * t442 + t444 * t461;
t454 = -t393 * t429 + t395 * t430;
t414 = Icges(4,2) * t430 + t488;
t415 = Icges(4,1) * t429 + t487;
t453 = -t414 * t429 + t415 * t430;
t434 = qJD(2) * t442;
t452 = -t444 * t468 + t434;
t400 = t466 * t442;
t401 = t466 * t444;
t451 = t420 * t400 - t401 * t421 + t475;
t450 = (Icges(5,5) * t425 + Icges(5,6) * t426) * qJD(1) + (-Icges(5,3) * t444 + t442 * t456) * t421 + (Icges(5,3) * t442 + t444 * t456) * t420;
t449 = qJD(1) * t367 - t442 * t468 + t465;
t412 = pkin(4) * t425 - pkin(9) * t426;
t448 = t421 * t412 + (-t400 + t467) * qJD(1) + t452;
t447 = qJD(1) * t401 - t412 * t420 + t449;
t383 = -Icges(5,6) * t444 + t442 * t458;
t384 = Icges(5,6) * t442 + t444 * t458;
t385 = -Icges(5,5) * t444 + t442 * t460;
t386 = Icges(5,5) * t442 + t444 * t460;
t409 = Icges(5,2) * t426 + t486;
t410 = Icges(5,1) * t425 + t485;
t446 = (-t384 * t425 + t386 * t426) * t420 + (-t383 * t425 + t385 * t426) * t421 + (-t409 * t425 + t410 * t426) * qJD(1);
t424 = rSges(2,1) * t444 - rSges(2,2) * t442;
t423 = rSges(2,1) * t442 + rSges(2,2) * t444;
t419 = -qJD(5) * t426 + qJD(1);
t416 = rSges(4,1) * t429 + rSges(4,2) * t430;
t413 = Icges(4,5) * t429 + Icges(4,6) * t430;
t411 = rSges(5,1) * t425 + rSges(5,2) * t426;
t403 = t442 * t469 + t421;
t402 = t444 * t469 + t420;
t397 = rSges(4,3) * t442 + t444 * t463;
t396 = -rSges(4,3) * t444 + t442 * t463;
t391 = Icges(4,3) * t442 + t444 * t457;
t390 = -Icges(4,3) * t444 + t442 * t457;
t388 = rSges(5,3) * t442 + t444 * t462;
t387 = -rSges(5,3) * t444 + t442 * t462;
t378 = -rSges(6,3) * t426 + (rSges(6,1) * t443 - rSges(6,2) * t441) * t425;
t369 = qJD(1) * t442 * rSges(3,3) + t418 + (qJD(1) * t464 - qJD(2)) * t444;
t368 = t434 + (t444 * rSges(3,3) - t442 * t464 - t422) * qJD(1);
t360 = rSges(6,1) * t407 - rSges(6,2) * t406 + rSges(6,3) * t483;
t358 = rSges(6,1) * t405 - rSges(6,2) * t404 + rSges(6,3) * t484;
t344 = (t396 * t442 + t397 * t444) * qJD(3);
t343 = qJD(1) * t397 - t416 * t433 + t465;
t342 = -t416 * t470 + t434 + (-t396 + t473) * qJD(1);
t341 = qJD(1) * t388 - t411 * t420 + t449;
t340 = t411 * t421 + (-t387 + t467) * qJD(1) + t452;
t339 = t387 * t420 - t388 * t421 + t475;
t338 = t360 * t419 - t378 * t402 + t447;
t337 = -t358 * t419 + t378 * t403 + t448;
t336 = t358 * t402 - t360 * t403 + t451;
t335 = qJD(6) * t404 - t402 * t474 + t419 * t476 + t447;
t334 = qJD(6) * t406 + t403 * t474 - t419 * t477 + t448;
t333 = qJD(6) * t425 * t441 + t402 * t477 - t403 * t476 + t451;
t1 = m(4) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(3) * (t368 ^ 2 + t369 ^ 2) / 0.2e1 + t421 * (t442 * t446 - t450 * t444) / 0.2e1 + t420 * (t442 * t450 + t444 * t446) / 0.2e1 - ((-t444 * t413 + t442 * t453) * qJD(1) + (t444 ^ 2 * t390 + (t454 * t442 + (-t391 + t455) * t444) * t442) * qJD(3)) * t470 / 0.2e1 + m(7) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(6) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(5) * (t339 ^ 2 + t340 ^ 2 + t341 ^ 2) / 0.2e1 + ((t442 * t413 + t444 * t453) * qJD(1) + (t442 ^ 2 * t391 + (t455 * t444 + (-t390 + t454) * t442) * t444) * qJD(3)) * t433 / 0.2e1 + ((t497 * t406 + t495 * t407 + t496 * t483) * t419 + (t503 * t406 + t499 * t407 + t501 * t483) * t403 + (t502 * t406 + t498 * t407 + t500 * t483) * t402) * t402 / 0.2e1 + ((t497 * t404 + t495 * t405 + t496 * t484) * t419 + (t503 * t404 + t499 * t405 + t501 * t484) * t403 + (t502 * t404 + t498 * t405 + t500 * t484) * t402) * t403 / 0.2e1 + ((-t500 * t402 - t501 * t403 - t496 * t419) * t426 + ((t497 * t441 + t495 * t443) * t419 + (t503 * t441 + t499 * t443) * t403 + (t502 * t441 + t498 * t443) * t402) * t425) * t419 / 0.2e1 + ((t384 * t426 + t386 * t425) * t420 + (t383 * t426 + t385 * t425) * t421 + ((t393 * t430 + t395 * t429) * t442 - (t392 * t430 + t394 * t429) * t444) * qJD(3) + (t426 * t409 + t425 * t410 + t430 * t414 + t429 * t415) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t423 ^ 2 + t424 ^ 2) + Icges(3,2) * t439 ^ 2 + (Icges(3,1) * t438 + 0.2e1 * Icges(3,4) * t439) * t438 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
