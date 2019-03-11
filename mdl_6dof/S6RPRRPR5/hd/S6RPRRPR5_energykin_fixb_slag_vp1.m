% Calculate kinetic energy for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:12:10
% EndTime: 2019-03-09 05:12:12
% DurationCPUTime: 2.59s
% Computational Cost: add. (1487->232), mult. (1390->360), div. (0->0), fcn. (1268->10), ass. (0->134)
t519 = Icges(5,4) + Icges(6,6);
t518 = Icges(5,1) + Icges(6,2);
t517 = -Icges(5,2) - Icges(6,3);
t436 = pkin(10) + qJ(3);
t430 = qJ(4) + t436;
t426 = cos(t430);
t516 = t519 * t426;
t425 = sin(t430);
t515 = t519 * t425;
t514 = Icges(6,4) - Icges(5,5);
t513 = Icges(6,5) - Icges(5,6);
t512 = t517 * t425 + t516;
t511 = t518 * t426 - t515;
t510 = Icges(6,1) + Icges(5,3);
t441 = sin(qJ(1));
t443 = cos(qJ(1));
t509 = t512 * t441 + t513 * t443;
t508 = -t513 * t441 + t512 * t443;
t507 = t511 * t441 + t514 * t443;
t506 = -t514 * t441 + t511 * t443;
t505 = t517 * t426 - t515;
t504 = t518 * t425 + t516;
t503 = t513 * t425 - t514 * t426;
t432 = qJD(3) * t441;
t420 = qJD(4) * t441 + t432;
t421 = (-qJD(3) - qJD(4)) * t443;
t502 = (-t509 * t425 + t507 * t426) * t421 + (-t508 * t425 + t506 * t426) * t420 + (t505 * t425 + t504 * t426) * qJD(1);
t501 = (t503 * t441 - t510 * t443) * t421 + (t510 * t441 + t503 * t443) * t420 + (-t514 * t425 - t513 * t426) * qJD(1);
t497 = pkin(9) * t425;
t438 = cos(pkin(10));
t495 = t438 * pkin(2);
t428 = sin(t436);
t494 = Icges(4,4) * t428;
t429 = cos(t436);
t493 = Icges(4,4) * t429;
t488 = t426 * t441;
t487 = t426 * t443;
t440 = sin(qJ(6));
t486 = t440 * t441;
t485 = t440 * t443;
t442 = cos(qJ(6));
t484 = t441 * t442;
t483 = t442 * t443;
t479 = pkin(3) * t429;
t355 = -pkin(8) * t443 + t441 * t479;
t356 = pkin(8) * t441 + t443 * t479;
t477 = qJD(3) * t443;
t481 = t355 * t432 + t356 * t477;
t422 = pkin(1) * t441 - qJ(2) * t443;
t480 = pkin(7) * t443 - t441 * t495 - t422;
t476 = qJD(5) * t425;
t475 = qJD(6) * t426;
t474 = pkin(3) * qJD(3) * t428;
t473 = -t355 + t480;
t466 = pkin(4) * t426 + qJ(5) * t425;
t392 = t466 * t441;
t472 = -t392 + t473;
t416 = qJD(1) * (pkin(1) * t443 + qJ(2) * t441);
t471 = -qJD(2) * t443 + qJD(1) * (pkin(7) * t441 + t443 * t495) + t416;
t437 = sin(pkin(10));
t470 = rSges(3,1) * t438 - rSges(3,2) * t437;
t469 = rSges(4,1) * t429 - rSges(4,2) * t428;
t468 = rSges(5,1) * t426 - rSges(5,2) * t425;
t467 = -rSges(6,2) * t426 + rSges(6,3) * t425;
t465 = Icges(4,1) * t429 - t494;
t463 = -Icges(4,2) * t428 + t493;
t460 = Icges(4,5) * t429 - Icges(4,6) * t428;
t385 = -Icges(4,6) * t443 + t441 * t463;
t387 = -Icges(4,5) * t443 + t441 * t465;
t456 = t385 * t428 - t387 * t429;
t386 = Icges(4,6) * t441 + t443 * t463;
t388 = Icges(4,5) * t441 + t443 * t465;
t455 = -t386 * t428 + t388 * t429;
t412 = Icges(4,2) * t429 + t494;
t413 = Icges(4,1) * t428 + t493;
t454 = -t412 * t428 + t413 * t429;
t433 = qJD(2) * t441;
t453 = -t443 * t474 + t433;
t452 = -qJD(5) * t426 + t420 * t392 + t481;
t406 = pkin(4) * t425 - qJ(5) * t426;
t451 = t421 * t406 + t443 * t476 + t453;
t448 = qJD(1) * t356 - t441 * t474 + t471;
t393 = t466 * t443;
t447 = qJD(1) * t393 + t441 * t476 + t448;
t424 = rSges(2,1) * t443 - rSges(2,2) * t441;
t423 = rSges(2,1) * t441 + rSges(2,2) * t443;
t419 = qJD(6) * t425 + qJD(1);
t414 = rSges(4,1) * t428 + rSges(4,2) * t429;
t411 = Icges(4,5) * t428 + Icges(4,6) * t429;
t410 = -pkin(5) * t443 + pkin(9) * t488;
t409 = pkin(5) * t441 + pkin(9) * t487;
t408 = rSges(5,1) * t425 + rSges(5,2) * t426;
t407 = -rSges(6,2) * t425 - rSges(6,3) * t426;
t399 = t425 * t486 - t483;
t398 = t425 * t484 + t485;
t397 = t425 * t485 + t484;
t396 = t425 * t483 - t486;
t395 = t441 * t475 + t421;
t394 = t443 * t475 + t420;
t391 = rSges(4,3) * t441 + t443 * t469;
t390 = -rSges(4,3) * t443 + t441 * t469;
t384 = Icges(4,3) * t441 + t443 * t460;
t383 = -Icges(4,3) * t443 + t441 * t460;
t381 = -rSges(6,1) * t443 + t441 * t467;
t380 = rSges(6,1) * t441 + t443 * t467;
t379 = rSges(5,3) * t441 + t443 * t468;
t378 = -rSges(5,3) * t443 + t441 * t468;
t363 = rSges(7,3) * t425 + (-rSges(7,1) * t440 - rSges(7,2) * t442) * t426;
t362 = Icges(7,5) * t425 + (-Icges(7,1) * t440 - Icges(7,4) * t442) * t426;
t361 = Icges(7,6) * t425 + (-Icges(7,4) * t440 - Icges(7,2) * t442) * t426;
t360 = Icges(7,3) * t425 + (-Icges(7,5) * t440 - Icges(7,6) * t442) * t426;
t358 = qJD(1) * t441 * rSges(3,3) + t416 + (qJD(1) * t470 - qJD(2)) * t443;
t357 = t433 + (t443 * rSges(3,3) - t441 * t470 - t422) * qJD(1);
t351 = rSges(7,1) * t399 + rSges(7,2) * t398 + rSges(7,3) * t488;
t350 = rSges(7,1) * t397 + rSges(7,2) * t396 + rSges(7,3) * t487;
t349 = Icges(7,1) * t399 + Icges(7,4) * t398 + Icges(7,5) * t488;
t348 = Icges(7,1) * t397 + Icges(7,4) * t396 + Icges(7,5) * t487;
t347 = Icges(7,4) * t399 + Icges(7,2) * t398 + Icges(7,6) * t488;
t346 = Icges(7,4) * t397 + Icges(7,2) * t396 + Icges(7,6) * t487;
t345 = Icges(7,5) * t399 + Icges(7,6) * t398 + Icges(7,3) * t488;
t344 = Icges(7,5) * t397 + Icges(7,6) * t396 + Icges(7,3) * t487;
t343 = (t390 * t441 + t391 * t443) * qJD(3);
t342 = qJD(1) * t391 - t414 * t432 + t471;
t341 = -t414 * t477 + t433 + (-t390 + t480) * qJD(1);
t340 = qJD(1) * t379 - t408 * t420 + t448;
t339 = t408 * t421 + (-t378 + t473) * qJD(1) + t453;
t338 = t378 * t420 - t379 * t421 + t481;
t337 = qJD(1) * t380 + (-t406 - t407) * t420 + t447;
t336 = t407 * t421 + (-t381 + t472) * qJD(1) + t451;
t335 = t381 * t420 + (-t380 - t393) * t421 + t452;
t334 = qJD(1) * t409 + t350 * t419 - t363 * t394 + (-t406 - t497) * t420 + t447;
t333 = t421 * t497 - t351 * t419 + t363 * t395 + (-t410 + t472) * qJD(1) + t451;
t332 = -t350 * t395 + t351 * t394 + t410 * t420 + (-t393 - t409) * t421 + t452;
t1 = t419 * ((t344 * t394 + t345 * t395 + t360 * t419) * t425 + ((-t346 * t442 - t348 * t440) * t394 + (-t347 * t442 - t349 * t440) * t395 + (-t361 * t442 - t362 * t440) * t419) * t426) / 0.2e1 + t394 * ((t344 * t487 + t396 * t346 + t397 * t348) * t394 + (t345 * t487 + t347 * t396 + t349 * t397) * t395 + (t360 * t487 + t361 * t396 + t362 * t397) * t419) / 0.2e1 + t395 * ((t344 * t488 + t346 * t398 + t348 * t399) * t394 + (t345 * t488 + t398 * t347 + t399 * t349) * t395 + (t360 * t488 + t361 * t398 + t362 * t399) * t419) / 0.2e1 - ((-t443 * t411 + t441 * t454) * qJD(1) + (t443 ^ 2 * t383 + (t455 * t441 + (-t384 + t456) * t443) * t441) * qJD(3)) * t477 / 0.2e1 + m(7) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(6) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(5) * (t338 ^ 2 + t339 ^ 2 + t340 ^ 2) / 0.2e1 + m(4) * (t341 ^ 2 + t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(3) * (t357 ^ 2 + t358 ^ 2) / 0.2e1 + ((t441 * t411 + t443 * t454) * qJD(1) + (t441 ^ 2 * t384 + (t456 * t443 + (-t383 + t455) * t441) * t443) * qJD(3)) * t432 / 0.2e1 + (t501 * t441 + t502 * t443) * t420 / 0.2e1 + (t502 * t441 - t501 * t443) * t421 / 0.2e1 + (m(2) * (t423 ^ 2 + t424 ^ 2) + Icges(3,2) * t438 ^ 2 + (Icges(3,1) * t437 + 0.2e1 * Icges(3,4) * t438) * t437 + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t386 * t429 + t388 * t428) * t441 - (t385 * t429 + t387 * t428) * t443) * qJD(3) + (t507 * t425 + t509 * t426) * t421 + (t506 * t425 + t508 * t426) * t420 + (t429 * t412 + t428 * t413 + t504 * t425 - t505 * t426) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
