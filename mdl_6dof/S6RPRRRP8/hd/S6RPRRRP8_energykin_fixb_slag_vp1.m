% Calculate kinetic energy for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:49
% EndTime: 2019-03-09 06:22:51
% DurationCPUTime: 1.95s
% Computational Cost: add. (1121->220), mult. (1570->343), div. (0->0), fcn. (1514->8), ass. (0->125)
t498 = Icges(6,1) + Icges(7,1);
t497 = Icges(6,4) - Icges(7,5);
t496 = Icges(7,4) + Icges(6,5);
t495 = Icges(6,2) + Icges(7,3);
t494 = Icges(7,6) - Icges(6,6);
t493 = Icges(6,3) + Icges(7,2);
t492 = rSges(7,1) + pkin(5);
t491 = rSges(7,3) + qJ(6);
t425 = qJ(3) + qJ(4);
t423 = sin(t425);
t429 = cos(qJ(5));
t431 = cos(qJ(1));
t464 = t431 * t429;
t426 = sin(qJ(5));
t428 = sin(qJ(1));
t467 = t426 * t428;
t395 = t423 * t467 - t464;
t465 = t428 * t429;
t466 = t426 * t431;
t396 = t423 * t465 + t466;
t424 = cos(t425);
t469 = t424 * t428;
t490 = t495 * t395 - t497 * t396 - t494 * t469;
t397 = t423 * t466 + t465;
t398 = -t423 * t464 + t467;
t468 = t424 * t431;
t489 = -t495 * t397 - t497 * t398 + t494 * t468;
t488 = t494 * t395 + t496 * t396 - t493 * t469;
t487 = -t494 * t397 + t496 * t398 + t493 * t468;
t486 = -t497 * t395 + t498 * t396 - t496 * t469;
t485 = t497 * t397 + t498 * t398 + t496 * t468;
t484 = (t495 * t426 - t497 * t429) * t424 + t494 * t423;
t483 = (t494 * t426 + t496 * t429) * t424 + t493 * t423;
t482 = (-t497 * t426 + t498 * t429) * t424 + t496 * t423;
t471 = Icges(5,4) * t423;
t447 = Icges(5,2) * t424 + t471;
t371 = Icges(5,6) * t431 + t428 * t447;
t372 = Icges(5,6) * t428 - t431 * t447;
t470 = Icges(5,4) * t424;
t449 = Icges(5,1) * t423 + t470;
t373 = Icges(5,5) * t431 + t428 * t449;
t374 = Icges(5,5) * t428 - t431 * t449;
t400 = -Icges(5,2) * t423 + t470;
t401 = Icges(5,1) * t424 - t471;
t420 = qJD(3) * t428;
t405 = qJD(4) * t428 + t420;
t421 = qJD(3) * t431;
t406 = qJD(4) * t431 + t421;
t481 = (t371 * t424 + t373 * t423) * t406 + (t372 * t424 + t374 * t423) * t405 + (t400 * t424 + t401 * t423) * qJD(1);
t427 = sin(qJ(3));
t476 = pkin(3) * t427;
t473 = Icges(4,4) * t427;
t430 = cos(qJ(3));
t472 = Icges(4,4) * t430;
t463 = -rSges(7,2) * t469 + t491 * t395 + t492 * t396;
t462 = rSges(7,2) * t468 - t491 * t397 + t492 * t398;
t461 = rSges(7,2) * t423 + (t491 * t426 + t492 * t429) * t424;
t404 = qJD(1) * (pkin(1) * t431 + qJ(2) * t428);
t460 = qJD(1) * t431 * pkin(7) + t404;
t422 = qJD(2) * t428;
t457 = pkin(3) * qJD(3) * t430;
t459 = t428 * t457 + t422;
t458 = qJD(5) * t424;
t410 = pkin(1) * t428 - qJ(2) * t431;
t456 = -pkin(7) * t428 - t410;
t393 = pkin(8) * t428 - t431 * t476;
t455 = -t393 + t456;
t454 = pkin(4) * t423 - pkin(9) * t424;
t394 = pkin(8) * t431 + t428 * t476;
t453 = t393 * t421 - t394 * t420;
t452 = rSges(4,1) * t427 + rSges(4,2) * t430;
t451 = rSges(5,1) * t423 + rSges(5,2) * t424;
t450 = Icges(4,1) * t427 + t472;
t448 = Icges(4,2) * t430 + t473;
t446 = Icges(4,5) * t427 + Icges(4,6) * t430;
t445 = Icges(5,5) * t423 + Icges(5,6) * t424;
t381 = Icges(4,6) * t431 + t428 * t448;
t383 = Icges(4,5) * t431 + t428 * t450;
t442 = -t381 * t430 - t383 * t427;
t382 = Icges(4,6) * t428 - t431 * t448;
t384 = Icges(4,5) * t428 - t431 * t450;
t441 = t382 * t430 + t384 * t427;
t408 = -Icges(4,2) * t427 + t472;
t409 = Icges(4,1) * t430 - t473;
t439 = t408 * t430 + t409 * t427;
t438 = (Icges(5,5) * t424 - Icges(5,6) * t423) * qJD(1) + (Icges(5,3) * t431 + t428 * t445) * t406 + (Icges(5,3) * t428 - t431 * t445) * t405;
t391 = t454 * t428;
t392 = t454 * t431;
t437 = -t391 * t405 - t406 * t392 + t453;
t436 = qJD(1) * t394 + (-qJD(2) - t457) * t431 + t460;
t403 = pkin(4) * t424 + pkin(9) * t423;
t435 = t405 * t403 + (t392 + t455) * qJD(1) + t459;
t434 = qJD(1) * t391 - t403 * t406 + t436;
t414 = qJD(5) * t423 + qJD(1);
t413 = rSges(2,1) * t431 - rSges(2,2) * t428;
t412 = rSges(4,1) * t430 - rSges(4,2) * t427;
t411 = rSges(2,1) * t428 + rSges(2,2) * t431;
t407 = Icges(4,5) * t430 - Icges(4,6) * t427;
t402 = rSges(5,1) * t424 - rSges(5,2) * t423;
t389 = rSges(4,3) * t428 - t431 * t452;
t388 = rSges(4,3) * t431 + t428 * t452;
t387 = -t428 * t458 + t406;
t386 = t431 * t458 + t405;
t380 = Icges(4,3) * t428 - t431 * t446;
t379 = Icges(4,3) * t431 + t428 * t446;
t376 = rSges(5,3) * t428 - t431 * t451;
t375 = rSges(5,3) * t431 + t428 * t451;
t367 = rSges(6,3) * t423 + (rSges(6,1) * t429 - rSges(6,2) * t426) * t424;
t359 = t404 - qJD(2) * t431 + qJD(1) * (-rSges(3,2) * t431 + rSges(3,3) * t428);
t358 = t422 + (rSges(3,2) * t428 + rSges(3,3) * t431 - t410) * qJD(1);
t354 = rSges(6,1) * t398 + rSges(6,2) * t397 + rSges(6,3) * t468;
t352 = rSges(6,1) * t396 - rSges(6,2) * t395 - rSges(6,3) * t469;
t338 = (-t388 * t428 + t389 * t431) * qJD(3);
t337 = qJD(1) * t388 + (-qJD(3) * t412 - qJD(2)) * t431 + t460;
t336 = t412 * t420 + t422 + (-t389 + t456) * qJD(1);
t335 = qJD(1) * t375 - t402 * t406 + t436;
t334 = t402 * t405 + (-t376 + t455) * qJD(1) + t459;
t333 = -t375 * t405 + t376 * t406 + t453;
t332 = t352 * t414 - t367 * t387 + t434;
t331 = -t354 * t414 + t367 * t386 + t435;
t330 = -t352 * t386 + t354 * t387 + t437;
t329 = -qJD(6) * t397 - t387 * t461 + t414 * t463 + t434;
t328 = qJD(6) * t395 + t386 * t461 - t414 * t462 + t435;
t327 = qJD(6) * t424 * t426 - t386 * t463 + t387 * t462 + t437;
t1 = m(6) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(7) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + ((t428 * t407 - t431 * t439) * qJD(1) + (t428 ^ 2 * t380 + (t442 * t431 + (t379 - t441) * t428) * t431) * qJD(3)) * t420 / 0.2e1 + t406 * (t481 * t428 + t438 * t431) / 0.2e1 + t405 * (t438 * t428 - t481 * t431) / 0.2e1 + m(3) * (t358 ^ 2 + t359 ^ 2) / 0.2e1 + m(4) * (t336 ^ 2 + t337 ^ 2 + t338 ^ 2) / 0.2e1 + m(5) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + ((t431 * t407 + t428 * t439) * qJD(1) + (t431 ^ 2 * t379 + (t441 * t428 + (t380 - t442) * t431) * t428) * qJD(3)) * t421 / 0.2e1 + ((-t397 * t484 + t398 * t482 + t468 * t483) * t414 + (-t397 * t490 + t486 * t398 + t488 * t468) * t387 + (-t489 * t397 + t485 * t398 + t487 * t468) * t386) * t386 / 0.2e1 + ((t395 * t484 + t396 * t482 - t469 * t483) * t414 + (t490 * t395 + t486 * t396 - t488 * t469) * t387 + (t395 * t489 + t396 * t485 - t469 * t487) * t386) * t387 / 0.2e1 + (((t426 * t484 + t429 * t482) * t414 + (t426 * t490 + t486 * t429) * t387 + (t426 * t489 + t429 * t485) * t386) * t424 + (t386 * t487 + t387 * t488 + t414 * t483) * t423) * t414 / 0.2e1 + (((-t381 * t427 + t383 * t430) * t431 + (-t382 * t427 + t384 * t430) * t428) * qJD(3) + (-t371 * t423 + t373 * t424) * t406 + (-t372 * t423 + t374 * t424) * t405 + (-t400 * t423 + t401 * t424 - t408 * t427 + t409 * t430) * qJD(1)) * qJD(1) / 0.2e1 + (Icges(2,3) + Icges(3,1) + m(2) * (t411 ^ 2 + t413 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
