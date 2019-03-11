% Calculate kinetic energy for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:56:11
% EndTime: 2019-03-09 02:56:13
% DurationCPUTime: 2.25s
% Computational Cost: add. (896->227), mult. (1268->335), div. (0->0), fcn. (1146->8), ass. (0->125)
t507 = -Icges(5,4) - Icges(6,6);
t506 = Icges(5,1) + Icges(6,2);
t505 = Icges(5,2) + Icges(6,3);
t418 = qJ(3) + pkin(9);
t414 = cos(t418);
t504 = t507 * t414;
t413 = sin(t418);
t503 = t507 * t413;
t502 = -Icges(6,4) + Icges(5,5);
t501 = Icges(6,5) - Icges(5,6);
t500 = -t505 * t414 + t503;
t499 = t506 * t413 - t504;
t422 = sin(qJ(1));
t425 = cos(qJ(1));
t498 = t500 * t422 + t501 * t425;
t497 = -t501 * t422 + t500 * t425;
t496 = t499 * t422 + t502 * t425;
t495 = t502 * t422 - t499 * t425;
t494 = t505 * t413 + t504;
t493 = t506 * t414 + t503;
t492 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t421 = sin(qJ(3));
t424 = cos(qJ(3));
t491 = Icges(4,5) * t421 + Icges(4,6) * t424 + t502 * t413 - t501 * t414;
t490 = t491 * t422 + t492 * t425;
t489 = t492 * t422 - t491 * t425;
t488 = Icges(4,5) * t424 - Icges(4,6) * t421 + t501 * t413 + t502 * t414;
t476 = Icges(4,4) * t424;
t402 = -Icges(4,2) * t421 + t476;
t477 = Icges(4,4) * t421;
t403 = Icges(4,1) * t424 - t477;
t487 = t402 * t424 + t403 * t421 + t493 * t413 - t494 * t414;
t444 = Icges(4,2) * t424 + t477;
t372 = Icges(4,6) * t425 + t444 * t422;
t446 = Icges(4,1) * t421 + t476;
t374 = Icges(4,5) * t425 + t446 * t422;
t486 = -t372 * t424 - t374 * t421 - t496 * t413 + t498 * t414;
t373 = Icges(4,6) * t422 - t444 * t425;
t375 = Icges(4,5) * t422 - t446 * t425;
t485 = t373 * t424 + t375 * t421 + t495 * t413 + t497 * t414;
t481 = pkin(3) * t421;
t480 = pkin(3) * t424;
t471 = t413 * t422;
t470 = t413 * t425;
t420 = sin(qJ(6));
t469 = t420 * t422;
t468 = t420 * t425;
t423 = cos(qJ(6));
t467 = t422 * t423;
t466 = t423 * t425;
t447 = pkin(4) * t413 - qJ(5) * t414;
t376 = t447 * t422;
t382 = qJ(4) * t425 + t422 * t481;
t465 = -t376 - t382;
t400 = qJD(1) * (pkin(1) * t425 + qJ(2) * t422);
t464 = qJD(1) * t425 * pkin(7) + t400;
t463 = qJD(3) * t422;
t462 = qJD(3) * t425;
t461 = qJD(5) * t414;
t460 = qJD(6) * t413;
t381 = qJ(4) * t422 - t425 * t481;
t367 = t381 * t462;
t377 = t447 * t425;
t459 = qJD(5) * t413 - t377 * t462 + t367;
t417 = qJD(2) * t422;
t458 = qJD(4) * t425 + t463 * t480 + t417;
t396 = pkin(4) * t414 + qJ(5) * t413;
t455 = -t396 - t480;
t405 = pkin(1) * t422 - qJ(2) * t425;
t454 = -pkin(7) * t422 - t405;
t453 = qJD(1) * t382 + qJD(4) * t422 + t464;
t452 = t396 * t463 + t458;
t451 = -t381 + t454;
t450 = rSges(4,1) * t421 + rSges(4,2) * t424;
t449 = rSges(5,1) * t413 + rSges(5,2) * t414;
t448 = rSges(6,2) * t413 + rSges(6,3) * t414;
t428 = t377 + t451;
t427 = qJD(1) * t376 + t425 * t461 + t453;
t409 = qJD(6) * t414 + qJD(1);
t408 = rSges(2,1) * t425 - rSges(2,2) * t422;
t407 = rSges(4,1) * t424 - rSges(4,2) * t421;
t406 = rSges(2,1) * t422 + rSges(2,2) * t425;
t399 = pkin(5) * t425 + pkin(8) * t471;
t398 = rSges(5,1) * t414 - rSges(5,2) * t413;
t397 = -rSges(6,2) * t414 + rSges(6,3) * t413;
t395 = pkin(5) * t422 - pkin(8) * t470;
t388 = t422 * t460 + t462;
t387 = -t425 * t460 + t463;
t386 = -t414 * t469 + t466;
t385 = -t414 * t467 - t468;
t384 = t414 * t468 + t467;
t383 = t414 * t466 - t469;
t379 = rSges(4,3) * t422 - t450 * t425;
t378 = rSges(4,3) * t425 + t450 * t422;
t366 = rSges(6,1) * t425 - t448 * t422;
t365 = rSges(6,1) * t422 + t448 * t425;
t364 = rSges(5,3) * t422 - t449 * t425;
t363 = rSges(5,3) * t425 + t449 * t422;
t349 = rSges(7,3) * t414 + (rSges(7,1) * t420 + rSges(7,2) * t423) * t413;
t348 = Icges(7,5) * t414 + (Icges(7,1) * t420 + Icges(7,4) * t423) * t413;
t347 = Icges(7,6) * t414 + (Icges(7,4) * t420 + Icges(7,2) * t423) * t413;
t346 = Icges(7,3) * t414 + (Icges(7,5) * t420 + Icges(7,6) * t423) * t413;
t345 = t400 - qJD(2) * t425 + qJD(1) * (-rSges(3,2) * t425 + rSges(3,3) * t422);
t344 = t417 + (rSges(3,2) * t422 + rSges(3,3) * t425 - t405) * qJD(1);
t343 = (-t378 * t422 + t379 * t425) * qJD(3);
t342 = rSges(7,1) * t386 + rSges(7,2) * t385 + rSges(7,3) * t471;
t341 = rSges(7,1) * t384 + rSges(7,2) * t383 - rSges(7,3) * t470;
t340 = Icges(7,1) * t386 + Icges(7,4) * t385 + Icges(7,5) * t471;
t339 = Icges(7,1) * t384 + Icges(7,4) * t383 - Icges(7,5) * t470;
t338 = Icges(7,4) * t386 + Icges(7,2) * t385 + Icges(7,6) * t471;
t337 = Icges(7,4) * t384 + Icges(7,2) * t383 - Icges(7,6) * t470;
t336 = Icges(7,5) * t386 + Icges(7,6) * t385 + Icges(7,3) * t471;
t335 = Icges(7,5) * t384 + Icges(7,6) * t383 - Icges(7,3) * t470;
t334 = qJD(1) * t378 + (-qJD(3) * t407 - qJD(2)) * t425 + t464;
t333 = t407 * t463 + t417 + (-t379 + t454) * qJD(1);
t332 = qJD(1) * t363 + (-qJD(2) + (-t398 - t480) * qJD(3)) * t425 + t453;
t331 = t398 * t463 + (-t364 + t451) * qJD(1) + t458;
t330 = t367 + (t364 * t425 + (-t363 - t382) * t422) * qJD(3);
t329 = qJD(1) * t366 + (-qJD(2) + (-t397 + t455) * qJD(3)) * t425 + t427;
t328 = (qJD(3) * t397 - t461) * t422 + (-t365 + t428) * qJD(1) + t452;
t327 = (t365 * t425 + (-t366 + t465) * t422) * qJD(3) + t459;
t326 = qJD(1) * t399 + t342 * t409 - t349 * t388 + (-qJD(2) + (-pkin(8) * t414 + t455) * qJD(3)) * t425 + t427;
t325 = -t341 * t409 + t349 * t387 + (pkin(8) * qJD(3) - qJD(5)) * t422 * t414 + (-t395 + t428) * qJD(1) + t452;
t324 = t341 * t388 - t342 * t387 + (t395 * t425 + (-t399 + t465) * t422) * qJD(3) + t459;
t1 = t388 * ((t336 * t471 + t385 * t338 + t386 * t340) * t388 + (t335 * t471 + t337 * t385 + t339 * t386) * t387 + (t346 * t471 + t347 * t385 + t348 * t386) * t409) / 0.2e1 + t387 * ((-t336 * t470 + t338 * t383 + t340 * t384) * t388 + (-t335 * t470 + t383 * t337 + t384 * t339) * t387 + (-t346 * t470 + t347 * t383 + t348 * t384) * t409) / 0.2e1 + t409 * ((t335 * t387 + t336 * t388 + t346 * t409) * t414 + ((t338 * t423 + t340 * t420) * t388 + (t337 * t423 + t339 * t420) * t387 + (t347 * t423 + t348 * t420) * t409) * t413) / 0.2e1 + m(5) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(6) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + m(7) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(4) * (t333 ^ 2 + t334 ^ 2 + t343 ^ 2) / 0.2e1 + m(3) * (t344 ^ 2 + t345 ^ 2) / 0.2e1 + (Icges(3,1) + Icges(2,3) + m(2) * (t406 ^ 2 + t408 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t372 * t421 + t374 * t424 + t498 * t413 + t496 * t414) * t425 + (-t373 * t421 + t375 * t424 - t497 * t413 + t495 * t414) * t422) * qJD(3) + (-t421 * t402 + t424 * t403 + t494 * t413 + t493 * t414) * qJD(1)) * qJD(1) / 0.2e1 + ((t489 * t422 ^ 2 + (t486 * t425 + (-t485 + t490) * t422) * t425) * qJD(3) + (t488 * t422 - t487 * t425) * qJD(1)) * t463 / 0.2e1 + ((t490 * t425 ^ 2 + (t485 * t422 + (-t486 + t489) * t425) * t422) * qJD(3) + (t487 * t422 + t488 * t425) * qJD(1)) * t462 / 0.2e1;
T  = t1;
