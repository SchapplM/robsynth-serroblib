% Calculate kinetic energy for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:53:14
% EndTime: 2019-03-09 02:53:16
% DurationCPUTime: 2.03s
% Computational Cost: add. (1207->276), mult. (1502->402), div. (0->0), fcn. (1436->10), ass. (0->142)
t492 = Icges(4,3) + Icges(5,3);
t419 = qJ(3) + pkin(9);
t412 = sin(t419);
t414 = cos(t419);
t424 = sin(qJ(3));
t426 = cos(qJ(3));
t491 = Icges(4,5) * t424 + Icges(5,5) * t412 + Icges(4,6) * t426 + Icges(5,6) * t414;
t425 = sin(qJ(1));
t427 = cos(qJ(1));
t490 = t491 * t425 + t492 * t427;
t489 = t492 * t425 - t491 * t427;
t488 = Icges(4,5) * t426 + Icges(5,5) * t414 - Icges(4,6) * t424 - Icges(5,6) * t412;
t472 = Icges(5,4) * t414;
t393 = -Icges(5,2) * t412 + t472;
t473 = Icges(5,4) * t412;
t394 = Icges(5,1) * t414 - t473;
t474 = Icges(4,4) * t426;
t399 = -Icges(4,2) * t424 + t474;
t475 = Icges(4,4) * t424;
t400 = Icges(4,1) * t426 - t475;
t487 = t393 * t414 + t394 * t412 + t399 * t426 + t400 * t424;
t439 = Icges(5,2) * t414 + t473;
t360 = Icges(5,6) * t425 - t427 * t439;
t441 = Icges(5,1) * t412 + t472;
t362 = Icges(5,5) * t425 - t427 * t441;
t440 = Icges(4,2) * t426 + t475;
t376 = Icges(4,6) * t425 - t427 * t440;
t442 = Icges(4,1) * t424 + t474;
t378 = Icges(4,5) * t425 - t427 * t442;
t486 = t360 * t414 + t362 * t412 + t376 * t426 + t378 * t424;
t359 = Icges(5,6) * t427 + t425 * t439;
t361 = Icges(5,5) * t427 + t425 * t441;
t375 = Icges(4,6) * t427 + t425 * t440;
t377 = Icges(4,5) * t427 + t425 * t442;
t485 = -t359 * t414 - t361 * t412 - t375 * t426 - t377 * t424;
t421 = cos(pkin(10));
t477 = pkin(5) * t421;
t484 = -pkin(8) * t414 + t412 * t477;
t480 = pkin(3) * t424;
t479 = pkin(3) * t426;
t418 = pkin(10) + qJ(6);
t411 = sin(t418);
t471 = t411 * t425;
t470 = t411 * t427;
t413 = cos(t418);
t469 = t413 * t425;
t468 = t413 * t427;
t467 = t414 * t425;
t466 = t414 * t427;
t420 = sin(pkin(10));
t465 = t420 * t425;
t464 = t420 * t427;
t463 = t421 * t425;
t462 = t421 * t427;
t443 = pkin(4) * t412 - qJ(5) * t414;
t379 = t443 * t425;
t389 = qJ(4) * t427 + t425 * t480;
t460 = -t379 - t389;
t397 = qJD(1) * (pkin(1) * t427 + qJ(2) * t425);
t459 = qJD(1) * t427 * pkin(7) + t397;
t458 = qJD(3) * t425;
t457 = qJD(3) * t427;
t456 = qJD(5) * t414;
t455 = qJD(6) * t414;
t388 = qJ(4) * t425 - t427 * t480;
t366 = t388 * t457;
t380 = t443 * t427;
t454 = qJD(5) * t412 - t380 * t457 + t366;
t417 = qJD(2) * t425;
t453 = qJD(4) * t427 + t458 * t479 + t417;
t395 = pkin(4) * t414 + qJ(5) * t412;
t450 = -t395 - t479;
t402 = pkin(1) * t425 - qJ(2) * t427;
t449 = -pkin(7) * t425 - t402;
t448 = qJD(1) * t389 + qJD(4) * t425 + t459;
t447 = t395 * t458 + t453;
t446 = -t388 + t449;
t445 = rSges(4,1) * t424 + rSges(4,2) * t426;
t444 = rSges(5,1) * t412 + rSges(5,2) * t414;
t430 = t380 + t446;
t429 = qJD(1) * t379 + t427 * t456 + t448;
t406 = qJD(6) * t412 + qJD(1);
t405 = rSges(2,1) * t427 - rSges(2,2) * t425;
t404 = rSges(4,1) * t426 - rSges(4,2) * t424;
t403 = rSges(2,1) * t425 + rSges(2,2) * t427;
t396 = rSges(5,1) * t414 - rSges(5,2) * t412;
t391 = -t425 * t455 + t457;
t390 = t427 * t455 + t458;
t387 = -t412 * t462 + t465;
t386 = t412 * t464 + t463;
t385 = t412 * t463 + t464;
t384 = -t412 * t465 + t462;
t382 = t425 * rSges(4,3) - t427 * t445;
t381 = t427 * rSges(4,3) + t425 * t445;
t370 = -t412 * t468 + t471;
t369 = t412 * t470 + t469;
t368 = t412 * t469 + t470;
t367 = -t412 * t471 + t468;
t365 = t425 * rSges(5,3) - t427 * t444;
t364 = t427 * rSges(5,3) + t425 * t444;
t356 = t412 * rSges(6,3) + (rSges(6,1) * t421 - rSges(6,2) * t420) * t414;
t355 = Icges(6,5) * t412 + (Icges(6,1) * t421 - Icges(6,4) * t420) * t414;
t354 = Icges(6,6) * t412 + (Icges(6,4) * t421 - Icges(6,2) * t420) * t414;
t353 = Icges(6,3) * t412 + (Icges(6,5) * t421 - Icges(6,6) * t420) * t414;
t352 = t397 - qJD(2) * t427 + qJD(1) * (-rSges(3,2) * t427 + rSges(3,3) * t425);
t351 = t417 + (rSges(3,2) * t425 + rSges(3,3) * t427 - t402) * qJD(1);
t350 = t412 * rSges(7,3) + (rSges(7,1) * t413 - rSges(7,2) * t411) * t414;
t349 = Icges(7,5) * t412 + (Icges(7,1) * t413 - Icges(7,4) * t411) * t414;
t348 = Icges(7,6) * t412 + (Icges(7,4) * t413 - Icges(7,2) * t411) * t414;
t347 = Icges(7,3) * t412 + (Icges(7,5) * t413 - Icges(7,6) * t411) * t414;
t346 = pkin(8) * t412 + t414 * t477;
t345 = (-t381 * t425 + t382 * t427) * qJD(3);
t344 = rSges(6,1) * t387 + rSges(6,2) * t386 + rSges(6,3) * t466;
t343 = rSges(6,1) * t385 + rSges(6,2) * t384 - rSges(6,3) * t467;
t342 = Icges(6,1) * t387 + Icges(6,4) * t386 + Icges(6,5) * t466;
t341 = Icges(6,1) * t385 + Icges(6,4) * t384 - Icges(6,5) * t467;
t340 = Icges(6,4) * t387 + Icges(6,2) * t386 + Icges(6,6) * t466;
t339 = Icges(6,4) * t385 + Icges(6,2) * t384 - Icges(6,6) * t467;
t338 = Icges(6,5) * t387 + Icges(6,6) * t386 + Icges(6,3) * t466;
t337 = Icges(6,5) * t385 + Icges(6,6) * t384 - Icges(6,3) * t467;
t336 = pkin(5) * t465 - t484 * t427;
t335 = pkin(5) * t464 + t484 * t425;
t334 = rSges(7,1) * t370 + rSges(7,2) * t369 + rSges(7,3) * t466;
t333 = rSges(7,1) * t368 + rSges(7,2) * t367 - rSges(7,3) * t467;
t332 = qJD(1) * t381 + (-qJD(3) * t404 - qJD(2)) * t427 + t459;
t331 = t404 * t458 + t417 + (-t382 + t449) * qJD(1);
t330 = Icges(7,1) * t370 + Icges(7,4) * t369 + Icges(7,5) * t466;
t329 = Icges(7,1) * t368 + Icges(7,4) * t367 - Icges(7,5) * t467;
t328 = Icges(7,4) * t370 + Icges(7,2) * t369 + Icges(7,6) * t466;
t327 = Icges(7,4) * t368 + Icges(7,2) * t367 - Icges(7,6) * t467;
t326 = Icges(7,5) * t370 + Icges(7,6) * t369 + Icges(7,3) * t466;
t325 = Icges(7,5) * t368 + Icges(7,6) * t367 - Icges(7,3) * t467;
t324 = qJD(1) * t364 + (-qJD(2) + (-t396 - t479) * qJD(3)) * t427 + t448;
t323 = t396 * t458 + (-t365 + t446) * qJD(1) + t453;
t322 = t366 + (t365 * t427 + (-t364 - t389) * t425) * qJD(3);
t321 = qJD(1) * t343 + (-qJD(2) + (-t356 + t450) * qJD(3)) * t427 + t429;
t320 = (qJD(3) * t356 - t456) * t425 + (-t344 + t430) * qJD(1) + t447;
t319 = (t344 * t427 + (-t343 + t460) * t425) * qJD(3) + t454;
t318 = qJD(1) * t335 + t333 * t406 - t350 * t391 + (-qJD(2) + (-t346 + t450) * qJD(3)) * t427 + t429;
t317 = -t334 * t406 + t350 * t390 + (qJD(3) * t346 - t456) * t425 + (-t336 + t430) * qJD(1) + t447;
t316 = -t333 * t390 + t334 * t391 + (t336 * t427 + (-t335 + t460) * t425) * qJD(3) + t454;
t1 = m(7) * (t316 ^ 2 + t317 ^ 2 + t318 ^ 2) / 0.2e1 + m(5) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + m(6) * (t319 ^ 2 + t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(4) * (t331 ^ 2 + t332 ^ 2 + t345 ^ 2) / 0.2e1 + m(3) * (t351 ^ 2 + t352 ^ 2) / 0.2e1 + t406 * ((t325 * t391 + t326 * t390 + t347 * t406) * t412 + ((-t327 * t411 + t329 * t413) * t391 + (-t328 * t411 + t330 * t413) * t390 + (-t348 * t411 + t349 * t413) * t406) * t414) / 0.2e1 + t391 * ((-t325 * t467 + t367 * t327 + t368 * t329) * t391 + (-t326 * t467 + t328 * t367 + t330 * t368) * t390 + (-t347 * t467 + t348 * t367 + t349 * t368) * t406) / 0.2e1 + t390 * ((t325 * t466 + t327 * t369 + t329 * t370) * t391 + (t326 * t466 + t328 * t369 + t330 * t370) * t390 + (t347 * t466 + t348 * t369 + t349 * t370) * t406) / 0.2e1 + (m(2) * (t403 ^ 2 + t405 ^ 2) + Icges(3,1) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t424 * t375 + t426 * t377 + (-t339 * t420 + t341 * t421 + t361) * t414 + (-t359 + t337) * t412) * t427 + (-t424 * t376 + t426 * t378 + (-t340 * t420 + t342 * t421 + t362) * t414 + (-t360 + t338) * t412) * t425) * qJD(3) + (-t424 * t399 + t426 * t400 + (-t354 * t420 + t355 * t421 + t394) * t414 + (-t393 + t353) * t412) * qJD(1)) * qJD(1) / 0.2e1 + (((t337 * t466 + t339 * t386 + t341 * t387 + t485 * t427) * t427 + (t338 * t466 + t386 * t340 + t387 * t342 + (-t486 + t490) * t427 + t489 * t425) * t425) * qJD(3) + (t353 * t466 + t386 * t354 + t387 * t355 + t488 * t425 - t487 * t427) * qJD(1)) * t458 / 0.2e1 + (((-t338 * t467 + t340 * t384 + t342 * t385 + t486 * t425) * t425 + (-t337 * t467 + t384 * t339 + t385 * t341 + (-t485 + t489) * t425 + t490 * t427) * t427) * qJD(3) + (-t353 * t467 + t384 * t354 + t385 * t355 + t487 * t425 + t488 * t427) * qJD(1)) * t457 / 0.2e1;
T  = t1;
