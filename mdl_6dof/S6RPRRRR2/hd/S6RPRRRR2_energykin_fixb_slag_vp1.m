% Calculate kinetic energy for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:56:59
% EndTime: 2019-03-09 06:57:01
% DurationCPUTime: 1.79s
% Computational Cost: add. (2038->271), mult. (1605->441), div. (0->0), fcn. (1538->12), ass. (0->146)
t423 = sin(qJ(1));
t475 = pkin(1) * t423;
t425 = cos(qJ(3));
t473 = pkin(3) * t425;
t424 = cos(qJ(5));
t472 = pkin(5) * t424;
t422 = sin(qJ(3));
t469 = Icges(4,4) * t422;
t468 = Icges(4,4) * t425;
t420 = qJ(3) + qJ(4);
t415 = sin(t420);
t467 = Icges(5,4) * t415;
t417 = cos(t420);
t466 = Icges(5,4) * t417;
t418 = qJ(1) + pkin(11);
t412 = sin(t418);
t465 = t412 * t415;
t421 = sin(qJ(5));
t464 = t412 * t421;
t413 = cos(t418);
t463 = t413 * t415;
t462 = t413 * t421;
t419 = qJ(5) + qJ(6);
t414 = sin(t419);
t461 = t414 * t417;
t416 = cos(t419);
t460 = t416 * t417;
t459 = t417 * t421;
t458 = t417 * t424;
t426 = cos(qJ(1));
t411 = qJD(1) * t426 * pkin(1);
t457 = qJD(1) * (pkin(2) * t413 + pkin(7) * t412) + t411;
t408 = qJD(3) * t412;
t390 = qJD(4) * t412 + t408;
t456 = qJD(3) * t413;
t455 = qJD(5) * t415;
t454 = qJD(6) * t415;
t453 = pkin(3) * qJD(3) * t422;
t367 = t413 * t455 + t390;
t345 = -pkin(8) * t413 + t412 * t473;
t346 = pkin(8) * t412 + t413 * t473;
t452 = t345 * t408 + t346 * t456 + qJD(2);
t451 = -pkin(2) * t412 + pkin(7) * t413 - t475;
t391 = (-qJD(3) - qJD(4)) * t413;
t450 = t413 * t453;
t449 = -t345 + t451;
t448 = pkin(4) * t417 + pkin(9) * t415;
t447 = rSges(4,1) * t425 - rSges(4,2) * t422;
t446 = rSges(5,1) * t417 - rSges(5,2) * t415;
t445 = Icges(4,1) * t425 - t469;
t444 = Icges(5,1) * t417 - t467;
t443 = -Icges(4,2) * t422 + t468;
t442 = -Icges(5,2) * t415 + t466;
t441 = Icges(4,5) * t425 - Icges(4,6) * t422;
t440 = Icges(5,5) * t417 - Icges(5,6) * t415;
t363 = -Icges(4,6) * t413 + t412 * t443;
t365 = -Icges(4,5) * t413 + t412 * t445;
t439 = t363 * t422 - t365 * t425;
t364 = Icges(4,6) * t412 + t413 * t443;
t366 = Icges(4,5) * t412 + t413 * t445;
t438 = -t364 * t422 + t366 * t425;
t401 = Icges(4,2) * t425 + t469;
t402 = Icges(4,1) * t422 + t468;
t437 = -t401 * t422 + t402 * t425;
t368 = t412 * t455 + t391;
t380 = t448 * t412;
t381 = t448 * t413;
t436 = t390 * t380 - t381 * t391 + t452;
t435 = qJD(1) * t346 - t412 * t453 + t457;
t434 = pkin(10) * t415 + t417 * t472;
t433 = qJD(1) * (Icges(5,5) * t415 + Icges(5,6) * t417) + (-Icges(5,3) * t413 + t412 * t440) * t391 + (Icges(5,3) * t412 + t413 * t440) * t390;
t397 = pkin(4) * t415 - pkin(9) * t417;
t432 = qJD(1) * t381 - t390 * t397 + t435;
t431 = t391 * t397 + (-t380 + t449) * qJD(1) - t450;
t350 = -Icges(5,6) * t413 + t412 * t442;
t351 = Icges(5,6) * t412 + t413 * t442;
t352 = -Icges(5,5) * t413 + t412 * t444;
t353 = Icges(5,5) * t412 + t413 * t444;
t394 = Icges(5,2) * t417 + t467;
t395 = Icges(5,1) * t415 + t466;
t430 = (-t351 * t415 + t353 * t417) * t390 + (-t350 * t415 + t352 * t417) * t391 + (-t394 * t415 + t395 * t417) * qJD(1);
t406 = -qJD(5) * t417 + qJD(1);
t405 = rSges(2,1) * t426 - rSges(2,2) * t423;
t404 = rSges(2,1) * t423 + rSges(2,2) * t426;
t403 = rSges(4,1) * t422 + rSges(4,2) * t425;
t400 = Icges(4,5) * t422 + Icges(4,6) * t425;
t396 = rSges(5,1) * t415 + rSges(5,2) * t417;
t389 = qJD(1) + (-qJD(5) - qJD(6)) * t417;
t387 = t413 * t458 + t464;
t386 = t412 * t424 - t413 * t459;
t385 = t412 * t458 - t462;
t384 = -t412 * t459 - t413 * t424;
t383 = t411 + qJD(1) * (rSges(3,1) * t413 - rSges(3,2) * t412);
t382 = (-rSges(3,1) * t412 - rSges(3,2) * t413 - t475) * qJD(1);
t378 = -rSges(6,3) * t417 + (rSges(6,1) * t424 - rSges(6,2) * t421) * t415;
t377 = -Icges(6,5) * t417 + (Icges(6,1) * t424 - Icges(6,4) * t421) * t415;
t376 = -Icges(6,6) * t417 + (Icges(6,4) * t424 - Icges(6,2) * t421) * t415;
t375 = -Icges(6,3) * t417 + (Icges(6,5) * t424 - Icges(6,6) * t421) * t415;
t374 = t412 * t414 + t413 * t460;
t373 = t412 * t416 - t413 * t461;
t372 = t412 * t460 - t413 * t414;
t371 = -t412 * t461 - t413 * t416;
t370 = rSges(4,3) * t412 + t413 * t447;
t369 = -rSges(4,3) * t413 + t412 * t447;
t362 = Icges(4,3) * t412 + t413 * t441;
t361 = -Icges(4,3) * t413 + t412 * t441;
t359 = -rSges(7,3) * t417 + (rSges(7,1) * t416 - rSges(7,2) * t414) * t415;
t358 = -Icges(7,5) * t417 + (Icges(7,1) * t416 - Icges(7,4) * t414) * t415;
t357 = -Icges(7,6) * t417 + (Icges(7,4) * t416 - Icges(7,2) * t414) * t415;
t356 = -Icges(7,3) * t417 + (Icges(7,5) * t416 - Icges(7,6) * t414) * t415;
t355 = rSges(5,3) * t412 + t413 * t446;
t354 = -rSges(5,3) * t413 + t412 * t446;
t347 = -pkin(10) * t417 + t415 * t472;
t343 = t412 * t454 + t368;
t342 = t413 * t454 + t367;
t338 = rSges(6,1) * t387 + rSges(6,2) * t386 + rSges(6,3) * t463;
t337 = rSges(6,1) * t385 + rSges(6,2) * t384 + rSges(6,3) * t465;
t336 = Icges(6,1) * t387 + Icges(6,4) * t386 + Icges(6,5) * t463;
t335 = Icges(6,1) * t385 + Icges(6,4) * t384 + Icges(6,5) * t465;
t334 = Icges(6,4) * t387 + Icges(6,2) * t386 + Icges(6,6) * t463;
t333 = Icges(6,4) * t385 + Icges(6,2) * t384 + Icges(6,6) * t465;
t332 = Icges(6,5) * t387 + Icges(6,6) * t386 + Icges(6,3) * t463;
t331 = Icges(6,5) * t385 + Icges(6,6) * t384 + Icges(6,3) * t465;
t330 = pkin(5) * t464 + t413 * t434;
t329 = -pkin(5) * t462 + t412 * t434;
t328 = rSges(7,1) * t374 + rSges(7,2) * t373 + rSges(7,3) * t463;
t327 = rSges(7,1) * t372 + rSges(7,2) * t371 + rSges(7,3) * t465;
t326 = qJD(1) * t370 - t403 * t408 + t457;
t325 = -t403 * t456 + (-t369 + t451) * qJD(1);
t324 = Icges(7,1) * t374 + Icges(7,4) * t373 + Icges(7,5) * t463;
t323 = Icges(7,1) * t372 + Icges(7,4) * t371 + Icges(7,5) * t465;
t322 = Icges(7,4) * t374 + Icges(7,2) * t373 + Icges(7,6) * t463;
t321 = Icges(7,4) * t372 + Icges(7,2) * t371 + Icges(7,6) * t465;
t320 = Icges(7,5) * t374 + Icges(7,6) * t373 + Icges(7,3) * t463;
t319 = Icges(7,5) * t372 + Icges(7,6) * t371 + Icges(7,3) * t465;
t318 = qJD(2) + (t369 * t412 + t370 * t413) * qJD(3);
t317 = qJD(1) * t355 - t390 * t396 + t435;
t316 = -t450 + t391 * t396 + (-t354 + t449) * qJD(1);
t315 = t354 * t390 - t355 * t391 + t452;
t314 = t338 * t406 - t367 * t378 + t432;
t313 = -t337 * t406 + t368 * t378 + t431;
t312 = t337 * t367 - t338 * t368 + t436;
t311 = t328 * t389 + t330 * t406 - t342 * t359 - t347 * t367 + t432;
t310 = -t327 * t389 - t329 * t406 + t343 * t359 + t347 * t368 + t431;
t309 = t327 * t342 - t328 * t343 + t329 * t367 - t330 * t368 + t436;
t1 = ((t412 * t400 + t413 * t437) * qJD(1) + (t412 ^ 2 * t362 + (t439 * t413 + (-t361 + t438) * t412) * t413) * qJD(3)) * t408 / 0.2e1 + m(7) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + m(6) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(5) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(4) * (t318 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + t342 * ((t320 * t463 + t373 * t322 + t374 * t324) * t342 + (t319 * t463 + t321 * t373 + t323 * t374) * t343 + (t356 * t463 + t357 * t373 + t358 * t374) * t389) / 0.2e1 + t343 * ((t320 * t465 + t322 * t371 + t324 * t372) * t342 + (t319 * t465 + t371 * t321 + t372 * t323) * t343 + (t356 * t465 + t357 * t371 + t358 * t372) * t389) / 0.2e1 + t389 * ((-t319 * t343 - t320 * t342 - t356 * t389) * t417 + ((-t322 * t414 + t324 * t416) * t342 + (-t321 * t414 + t323 * t416) * t343 + (-t357 * t414 + t358 * t416) * t389) * t415) / 0.2e1 + t367 * ((t332 * t463 + t334 * t386 + t336 * t387) * t367 + (t331 * t463 + t333 * t386 + t335 * t387) * t368 + (t375 * t463 + t376 * t386 + t377 * t387) * t406) / 0.2e1 + t368 * ((t332 * t465 + t334 * t384 + t336 * t385) * t367 + (t331 * t465 + t333 * t384 + t335 * t385) * t368 + (t375 * t465 + t376 * t384 + t377 * t385) * t406) / 0.2e1 + t406 * ((-t331 * t368 - t332 * t367 - t375 * t406) * t417 + ((-t334 * t421 + t336 * t424) * t367 + (-t333 * t421 + t335 * t424) * t368 + (-t376 * t421 + t377 * t424) * t406) * t415) / 0.2e1 + t390 * (t433 * t412 + t430 * t413) / 0.2e1 + t391 * (t430 * t412 - t433 * t413) / 0.2e1 - ((-t413 * t400 + t412 * t437) * qJD(1) + (t413 ^ 2 * t361 + (t438 * t412 + (-t362 + t439) * t413) * t412) * qJD(3)) * t456 / 0.2e1 + ((t351 * t417 + t353 * t415) * t390 + (t350 * t417 + t352 * t415) * t391 + ((t364 * t425 + t366 * t422) * t412 - (t363 * t425 + t365 * t422) * t413) * qJD(3) + (t417 * t394 + t415 * t395 + t425 * t401 + t422 * t402) * qJD(1)) * qJD(1) / 0.2e1 + (Icges(2,3) + Icges(3,3) + m(2) * (t404 ^ 2 + t405 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
