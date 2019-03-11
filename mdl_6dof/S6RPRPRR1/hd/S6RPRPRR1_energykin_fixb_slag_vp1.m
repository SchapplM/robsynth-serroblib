% Calculate kinetic energy for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:33:50
% EndTime: 2019-03-09 03:33:52
% DurationCPUTime: 1.74s
% Computational Cost: add. (1819->242), mult. (1313->378), div. (0->0), fcn. (1190->12), ass. (0->139)
t484 = Icges(4,3) + Icges(5,3);
t409 = qJ(3) + pkin(11);
t402 = sin(t409);
t404 = cos(t409);
t413 = sin(qJ(3));
t416 = cos(qJ(3));
t483 = Icges(4,5) * t416 + Icges(5,5) * t404 - Icges(4,6) * t413 - Icges(5,6) * t402;
t410 = qJ(1) + pkin(10);
t403 = sin(t410);
t405 = cos(t410);
t482 = t483 * t403 - t405 * t484;
t481 = t403 * t484 + t483 * t405;
t480 = Icges(4,5) * t413 + Icges(5,5) * t402 + Icges(4,6) * t416 + Icges(5,6) * t404;
t466 = Icges(5,4) * t402;
t383 = Icges(5,2) * t404 + t466;
t465 = Icges(5,4) * t404;
t384 = Icges(5,1) * t402 + t465;
t468 = Icges(4,4) * t413;
t390 = Icges(4,2) * t416 + t468;
t467 = Icges(4,4) * t416;
t391 = Icges(4,1) * t413 + t467;
t479 = -t383 * t402 + t384 * t404 - t390 * t413 + t391 * t416;
t436 = -Icges(5,2) * t402 + t465;
t347 = Icges(5,6) * t403 + t405 * t436;
t439 = Icges(5,1) * t404 - t466;
t349 = Icges(5,5) * t403 + t405 * t439;
t437 = -Icges(4,2) * t413 + t467;
t361 = Icges(4,6) * t403 + t405 * t437;
t440 = Icges(4,1) * t416 - t468;
t363 = Icges(4,5) * t403 + t405 * t440;
t478 = -t347 * t402 + t349 * t404 - t361 * t413 + t363 * t416;
t346 = -Icges(5,6) * t405 + t403 * t436;
t348 = -Icges(5,5) * t405 + t403 * t439;
t360 = -Icges(4,6) * t405 + t403 * t437;
t362 = -Icges(4,5) * t405 + t403 * t440;
t477 = t346 * t402 - t348 * t404 + t360 * t413 - t362 * t416;
t414 = sin(qJ(1));
t473 = pkin(1) * t414;
t472 = pkin(3) * t413;
t470 = t416 * pkin(3);
t406 = qJ(5) + t409;
t398 = sin(t406);
t464 = Icges(6,4) * t398;
t399 = cos(t406);
t463 = Icges(6,4) * t399;
t462 = t398 * t403;
t461 = t398 * t405;
t412 = sin(qJ(6));
t460 = t403 * t412;
t415 = cos(qJ(6));
t459 = t403 * t415;
t458 = t405 * t412;
t457 = t405 * t415;
t417 = cos(qJ(1));
t401 = qJD(1) * t417 * pkin(1);
t456 = qJD(1) * (pkin(2) * t405 + pkin(7) * t403) + t401;
t455 = pkin(4) * t404;
t397 = qJD(3) * t403;
t380 = qJD(5) * t403 + t397;
t453 = qJD(3) * t405;
t452 = qJD(6) * t398;
t334 = -qJ(4) * t405 + t403 * t470;
t335 = qJ(4) * t403 + t405 * t470;
t451 = t334 * t397 + t335 * t453 + qJD(2);
t448 = -pkin(2) * t403 + pkin(7) * t405 - t473;
t381 = (-qJD(3) - qJD(5)) * t405;
t447 = -t334 + t448;
t446 = pkin(5) * t399 + pkin(9) * t398;
t445 = rSges(4,1) * t416 - rSges(4,2) * t413;
t444 = rSges(5,1) * t404 - rSges(5,2) * t402;
t443 = rSges(6,1) * t399 - rSges(6,2) * t398;
t442 = qJD(3) * (-rSges(5,1) * t402 - rSges(5,2) * t404 - t472);
t329 = -pkin(8) * t405 + t403 * t455;
t330 = pkin(8) * t403 + t405 * t455;
t441 = t329 * t397 + t330 * t453 + t451;
t438 = Icges(6,1) * t399 - t464;
t435 = -Icges(6,2) * t398 + t463;
t432 = Icges(6,5) * t399 - Icges(6,6) * t398;
t425 = -t329 + t447;
t424 = qJD(1) * t335 - qJD(4) * t405 + t456;
t423 = qJD(3) * (-pkin(4) * t402 - t472);
t422 = qJD(1) * (Icges(6,5) * t398 + Icges(6,6) * t399) + (-Icges(6,3) * t405 + t403 * t432) * t381 + (Icges(6,3) * t403 + t405 * t432) * t380;
t396 = qJD(4) * t403;
t421 = t405 * t423 + t396;
t420 = qJD(1) * t330 + t403 * t423 + t424;
t338 = -Icges(6,6) * t405 + t403 * t435;
t339 = Icges(6,6) * t403 + t405 * t435;
t340 = -Icges(6,5) * t405 + t403 * t438;
t341 = Icges(6,5) * t403 + t405 * t438;
t375 = Icges(6,2) * t399 + t464;
t376 = Icges(6,1) * t398 + t463;
t419 = (-t339 * t398 + t341 * t399) * t380 + (-t338 * t398 + t340 * t399) * t381 + (-t375 * t398 + t376 * t399) * qJD(1);
t394 = rSges(2,1) * t417 - rSges(2,2) * t414;
t393 = rSges(2,1) * t414 + rSges(2,2) * t417;
t392 = rSges(4,1) * t413 + rSges(4,2) * t416;
t388 = -qJD(6) * t399 + qJD(1);
t378 = pkin(5) * t398 - pkin(9) * t399;
t377 = rSges(6,1) * t398 + rSges(6,2) * t399;
t373 = t401 + qJD(1) * (rSges(3,1) * t405 - rSges(3,2) * t403);
t372 = (-rSges(3,1) * t403 - rSges(3,2) * t405 - t473) * qJD(1);
t371 = t399 * t457 + t460;
t370 = -t399 * t458 + t459;
t369 = t399 * t459 - t458;
t368 = -t399 * t460 - t457;
t367 = t446 * t405;
t366 = t446 * t403;
t365 = rSges(4,3) * t403 + t405 * t445;
t364 = -rSges(4,3) * t405 + t403 * t445;
t357 = t403 * t452 + t381;
t356 = t405 * t452 + t380;
t355 = -rSges(7,3) * t399 + (rSges(7,1) * t415 - rSges(7,2) * t412) * t398;
t354 = -Icges(7,5) * t399 + (Icges(7,1) * t415 - Icges(7,4) * t412) * t398;
t353 = -Icges(7,6) * t399 + (Icges(7,4) * t415 - Icges(7,2) * t412) * t398;
t352 = -Icges(7,3) * t399 + (Icges(7,5) * t415 - Icges(7,6) * t412) * t398;
t351 = rSges(5,3) * t403 + t405 * t444;
t350 = -rSges(5,3) * t405 + t403 * t444;
t343 = rSges(6,3) * t403 + t405 * t443;
t342 = -rSges(6,3) * t405 + t403 * t443;
t325 = rSges(7,1) * t371 + rSges(7,2) * t370 + rSges(7,3) * t461;
t324 = rSges(7,1) * t369 + rSges(7,2) * t368 + rSges(7,3) * t462;
t323 = Icges(7,1) * t371 + Icges(7,4) * t370 + Icges(7,5) * t461;
t322 = Icges(7,1) * t369 + Icges(7,4) * t368 + Icges(7,5) * t462;
t321 = Icges(7,4) * t371 + Icges(7,2) * t370 + Icges(7,6) * t461;
t320 = Icges(7,4) * t369 + Icges(7,2) * t368 + Icges(7,6) * t462;
t319 = Icges(7,5) * t371 + Icges(7,6) * t370 + Icges(7,3) * t461;
t318 = Icges(7,5) * t369 + Icges(7,6) * t368 + Icges(7,3) * t462;
t317 = qJD(1) * t365 - t392 * t397 + t456;
t316 = -t392 * t453 + (-t364 + t448) * qJD(1);
t315 = qJD(2) + (t364 * t403 + t365 * t405) * qJD(3);
t314 = qJD(1) * t351 + t403 * t442 + t424;
t313 = t396 + t405 * t442 + (-t350 + t447) * qJD(1);
t312 = (t350 * t403 + t351 * t405) * qJD(3) + t451;
t311 = qJD(1) * t343 - t377 * t380 + t420;
t310 = t377 * t381 + (-t342 + t425) * qJD(1) + t421;
t309 = t342 * t380 - t343 * t381 + t441;
t308 = qJD(1) * t367 + t325 * t388 - t355 * t356 - t378 * t380 + t420;
t307 = -t324 * t388 + t355 * t357 + t378 * t381 + (-t366 + t425) * qJD(1) + t421;
t306 = t324 * t356 - t325 * t357 + t366 * t380 - t367 * t381 + t441;
t1 = m(5) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(6) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + m(7) * (t306 ^ 2 + t307 ^ 2 + t308 ^ 2) / 0.2e1 + t380 * (t422 * t403 + t419 * t405) / 0.2e1 + t381 * (t419 * t403 - t422 * t405) / 0.2e1 + t388 * ((-t318 * t357 - t319 * t356 - t352 * t388) * t399 + ((-t321 * t412 + t323 * t415) * t356 + (-t320 * t412 + t322 * t415) * t357 + (-t353 * t412 + t354 * t415) * t388) * t398) / 0.2e1 + t356 * ((t319 * t461 + t370 * t321 + t371 * t323) * t356 + (t318 * t461 + t320 * t370 + t322 * t371) * t357 + (t352 * t461 + t353 * t370 + t354 * t371) * t388) / 0.2e1 + t357 * ((t319 * t462 + t321 * t368 + t323 * t369) * t356 + (t318 * t462 + t368 * t320 + t369 * t322) * t357 + (t352 * t462 + t353 * t368 + t354 * t369) * t388) / 0.2e1 + m(4) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + ((t481 * t403 ^ 2 + (t477 * t405 + (t478 - t482) * t403) * t405) * qJD(3) + (t403 * t480 + t405 * t479) * qJD(1)) * t397 / 0.2e1 - ((t482 * t405 ^ 2 + (t478 * t403 + (t477 - t481) * t405) * t403) * qJD(3) + (t403 * t479 - t405 * t480) * qJD(1)) * t453 / 0.2e1 + (Icges(2,3) + Icges(3,3) + m(2) * (t393 ^ 2 + t394 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t339 * t399 + t341 * t398) * t380 + (t338 * t399 + t340 * t398) * t381 + ((-t346 * t404 - t348 * t402 - t360 * t416 - t362 * t413) * t405 + (t347 * t404 + t349 * t402 + t361 * t416 + t363 * t413) * t403) * qJD(3) + (t399 * t375 + t398 * t376 + t404 * t383 + t402 * t384 + t416 * t390 + t413 * t391) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
