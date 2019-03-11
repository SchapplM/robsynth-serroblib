% Calculate kinetic energy for
% S6RPRRRR1
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:54:08
% EndTime: 2019-03-09 06:54:10
% DurationCPUTime: 1.51s
% Computational Cost: add. (1903->246), mult. (1355->397), div. (0->0), fcn. (1232->12), ass. (0->144)
t412 = sin(qJ(1));
t470 = pkin(1) * t412;
t409 = qJ(3) + qJ(4);
t403 = sin(t409);
t469 = pkin(4) * t403;
t414 = cos(qJ(3));
t467 = t414 * pkin(3);
t411 = sin(qJ(3));
t465 = Icges(4,4) * t411;
t464 = Icges(4,4) * t414;
t463 = Icges(5,4) * t403;
t404 = cos(t409);
t462 = Icges(5,4) * t404;
t405 = qJ(5) + t409;
t397 = sin(t405);
t461 = Icges(6,4) * t397;
t398 = cos(t405);
t460 = Icges(6,4) * t398;
t407 = qJ(1) + pkin(11);
t401 = sin(t407);
t459 = t397 * t401;
t402 = cos(t407);
t458 = t397 * t402;
t410 = sin(qJ(6));
t457 = t401 * t410;
t413 = cos(qJ(6));
t456 = t401 * t413;
t455 = t402 * t410;
t454 = t402 * t413;
t415 = cos(qJ(1));
t400 = qJD(1) * t415 * pkin(1);
t453 = qJD(1) * (pkin(2) * t402 + pkin(7) * t401) + t400;
t452 = pkin(4) * t404;
t396 = qJD(3) * t401;
t379 = qJD(4) * t401 + t396;
t450 = qJD(3) * t402;
t449 = qJD(6) * t397;
t448 = -qJD(3) - qJD(4);
t447 = pkin(3) * qJD(3) * t411;
t370 = qJD(5) * t401 + t379;
t332 = -pkin(8) * t402 + t467 * t401;
t333 = pkin(8) * t401 + t467 * t402;
t446 = t332 * t396 + t333 * t450 + qJD(2);
t445 = -pkin(2) * t401 + pkin(7) * t402 - t470;
t444 = t402 * t447;
t443 = -t332 + t445;
t442 = pkin(5) * t398 + pkin(10) * t397;
t441 = rSges(4,1) * t414 - rSges(4,2) * t411;
t440 = rSges(5,1) * t404 - rSges(5,2) * t403;
t439 = rSges(6,1) * t398 - rSges(6,2) * t397;
t371 = (-qJD(5) + t448) * t402;
t438 = Icges(4,1) * t414 - t465;
t437 = Icges(5,1) * t404 - t463;
t436 = Icges(6,1) * t398 - t461;
t435 = -Icges(4,2) * t411 + t464;
t434 = -Icges(5,2) * t403 + t462;
t433 = -Icges(6,2) * t397 + t460;
t432 = Icges(4,5) * t414 - Icges(4,6) * t411;
t431 = Icges(5,5) * t404 - Icges(5,6) * t403;
t430 = Icges(6,5) * t398 - Icges(6,6) * t397;
t356 = -Icges(4,6) * t402 + t435 * t401;
t358 = -Icges(4,5) * t402 + t438 * t401;
t429 = t356 * t411 - t358 * t414;
t357 = Icges(4,6) * t401 + t435 * t402;
t359 = Icges(4,5) * t401 + t438 * t402;
t428 = -t357 * t411 + t359 * t414;
t389 = Icges(4,2) * t414 + t465;
t390 = Icges(4,1) * t411 + t464;
t427 = -t389 * t411 + t390 * t414;
t325 = -pkin(9) * t402 + t452 * t401;
t426 = -t325 + t443;
t380 = t448 * t402;
t425 = t380 * t469 - t444;
t326 = pkin(9) * t401 + t452 * t402;
t424 = t379 * t325 - t326 * t380 + t446;
t423 = qJD(1) * t333 - t401 * t447 + t453;
t422 = (Icges(6,5) * t397 + Icges(6,6) * t398) * qJD(1) + (-Icges(6,3) * t402 + t430 * t401) * t371 + (Icges(6,3) * t401 + t430 * t402) * t370;
t421 = (Icges(5,5) * t403 + Icges(5,6) * t404) * qJD(1) + (-Icges(5,3) * t402 + t431 * t401) * t380 + (Icges(5,3) * t401 + t431 * t402) * t379;
t420 = qJD(1) * t326 - t379 * t469 + t423;
t336 = -Icges(6,6) * t402 + t433 * t401;
t337 = Icges(6,6) * t401 + t433 * t402;
t338 = -Icges(6,5) * t402 + t436 * t401;
t339 = Icges(6,5) * t401 + t436 * t402;
t374 = Icges(6,2) * t398 + t461;
t375 = Icges(6,1) * t397 + t460;
t419 = (-t337 * t397 + t339 * t398) * t370 + (-t336 * t397 + t338 * t398) * t371 + (-t374 * t397 + t375 * t398) * qJD(1);
t344 = -Icges(5,6) * t402 + t434 * t401;
t345 = Icges(5,6) * t401 + t434 * t402;
t346 = -Icges(5,5) * t402 + t437 * t401;
t347 = Icges(5,5) * t401 + t437 * t402;
t383 = Icges(5,2) * t404 + t463;
t384 = Icges(5,1) * t403 + t462;
t418 = (-t345 * t403 + t347 * t404) * t379 + (-t344 * t403 + t346 * t404) * t380 + (-t383 * t403 + t384 * t404) * qJD(1);
t393 = rSges(2,1) * t415 - rSges(2,2) * t412;
t392 = rSges(2,1) * t412 + rSges(2,2) * t415;
t391 = rSges(4,1) * t411 + rSges(4,2) * t414;
t388 = Icges(4,5) * t411 + Icges(4,6) * t414;
t387 = -qJD(6) * t398 + qJD(1);
t385 = rSges(5,1) * t403 + rSges(5,2) * t404;
t378 = pkin(5) * t397 - pkin(10) * t398;
t376 = rSges(6,1) * t397 + rSges(6,2) * t398;
t369 = t400 + qJD(1) * (rSges(3,1) * t402 - rSges(3,2) * t401);
t368 = (-rSges(3,1) * t401 - rSges(3,2) * t402 - t470) * qJD(1);
t367 = t398 * t454 + t457;
t366 = -t398 * t455 + t456;
t365 = t398 * t456 - t455;
t364 = -t398 * t457 - t454;
t363 = t442 * t402;
t362 = t442 * t401;
t361 = rSges(4,3) * t401 + t441 * t402;
t360 = -rSges(4,3) * t402 + t441 * t401;
t355 = Icges(4,3) * t401 + t432 * t402;
t354 = -Icges(4,3) * t402 + t432 * t401;
t353 = -rSges(7,3) * t398 + (rSges(7,1) * t413 - rSges(7,2) * t410) * t397;
t352 = -Icges(7,5) * t398 + (Icges(7,1) * t413 - Icges(7,4) * t410) * t397;
t351 = -Icges(7,6) * t398 + (Icges(7,4) * t413 - Icges(7,2) * t410) * t397;
t350 = -Icges(7,3) * t398 + (Icges(7,5) * t413 - Icges(7,6) * t410) * t397;
t349 = rSges(5,3) * t401 + t440 * t402;
t348 = -rSges(5,3) * t402 + t440 * t401;
t341 = rSges(6,3) * t401 + t439 * t402;
t340 = -rSges(6,3) * t402 + t439 * t401;
t331 = t401 * t449 + t371;
t330 = t402 * t449 + t370;
t322 = rSges(7,1) * t367 + rSges(7,2) * t366 + rSges(7,3) * t458;
t321 = rSges(7,1) * t365 + rSges(7,2) * t364 + rSges(7,3) * t459;
t320 = Icges(7,1) * t367 + Icges(7,4) * t366 + Icges(7,5) * t458;
t319 = Icges(7,1) * t365 + Icges(7,4) * t364 + Icges(7,5) * t459;
t318 = Icges(7,4) * t367 + Icges(7,2) * t366 + Icges(7,6) * t458;
t317 = Icges(7,4) * t365 + Icges(7,2) * t364 + Icges(7,6) * t459;
t316 = Icges(7,5) * t367 + Icges(7,6) * t366 + Icges(7,3) * t458;
t315 = Icges(7,5) * t365 + Icges(7,6) * t364 + Icges(7,3) * t459;
t314 = qJD(1) * t361 - t391 * t396 + t453;
t313 = -t391 * t450 + (-t360 + t445) * qJD(1);
t312 = qJD(2) + (t360 * t401 + t361 * t402) * qJD(3);
t311 = qJD(1) * t349 - t379 * t385 + t423;
t310 = -t444 + t380 * t385 + (-t348 + t443) * qJD(1);
t309 = t348 * t379 - t349 * t380 + t446;
t308 = qJD(1) * t341 - t370 * t376 + t420;
t307 = t371 * t376 + (-t340 + t426) * qJD(1) + t425;
t306 = t340 * t370 - t341 * t371 + t424;
t305 = qJD(1) * t363 + t322 * t387 - t330 * t353 - t370 * t378 + t420;
t304 = -t321 * t387 + t331 * t353 + t371 * t378 + (-t362 + t426) * qJD(1) + t425;
t303 = t321 * t330 - t322 * t331 + t362 * t370 - t363 * t371 + t424;
t1 = ((t401 * t388 + t427 * t402) * qJD(1) + (t401 ^ 2 * t355 + (t429 * t402 + (-t354 + t428) * t401) * t402) * qJD(3)) * t396 / 0.2e1 + t330 * ((t316 * t458 + t366 * t318 + t367 * t320) * t330 + (t315 * t458 + t317 * t366 + t319 * t367) * t331 + (t350 * t458 + t351 * t366 + t352 * t367) * t387) / 0.2e1 + t331 * ((t316 * t459 + t318 * t364 + t320 * t365) * t330 + (t315 * t459 + t364 * t317 + t365 * t319) * t331 + (t350 * t459 + t351 * t364 + t352 * t365) * t387) / 0.2e1 + t387 * ((-t315 * t331 - t316 * t330 - t350 * t387) * t398 + ((-t318 * t410 + t320 * t413) * t330 + (-t317 * t410 + t319 * t413) * t331 + (-t351 * t410 + t352 * t413) * t387) * t397) / 0.2e1 + t370 * (t422 * t401 + t419 * t402) / 0.2e1 + t371 * (t419 * t401 - t422 * t402) / 0.2e1 + t379 * (t421 * t401 + t418 * t402) / 0.2e1 + t380 * (t418 * t401 - t421 * t402) / 0.2e1 - ((-t402 * t388 + t427 * t401) * qJD(1) + (t402 ^ 2 * t354 + (t428 * t401 + (-t355 + t429) * t402) * t401) * qJD(3)) * t450 / 0.2e1 + m(7) * (t303 ^ 2 + t304 ^ 2 + t305 ^ 2) / 0.2e1 + m(6) * (t306 ^ 2 + t307 ^ 2 + t308 ^ 2) / 0.2e1 + m(5) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + m(4) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + (Icges(3,3) + Icges(2,3) + m(2) * (t392 ^ 2 + t393 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t337 * t398 + t339 * t397) * t370 + (t336 * t398 + t338 * t397) * t371 + (t345 * t404 + t347 * t403) * t379 + (t344 * t404 + t346 * t403) * t380 + ((t357 * t414 + t359 * t411) * t401 - (t356 * t414 + t358 * t411) * t402) * qJD(3) + (t398 * t374 + t397 * t375 + t404 * t383 + t403 * t384 + t414 * t389 + t411 * t390) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
