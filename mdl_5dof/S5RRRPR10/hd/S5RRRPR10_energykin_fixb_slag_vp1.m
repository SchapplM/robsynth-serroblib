% Calculate kinetic energy for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR10_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR10_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:36
% EndTime: 2019-12-31 21:26:38
% DurationCPUTime: 2.30s
% Computational Cost: add. (2060->301), mult. (3889->457), div. (0->0), fcn. (4630->12), ass. (0->135)
t469 = Icges(4,3) + Icges(5,3);
t429 = sin(qJ(2));
t430 = sin(qJ(1));
t433 = cos(qJ(2));
t434 = cos(qJ(1));
t458 = cos(pkin(5));
t445 = t434 * t458;
t407 = t429 * t430 - t433 * t445;
t408 = t429 * t445 + t430 * t433;
t425 = sin(pkin(5));
t454 = t434 * t425;
t369 = Icges(3,5) * t408 - Icges(3,6) * t407 - Icges(3,3) * t454;
t446 = t430 * t458;
t409 = t434 * t429 + t433 * t446;
t410 = -t429 * t446 + t434 * t433;
t455 = t430 * t425;
t370 = Icges(3,5) * t410 - Icges(3,6) * t409 + Icges(3,3) * t455;
t468 = (t369 * t434 - t370 * t430) * t425;
t450 = qJ(3) + pkin(10);
t424 = sin(t450);
t443 = cos(t450);
t442 = t425 * t443;
t385 = t408 * t424 + t434 * t442;
t386 = t408 * t443 - t424 * t454;
t428 = sin(qJ(3));
t432 = cos(qJ(3));
t389 = -t408 * t428 - t432 * t454;
t448 = t428 * t454;
t390 = t408 * t432 - t448;
t467 = Icges(4,5) * t390 + Icges(5,5) * t386 + Icges(4,6) * t389 - Icges(5,6) * t385 + t469 * t407;
t387 = t410 * t424 - t430 * t442;
t388 = t410 * t443 + t424 * t455;
t391 = -t410 * t428 + t432 * t455;
t449 = t428 * t455;
t392 = t410 * t432 + t449;
t466 = Icges(4,5) * t392 + Icges(5,5) * t388 + Icges(4,6) * t391 - Icges(5,6) * t387 + t469 * t409;
t457 = t425 * t429;
t398 = t424 * t457 - t443 * t458;
t399 = t424 * t458 + t429 * t442;
t405 = -t428 * t457 + t432 * t458;
t444 = t458 * t428;
t406 = t432 * t457 + t444;
t456 = t425 * t433;
t465 = Icges(4,5) * t406 + Icges(5,5) * t399 + Icges(4,6) * t405 - Icges(5,6) * t398 - t469 * t456;
t460 = pkin(3) * t432;
t381 = pkin(2) * t408 + pkin(8) * t407;
t382 = pkin(2) * t410 + pkin(8) * t409;
t451 = qJD(2) * t425;
t420 = t430 * t451;
t447 = t434 * t451;
t453 = t381 * t420 + t382 * t447;
t393 = qJD(3) * t409 + t420;
t452 = qJD(1) * (pkin(1) * t430 - pkin(7) * t454);
t421 = qJD(2) * t458 + qJD(1);
t394 = qJD(3) * t407 - t447;
t412 = -qJD(3) * t456 + t421;
t411 = (pkin(2) * t429 - pkin(8) * t433) * t425;
t413 = qJD(1) * (pkin(1) * t434 + pkin(7) * t455);
t440 = t421 * t382 - t411 * t420 + t413;
t340 = -pkin(3) * t448 + qJ(4) * t407 + t408 * t460;
t439 = -qJD(4) * t456 + t393 * t340 + t453;
t341 = pkin(3) * t449 + qJ(4) * t409 + t410 * t460;
t438 = qJD(4) * t407 + t412 * t341 + t440;
t437 = -t381 * t421 - t411 * t447 - t452;
t367 = pkin(3) * t444 + (-qJ(4) * t433 + t429 * t460) * t425;
t436 = qJD(4) * t409 + t394 * t367 + t437;
t431 = cos(qJ(5));
t427 = sin(qJ(5));
t417 = rSges(2,1) * t434 - rSges(2,2) * t430;
t416 = rSges(2,1) * t430 + rSges(2,2) * t434;
t400 = t458 * rSges(3,3) + (rSges(3,1) * t429 + rSges(3,2) * t433) * t425;
t397 = Icges(3,5) * t458 + (Icges(3,1) * t429 + Icges(3,4) * t433) * t425;
t396 = Icges(3,6) * t458 + (Icges(3,4) * t429 + Icges(3,2) * t433) * t425;
t395 = Icges(3,3) * t458 + (Icges(3,5) * t429 + Icges(3,6) * t433) * t425;
t384 = t399 * t431 - t427 * t456;
t383 = -t399 * t427 - t431 * t456;
t380 = qJD(5) * t398 + t412;
t377 = rSges(3,1) * t410 - rSges(3,2) * t409 + rSges(3,3) * t455;
t376 = rSges(3,1) * t408 - rSges(3,2) * t407 - rSges(3,3) * t454;
t374 = Icges(3,1) * t410 - Icges(3,4) * t409 + Icges(3,5) * t455;
t373 = Icges(3,1) * t408 - Icges(3,4) * t407 - Icges(3,5) * t454;
t372 = Icges(3,4) * t410 - Icges(3,2) * t409 + Icges(3,6) * t455;
t371 = Icges(3,4) * t408 - Icges(3,2) * t407 - Icges(3,6) * t454;
t368 = rSges(4,1) * t406 + rSges(4,2) * t405 - rSges(4,3) * t456;
t366 = Icges(4,1) * t406 + Icges(4,4) * t405 - Icges(4,5) * t456;
t365 = Icges(4,4) * t406 + Icges(4,2) * t405 - Icges(4,6) * t456;
t363 = pkin(4) * t399 + pkin(9) * t398;
t362 = rSges(5,1) * t399 - rSges(5,2) * t398 - rSges(5,3) * t456;
t361 = Icges(5,1) * t399 - Icges(5,4) * t398 - Icges(5,5) * t456;
t360 = Icges(5,4) * t399 - Icges(5,2) * t398 - Icges(5,6) * t456;
t358 = t388 * t431 + t409 * t427;
t357 = -t388 * t427 + t409 * t431;
t356 = t386 * t431 + t407 * t427;
t355 = -t386 * t427 + t407 * t431;
t354 = qJD(5) * t385 + t394;
t353 = qJD(5) * t387 + t393;
t352 = pkin(4) * t388 + pkin(9) * t387;
t351 = pkin(4) * t386 + pkin(9) * t385;
t349 = rSges(4,1) * t392 + rSges(4,2) * t391 + rSges(4,3) * t409;
t348 = rSges(4,1) * t390 + rSges(4,2) * t389 + rSges(4,3) * t407;
t347 = Icges(4,1) * t392 + Icges(4,4) * t391 + Icges(4,5) * t409;
t346 = Icges(4,1) * t390 + Icges(4,4) * t389 + Icges(4,5) * t407;
t345 = Icges(4,4) * t392 + Icges(4,2) * t391 + Icges(4,6) * t409;
t344 = Icges(4,4) * t390 + Icges(4,2) * t389 + Icges(4,6) * t407;
t339 = rSges(5,1) * t388 - rSges(5,2) * t387 + rSges(5,3) * t409;
t338 = rSges(5,1) * t386 - rSges(5,2) * t385 + rSges(5,3) * t407;
t337 = Icges(5,1) * t388 - Icges(5,4) * t387 + Icges(5,5) * t409;
t336 = Icges(5,1) * t386 - Icges(5,4) * t385 + Icges(5,5) * t407;
t335 = Icges(5,4) * t388 - Icges(5,2) * t387 + Icges(5,6) * t409;
t334 = Icges(5,4) * t386 - Icges(5,2) * t385 + Icges(5,6) * t407;
t331 = t377 * t421 - t400 * t420 + t413;
t330 = -t376 * t421 - t400 * t447 - t452;
t329 = rSges(6,1) * t384 + rSges(6,2) * t383 + rSges(6,3) * t398;
t328 = Icges(6,1) * t384 + Icges(6,4) * t383 + Icges(6,5) * t398;
t327 = Icges(6,4) * t384 + Icges(6,2) * t383 + Icges(6,6) * t398;
t326 = Icges(6,5) * t384 + Icges(6,6) * t383 + Icges(6,3) * t398;
t324 = (t376 * t430 + t377 * t434) * t451;
t322 = rSges(6,1) * t358 + rSges(6,2) * t357 + rSges(6,3) * t387;
t321 = rSges(6,1) * t356 + rSges(6,2) * t355 + rSges(6,3) * t385;
t320 = Icges(6,1) * t358 + Icges(6,4) * t357 + Icges(6,5) * t387;
t319 = Icges(6,1) * t356 + Icges(6,4) * t355 + Icges(6,5) * t385;
t318 = Icges(6,4) * t358 + Icges(6,2) * t357 + Icges(6,6) * t387;
t317 = Icges(6,4) * t356 + Icges(6,2) * t355 + Icges(6,6) * t385;
t316 = Icges(6,5) * t358 + Icges(6,6) * t357 + Icges(6,3) * t387;
t315 = Icges(6,5) * t356 + Icges(6,6) * t355 + Icges(6,3) * t385;
t314 = t349 * t412 - t368 * t393 + t440;
t313 = -t348 * t412 + t368 * t394 + t437;
t312 = t348 * t393 - t349 * t394 + t453;
t311 = t339 * t412 + (-t362 - t367) * t393 + t438;
t310 = t362 * t394 + (-t338 - t340) * t412 + t436;
t309 = t338 * t393 + (-t339 - t341) * t394 + t439;
t308 = t322 * t380 - t329 * t353 + t352 * t412 + (-t363 - t367) * t393 + t438;
t307 = -t321 * t380 + t329 * t354 + t363 * t394 + (-t340 - t351) * t412 + t436;
t306 = t321 * t353 - t322 * t354 + t351 * t393 + (-t341 - t352) * t394 + t439;
t1 = m(3) * (t324 ^ 2 + t330 ^ 2 + t331 ^ 2) / 0.2e1 + ((t395 * t455 - t396 * t409 + t397 * t410) * t421 + (-(-t371 * t409 + t373 * t410) * t434 + (-t409 * t372 + t410 * t374 - t468) * t430) * t451) * t420 / 0.2e1 - ((-t395 * t454 - t396 * t407 + t397 * t408) * t421 + ((-t372 * t407 + t374 * t408) * t430 + (t407 * t371 - t408 * t373 + t468) * t434) * t451) * t447 / 0.2e1 + t421 * ((t458 * t370 + (t372 * t433 + t374 * t429) * t425) * t420 - (t458 * t369 + (t371 * t433 + t373 * t429) * t425) * t447 + (t458 * t395 + (t396 * t433 + t397 * t429) * t425) * t421) / 0.2e1 + m(4) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(5) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + m(6) * (t306 ^ 2 + t307 ^ 2 + t308 ^ 2) / 0.2e1 + t353 * ((t387 * t316 + t357 * t318 + t358 * t320) * t353 + (t315 * t387 + t317 * t357 + t319 * t358) * t354 + (t326 * t387 + t327 * t357 + t328 * t358) * t380) / 0.2e1 + t354 * ((t316 * t385 + t318 * t355 + t320 * t356) * t353 + (t385 * t315 + t355 * t317 + t356 * t319) * t354 + (t326 * t385 + t327 * t355 + t328 * t356) * t380) / 0.2e1 + t380 * ((t316 * t398 + t318 * t383 + t320 * t384) * t353 + (t315 * t398 + t317 * t383 + t319 * t384) * t354 + (t398 * t326 + t383 * t327 + t384 * t328) * t380) / 0.2e1 + ((-t360 * t387 + t361 * t388 + t365 * t391 + t366 * t392 + t465 * t409) * t412 + (-t334 * t387 + t336 * t388 + t344 * t391 + t346 * t392 + t467 * t409) * t394 + (-t387 * t335 + t388 * t337 + t391 * t345 + t392 * t347 + t466 * t409) * t393) * t393 / 0.2e1 + ((-t360 * t385 + t361 * t386 + t365 * t389 + t366 * t390 + t465 * t407) * t412 + (-t385 * t334 + t386 * t336 + t389 * t344 + t390 * t346 + t467 * t407) * t394 + (-t335 * t385 + t337 * t386 + t345 * t389 + t347 * t390 + t466 * t407) * t393) * t394 / 0.2e1 + ((-t360 * t398 + t361 * t399 + t365 * t405 + t366 * t406 - t465 * t456) * t412 + (-t334 * t398 + t336 * t399 + t344 * t405 + t346 * t406 - t467 * t456) * t394 + (-t335 * t398 + t337 * t399 + t345 * t405 + t347 * t406 - t466 * t456) * t393) * t412 / 0.2e1 + (m(2) * (t416 ^ 2 + t417 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
