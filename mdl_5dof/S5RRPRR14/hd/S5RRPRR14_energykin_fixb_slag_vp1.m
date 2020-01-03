% Calculate kinetic energy for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR14_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR14_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR14_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:45
% EndTime: 2019-12-31 20:35:47
% DurationCPUTime: 2.00s
% Computational Cost: add. (2015->305), mult. (3799->472), div. (0->0), fcn. (4522->12), ass. (0->141)
t428 = sin(qJ(2));
t429 = sin(qJ(1));
t431 = cos(qJ(2));
t432 = cos(qJ(1));
t462 = cos(pkin(5));
t446 = t432 * t462;
t405 = t428 * t429 - t431 * t446;
t406 = t428 * t446 + t429 * t431;
t424 = sin(pkin(5));
t458 = t424 * t432;
t367 = Icges(3,5) * t406 - Icges(3,6) * t405 - Icges(3,3) * t458;
t447 = t429 * t462;
t407 = t432 * t428 + t431 * t447;
t408 = -t428 * t447 + t432 * t431;
t460 = t424 * t429;
t368 = Icges(3,5) * t408 - Icges(3,6) * t407 + Icges(3,3) * t460;
t469 = (t367 * t432 - t368 * t429) * t424;
t423 = sin(pkin(10));
t425 = cos(pkin(10));
t387 = -t406 * t423 - t425 * t458;
t450 = t423 * t458;
t388 = t406 * t425 - t450;
t341 = Icges(4,5) * t388 + Icges(4,6) * t387 + Icges(4,3) * t405;
t369 = Icges(3,4) * t406 - Icges(3,2) * t405 - Icges(3,6) * t458;
t468 = -t341 + t369;
t389 = -t408 * t423 + t425 * t460;
t451 = t423 * t460;
t390 = t408 * t425 + t451;
t342 = Icges(4,5) * t390 + Icges(4,6) * t389 + Icges(4,3) * t407;
t370 = Icges(3,4) * t408 - Icges(3,2) * t407 + Icges(3,6) * t460;
t467 = t342 - t370;
t461 = t424 * t428;
t403 = -t423 * t461 + t425 * t462;
t445 = t462 * t423;
t404 = t425 * t461 + t445;
t459 = t424 * t431;
t361 = Icges(4,5) * t404 + Icges(4,6) * t403 - Icges(4,3) * t459;
t394 = Icges(3,6) * t462 + (Icges(3,4) * t428 + Icges(3,2) * t431) * t424;
t466 = t361 - t394;
t463 = pkin(3) * t425;
t379 = pkin(2) * t406 + qJ(3) * t405;
t380 = pkin(2) * t408 + qJ(3) * t407;
t454 = qJD(2) * t424;
t418 = t429 * t454;
t448 = t432 * t454;
t456 = t379 * t418 + t380 * t448;
t391 = qJD(4) * t407 + t418;
t455 = qJD(1) * (pkin(1) * t429 - pkin(7) * t458);
t453 = qJD(3) * t431;
t419 = qJD(2) * t462 + qJD(1);
t452 = pkin(10) + qJ(4);
t411 = qJD(1) * (pkin(1) * t432 + pkin(7) * t460);
t449 = qJD(3) * t405 + t419 * t380 + t411;
t444 = cos(t452);
t443 = qJD(3) * t407 - t455;
t409 = (pkin(2) * t428 - qJ(3) * t431) * t424;
t440 = (-rSges(4,1) * t404 - rSges(4,2) * t403 + rSges(4,3) * t459 - t409) * t454;
t439 = (-pkin(3) * t445 - (-pkin(8) * t431 + t428 * t463) * t424 - t409) * t454;
t392 = qJD(4) * t405 - t448;
t438 = t424 * t444;
t410 = -qJD(4) * t459 + t419;
t339 = -pkin(3) * t450 + pkin(8) * t405 + t406 * t463;
t340 = pkin(3) * t451 + pkin(8) * t407 + t408 * t463;
t436 = t339 * t418 + t340 * t448 - t424 * t453 + t456;
t435 = t419 * t340 + t429 * t439 + t449;
t434 = (-t339 - t379) * t419 + t432 * t439 + t443;
t430 = cos(qJ(5));
t427 = sin(qJ(5));
t422 = sin(t452);
t415 = rSges(2,1) * t432 - rSges(2,2) * t429;
t414 = rSges(2,1) * t429 + rSges(2,2) * t432;
t398 = t462 * rSges(3,3) + (rSges(3,1) * t428 + rSges(3,2) * t431) * t424;
t397 = t422 * t462 + t428 * t438;
t396 = t422 * t461 - t444 * t462;
t395 = Icges(3,5) * t462 + (Icges(3,1) * t428 + Icges(3,4) * t431) * t424;
t393 = Icges(3,3) * t462 + (Icges(3,5) * t428 + Icges(3,6) * t431) * t424;
t386 = t408 * t444 + t422 * t460;
t385 = t408 * t422 - t429 * t438;
t384 = t406 * t444 - t422 * t458;
t383 = t406 * t422 + t432 * t438;
t382 = t397 * t430 - t427 * t459;
t381 = -t397 * t427 - t430 * t459;
t378 = qJD(5) * t396 + t410;
t377 = rSges(3,1) * t408 - rSges(3,2) * t407 + rSges(3,3) * t460;
t376 = rSges(3,1) * t406 - rSges(3,2) * t405 - rSges(3,3) * t458;
t372 = Icges(3,1) * t408 - Icges(3,4) * t407 + Icges(3,5) * t460;
t371 = Icges(3,1) * t406 - Icges(3,4) * t405 - Icges(3,5) * t458;
t365 = pkin(4) * t397 + pkin(9) * t396;
t363 = Icges(4,1) * t404 + Icges(4,4) * t403 - Icges(4,5) * t459;
t362 = Icges(4,4) * t404 + Icges(4,2) * t403 - Icges(4,6) * t459;
t360 = rSges(5,1) * t397 - rSges(5,2) * t396 - rSges(5,3) * t459;
t359 = Icges(5,1) * t397 - Icges(5,4) * t396 - Icges(5,5) * t459;
t358 = Icges(5,4) * t397 - Icges(5,2) * t396 - Icges(5,6) * t459;
t357 = Icges(5,5) * t397 - Icges(5,6) * t396 - Icges(5,3) * t459;
t356 = t386 * t430 + t407 * t427;
t355 = -t386 * t427 + t407 * t430;
t354 = t384 * t430 + t405 * t427;
t353 = -t384 * t427 + t405 * t430;
t352 = qJD(5) * t383 + t392;
t351 = qJD(5) * t385 + t391;
t350 = pkin(4) * t386 + pkin(9) * t385;
t349 = pkin(4) * t384 + pkin(9) * t383;
t348 = rSges(4,1) * t390 + rSges(4,2) * t389 + rSges(4,3) * t407;
t347 = rSges(4,1) * t388 + rSges(4,2) * t387 + rSges(4,3) * t405;
t346 = Icges(4,1) * t390 + Icges(4,4) * t389 + Icges(4,5) * t407;
t345 = Icges(4,1) * t388 + Icges(4,4) * t387 + Icges(4,5) * t405;
t344 = Icges(4,4) * t390 + Icges(4,2) * t389 + Icges(4,6) * t407;
t343 = Icges(4,4) * t388 + Icges(4,2) * t387 + Icges(4,6) * t405;
t338 = rSges(5,1) * t386 - rSges(5,2) * t385 + rSges(5,3) * t407;
t337 = rSges(5,1) * t384 - rSges(5,2) * t383 + rSges(5,3) * t405;
t336 = Icges(5,1) * t386 - Icges(5,4) * t385 + Icges(5,5) * t407;
t335 = Icges(5,1) * t384 - Icges(5,4) * t383 + Icges(5,5) * t405;
t334 = Icges(5,4) * t386 - Icges(5,2) * t385 + Icges(5,6) * t407;
t333 = Icges(5,4) * t384 - Icges(5,2) * t383 + Icges(5,6) * t405;
t332 = Icges(5,5) * t386 - Icges(5,6) * t385 + Icges(5,3) * t407;
t331 = Icges(5,5) * t384 - Icges(5,6) * t383 + Icges(5,3) * t405;
t330 = t377 * t419 - t398 * t418 + t411;
t329 = -t376 * t419 - t398 * t448 - t455;
t328 = rSges(6,1) * t382 + rSges(6,2) * t381 + rSges(6,3) * t396;
t327 = Icges(6,1) * t382 + Icges(6,4) * t381 + Icges(6,5) * t396;
t326 = Icges(6,4) * t382 + Icges(6,2) * t381 + Icges(6,6) * t396;
t325 = Icges(6,5) * t382 + Icges(6,6) * t381 + Icges(6,3) * t396;
t321 = (t376 * t429 + t377 * t432) * t454;
t320 = rSges(6,1) * t356 + rSges(6,2) * t355 + rSges(6,3) * t385;
t319 = rSges(6,1) * t354 + rSges(6,2) * t353 + rSges(6,3) * t383;
t318 = Icges(6,1) * t356 + Icges(6,4) * t355 + Icges(6,5) * t385;
t317 = Icges(6,1) * t354 + Icges(6,4) * t353 + Icges(6,5) * t383;
t316 = Icges(6,4) * t356 + Icges(6,2) * t355 + Icges(6,6) * t385;
t315 = Icges(6,4) * t354 + Icges(6,2) * t353 + Icges(6,6) * t383;
t314 = Icges(6,5) * t356 + Icges(6,6) * t355 + Icges(6,3) * t385;
t313 = Icges(6,5) * t354 + Icges(6,6) * t353 + Icges(6,3) * t383;
t312 = t348 * t419 + t429 * t440 + t449;
t311 = (-t347 - t379) * t419 + t432 * t440 + t443;
t310 = (-t453 + (t347 * t429 + t348 * t432) * qJD(2)) * t424 + t456;
t309 = t338 * t410 - t360 * t391 + t435;
t308 = -t337 * t410 + t360 * t392 + t434;
t307 = t337 * t391 - t338 * t392 + t436;
t306 = t320 * t378 - t328 * t351 + t350 * t410 - t365 * t391 + t435;
t305 = -t319 * t378 + t328 * t352 - t349 * t410 + t365 * t392 + t434;
t304 = t319 * t351 - t320 * t352 + t349 * t391 - t350 * t392 + t436;
t1 = m(3) * (t321 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(4) * (t310 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + m(5) * (t307 ^ 2 + t308 ^ 2 + t309 ^ 2) / 0.2e1 + t391 * ((t407 * t332 - t385 * t334 + t386 * t336) * t391 + (t331 * t407 - t333 * t385 + t335 * t386) * t392 + (t357 * t407 - t358 * t385 + t359 * t386) * t410) / 0.2e1 + t392 * ((t332 * t405 - t334 * t383 + t336 * t384) * t391 + (t405 * t331 - t383 * t333 + t384 * t335) * t392 + (t357 * t405 - t358 * t383 + t359 * t384) * t410) / 0.2e1 + t410 * ((-t332 * t459 - t334 * t396 + t336 * t397) * t391 + (-t331 * t459 - t333 * t396 + t335 * t397) * t392 + (-t357 * t459 - t396 * t358 + t397 * t359) * t410) / 0.2e1 + m(6) * (t304 ^ 2 + t305 ^ 2 + t306 ^ 2) / 0.2e1 + t351 * ((t385 * t314 + t355 * t316 + t356 * t318) * t351 + (t313 * t385 + t315 * t355 + t317 * t356) * t352 + (t325 * t385 + t326 * t355 + t327 * t356) * t378) / 0.2e1 + t352 * ((t314 * t383 + t316 * t353 + t318 * t354) * t351 + (t383 * t313 + t353 * t315 + t354 * t317) * t352 + (t325 * t383 + t326 * t353 + t327 * t354) * t378) / 0.2e1 + t378 * ((t314 * t396 + t316 * t381 + t318 * t382) * t351 + (t313 * t396 + t315 * t381 + t317 * t382) * t352 + (t396 * t325 + t381 * t326 + t382 * t327) * t378) / 0.2e1 + ((t462 * t368 + (t370 * t431 + t372 * t428) * t424) * t418 - (t462 * t367 + (t369 * t431 + t371 * t428) * t424) * t448 + ((t344 * t403 + t346 * t404) * t429 - (t343 * t403 + t345 * t404) * t432 + (t341 * t432 - t342 * t429) * t459) * t454 + (t462 * t393 + (t394 * t431 + t395 * t428) * t424 - t361 * t459 + t403 * t362 + t404 * t363) * t419) * t419 / 0.2e1 + (m(2) * (t414 ^ 2 + t415 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t343 * t389 - t345 * t390 - t371 * t408 + t468 * t407) * t432 + (t344 * t389 + t346 * t390 + t408 * t372 + t467 * t407 - t469) * t429) * t454 + (t362 * t389 + t363 * t390 + t393 * t460 + t395 * t408 + t466 * t407) * t419) * t418 / 0.2e1 - (((-t343 * t387 - t345 * t388 - t406 * t371 + t468 * t405 + t469) * t432 + (t344 * t387 + t346 * t388 + t372 * t406 + t467 * t405) * t429) * t454 + (t362 * t387 + t363 * t388 - t393 * t458 + t406 * t395 + t466 * t405) * t419) * t448 / 0.2e1;
T = t1;
