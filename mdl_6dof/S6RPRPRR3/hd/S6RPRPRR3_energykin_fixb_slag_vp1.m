% Calculate kinetic energy for
% S6RPRPRR3
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
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:31
% EndTime: 2019-03-09 03:40:33
% DurationCPUTime: 2.11s
% Computational Cost: add. (1996->296), mult. (1852->466), div. (0->0), fcn. (1841->12), ass. (0->144)
t432 = sin(qJ(1));
t477 = pkin(1) * t432;
t429 = cos(pkin(11));
t474 = t429 * pkin(4);
t431 = sin(qJ(3));
t473 = Icges(4,4) * t431;
t433 = cos(qJ(3));
t472 = Icges(4,4) * t433;
t427 = qJ(1) + pkin(10);
t420 = sin(t427);
t428 = sin(pkin(11));
t471 = t420 * t428;
t470 = t420 * t431;
t422 = cos(t427);
t469 = t422 * t428;
t468 = t422 * t431;
t467 = t433 * t420;
t466 = t433 * t422;
t434 = cos(qJ(1));
t418 = qJD(1) * t434 * pkin(1);
t464 = qJD(1) * (pkin(2) * t422 + pkin(7) * t420) + t418;
t426 = pkin(11) + qJ(5);
t421 = cos(t426);
t463 = pkin(5) * t421;
t413 = qJD(3) * t420;
t459 = qJD(5) * t431;
t395 = t422 * t459 + t413;
t461 = qJD(3) * t422;
t460 = qJD(4) * t431;
t458 = qJD(6) * t431;
t455 = -pkin(2) * t420 + pkin(7) * t422 - t477;
t419 = sin(t426);
t454 = pkin(5) * t419;
t409 = pkin(3) * t431 - qJ(4) * t433;
t453 = qJD(3) * (pkin(8) * t433 - t431 * t474 - t409);
t452 = qJD(3) * (rSges(5,3) * t433 - (rSges(5,1) * t429 - rSges(5,2) * t428) * t431 - t409);
t396 = t420 * t459 - t461;
t448 = pkin(3) * t433 + qJ(4) * t431;
t394 = t448 * t422;
t451 = qJD(1) * t394 + t420 * t460 + t464;
t393 = t448 * t420;
t450 = -t393 + t455;
t449 = rSges(4,1) * t433 - rSges(4,2) * t431;
t447 = Icges(4,1) * t433 - t473;
t446 = -Icges(4,2) * t431 + t472;
t445 = Icges(4,5) * t433 - Icges(4,6) * t431;
t364 = -Icges(4,6) * t422 + t420 * t446;
t366 = -Icges(4,5) * t422 + t420 * t447;
t444 = t364 * t431 - t366 * t433;
t365 = Icges(4,6) * t420 + t422 * t446;
t367 = Icges(4,5) * t420 + t422 * t447;
t443 = -t365 * t431 + t367 * t433;
t403 = Icges(4,2) * t433 + t473;
t404 = Icges(4,1) * t431 + t472;
t442 = -t403 * t431 + t404 * t433;
t441 = -qJD(4) * t433 + t393 * t413 + t394 * t461 + qJD(2);
t439 = pkin(8) * t431 + t433 * t474;
t350 = -pkin(4) * t469 + t420 * t439;
t351 = pkin(4) * t471 + t422 * t439;
t440 = t350 * t413 + t351 * t461 + t441;
t438 = pkin(9) * t431 + t433 * t463;
t437 = qJD(1) * t351 + t420 * t453 + t451;
t408 = t422 * t460;
t436 = t408 + (-t350 + t450) * qJD(1) + t422 * t453;
t423 = qJ(6) + t426;
t416 = cos(t423);
t415 = sin(t423);
t414 = -qJD(5) * t433 + qJD(1);
t412 = rSges(2,1) * t434 - rSges(2,2) * t432;
t411 = rSges(2,1) * t432 + rSges(2,2) * t434;
t410 = rSges(4,1) * t431 + rSges(4,2) * t433;
t402 = Icges(4,5) * t431 + Icges(4,6) * t433;
t401 = qJD(1) + (-qJD(5) - qJD(6)) * t433;
t391 = -Icges(5,5) * t433 + (Icges(5,1) * t429 - Icges(5,4) * t428) * t431;
t390 = -Icges(5,6) * t433 + (Icges(5,4) * t429 - Icges(5,2) * t428) * t431;
t389 = -Icges(5,3) * t433 + (Icges(5,5) * t429 - Icges(5,6) * t428) * t431;
t388 = t429 * t466 + t471;
t387 = t420 * t429 - t428 * t466;
t386 = t429 * t467 - t469;
t385 = -t422 * t429 - t428 * t467;
t383 = t418 + qJD(1) * (rSges(3,1) * t422 - rSges(3,2) * t420);
t382 = (-rSges(3,1) * t420 - rSges(3,2) * t422 - t477) * qJD(1);
t381 = -rSges(6,3) * t433 + (rSges(6,1) * t421 - rSges(6,2) * t419) * t431;
t380 = -Icges(6,5) * t433 + (Icges(6,1) * t421 - Icges(6,4) * t419) * t431;
t379 = -Icges(6,6) * t433 + (Icges(6,4) * t421 - Icges(6,2) * t419) * t431;
t378 = -Icges(6,3) * t433 + (Icges(6,5) * t421 - Icges(6,6) * t419) * t431;
t377 = t419 * t420 + t421 * t466;
t376 = -t419 * t466 + t420 * t421;
t375 = -t419 * t422 + t421 * t467;
t374 = -t419 * t467 - t421 * t422;
t372 = rSges(4,3) * t420 + t422 * t449;
t371 = -rSges(4,3) * t422 + t420 * t449;
t370 = -rSges(7,3) * t433 + (rSges(7,1) * t416 - rSges(7,2) * t415) * t431;
t363 = Icges(4,3) * t420 + t422 * t445;
t362 = -Icges(4,3) * t422 + t420 * t445;
t361 = -Icges(7,5) * t433 + (Icges(7,1) * t416 - Icges(7,4) * t415) * t431;
t360 = -Icges(7,6) * t433 + (Icges(7,4) * t416 - Icges(7,2) * t415) * t431;
t359 = -Icges(7,3) * t433 + (Icges(7,5) * t416 - Icges(7,6) * t415) * t431;
t358 = t415 * t420 + t416 * t466;
t357 = -t415 * t466 + t416 * t420;
t356 = -t415 * t422 + t416 * t467;
t355 = -t415 * t467 - t416 * t422;
t354 = t420 * t458 + t396;
t353 = t422 * t458 + t395;
t352 = -pkin(9) * t433 + t431 * t463;
t349 = rSges(5,1) * t388 + rSges(5,2) * t387 + rSges(5,3) * t468;
t348 = rSges(5,1) * t386 + rSges(5,2) * t385 + rSges(5,3) * t470;
t347 = Icges(5,1) * t388 + Icges(5,4) * t387 + Icges(5,5) * t468;
t346 = Icges(5,1) * t386 + Icges(5,4) * t385 + Icges(5,5) * t470;
t345 = Icges(5,4) * t388 + Icges(5,2) * t387 + Icges(5,6) * t468;
t344 = Icges(5,4) * t386 + Icges(5,2) * t385 + Icges(5,6) * t470;
t343 = Icges(5,5) * t388 + Icges(5,6) * t387 + Icges(5,3) * t468;
t342 = Icges(5,5) * t386 + Icges(5,6) * t385 + Icges(5,3) * t470;
t338 = rSges(6,1) * t377 + rSges(6,2) * t376 + rSges(6,3) * t468;
t337 = rSges(6,1) * t375 + rSges(6,2) * t374 + rSges(6,3) * t470;
t336 = Icges(6,1) * t377 + Icges(6,4) * t376 + Icges(6,5) * t468;
t335 = Icges(6,1) * t375 + Icges(6,4) * t374 + Icges(6,5) * t470;
t334 = Icges(6,4) * t377 + Icges(6,2) * t376 + Icges(6,6) * t468;
t333 = Icges(6,4) * t375 + Icges(6,2) * t374 + Icges(6,6) * t470;
t332 = Icges(6,5) * t377 + Icges(6,6) * t376 + Icges(6,3) * t468;
t331 = Icges(6,5) * t375 + Icges(6,6) * t374 + Icges(6,3) * t470;
t330 = qJD(1) * t372 - t410 * t413 + t464;
t329 = -t410 * t461 + (-t371 + t455) * qJD(1);
t328 = rSges(7,1) * t358 + rSges(7,2) * t357 + rSges(7,3) * t468;
t327 = rSges(7,1) * t356 + rSges(7,2) * t355 + rSges(7,3) * t470;
t326 = Icges(7,1) * t358 + Icges(7,4) * t357 + Icges(7,5) * t468;
t325 = Icges(7,1) * t356 + Icges(7,4) * t355 + Icges(7,5) * t470;
t324 = Icges(7,4) * t358 + Icges(7,2) * t357 + Icges(7,6) * t468;
t323 = Icges(7,4) * t356 + Icges(7,2) * t355 + Icges(7,6) * t470;
t322 = Icges(7,5) * t358 + Icges(7,6) * t357 + Icges(7,3) * t468;
t321 = Icges(7,5) * t356 + Icges(7,6) * t355 + Icges(7,3) * t470;
t320 = qJD(2) + (t371 * t420 + t372 * t422) * qJD(3);
t319 = t420 * t454 + t422 * t438;
t318 = t420 * t438 - t422 * t454;
t317 = qJD(1) * t349 + t420 * t452 + t451;
t316 = t408 + t422 * t452 + (-t348 + t450) * qJD(1);
t315 = (t348 * t420 + t349 * t422) * qJD(3) + t441;
t314 = t338 * t414 - t381 * t395 + t437;
t313 = -t337 * t414 + t381 * t396 + t436;
t312 = t337 * t395 - t338 * t396 + t440;
t311 = t319 * t414 + t328 * t401 - t352 * t395 - t353 * t370 + t437;
t310 = -t318 * t414 - t327 * t401 + t352 * t396 + t354 * t370 + t436;
t309 = t318 * t395 - t319 * t396 + t327 * t353 - t328 * t354 + t440;
t1 = t401 * ((-t321 * t354 - t322 * t353 - t359 * t401) * t433 + ((-t324 * t415 + t326 * t416) * t353 + (-t323 * t415 + t325 * t416) * t354 + (-t360 * t415 + t361 * t416) * t401) * t431) / 0.2e1 + t353 * ((t322 * t468 + t357 * t324 + t358 * t326) * t353 + (t321 * t468 + t323 * t357 + t325 * t358) * t354 + (t357 * t360 + t358 * t361 + t359 * t468) * t401) / 0.2e1 + t354 * ((t322 * t470 + t324 * t355 + t326 * t356) * t353 + (t321 * t470 + t355 * t323 + t356 * t325) * t354 + (t355 * t360 + t356 * t361 + t359 * t470) * t401) / 0.2e1 + t414 * ((-t331 * t396 - t332 * t395 - t378 * t414) * t433 + ((-t334 * t419 + t336 * t421) * t395 + (-t333 * t419 + t335 * t421) * t396 + (-t379 * t419 + t380 * t421) * t414) * t431) / 0.2e1 + t396 * ((t332 * t470 + t334 * t374 + t336 * t375) * t395 + (t331 * t470 + t374 * t333 + t375 * t335) * t396 + (t374 * t379 + t375 * t380 + t378 * t470) * t414) / 0.2e1 + t395 * ((t332 * t468 + t376 * t334 + t377 * t336) * t395 + (t331 * t468 + t333 * t376 + t335 * t377) * t396 + (t376 * t379 + t377 * t380 + t378 * t468) * t414) / 0.2e1 + m(7) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + m(6) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(5) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(4) * (t320 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + (((t365 * t433 + t367 * t431) * t420 - (t364 * t433 + t366 * t431) * t422 + (t342 * t422 - t343 * t420) * t433 + ((-t345 * t428 + t347 * t429) * t420 - (-t344 * t428 + t346 * t429) * t422) * t431) * qJD(3) + ((t403 - t389) * t433 + (-t390 * t428 + t391 * t429 + t404) * t431) * qJD(1)) * qJD(1) / 0.2e1 + (((-t342 * t468 - t344 * t387 - t346 * t388 + t444 * t422) * t422 + ((-t362 + t443) * t422 + t343 * t468 + t345 * t387 + t347 * t388 + t363 * t420) * t420) * qJD(3) + (t387 * t390 + t388 * t391 + t389 * t468 + t420 * t402 + t442 * t422) * qJD(1)) * t413 / 0.2e1 - (((-t342 * t470 - t344 * t385 - t346 * t386 + t362 * t422) * t422 + ((-t363 + t444) * t422 + t343 * t470 + t345 * t385 + t347 * t386 + t443 * t420) * t420) * qJD(3) + (t385 * t390 + t386 * t391 + t389 * t470 - t422 * t402 + t420 * t442) * qJD(1)) * t461 / 0.2e1 + (Icges(3,3) + Icges(2,3) + m(2) * (t411 ^ 2 + t412 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
