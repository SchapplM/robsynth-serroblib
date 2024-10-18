% Calculate kinetic energy for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR14_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR14_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR14_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:19
% EndTime: 2024-09-27 18:42:19
% DurationCPUTime: 0.36s
% Computational Cost: add. (2509->256), mult. (1778->394), div. (0->0), fcn. (1648->22), ass. (0->136)
t423 = sin(pkin(5));
t422 = qJ(1) + qJ(2);
t414 = sin(t422);
t416 = cos(t422);
t426 = sin(qJ(3));
t424 = cos(pkin(5));
t428 = cos(qJ(3));
t446 = t424 * t428;
t371 = -t414 * t426 + t416 * t446;
t447 = t424 * t426;
t372 = t414 * t428 + t416 * t447;
t373 = -t414 * t446 - t416 * t426;
t374 = -t414 * t447 + t416 * t428;
t448 = t416 * t423;
t449 = t414 * t423;
t435 = (Icges(4,5) * t372 + Icges(4,6) * t371 - Icges(4,3) * t448) * t416 - (Icges(4,5) * t374 + Icges(4,6) * t373 + Icges(4,3) * t449) * t414;
t453 = t435 * t423;
t430 = pkin(8) + pkin(9);
t451 = t428 * pkin(3);
t450 = pkin(1) * qJD(1);
t421 = qJ(3) + qJ(4);
t386 = pkin(3) * t447 - t423 * t430;
t436 = pkin(8) * t423 + t386;
t340 = t451 * t414 + t436 * t416;
t341 = -t436 * t414 + t451 * t416;
t441 = qJD(3) * t423;
t390 = t414 * t441;
t437 = t416 * t441;
t445 = t340 * t390 + t341 * t437;
t406 = cos(qJ(4)) * pkin(4) + pkin(3);
t420 = pkin(10) + t430;
t439 = pkin(4) * sin(qJ(4)) * t428;
t444 = -t420 * t423 + (t406 * t426 + t439) * t424 - t386;
t429 = cos(qJ(1));
t408 = t429 * t450;
t419 = qJD(1) + qJD(2);
t443 = t419 * (pkin(2) * t416 + pkin(8) * t449) + t408;
t415 = cos(t421);
t442 = pkin(4) * t415;
t369 = qJD(4) * t449 + t390;
t440 = -qJD(3) - qJD(4);
t427 = sin(qJ(1));
t438 = t427 * t450;
t412 = pkin(5) - t421;
t411 = pkin(5) + t421;
t395 = qJD(3) * t424 + t419;
t387 = qJD(4) * t424 + t395;
t434 = -(pkin(2) * t414 - pkin(8) * t448) * t419 - t438;
t376 = pkin(3) * t423 * t426 + (-pkin(8) + t430) * t424;
t433 = t395 * t341 - t376 * t390 + t443;
t432 = -t340 * t395 - t376 * t437 + t434;
t417 = qJ(5) + t421;
t413 = sin(t421);
t405 = cos(t417);
t404 = sin(t417);
t403 = -qJ(5) + t412;
t402 = qJ(5) + t411;
t401 = cos(t411);
t400 = sin(t412);
t399 = cos(t402);
t398 = sin(t403);
t397 = cos(t412) / 0.2e1;
t396 = sin(t411) / 0.2e1;
t394 = cos(t403) / 0.2e1;
t393 = sin(t402) / 0.2e1;
t392 = rSges(2,1) * t429 - rSges(2,2) * t427;
t391 = rSges(2,1) * t427 + rSges(2,2) * t429;
t385 = t397 - t401 / 0.2e1;
t384 = t397 + t401 / 0.2e1;
t383 = t396 - t400 / 0.2e1;
t382 = t396 + t400 / 0.2e1;
t380 = t394 - t399 / 0.2e1;
t379 = t394 + t399 / 0.2e1;
t378 = t393 - t398 / 0.2e1;
t377 = t393 + t398 / 0.2e1;
t375 = qJD(5) * t424 + t387;
t370 = t440 * t448;
t367 = rSges(4,3) * t424 + (rSges(4,1) * t426 + rSges(4,2) * t428) * t423;
t366 = Icges(4,5) * t424 + (Icges(4,1) * t426 + Icges(4,4) * t428) * t423;
t365 = Icges(4,6) * t424 + (Icges(4,4) * t426 + Icges(4,2) * t428) * t423;
t364 = Icges(4,3) * t424 + (Icges(4,5) * t426 + Icges(4,6) * t428) * t423;
t363 = t408 + t419 * (rSges(3,1) * t416 - rSges(3,2) * t414);
t362 = -t438 - t419 * (rSges(3,1) * t414 + rSges(3,2) * t416);
t361 = (-qJD(5) + t440) * t448;
t360 = qJD(5) * t449 + t369;
t358 = -t383 * t414 + t415 * t416;
t357 = -t384 * t414 - t413 * t416;
t356 = t383 * t416 + t414 * t415;
t355 = t384 * t416 - t413 * t414;
t354 = -t378 * t414 + t405 * t416;
t353 = -t379 * t414 - t404 * t416;
t352 = t378 * t416 + t405 * t414;
t351 = t379 * t416 - t404 * t414;
t350 = rSges(5,1) * t385 + rSges(5,2) * t382 + rSges(5,3) * t424;
t349 = Icges(5,1) * t385 + Icges(5,4) * t382 + Icges(5,5) * t424;
t348 = Icges(5,4) * t385 + Icges(5,2) * t382 + Icges(5,6) * t424;
t347 = Icges(5,5) * t385 + Icges(5,6) * t382 + Icges(5,3) * t424;
t346 = rSges(6,1) * t380 + rSges(6,2) * t377 + rSges(6,3) * t424;
t345 = (t420 - t430) * t424 + (t439 + (-pkin(3) + t406) * t426) * t423;
t344 = Icges(6,1) * t380 + Icges(6,4) * t377 + Icges(6,5) * t424;
t343 = Icges(6,4) * t380 + Icges(6,2) * t377 + Icges(6,6) * t424;
t342 = Icges(6,5) * t380 + Icges(6,6) * t377 + Icges(6,3) * t424;
t339 = rSges(4,1) * t374 + rSges(4,2) * t373 + rSges(4,3) * t449;
t338 = rSges(4,1) * t372 + rSges(4,2) * t371 - rSges(4,3) * t448;
t337 = Icges(4,1) * t374 + Icges(4,4) * t373 + Icges(4,5) * t449;
t336 = Icges(4,1) * t372 + Icges(4,4) * t371 - Icges(4,5) * t448;
t335 = Icges(4,4) * t374 + Icges(4,2) * t373 + Icges(4,6) * t449;
t334 = Icges(4,4) * t372 + Icges(4,2) * t371 - Icges(4,6) * t448;
t328 = rSges(5,1) * t358 + rSges(5,2) * t357 + rSges(5,3) * t449;
t327 = rSges(5,1) * t356 + rSges(5,2) * t355 - rSges(5,3) * t448;
t326 = Icges(5,1) * t358 + Icges(5,4) * t357 + Icges(5,5) * t449;
t325 = Icges(5,1) * t356 + Icges(5,4) * t355 - Icges(5,5) * t448;
t324 = Icges(5,4) * t358 + Icges(5,2) * t357 + Icges(5,6) * t449;
t323 = Icges(5,4) * t356 + Icges(5,2) * t355 - Icges(5,6) * t448;
t322 = Icges(5,5) * t358 + Icges(5,6) * t357 + Icges(5,3) * t449;
t321 = Icges(5,5) * t356 + Icges(5,6) * t355 - Icges(5,3) * t448;
t320 = rSges(6,1) * t354 + rSges(6,2) * t353 + rSges(6,3) * t449;
t319 = rSges(6,1) * t352 + rSges(6,2) * t351 - rSges(6,3) * t448;
t318 = Icges(6,1) * t354 + Icges(6,4) * t353 + Icges(6,5) * t449;
t317 = Icges(6,1) * t352 + Icges(6,4) * t351 - Icges(6,5) * t448;
t316 = Icges(6,4) * t354 + Icges(6,2) * t353 + Icges(6,6) * t449;
t315 = Icges(6,4) * t352 + Icges(6,2) * t351 - Icges(6,6) * t448;
t314 = Icges(6,5) * t354 + Icges(6,6) * t353 + Icges(6,3) * t449;
t313 = Icges(6,5) * t352 + Icges(6,6) * t351 - Icges(6,3) * t448;
t312 = -t444 * t414 + t442 * t416;
t311 = t442 * t414 + t444 * t416;
t310 = t339 * t395 - t367 * t390 + t443;
t309 = -t338 * t395 - t367 * t437 + t434;
t308 = (t338 * t414 + t339 * t416) * t441;
t307 = t328 * t387 - t350 * t369 + t433;
t306 = -t327 * t387 + t350 * t370 + t432;
t305 = t327 * t369 - t328 * t370 + t445;
t304 = t312 * t387 + t320 * t375 - t345 * t369 - t346 * t360 + t433;
t303 = -t311 * t387 - t319 * t375 + t345 * t370 + t346 * t361 + t432;
t302 = t311 * t369 - t312 * t370 + t319 * t360 - t320 * t361 + t445;
t1 = m(3) * (t362 ^ 2 + t363 ^ 2) / 0.2e1 + t419 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t308 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + ((t364 * t449 + t365 * t373 + t366 * t374) * t395 + (-(t334 * t373 + t336 * t374) * t416 + (t373 * t335 + t374 * t337 - t453) * t414) * t441) * t390 / 0.2e1 - ((-t364 * t448 + t365 * t371 + t366 * t372) * t395 + ((t335 * t371 + t337 * t372) * t414 + (-t371 * t334 - t372 * t336 + t453) * t416) * t441) * t437 / 0.2e1 + t395 * ((t424 * t364 + (t365 * t428 + t366 * t426) * t423) * t395 + (((t335 * t428 + t337 * t426) * t414 - (t334 * t428 + t336 * t426) * t416) * t423 - t435 * t424) * t441) / 0.2e1 + m(5) * (t305 ^ 2 + t306 ^ 2 + t307 ^ 2) / 0.2e1 + t369 * ((t322 * t449 + t357 * t324 + t358 * t326) * t369 + (t321 * t449 + t323 * t357 + t325 * t358) * t370 + (t347 * t449 + t348 * t357 + t349 * t358) * t387) / 0.2e1 + t370 * ((-t322 * t448 + t324 * t355 + t326 * t356) * t369 + (-t321 * t448 + t323 * t355 + t325 * t356) * t370 + (-t347 * t448 + t348 * t355 + t349 * t356) * t387) / 0.2e1 + t387 * ((t322 * t424 + t324 * t382 + t326 * t385) * t369 + (t321 * t424 + t323 * t382 + t325 * t385) * t370 + (t347 * t424 + t348 * t382 + t349 * t385) * t387) / 0.2e1 + m(6) * (t302 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + t360 * ((t314 * t449 + t353 * t316 + t354 * t318) * t360 + (t313 * t449 + t315 * t353 + t317 * t354) * t361 + (t342 * t449 + t343 * t353 + t344 * t354) * t375) / 0.2e1 + t361 * ((-t314 * t448 + t316 * t351 + t318 * t352) * t360 + (-t313 * t448 + t351 * t315 + t352 * t317) * t361 + (-t342 * t448 + t343 * t351 + t344 * t352) * t375) / 0.2e1 + t375 * ((t314 * t424 + t316 * t377 + t318 * t380) * t360 + (t313 * t424 + t315 * t377 + t317 * t380) * t361 + (t342 * t424 + t343 * t377 + t344 * t380) * t375) / 0.2e1 + (m(2) * (t391 ^ 2 + t392 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
