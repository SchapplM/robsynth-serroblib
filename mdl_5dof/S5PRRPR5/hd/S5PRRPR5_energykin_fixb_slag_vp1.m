% Calculate kinetic energy for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:16
% EndTime: 2019-12-05 16:25:19
% DurationCPUTime: 2.24s
% Computational Cost: add. (2008->293), mult. (3845->449), div. (0->0), fcn. (4602->12), ass. (0->132)
t464 = Icges(4,3) + Icges(5,3);
t422 = sin(pkin(9));
t424 = cos(pkin(9));
t432 = cos(qJ(2));
t425 = cos(pkin(5));
t429 = sin(qJ(2));
t448 = t425 * t429;
t408 = t422 * t432 + t424 * t448;
t444 = qJ(3) + pkin(10);
t421 = sin(t444);
t423 = sin(pkin(5));
t439 = cos(t444);
t438 = t423 * t439;
t385 = t408 * t421 + t424 * t438;
t450 = t424 * t423;
t386 = t408 * t439 - t421 * t450;
t428 = sin(qJ(3));
t431 = cos(qJ(3));
t391 = -t408 * t428 - t431 * t450;
t442 = t428 * t450;
t392 = t408 * t431 - t442;
t447 = t425 * t432;
t407 = t422 * t429 - t424 * t447;
t463 = Icges(4,5) * t392 + Icges(5,5) * t386 + Icges(4,6) * t391 - Icges(5,6) * t385 + t464 * t407;
t410 = -t422 * t448 + t424 * t432;
t387 = t410 * t421 - t422 * t438;
t454 = t422 * t423;
t388 = t410 * t439 + t421 * t454;
t452 = t423 * t431;
t393 = -t410 * t428 + t422 * t452;
t443 = t428 * t454;
t394 = t410 * t431 + t443;
t409 = t422 * t447 + t424 * t429;
t462 = Icges(4,5) * t394 + Icges(5,5) * t388 + Icges(4,6) * t393 - Icges(5,6) * t387 + t464 * t409;
t453 = t423 * t429;
t400 = t421 * t453 - t425 * t439;
t401 = t425 * t421 + t429 * t438;
t411 = t425 * t431 - t428 * t453;
t449 = t425 * t428;
t412 = t429 * t452 + t449;
t451 = t423 * t432;
t461 = Icges(4,5) * t412 + Icges(5,5) * t401 + Icges(4,6) * t411 - Icges(5,6) * t400 - t464 * t451;
t460 = qJD(2) ^ 2;
t456 = t431 * pkin(3);
t446 = qJD(2) * t423;
t418 = t422 * t446;
t395 = qJD(3) * t409 + t418;
t420 = qJD(2) * t425;
t441 = t424 * t446;
t382 = t408 * pkin(2) + t407 * pkin(7);
t383 = t410 * pkin(2) + t409 * pkin(7);
t440 = t382 * t418 + t383 * t441 + qJD(1);
t396 = qJD(3) * t407 - t441;
t414 = -qJD(3) * t451 + t420;
t413 = (pkin(2) * t429 - pkin(7) * t432) * t423;
t437 = t383 * t420 - t413 * t418;
t341 = pkin(3) * t443 + qJ(4) * t409 + t410 * t456;
t436 = qJD(4) * t407 + t414 * t341 + t437;
t435 = (-t382 * t425 - t413 * t450) * qJD(2);
t340 = -pkin(3) * t442 + qJ(4) * t407 + t408 * t456;
t434 = -qJD(4) * t451 + t395 * t340 + t440;
t379 = pkin(3) * t449 + (-qJ(4) * t432 + t429 * t456) * t423;
t433 = qJD(4) * t409 + t396 * t379 + t435;
t430 = cos(qJ(5));
t427 = sin(qJ(5));
t402 = t425 * rSges(3,3) + (rSges(3,1) * t429 + rSges(3,2) * t432) * t423;
t399 = Icges(3,5) * t425 + (Icges(3,1) * t429 + Icges(3,4) * t432) * t423;
t398 = Icges(3,6) * t425 + (Icges(3,4) * t429 + Icges(3,2) * t432) * t423;
t397 = Icges(3,3) * t425 + (Icges(3,5) * t429 + Icges(3,6) * t432) * t423;
t390 = t401 * t430 - t427 * t451;
t389 = -t401 * t427 - t430 * t451;
t384 = qJD(5) * t400 + t414;
t380 = t412 * rSges(4,1) + t411 * rSges(4,2) - rSges(4,3) * t451;
t378 = Icges(4,1) * t412 + Icges(4,4) * t411 - Icges(4,5) * t451;
t377 = Icges(4,4) * t412 + Icges(4,2) * t411 - Icges(4,6) * t451;
t373 = t410 * rSges(3,1) - t409 * rSges(3,2) + rSges(3,3) * t454;
t372 = t408 * rSges(3,1) - t407 * rSges(3,2) - rSges(3,3) * t450;
t371 = t401 * pkin(4) + t400 * pkin(8);
t370 = Icges(3,1) * t410 - Icges(3,4) * t409 + Icges(3,5) * t454;
t369 = Icges(3,1) * t408 - Icges(3,4) * t407 - Icges(3,5) * t450;
t368 = Icges(3,4) * t410 - Icges(3,2) * t409 + Icges(3,6) * t454;
t367 = Icges(3,4) * t408 - Icges(3,2) * t407 - Icges(3,6) * t450;
t366 = Icges(3,5) * t410 - Icges(3,6) * t409 + Icges(3,3) * t454;
t365 = Icges(3,5) * t408 - Icges(3,6) * t407 - Icges(3,3) * t450;
t364 = t401 * rSges(5,1) - t400 * rSges(5,2) - rSges(5,3) * t451;
t363 = Icges(5,1) * t401 - Icges(5,4) * t400 - Icges(5,5) * t451;
t362 = Icges(5,4) * t401 - Icges(5,2) * t400 - Icges(5,6) * t451;
t360 = t388 * t430 + t409 * t427;
t359 = -t388 * t427 + t409 * t430;
t358 = t386 * t430 + t407 * t427;
t357 = -t386 * t427 + t407 * t430;
t356 = qJD(5) * t385 + t396;
t355 = qJD(5) * t387 + t395;
t353 = t388 * pkin(4) + t387 * pkin(8);
t352 = t386 * pkin(4) + t385 * pkin(8);
t351 = (-t372 * t425 - t402 * t450) * qJD(2);
t350 = (t373 * t425 - t402 * t454) * qJD(2);
t349 = t394 * rSges(4,1) + t393 * rSges(4,2) + t409 * rSges(4,3);
t348 = t392 * rSges(4,1) + t391 * rSges(4,2) + t407 * rSges(4,3);
t347 = Icges(4,1) * t394 + Icges(4,4) * t393 + Icges(4,5) * t409;
t346 = Icges(4,1) * t392 + Icges(4,4) * t391 + Icges(4,5) * t407;
t345 = Icges(4,4) * t394 + Icges(4,2) * t393 + Icges(4,6) * t409;
t344 = Icges(4,4) * t392 + Icges(4,2) * t391 + Icges(4,6) * t407;
t339 = t388 * rSges(5,1) - t387 * rSges(5,2) + t409 * rSges(5,3);
t338 = t386 * rSges(5,1) - t385 * rSges(5,2) + t407 * rSges(5,3);
t337 = Icges(5,1) * t388 - Icges(5,4) * t387 + Icges(5,5) * t409;
t336 = Icges(5,1) * t386 - Icges(5,4) * t385 + Icges(5,5) * t407;
t335 = Icges(5,4) * t388 - Icges(5,2) * t387 + Icges(5,6) * t409;
t334 = Icges(5,4) * t386 - Icges(5,2) * t385 + Icges(5,6) * t407;
t331 = t390 * rSges(6,1) + t389 * rSges(6,2) + t400 * rSges(6,3);
t330 = Icges(6,1) * t390 + Icges(6,4) * t389 + Icges(6,5) * t400;
t329 = Icges(6,4) * t390 + Icges(6,2) * t389 + Icges(6,6) * t400;
t328 = Icges(6,5) * t390 + Icges(6,6) * t389 + Icges(6,3) * t400;
t326 = qJD(1) + (t372 * t422 + t373 * t424) * t446;
t324 = t360 * rSges(6,1) + t359 * rSges(6,2) + t387 * rSges(6,3);
t323 = t358 * rSges(6,1) + t357 * rSges(6,2) + t385 * rSges(6,3);
t322 = Icges(6,1) * t360 + Icges(6,4) * t359 + Icges(6,5) * t387;
t321 = Icges(6,1) * t358 + Icges(6,4) * t357 + Icges(6,5) * t385;
t320 = Icges(6,4) * t360 + Icges(6,2) * t359 + Icges(6,6) * t387;
t319 = Icges(6,4) * t358 + Icges(6,2) * t357 + Icges(6,6) * t385;
t318 = Icges(6,5) * t360 + Icges(6,6) * t359 + Icges(6,3) * t387;
t317 = Icges(6,5) * t358 + Icges(6,6) * t357 + Icges(6,3) * t385;
t316 = -t414 * t348 + t396 * t380 + t435;
t315 = t414 * t349 - t395 * t380 + t437;
t314 = t395 * t348 - t396 * t349 + t440;
t313 = t396 * t364 + (-t338 - t340) * t414 + t433;
t312 = t414 * t339 + (-t364 - t379) * t395 + t436;
t311 = t395 * t338 + (-t339 - t341) * t396 + t434;
t310 = -t384 * t323 + t356 * t331 + t396 * t371 + (-t340 - t352) * t414 + t433;
t309 = t384 * t324 - t355 * t331 + t414 * t353 + (-t371 - t379) * t395 + t436;
t308 = t355 * t323 - t356 * t324 + t395 * t352 + (-t341 - t353) * t396 + t434;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t326 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 - t460 * ((-t366 * t450 - t407 * t368 + t408 * t370) * t454 - (-t365 * t450 - t407 * t367 + t408 * t369) * t450 + (-t397 * t450 - t407 * t398 + t408 * t399) * t425) * t450 / 0.2e1 + m(4) * (t314 ^ 2 + t315 ^ 2 + t316 ^ 2) / 0.2e1 + m(5) * (t311 ^ 2 + t312 ^ 2 + t313 ^ 2) / 0.2e1 + m(6) * (t308 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + t355 * ((t387 * t318 + t359 * t320 + t360 * t322) * t355 + (t387 * t317 + t359 * t319 + t360 * t321) * t356 + (t387 * t328 + t359 * t329 + t360 * t330) * t384) / 0.2e1 + t356 * ((t385 * t318 + t357 * t320 + t358 * t322) * t355 + (t385 * t317 + t357 * t319 + t358 * t321) * t356 + (t385 * t328 + t357 * t329 + t358 * t330) * t384) / 0.2e1 + t384 * ((t400 * t318 + t389 * t320 + t390 * t322) * t355 + (t400 * t317 + t389 * t319 + t390 * t321) * t356 + (t400 * t328 + t389 * t329 + t390 * t330) * t384) / 0.2e1 + ((-t387 * t362 + t388 * t363 + t393 * t377 + t394 * t378 + t461 * t409) * t414 + (-t387 * t334 + t388 * t336 + t393 * t344 + t394 * t346 + t463 * t409) * t396 + (-t387 * t335 + t388 * t337 + t393 * t345 + t394 * t347 + t462 * t409) * t395) * t395 / 0.2e1 + ((-t385 * t362 + t386 * t363 + t391 * t377 + t392 * t378 + t461 * t407) * t414 + (-t385 * t334 + t386 * t336 + t391 * t344 + t392 * t346 + t463 * t407) * t396 + (-t385 * t335 + t386 * t337 + t391 * t345 + t392 * t347 + t462 * t407) * t395) * t396 / 0.2e1 + ((-t400 * t362 + t401 * t363 + t411 * t377 + t412 * t378 - t461 * t451) * t414 + (-t400 * t334 + t401 * t336 + t411 * t344 + t412 * t346 - t463 * t451) * t396 + (-t400 * t335 + t401 * t337 + t411 * t345 + t412 * t347 - t462 * t451) * t395) * t414 / 0.2e1 + (((t366 * t454 - t409 * t368 + t410 * t370) * t454 - (t365 * t454 - t409 * t367 + t410 * t369) * t450 + (t397 * t454 - t409 * t398 + t410 * t399) * t425) * t454 + t425 * (t425 ^ 2 * t397 + (((t368 * t432 + t370 * t429) * t422 - (t367 * t432 + t369 * t429) * t424) * t423 + (-t365 * t424 + t366 * t422 + t398 * t432 + t399 * t429) * t425) * t423)) * t460 / 0.2e1;
T = t1;
