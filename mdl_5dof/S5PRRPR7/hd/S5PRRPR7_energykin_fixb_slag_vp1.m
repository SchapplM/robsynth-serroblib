% Calculate kinetic energy for
% S5PRRPR7
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:34:58
% EndTime: 2019-12-05 16:35:00
% DurationCPUTime: 2.28s
% Computational Cost: add. (2008->289), mult. (5177->443), div. (0->0), fcn. (6474->12), ass. (0->126)
t459 = Icges(4,2) + Icges(5,3);
t425 = sin(pkin(9));
t427 = cos(pkin(9));
t433 = cos(qJ(2));
t428 = cos(pkin(5));
t431 = sin(qJ(2));
t445 = t428 * t431;
t412 = t425 * t433 + t427 * t445;
t430 = sin(qJ(3));
t426 = sin(pkin(5));
t446 = t427 * t426;
t451 = cos(qJ(3));
t397 = t412 * t451 - t430 * t446;
t444 = t428 * t433;
t411 = t425 * t431 - t427 * t444;
t424 = sin(pkin(10));
t450 = cos(pkin(10));
t367 = t397 * t424 - t411 * t450;
t368 = t397 * t450 + t411 * t424;
t441 = t426 * t451;
t396 = t412 * t430 + t427 * t441;
t458 = -Icges(4,4) * t397 + Icges(5,5) * t368 - Icges(4,6) * t411 - Icges(5,6) * t367 + t459 * t396;
t414 = -t425 * t445 + t427 * t433;
t448 = t426 * t430;
t399 = t414 * t451 + t425 * t448;
t413 = t425 * t444 + t427 * t431;
t369 = t399 * t424 - t413 * t450;
t370 = t399 * t450 + t413 * t424;
t398 = t414 * t430 - t425 * t441;
t457 = -Icges(4,4) * t399 + Icges(5,5) * t370 - Icges(4,6) * t413 - Icges(5,6) * t369 + t459 * t398;
t416 = t428 * t430 + t431 * t441;
t447 = t426 * t433;
t394 = t416 * t424 + t447 * t450;
t395 = t416 * t450 - t424 * t447;
t415 = -t428 * t451 + t431 * t448;
t456 = -Icges(4,4) * t416 + Icges(5,5) * t395 + Icges(4,6) * t447 - Icges(5,6) * t394 + t459 * t415;
t455 = qJD(2) ^ 2;
t449 = t425 * t426;
t443 = qJD(2) * t426;
t421 = t425 * t443;
t400 = qJD(3) * t413 + t421;
t423 = qJD(2) * t428;
t440 = t427 * t443;
t389 = pkin(2) * t412 + pkin(7) * t411;
t390 = pkin(2) * t414 + pkin(7) * t413;
t439 = t389 * t421 + t390 * t440 + qJD(1);
t401 = qJD(3) * t411 - t440;
t418 = -qJD(3) * t447 + t423;
t364 = pkin(3) * t397 + qJ(4) * t396;
t438 = qJD(4) * t415 + t400 * t364 + t439;
t417 = (pkin(2) * t431 - pkin(7) * t433) * t426;
t437 = t390 * t423 - t417 * t421;
t365 = pkin(3) * t399 + qJ(4) * t398;
t436 = qJD(4) * t396 + t418 * t365 + t437;
t435 = (-t389 * t428 - t417 * t446) * qJD(2);
t391 = pkin(3) * t416 + qJ(4) * t415;
t434 = qJD(4) * t398 + t401 * t391 + t435;
t432 = cos(qJ(5));
t429 = sin(qJ(5));
t405 = t428 * rSges(3,3) + (rSges(3,1) * t431 + rSges(3,2) * t433) * t426;
t404 = Icges(3,5) * t428 + (Icges(3,1) * t431 + Icges(3,4) * t433) * t426;
t403 = Icges(3,6) * t428 + (Icges(3,4) * t431 + Icges(3,2) * t433) * t426;
t402 = Icges(3,3) * t428 + (Icges(3,5) * t431 + Icges(3,6) * t433) * t426;
t387 = t416 * rSges(4,1) - t415 * rSges(4,2) - rSges(4,3) * t447;
t386 = Icges(4,1) * t416 - Icges(4,4) * t415 - Icges(4,5) * t447;
t384 = Icges(4,5) * t416 - Icges(4,6) * t415 - Icges(4,3) * t447;
t381 = rSges(3,1) * t414 - rSges(3,2) * t413 + rSges(3,3) * t449;
t380 = rSges(3,1) * t412 - rSges(3,2) * t411 - rSges(3,3) * t446;
t379 = Icges(3,1) * t414 - Icges(3,4) * t413 + Icges(3,5) * t449;
t378 = Icges(3,1) * t412 - Icges(3,4) * t411 - Icges(3,5) * t446;
t377 = Icges(3,4) * t414 - Icges(3,2) * t413 + Icges(3,6) * t449;
t376 = Icges(3,4) * t412 - Icges(3,2) * t411 - Icges(3,6) * t446;
t375 = Icges(3,5) * t414 - Icges(3,6) * t413 + Icges(3,3) * t449;
t374 = Icges(3,5) * t412 - Icges(3,6) * t411 - Icges(3,3) * t446;
t373 = qJD(5) * t394 + t418;
t372 = t395 * t432 + t415 * t429;
t371 = -t395 * t429 + t415 * t432;
t366 = pkin(4) * t395 + pkin(8) * t394;
t361 = (-t380 * t428 - t405 * t446) * qJD(2);
t360 = (t381 * t428 - t405 * t449) * qJD(2);
t359 = rSges(5,1) * t395 - rSges(5,2) * t394 + rSges(5,3) * t415;
t358 = rSges(4,1) * t399 - rSges(4,2) * t398 + rSges(4,3) * t413;
t357 = rSges(4,1) * t397 - rSges(4,2) * t396 + rSges(4,3) * t411;
t356 = Icges(5,1) * t395 - Icges(5,4) * t394 + Icges(5,5) * t415;
t355 = Icges(5,4) * t395 - Icges(5,2) * t394 + Icges(5,6) * t415;
t353 = Icges(4,1) * t399 - Icges(4,4) * t398 + Icges(4,5) * t413;
t352 = Icges(4,1) * t397 - Icges(4,4) * t396 + Icges(4,5) * t411;
t349 = Icges(4,5) * t399 - Icges(4,6) * t398 + Icges(4,3) * t413;
t348 = Icges(4,5) * t397 - Icges(4,6) * t396 + Icges(4,3) * t411;
t347 = qJD(5) * t367 + t401;
t346 = qJD(5) * t369 + t400;
t345 = t370 * t432 + t398 * t429;
t344 = -t370 * t429 + t398 * t432;
t343 = t368 * t432 + t396 * t429;
t342 = -t368 * t429 + t396 * t432;
t340 = qJD(1) + (t380 * t425 + t381 * t427) * t443;
t339 = pkin(4) * t370 + pkin(8) * t369;
t338 = pkin(4) * t368 + pkin(8) * t367;
t337 = rSges(6,1) * t372 + rSges(6,2) * t371 + rSges(6,3) * t394;
t336 = Icges(6,1) * t372 + Icges(6,4) * t371 + Icges(6,5) * t394;
t335 = Icges(6,4) * t372 + Icges(6,2) * t371 + Icges(6,6) * t394;
t334 = Icges(6,5) * t372 + Icges(6,6) * t371 + Icges(6,3) * t394;
t333 = rSges(5,1) * t370 - rSges(5,2) * t369 + rSges(5,3) * t398;
t332 = rSges(5,1) * t368 - rSges(5,2) * t367 + rSges(5,3) * t396;
t331 = Icges(5,1) * t370 - Icges(5,4) * t369 + Icges(5,5) * t398;
t330 = Icges(5,1) * t368 - Icges(5,4) * t367 + Icges(5,5) * t396;
t329 = Icges(5,4) * t370 - Icges(5,2) * t369 + Icges(5,6) * t398;
t328 = Icges(5,4) * t368 - Icges(5,2) * t367 + Icges(5,6) * t396;
t325 = -t357 * t418 + t387 * t401 + t435;
t324 = t358 * t418 - t387 * t400 + t437;
t323 = rSges(6,1) * t345 + rSges(6,2) * t344 + rSges(6,3) * t369;
t322 = rSges(6,1) * t343 + rSges(6,2) * t342 + rSges(6,3) * t367;
t321 = Icges(6,1) * t345 + Icges(6,4) * t344 + Icges(6,5) * t369;
t320 = Icges(6,1) * t343 + Icges(6,4) * t342 + Icges(6,5) * t367;
t319 = Icges(6,4) * t345 + Icges(6,2) * t344 + Icges(6,6) * t369;
t318 = Icges(6,4) * t343 + Icges(6,2) * t342 + Icges(6,6) * t367;
t317 = Icges(6,5) * t345 + Icges(6,6) * t344 + Icges(6,3) * t369;
t316 = Icges(6,5) * t343 + Icges(6,6) * t342 + Icges(6,3) * t367;
t315 = t357 * t400 - t358 * t401 + t439;
t314 = t359 * t401 + (-t332 - t364) * t418 + t434;
t313 = t333 * t418 + (-t359 - t391) * t400 + t436;
t312 = t332 * t400 + (-t333 - t365) * t401 + t438;
t311 = -t322 * t373 + t337 * t347 + t366 * t401 + (-t338 - t364) * t418 + t434;
t310 = t323 * t373 - t337 * t346 + t339 * t418 + (-t366 - t391) * t400 + t436;
t309 = t322 * t346 - t323 * t347 + t338 * t400 + (-t339 - t365) * t401 + t438;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t340 ^ 2 + t360 ^ 2 + t361 ^ 2) / 0.2e1 - t455 * ((-t375 * t446 - t377 * t411 + t379 * t412) * t449 - (-t374 * t446 - t376 * t411 + t378 * t412) * t446 + (-t402 * t446 - t403 * t411 + t404 * t412) * t428) * t446 / 0.2e1 + m(4) * (t315 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + m(5) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + m(6) * (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + t346 * ((t369 * t317 + t344 * t319 + t345 * t321) * t346 + (t316 * t369 + t318 * t344 + t320 * t345) * t347 + (t334 * t369 + t335 * t344 + t336 * t345) * t373) / 0.2e1 + t347 * ((t317 * t367 + t319 * t342 + t321 * t343) * t346 + (t367 * t316 + t342 * t318 + t343 * t320) * t347 + (t334 * t367 + t335 * t342 + t336 * t343) * t373) / 0.2e1 + t373 * ((t317 * t394 + t319 * t371 + t321 * t372) * t346 + (t316 * t394 + t318 * t371 + t320 * t372) * t347 + (t394 * t334 + t371 * t335 + t372 * t336) * t373) / 0.2e1 + ((-t355 * t369 + t356 * t370 + t384 * t413 + t386 * t399 + t398 * t456) * t418 + (-t328 * t369 + t330 * t370 + t348 * t413 + t352 * t399 + t398 * t458) * t401 + (-t369 * t329 + t370 * t331 + t413 * t349 + t399 * t353 + t398 * t457) * t400) * t400 / 0.2e1 + ((-t355 * t367 + t356 * t368 + t384 * t411 + t386 * t397 + t396 * t456) * t418 + (-t367 * t328 + t368 * t330 + t411 * t348 + t397 * t352 + t396 * t458) * t401 + (-t329 * t367 + t331 * t368 + t349 * t411 + t353 * t397 + t396 * t457) * t400) * t401 / 0.2e1 + ((-t394 * t355 + t395 * t356 - t384 * t447 + t416 * t386 + t415 * t456) * t418 + (-t328 * t394 + t330 * t395 - t348 * t447 + t416 * t352 + t415 * t458) * t401 + (-t329 * t394 + t331 * t395 - t349 * t447 + t416 * t353 + t457 * t415) * t400) * t418 / 0.2e1 + (((t375 * t449 - t377 * t413 + t379 * t414) * t449 - (t374 * t449 - t376 * t413 + t378 * t414) * t446 + (t402 * t449 - t403 * t413 + t404 * t414) * t428) * t449 + t428 * (t428 ^ 2 * t402 + (((t377 * t433 + t379 * t431) * t425 - (t376 * t433 + t378 * t431) * t427) * t426 + (-t374 * t427 + t375 * t425 + t403 * t433 + t404 * t431) * t428) * t426)) * t455 / 0.2e1;
T = t1;
