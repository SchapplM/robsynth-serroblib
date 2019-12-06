% Calculate kinetic energy for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:11:57
% EndTime: 2019-12-05 16:11:59
% DurationCPUTime: 1.85s
% Computational Cost: add. (901->184), mult. (2280->275), div. (0->0), fcn. (2585->8), ass. (0->99)
t436 = Icges(5,1) + Icges(6,1);
t435 = Icges(5,4) - Icges(6,5);
t434 = Icges(6,4) + Icges(5,5);
t433 = Icges(5,2) + Icges(6,3);
t432 = Icges(6,2) + Icges(5,3);
t431 = Icges(5,6) - Icges(6,6);
t374 = cos(pkin(7));
t412 = t374 ^ 2;
t373 = sin(pkin(7));
t413 = t373 ^ 2;
t417 = t412 + t413;
t430 = qJD(2) * t417;
t429 = rSges(6,1) + pkin(4);
t428 = rSges(6,3) + qJ(5);
t427 = t373 * t374;
t375 = sin(qJ(3));
t377 = cos(qJ(3));
t378 = cos(qJ(2));
t403 = t377 * t378;
t357 = t373 * t403 - t374 * t375;
t372 = sin(pkin(8));
t376 = sin(qJ(2));
t408 = cos(pkin(8));
t390 = t376 * t408;
t337 = t357 * t372 - t373 * t390;
t407 = t373 * t376;
t338 = t357 * t408 + t372 * t407;
t404 = t375 * t378;
t356 = t373 * t404 + t374 * t377;
t426 = t433 * t337 - t435 * t338 - t431 * t356;
t361 = t373 * t375 + t374 * t403;
t339 = t361 * t372 - t374 * t390;
t406 = t374 * t376;
t340 = t361 * t408 + t372 * t406;
t360 = -t373 * t377 + t374 * t404;
t425 = t433 * t339 - t435 * t340 - t431 * t360;
t424 = -t431 * t337 + t434 * t338 + t432 * t356;
t423 = -t431 * t339 + t434 * t340 + t432 * t360;
t422 = -t435 * t337 + t436 * t338 + t434 * t356;
t421 = -t435 * t339 + t436 * t340 + t434 * t360;
t358 = t376 * t377 * t372 + t378 * t408;
t359 = -t378 * t372 + t377 * t390;
t405 = t375 * t376;
t420 = -t433 * t358 + t435 * t359 + t431 * t405;
t419 = t431 * t358 - t434 * t359 - t432 * t405;
t418 = t435 * t358 - t436 * t359 - t434 * t405;
t319 = Icges(4,4) * t357 - Icges(4,2) * t356 + Icges(4,6) * t407;
t416 = -t319 + t424;
t320 = Icges(4,4) * t361 - Icges(4,2) * t360 + Icges(4,6) * t406;
t415 = -t320 + t423;
t351 = -Icges(4,6) * t378 + (Icges(4,4) * t377 - Icges(4,2) * t375) * t376;
t414 = t351 + t419;
t411 = qJD(2) ^ 2;
t402 = rSges(6,2) * t356 + t428 * t337 + t429 * t338;
t336 = pkin(3) * t361 + qJ(4) * t360;
t401 = -rSges(5,1) * t340 + rSges(5,2) * t339 - rSges(5,3) * t360 - t336;
t400 = rSges(6,2) * t405 + t428 * t358 + t429 * t359;
t399 = qJD(2) * t373;
t398 = qJD(2) * t374;
t397 = qJD(3) * t376;
t396 = qJD(3) * t378;
t395 = -rSges(6,2) * t360 - t428 * t339 - t429 * t340 - t336;
t369 = pkin(2) * t376 - pkin(6) * t378;
t394 = t369 * t399;
t393 = t369 * t398;
t392 = qJD(1) + (pkin(2) * t378 + pkin(6) * t376) * t430;
t389 = qJD(4) * t356 - t394;
t334 = pkin(3) * t357 + qJ(4) * t356;
t363 = t374 * t397 + t399;
t388 = qJD(4) * t405 + t363 * t334 + t392;
t385 = Icges(3,5) * t378 - Icges(3,6) * t376;
t362 = (pkin(3) * t377 + qJ(4) * t375) * t376;
t364 = t373 * t397 - t398;
t380 = qJD(4) * t360 + t334 * t396 + t364 * t362 - t393;
t366 = rSges(3,1) * t376 + rSges(3,2) * t378;
t353 = -rSges(4,3) * t378 + (rSges(4,1) * t377 - rSges(4,2) * t375) * t376;
t352 = -Icges(4,5) * t378 + (Icges(4,1) * t377 - Icges(4,4) * t375) * t376;
t350 = -Icges(4,3) * t378 + (Icges(4,5) * t377 - Icges(4,6) * t375) * t376;
t343 = Icges(3,3) * t373 + t374 * t385;
t342 = -Icges(3,3) * t374 + t373 * t385;
t332 = rSges(5,1) * t359 - rSges(5,2) * t358 + rSges(5,3) * t405;
t330 = rSges(4,1) * t361 - rSges(4,2) * t360 + rSges(4,3) * t406;
t329 = rSges(4,1) * t357 - rSges(4,2) * t356 + rSges(4,3) * t407;
t322 = Icges(4,1) * t361 - Icges(4,4) * t360 + Icges(4,5) * t406;
t321 = Icges(4,1) * t357 - Icges(4,4) * t356 + Icges(4,5) * t407;
t318 = Icges(4,5) * t361 - Icges(4,6) * t360 + Icges(4,3) * t406;
t317 = Icges(4,5) * t357 - Icges(4,6) * t356 + Icges(4,3) * t407;
t315 = qJD(1) + (rSges(3,1) * t378 - rSges(3,2) * t376) * t430;
t310 = rSges(5,1) * t338 - rSges(5,2) * t337 + rSges(5,3) * t356;
t296 = t329 * t396 + t353 * t364 - t393;
t295 = -t330 * t396 - t353 * t363 - t394;
t294 = t329 * t363 - t330 * t364 + t392;
t293 = t310 * t396 + t332 * t364 + t380;
t292 = (-t332 - t362) * t363 + t401 * t396 + t389;
t291 = t310 * t363 + t364 * t401 + t388;
t290 = qJD(5) * t339 + t364 * t400 + t396 * t402 + t380;
t289 = qJD(5) * t337 + (-t362 - t400) * t363 + t395 * t396 + t389;
t288 = qJD(5) * t358 + t363 * t402 + t364 * t395 + t388;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t417 * t411 * t366 ^ 2 + t315 ^ 2) / 0.2e1 + t411 * t373 * (-t342 * t427 + t413 * t343) / 0.2e1 - t411 * t374 * (t412 * t342 - t343 * t427) / 0.2e1 + m(4) * (t294 ^ 2 + t295 ^ 2 + t296 ^ 2) / 0.2e1 + m(5) * (t291 ^ 2 + t292 ^ 2 + t293 ^ 2) / 0.2e1 + m(6) * (t288 ^ 2 + t289 ^ 2 + t290 ^ 2) / 0.2e1 + ((t420 * t339 + t418 * t340 - t350 * t406 - t352 * t361 + t414 * t360) * t396 + (t317 * t406 + t321 * t361 + t426 * t339 + t422 * t340 + t416 * t360) * t364 + (t318 * t406 + t361 * t322 + t425 * t339 + t421 * t340 + t415 * t360) * t363) * t363 / 0.2e1 + ((t420 * t337 + t418 * t338 - t350 * t407 - t352 * t357 + t414 * t356) * t396 + (t317 * t407 + t357 * t321 + t426 * t337 + t422 * t338 + t416 * t356) * t364 + (t318 * t407 + t322 * t357 + t425 * t337 + t421 * t338 + t415 * t356) * t363) * t364 / 0.2e1 - ((t350 * t378 - (-t351 * t375 + t352 * t377) * t376 + t419 * t405 + t418 * t359 + t420 * t358) * t396 + (-t317 * t378 + (-t319 * t375 + t321 * t377) * t376 + t424 * t405 + t422 * t359 + t426 * t358) * t364 + (-t318 * t378 + (-t320 * t375 + t322 * t377) * t376 + t423 * t405 + t421 * t359 + t425 * t358) * t363) * t396 / 0.2e1;
T = t1;
