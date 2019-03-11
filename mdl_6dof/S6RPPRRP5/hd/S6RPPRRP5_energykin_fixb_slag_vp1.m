% Calculate kinetic energy for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:40
% EndTime: 2019-03-09 02:07:41
% DurationCPUTime: 1.32s
% Computational Cost: add. (565->184), mult. (1286->284), div. (0->0), fcn. (1249->6), ass. (0->97)
t424 = Icges(6,1) + Icges(7,1);
t423 = Icges(6,4) + Icges(7,4);
t422 = Icges(7,5) + Icges(6,5);
t421 = Icges(6,2) + Icges(7,2);
t420 = Icges(7,6) + Icges(6,6);
t419 = Icges(7,3) + Icges(6,3);
t366 = sin(qJ(4));
t368 = cos(qJ(5));
t370 = cos(qJ(1));
t396 = t368 * t370;
t365 = sin(qJ(5));
t367 = sin(qJ(1));
t400 = t365 * t367;
t339 = -t366 * t400 + t396;
t398 = t367 * t368;
t399 = t365 * t370;
t340 = t366 * t398 + t399;
t369 = cos(qJ(4));
t397 = t367 * t369;
t418 = t420 * t339 + t422 * t340 - t419 * t397;
t341 = -t366 * t399 - t398;
t342 = t366 * t396 - t400;
t395 = t369 * t370;
t417 = t420 * t341 + t422 * t342 - t419 * t395;
t416 = t421 * t339 + t423 * t340 - t420 * t397;
t415 = t421 * t341 + t423 * t342 - t420 * t395;
t414 = t423 * t339 + t424 * t340 - t422 * t397;
t413 = t423 * t341 + t424 * t342 - t422 * t395;
t412 = (-t420 * t365 + t422 * t368) * t369 + t419 * t366;
t411 = (-t421 * t365 + t423 * t368) * t369 + t420 * t366;
t410 = (-t423 * t365 + t424 * t368) * t369 + t422 * t366;
t405 = pkin(5) * t368;
t403 = pkin(7) * qJD(1);
t402 = Icges(5,4) * t366;
t401 = Icges(5,4) * t369;
t374 = -qJ(6) * t369 + t366 * t405;
t394 = rSges(7,1) * t340 + rSges(7,2) * t339 - rSges(7,3) * t397 + pkin(5) * t399 + t367 * t374;
t393 = rSges(7,1) * t342 + rSges(7,2) * t341 - rSges(7,3) * t395 - pkin(5) * t400 + t370 * t374;
t392 = (rSges(7,1) * t368 - rSges(7,2) * t365 + t405) * t369 + (qJ(6) + rSges(7,3)) * t366;
t363 = qJD(2) * t367;
t391 = qJD(3) * t370 + t363;
t390 = qJD(4) * t367;
t389 = qJD(4) * t370;
t388 = qJD(5) * t369;
t387 = qJD(6) * t369;
t386 = -qJD(2) * t370 + qJD(1) * (pkin(1) * t370 + qJ(2) * t367);
t385 = pkin(4) * t366 - pkin(8) * t369;
t384 = rSges(5,1) * t366 + rSges(5,2) * t369;
t383 = Icges(5,1) * t366 + t401;
t382 = Icges(5,2) * t369 + t402;
t381 = Icges(5,5) * t366 + Icges(5,6) * t369;
t328 = Icges(5,6) * t370 + t367 * t382;
t332 = Icges(5,5) * t370 + t367 * t383;
t380 = t328 * t369 + t332 * t366;
t329 = -Icges(5,6) * t367 + t370 * t382;
t333 = -Icges(5,5) * t367 + t370 * t383;
t379 = -t329 * t369 - t333 * t366;
t351 = -Icges(5,2) * t366 + t401;
t352 = Icges(5,1) * t369 - t402;
t378 = t351 * t369 + t352 * t366;
t377 = qJD(1) * t370 * qJ(3) + qJD(3) * t367 + t386;
t353 = pkin(1) * t367 - qJ(2) * t370;
t376 = -pkin(7) * t370 - qJ(3) * t367 - t353;
t343 = t385 * t367;
t344 = t385 * t370;
t375 = (-t343 * t367 - t344 * t370) * qJD(4);
t357 = pkin(4) * t369 + pkin(8) * t366;
t373 = qJD(1) * t344 + t357 * t390 + t377;
t372 = t357 * t389 + (-t343 + t376) * qJD(1) + t391;
t358 = qJD(5) * t366 + qJD(1);
t356 = rSges(2,1) * t370 - rSges(2,2) * t367;
t355 = rSges(5,1) * t369 - rSges(5,2) * t366;
t354 = rSges(2,1) * t367 + rSges(2,2) * t370;
t350 = Icges(5,5) * t369 - Icges(5,6) * t366;
t348 = -t367 * t388 + t389;
t347 = -t370 * t388 - t390;
t337 = -rSges(5,3) * t367 + t370 * t384;
t336 = rSges(6,3) * t366 + (rSges(6,1) * t368 - rSges(6,2) * t365) * t369;
t334 = rSges(5,3) * t370 + t367 * t384;
t325 = -Icges(5,3) * t367 + t370 * t381;
t324 = Icges(5,3) * t370 + t367 * t381;
t320 = qJD(1) * (-rSges(3,2) * t370 + rSges(3,3) * t367) + t386;
t319 = t363 + (rSges(3,2) * t367 + rSges(3,3) * t370 - t353) * qJD(1);
t318 = qJD(1) * (rSges(4,2) * t367 + rSges(4,3) * t370) + t377;
t317 = (t370 * rSges(4,2) - t353 + (-rSges(4,3) - qJ(3)) * t367) * qJD(1) + t391;
t316 = rSges(6,1) * t342 + rSges(6,2) * t341 - rSges(6,3) * t395;
t314 = rSges(6,1) * t340 + rSges(6,2) * t339 - rSges(6,3) * t397;
t298 = (-t334 * t367 - t337 * t370) * qJD(4);
t297 = t355 * t390 + (-pkin(7) * t367 + t337) * qJD(1) + t377;
t296 = t355 * t389 + (-t334 + t376) * qJD(1) + t391;
t295 = t316 * t358 - t336 * t347 - t367 * t403 + t373;
t294 = -t314 * t358 + t336 * t348 + t372;
t293 = t314 * t347 - t316 * t348 + t375;
t292 = (-t387 - t403) * t367 + t393 * t358 - t392 * t347 + t373;
t291 = t348 * t392 - t358 * t394 - t370 * t387 + t372;
t290 = qJD(6) * t366 + t347 * t394 - t348 * t393 + t375;
t1 = m(6) * (t293 ^ 2 + t294 ^ 2 + t295 ^ 2) / 0.2e1 + m(7) * (t290 ^ 2 + t291 ^ 2 + t292 ^ 2) / 0.2e1 + ((t370 * t350 + t367 * t378) * qJD(1) + (t370 ^ 2 * t324 + (t379 * t367 + (-t325 + t380) * t370) * t367) * qJD(4)) * t389 / 0.2e1 + qJD(1) * ((-t366 * t351 + t369 * t352) * qJD(1) + (-(-t329 * t366 + t333 * t369) * t367 + (-t328 * t366 + t332 * t369) * t370) * qJD(4)) / 0.2e1 + m(5) * (t296 ^ 2 + t297 ^ 2 + t298 ^ 2) / 0.2e1 + m(3) * (t319 ^ 2 + t320 ^ 2) / 0.2e1 + m(4) * (t317 ^ 2 + t318 ^ 2) / 0.2e1 - ((-t367 * t350 + t370 * t378) * qJD(1) + (t367 ^ 2 * t325 + (t380 * t370 + (-t324 + t379) * t367) * t370) * qJD(4)) * t390 / 0.2e1 + ((t411 * t341 + t410 * t342 - t412 * t395) * t358 + (t416 * t341 + t414 * t342 - t418 * t395) * t348 + (t415 * t341 + t413 * t342 - t417 * t395) * t347) * t347 / 0.2e1 + ((t411 * t339 + t410 * t340 - t412 * t397) * t358 + (t416 * t339 + t414 * t340 - t418 * t397) * t348 + (t415 * t339 + t413 * t340 - t417 * t397) * t347) * t348 / 0.2e1 + (((-t411 * t365 + t410 * t368) * t358 + (-t416 * t365 + t414 * t368) * t348 + (-t415 * t365 + t413 * t368) * t347) * t369 + (t417 * t347 + t418 * t348 + t412 * t358) * t366) * t358 / 0.2e1 + (Icges(4,1) + Icges(2,3) + Icges(3,1) + m(2) * (t354 ^ 2 + t356 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
