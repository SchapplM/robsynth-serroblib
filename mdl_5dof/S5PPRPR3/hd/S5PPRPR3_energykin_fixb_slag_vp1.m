% Calculate kinetic energy for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR3_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR3_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:44
% EndTime: 2019-12-05 15:04:46
% DurationCPUTime: 1.64s
% Computational Cost: add. (1125->214), mult. (1796->338), div. (0->0), fcn. (1990->10), ass. (0->103)
t411 = -Icges(5,3) - Icges(4,3);
t364 = qJ(3) + pkin(9);
t363 = cos(t364);
t366 = sin(pkin(7));
t367 = cos(pkin(8));
t362 = sin(t364);
t368 = cos(pkin(7));
t396 = t368 * t362;
t343 = -t366 * t363 + t367 * t396;
t395 = t368 * t363;
t344 = t366 * t362 + t367 * t395;
t371 = sin(qJ(3));
t393 = t368 * t371;
t373 = cos(qJ(3));
t397 = t366 * t373;
t354 = -t367 * t393 + t397;
t392 = t368 * t373;
t398 = t366 * t371;
t355 = t367 * t392 + t398;
t365 = sin(pkin(8));
t394 = t368 * t365;
t410 = Icges(4,5) * t355 + Icges(5,5) * t344 + Icges(4,6) * t354 - Icges(5,6) * t343 - t411 * t394;
t399 = t366 * t367;
t341 = t362 * t399 + t395;
t342 = t363 * t399 - t396;
t352 = -t367 * t398 - t392;
t353 = t367 * t397 - t393;
t400 = t366 * t365;
t409 = Icges(4,5) * t353 + Icges(5,5) * t342 + Icges(4,6) * t352 - Icges(5,6) * t341 - t411 * t400;
t408 = t411 * t367 + (Icges(4,5) * t373 + Icges(5,5) * t363 - Icges(4,6) * t371 - Icges(5,6) * t362) * t365;
t407 = -t408 * t367 + t410 * t394 + t409 * t400;
t404 = t373 * pkin(3);
t402 = t362 * t365;
t401 = t363 * t365;
t376 = qJ(4) * t365 + t404 * t367;
t326 = pkin(3) * t398 + t376 * t368;
t391 = -t344 * rSges(5,1) + t343 * rSges(5,2) - rSges(5,3) * t394 - t326;
t390 = -t344 * pkin(4) - t343 * pkin(6) - t326;
t389 = qJD(2) * t368;
t388 = qJD(3) * t365;
t387 = qJD(3) * t367;
t386 = qJD(4) * t365;
t385 = t366 * t388;
t384 = t368 * t388;
t381 = t366 * t386 - t389;
t325 = -pkin(3) * t393 + t376 * t366;
t336 = -qJ(4) * t367 + t404 * t365;
t361 = qJD(2) * t366;
t380 = t325 * t387 + t336 * t385 + t368 * t386 + t361;
t377 = -qJD(4) * t367 + t325 * t384 + qJD(1);
t375 = qJD(1) ^ 2;
t372 = cos(qJ(5));
t370 = sin(qJ(5));
t356 = qJD(5) * t402 - t387;
t351 = -t367 * t370 + t372 * t401;
t350 = -t367 * t372 - t370 * t401;
t349 = (pkin(4) * t363 + pkin(6) * t362) * t365;
t348 = -t367 * rSges(4,3) + (rSges(4,1) * t373 - rSges(4,2) * t371) * t365;
t347 = -Icges(4,5) * t367 + (Icges(4,1) * t373 - Icges(4,4) * t371) * t365;
t346 = -Icges(4,6) * t367 + (Icges(4,4) * t373 - Icges(4,2) * t371) * t365;
t340 = -t367 * rSges(5,3) + (rSges(5,1) * t363 - rSges(5,2) * t362) * t365;
t339 = -Icges(5,5) * t367 + (Icges(5,1) * t363 - Icges(5,4) * t362) * t365;
t338 = -Icges(5,6) * t367 + (Icges(5,4) * t363 - Icges(5,2) * t362) * t365;
t335 = qJD(5) * t343 + t384;
t334 = qJD(5) * t341 + t385;
t332 = t344 * t372 + t370 * t394;
t331 = -t344 * t370 + t372 * t394;
t330 = t342 * t372 + t370 * t400;
t329 = -t342 * t370 + t372 * t400;
t328 = t355 * rSges(4,1) + t354 * rSges(4,2) + rSges(4,3) * t394;
t327 = t353 * rSges(4,1) + t352 * rSges(4,2) + rSges(4,3) * t400;
t324 = Icges(4,1) * t355 + Icges(4,4) * t354 + Icges(4,5) * t394;
t323 = Icges(4,1) * t353 + Icges(4,4) * t352 + Icges(4,5) * t400;
t322 = Icges(4,4) * t355 + Icges(4,2) * t354 + Icges(4,6) * t394;
t321 = Icges(4,4) * t353 + Icges(4,2) * t352 + Icges(4,6) * t400;
t317 = t342 * pkin(4) + t341 * pkin(6);
t315 = t351 * rSges(6,1) + t350 * rSges(6,2) + rSges(6,3) * t402;
t314 = Icges(6,1) * t351 + Icges(6,4) * t350 + Icges(6,5) * t402;
t313 = Icges(6,4) * t351 + Icges(6,2) * t350 + Icges(6,6) * t402;
t312 = Icges(6,5) * t351 + Icges(6,6) * t350 + Icges(6,3) * t402;
t309 = t342 * rSges(5,1) - t341 * rSges(5,2) + rSges(5,3) * t400;
t308 = Icges(5,1) * t344 - Icges(5,4) * t343 + Icges(5,5) * t394;
t307 = Icges(5,1) * t342 - Icges(5,4) * t341 + Icges(5,5) * t400;
t306 = Icges(5,4) * t344 - Icges(5,2) * t343 + Icges(5,6) * t394;
t305 = Icges(5,4) * t342 - Icges(5,2) * t341 + Icges(5,6) * t400;
t302 = -t389 + (-t328 * t367 - t348 * t394) * qJD(3);
t301 = t361 + (t327 * t367 + t348 * t400) * qJD(3);
t300 = t332 * rSges(6,1) + t331 * rSges(6,2) + t343 * rSges(6,3);
t299 = t330 * rSges(6,1) + t329 * rSges(6,2) + t341 * rSges(6,3);
t298 = Icges(6,1) * t332 + Icges(6,4) * t331 + Icges(6,5) * t343;
t297 = Icges(6,1) * t330 + Icges(6,4) * t329 + Icges(6,5) * t341;
t296 = Icges(6,4) * t332 + Icges(6,2) * t331 + Icges(6,6) * t343;
t295 = Icges(6,4) * t330 + Icges(6,2) * t329 + Icges(6,6) * t341;
t294 = Icges(6,5) * t332 + Icges(6,6) * t331 + Icges(6,3) * t343;
t293 = Icges(6,5) * t330 + Icges(6,6) * t329 + Icges(6,3) * t341;
t292 = qJD(1) + (t327 * t368 - t328 * t366) * t388;
t291 = (t391 * t367 + (-t336 - t340) * t394) * qJD(3) + t381;
t290 = (t309 * t367 + t340 * t400) * qJD(3) + t380;
t289 = (t309 * t368 + t391 * t366) * t388 + t377;
t288 = t356 * t300 - t335 * t315 + (t390 * t367 + (-t336 - t349) * t394) * qJD(3) + t381;
t287 = -t356 * t299 + t334 * t315 + (t317 * t367 + t349 * t400) * qJD(3) + t380;
t286 = t335 * t299 - t334 * t300 + (t317 * t368 + t390 * t366) * t388 + t377;
t1 = m(2) * t375 / 0.2e1 + m(3) * (t375 + (t366 ^ 2 + t368 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t292 ^ 2 + t301 ^ 2 + t302 ^ 2) / 0.2e1 + m(5) * (t289 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + m(6) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + t335 * ((t343 * t294 + t331 * t296 + t332 * t298) * t335 + (t343 * t293 + t331 * t295 + t332 * t297) * t334 + (t343 * t312 + t331 * t313 + t332 * t314) * t356) / 0.2e1 + t334 * ((t341 * t294 + t329 * t296 + t330 * t298) * t335 + (t341 * t293 + t329 * t295 + t330 * t297) * t334 + (t341 * t312 + t329 * t313 + t330 * t314) * t356) / 0.2e1 + t356 * ((t294 * t402 + t350 * t296 + t351 * t298) * t335 + (t293 * t402 + t350 * t295 + t351 * t297) * t334 + (t312 * t402 + t350 * t313 + t351 * t314) * t356) / 0.2e1 + (-(t408 * t367 ^ 2 + (((-t306 * t362 + t308 * t363 - t322 * t371 + t324 * t373) * t368 + (-t305 * t362 + t307 * t363 - t321 * t371 + t323 * t373) * t366) * t365 + (t338 * t362 - t339 * t363 + t346 * t371 - t347 * t373 - t409 * t366 - t410 * t368) * t367) * t365) * t367 / 0.2e1 + (((-t341 * t306 + t342 * t308 + t352 * t322 + t353 * t324) * t394 + (t341 * t338 - t342 * t339 - t352 * t346 - t353 * t347) * t367 + (-t341 * t305 + t342 * t307 + t352 * t321 + t353 * t323 + t407) * t400) * t366 + ((-t343 * t305 + t344 * t307 + t354 * t321 + t355 * t323) * t400 + (t343 * t338 - t344 * t339 - t354 * t346 - t355 * t347) * t367 + (-t343 * t306 + t344 * t308 + t354 * t322 + t355 * t324 + t407) * t394) * t368) * t365 / 0.2e1) * qJD(3) ^ 2;
T = t1;
