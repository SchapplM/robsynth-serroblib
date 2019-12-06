% Calculate kinetic energy for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:02
% EndTime: 2019-12-05 18:27:03
% DurationCPUTime: 1.50s
% Computational Cost: add. (1106->207), mult. (1085->325), div. (0->0), fcn. (926->10), ass. (0->127)
t426 = Icges(3,3) + Icges(4,3);
t353 = qJ(2) + pkin(9);
t344 = sin(t353);
t345 = cos(t353);
t355 = sin(qJ(2));
t357 = cos(qJ(2));
t425 = Icges(3,5) * t357 + Icges(4,5) * t345 - Icges(3,6) * t355 - Icges(4,6) * t344;
t356 = sin(qJ(1));
t358 = cos(qJ(1));
t424 = t425 * t356 - t426 * t358;
t423 = t426 * t356 + t425 * t358;
t422 = Icges(3,5) * t355 + Icges(4,5) * t344 + Icges(3,6) * t357 + Icges(4,6) * t345;
t408 = Icges(4,4) * t344;
t323 = Icges(4,2) * t345 + t408;
t407 = Icges(4,4) * t345;
t324 = Icges(4,1) * t344 + t407;
t410 = Icges(3,4) * t355;
t331 = Icges(3,2) * t357 + t410;
t409 = Icges(3,4) * t357;
t332 = Icges(3,1) * t355 + t409;
t421 = -t323 * t344 + t324 * t345 - t331 * t355 + t332 * t357;
t379 = -Icges(4,2) * t344 + t407;
t298 = Icges(4,6) * t356 + t379 * t358;
t383 = Icges(4,1) * t345 - t408;
t300 = Icges(4,5) * t356 + t383 * t358;
t380 = -Icges(3,2) * t355 + t409;
t307 = Icges(3,6) * t356 + t380 * t358;
t384 = Icges(3,1) * t357 - t410;
t309 = Icges(3,5) * t356 + t384 * t358;
t420 = -t298 * t344 + t300 * t345 - t307 * t355 + t309 * t357;
t297 = -Icges(4,6) * t358 + t379 * t356;
t299 = -Icges(4,5) * t358 + t383 * t356;
t306 = -Icges(3,6) * t358 + t380 * t356;
t308 = -Icges(3,5) * t358 + t384 * t356;
t419 = t297 * t344 - t299 * t345 + t306 * t355 - t308 * t357;
t415 = pkin(2) * t355;
t346 = qJ(4) + t353;
t340 = sin(t346);
t414 = pkin(4) * t340;
t412 = t357 * pkin(2);
t406 = Icges(5,4) * t340;
t341 = cos(t346);
t405 = Icges(5,4) * t341;
t342 = qJ(5) + t346;
t337 = sin(t342);
t404 = Icges(6,4) * t337;
t338 = cos(t342);
t403 = Icges(6,4) * t338;
t293 = -qJ(3) * t358 + t412 * t356;
t294 = qJ(3) * t356 + t412 * t358;
t349 = qJD(2) * t356;
t396 = qJD(2) * t358;
t402 = t293 * t349 + t294 * t396;
t336 = pkin(1) * t356 - pkin(6) * t358;
t401 = -t293 - t336;
t400 = pkin(4) * t341;
t399 = pkin(3) * t345;
t328 = qJD(4) * t356 + t349;
t395 = -qJD(2) - qJD(4);
t272 = -pkin(7) * t358 + t399 * t356;
t394 = -t272 + t401;
t273 = pkin(7) * t356 + t399 * t358;
t391 = t272 * t349 + t273 * t396 + t402;
t327 = qJD(1) * (pkin(1) * t358 + pkin(6) * t356);
t390 = qJD(1) * t294 - qJD(3) * t358 + t327;
t389 = rSges(3,1) * t357 - rSges(3,2) * t355;
t388 = rSges(4,1) * t345 - rSges(4,2) * t344;
t387 = rSges(5,1) * t341 - rSges(5,2) * t340;
t386 = rSges(6,1) * t338 - rSges(6,2) * t337;
t385 = qJD(2) * (-rSges(4,1) * t344 - rSges(4,2) * t345 - t415);
t382 = Icges(5,1) * t341 - t406;
t381 = Icges(6,1) * t338 - t404;
t378 = -Icges(5,2) * t340 + t405;
t377 = -Icges(6,2) * t337 + t403;
t374 = Icges(5,5) * t341 - Icges(5,6) * t340;
t373 = Icges(6,5) * t338 - Icges(6,6) * t337;
t366 = qJD(2) * (-pkin(3) * t344 - t415);
t320 = qJD(5) * t356 + t328;
t321 = (-qJD(5) + t395) * t358;
t365 = (Icges(6,5) * t337 + Icges(6,6) * t338) * qJD(1) + (-Icges(6,3) * t358 + t373 * t356) * t321 + (Icges(6,3) * t356 + t373 * t358) * t320;
t329 = t395 * t358;
t364 = (Icges(5,5) * t340 + Icges(5,6) * t341) * qJD(1) + (-Icges(5,3) * t358 + t374 * t356) * t329 + (Icges(5,3) * t356 + t374 * t358) * t328;
t348 = qJD(3) * t356;
t363 = t358 * t366 + t348;
t362 = qJD(1) * t273 + t356 * t366 + t390;
t276 = -Icges(6,6) * t358 + t377 * t356;
t277 = Icges(6,6) * t356 + t377 * t358;
t278 = -Icges(6,5) * t358 + t381 * t356;
t279 = Icges(6,5) * t356 + t381 * t358;
t311 = Icges(6,2) * t338 + t404;
t312 = Icges(6,1) * t337 + t403;
t361 = (-t277 * t337 + t279 * t338) * t320 + (-t276 * t337 + t278 * t338) * t321 + (-t311 * t337 + t312 * t338) * qJD(1);
t286 = -Icges(5,6) * t358 + t378 * t356;
t287 = Icges(5,6) * t356 + t378 * t358;
t288 = -Icges(5,5) * t358 + t382 * t356;
t289 = Icges(5,5) * t356 + t382 * t358;
t317 = Icges(5,2) * t341 + t406;
t318 = Icges(5,1) * t340 + t405;
t360 = (-t287 * t340 + t289 * t341) * t328 + (-t286 * t340 + t288 * t341) * t329 + (-t317 * t340 + t318 * t341) * qJD(1);
t335 = rSges(2,1) * t358 - rSges(2,2) * t356;
t334 = rSges(2,1) * t356 + rSges(2,2) * t358;
t333 = rSges(3,1) * t355 + rSges(3,2) * t357;
t319 = rSges(5,1) * t340 + rSges(5,2) * t341;
t315 = rSges(6,1) * t337 + rSges(6,2) * t338;
t314 = rSges(3,3) * t356 + t389 * t358;
t313 = -rSges(3,3) * t358 + t389 * t356;
t302 = rSges(4,3) * t356 + t388 * t358;
t301 = -rSges(4,3) * t358 + t388 * t356;
t292 = rSges(5,3) * t356 + t387 * t358;
t291 = -rSges(5,3) * t358 + t387 * t356;
t283 = rSges(6,3) * t356 + t386 * t358;
t282 = -rSges(6,3) * t358 + t386 * t356;
t268 = qJD(1) * t314 - t333 * t349 + t327;
t267 = -t333 * t396 + (-t313 - t336) * qJD(1);
t266 = (t313 * t356 + t314 * t358) * qJD(2);
t265 = pkin(8) * t356 + t400 * t358;
t264 = -pkin(8) * t358 + t400 * t356;
t263 = qJD(1) * t302 + t356 * t385 + t390;
t262 = t348 + t358 * t385 + (-t301 + t401) * qJD(1);
t261 = (t301 * t356 + t302 * t358) * qJD(2) + t402;
t260 = qJD(1) * t292 - t319 * t328 + t362;
t259 = t319 * t329 + (-t291 + t394) * qJD(1) + t363;
t258 = t291 * t328 - t292 * t329 + t391;
t257 = -t328 * t414 - t315 * t320 + (t265 + t283) * qJD(1) + t362;
t256 = t329 * t414 + t315 * t321 + (-t264 - t282 + t394) * qJD(1) + t363;
t255 = t264 * t328 - t265 * t329 + t282 * t320 - t283 * t321 + t391;
t1 = m(3) * (t266 ^ 2 + t267 ^ 2 + t268 ^ 2) / 0.2e1 + m(4) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(5) * (t258 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + t328 * (t364 * t356 + t360 * t358) / 0.2e1 + t329 * (t360 * t356 - t364 * t358) / 0.2e1 + m(6) * (t255 ^ 2 + t256 ^ 2 + t257 ^ 2) / 0.2e1 + t320 * (t365 * t356 + t361 * t358) / 0.2e1 + t321 * (t361 * t356 - t365 * t358) / 0.2e1 + (m(2) * (t334 ^ 2 + t335 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t423 * t356 ^ 2 + (t419 * t358 + (t420 - t424) * t356) * t358) * qJD(2) + (t422 * t356 + t421 * t358) * qJD(1)) * t349 / 0.2e1 - ((t424 * t358 ^ 2 + (t420 * t356 + (t419 - t423) * t358) * t356) * qJD(2) + (t421 * t356 - t422 * t358) * qJD(1)) * t396 / 0.2e1 + ((t287 * t341 + t289 * t340) * t328 + (t286 * t341 + t288 * t340) * t329 + (t277 * t338 + t279 * t337) * t320 + (t276 * t338 + t278 * t337) * t321 + ((-t297 * t345 - t299 * t344 - t306 * t357 - t308 * t355) * t358 + (t298 * t345 + t300 * t344 + t307 * t357 + t309 * t355) * t356) * qJD(2) + (t338 * t311 + t337 * t312 + t341 * t317 + t340 * t318 + t345 * t323 + t344 * t324 + t357 * t331 + t355 * t332) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
