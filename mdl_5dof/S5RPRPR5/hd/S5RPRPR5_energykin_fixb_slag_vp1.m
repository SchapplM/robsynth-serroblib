% Calculate kinetic energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:23
% EndTime: 2020-01-03 11:41:25
% DurationCPUTime: 2.03s
% Computational Cost: add. (1145->235), mult. (1594->372), div. (0->0), fcn. (1632->10), ass. (0->115)
t425 = Icges(4,3) + Icges(5,3);
t371 = sin(pkin(8));
t372 = cos(pkin(8));
t370 = qJ(3) + pkin(9);
t366 = cos(t370);
t377 = cos(qJ(1));
t404 = t377 * t366;
t365 = sin(t370);
t375 = sin(qJ(1));
t412 = t375 * t365;
t341 = -t372 * t412 - t404;
t405 = t377 * t365;
t411 = t375 * t366;
t342 = t372 * t411 - t405;
t343 = t372 * t405 - t411;
t344 = -t372 * t404 - t412;
t376 = cos(qJ(3));
t401 = t377 * t376;
t374 = sin(qJ(3));
t409 = t375 * t374;
t348 = -t372 * t409 - t401;
t402 = t377 * t374;
t408 = t375 * t376;
t349 = t372 * t408 - t402;
t350 = t372 * t402 - t408;
t351 = -t372 * t401 - t409;
t403 = t377 * t371;
t410 = t375 * t371;
t422 = (-Icges(4,5) * t351 - Icges(5,5) * t344 - Icges(4,6) * t350 - Icges(5,6) * t343 + t425 * t403) * t377 + (Icges(4,5) * t349 + Icges(5,5) * t342 + Icges(4,6) * t348 + Icges(5,6) * t341 + t425 * t410) * t375;
t424 = t422 * t371;
t423 = t425 * t372 + (-Icges(4,5) * t376 - Icges(5,5) * t366 + Icges(4,6) * t374 + Icges(5,6) * t365) * t371;
t399 = pkin(4) * t366;
t421 = pkin(7) * t371 + t399 * t372;
t416 = t376 * pkin(3);
t420 = qJ(4) * t371 + t416 * t372;
t367 = qJ(5) + t370;
t362 = sin(t367);
t414 = t375 * t362;
t363 = cos(t367);
t413 = t375 * t363;
t407 = t377 * t362;
t406 = t377 * t363;
t355 = qJD(1) * (t375 * pkin(1) - t377 * qJ(2));
t386 = pkin(2) * t372 + pkin(6) * t371;
t400 = qJD(1) * t386 * t375 + t355;
t397 = qJD(3) * t371;
t396 = qJD(4) * t377;
t395 = qJD(3) + qJD(5);
t394 = t375 * t397;
t393 = t377 * t397;
t392 = pkin(4) * t365;
t328 = -qJ(4) * t372 + t416 * t371;
t391 = qJD(3) * (pkin(7) * t372 - t399 * t371 - t328);
t390 = qJD(3) * (-t328 + t372 * rSges(5,3) - (rSges(5,1) * t366 - rSges(5,2) * t365) * t371);
t358 = -t377 * pkin(1) - t375 * qJ(2);
t389 = (t386 * t377 - t358) * qJD(1);
t317 = -pkin(3) * t402 + t375 * t420;
t318 = -pkin(3) * t409 - t377 * t420;
t385 = -qJD(4) * t372 + t317 * t393 + t318 * t394;
t384 = rSges(3,1) * t372 - rSges(3,2) * t371;
t383 = -(-t372 * rSges(4,3) + (rSges(4,1) * t376 - rSges(4,2) * t374) * t371) * t397 - qJD(2);
t380 = qJD(4) * t410 + t389;
t361 = -qJD(3) * t372 + qJD(1);
t379 = -qJD(2) * t375 + t361 * t317 + t400;
t359 = -t377 * rSges(2,1) + t375 * rSges(2,2);
t357 = t375 * rSges(2,1) + t377 * rSges(2,2);
t354 = -t395 * t372 + qJD(1);
t347 = t395 * t403;
t346 = t395 * t410;
t339 = -Icges(4,5) * t372 + (Icges(4,1) * t376 - Icges(4,4) * t374) * t371;
t338 = -Icges(4,6) * t372 + (Icges(4,4) * t376 - Icges(4,2) * t374) * t371;
t336 = -t372 * t406 - t414;
t335 = t372 * t407 - t413;
t334 = t372 * t413 - t407;
t333 = -t372 * t414 - t406;
t331 = -Icges(5,5) * t372 + (Icges(5,1) * t366 - Icges(5,4) * t365) * t371;
t330 = -Icges(5,6) * t372 + (Icges(5,4) * t366 - Icges(5,2) * t365) * t371;
t327 = -t372 * rSges(6,3) + (rSges(6,1) * t363 - rSges(6,2) * t362) * t371;
t326 = -Icges(6,5) * t372 + (Icges(6,1) * t363 - Icges(6,4) * t362) * t371;
t325 = -Icges(6,6) * t372 + (Icges(6,4) * t363 - Icges(6,2) * t362) * t371;
t324 = -Icges(6,3) * t372 + (Icges(6,5) * t363 - Icges(6,6) * t362) * t371;
t323 = -qJD(2) * t377 + (t375 * rSges(3,3) + t384 * t377 - t358) * qJD(1);
t322 = -qJD(1) * t377 * rSges(3,3) + t355 + (qJD(1) * t384 - qJD(2)) * t375;
t320 = t351 * rSges(4,1) + t350 * rSges(4,2) - rSges(4,3) * t403;
t319 = t349 * rSges(4,1) + t348 * rSges(4,2) + rSges(4,3) * t410;
t316 = Icges(4,1) * t351 + Icges(4,4) * t350 - Icges(4,5) * t403;
t315 = Icges(4,1) * t349 + Icges(4,4) * t348 + Icges(4,5) * t410;
t314 = Icges(4,4) * t351 + Icges(4,2) * t350 - Icges(4,6) * t403;
t313 = Icges(4,4) * t349 + Icges(4,2) * t348 + Icges(4,6) * t410;
t307 = t344 * rSges(5,1) + t343 * rSges(5,2) - rSges(5,3) * t403;
t306 = t342 * rSges(5,1) + t341 * rSges(5,2) + rSges(5,3) * t410;
t305 = Icges(5,1) * t344 + Icges(5,4) * t343 - Icges(5,5) * t403;
t304 = Icges(5,1) * t342 + Icges(5,4) * t341 + Icges(5,5) * t410;
t303 = Icges(5,4) * t344 + Icges(5,2) * t343 - Icges(5,6) * t403;
t302 = Icges(5,4) * t342 + Icges(5,2) * t341 + Icges(5,6) * t410;
t299 = t336 * rSges(6,1) + t335 * rSges(6,2) - rSges(6,3) * t403;
t298 = t334 * rSges(6,1) + t333 * rSges(6,2) + rSges(6,3) * t410;
t297 = Icges(6,1) * t336 + Icges(6,4) * t335 - Icges(6,5) * t403;
t296 = Icges(6,1) * t334 + Icges(6,4) * t333 + Icges(6,5) * t410;
t295 = Icges(6,4) * t336 + Icges(6,2) * t335 - Icges(6,6) * t403;
t294 = Icges(6,4) * t334 + Icges(6,2) * t333 + Icges(6,6) * t410;
t293 = Icges(6,5) * t336 + Icges(6,6) * t335 - Icges(6,3) * t403;
t292 = Icges(6,5) * t334 + Icges(6,6) * t333 + Icges(6,3) * t410;
t291 = -t392 * t375 - t421 * t377;
t290 = t421 * t375 - t392 * t377;
t289 = (t319 * t377 + t320 * t375) * t397;
t288 = -t361 * t320 + t383 * t377 + t389;
t287 = t361 * t319 + t383 * t375 + t400;
t286 = (-t307 - t318) * t361 + (t371 * t390 - qJD(2)) * t377 + t380;
t285 = t361 * t306 + (t375 * t390 - t396) * t371 + t379;
t284 = (t306 * t377 + t307 * t375) * t397 + t385;
t283 = -t354 * t299 - t347 * t327 + (-t291 - t318) * t361 + (t371 * t391 - qJD(2)) * t377 + t380;
t282 = t361 * t290 + t354 * t298 - t346 * t327 + (t375 * t391 - t396) * t371 + t379;
t281 = t347 * t298 + t346 * t299 + (t290 * t377 + t291 * t375) * t397 + t385;
t1 = m(3) * (t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(4) * (t287 ^ 2 + t288 ^ 2 + t289 ^ 2) / 0.2e1 + m(5) * (t284 ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + m(6) * (t281 ^ 2 + t282 ^ 2 + t283 ^ 2) / 0.2e1 + t354 * ((-t292 * t346 + t293 * t347 - t324 * t354) * t372 + ((-t325 * t362 + t326 * t363) * t354 + (-t294 * t362 + t296 * t363) * t346 - (-t295 * t362 + t297 * t363) * t347) * t371) / 0.2e1 + t346 * ((t324 * t410 + t333 * t325 + t334 * t326) * t354 + (t292 * t410 + t333 * t294 + t334 * t296) * t346 - (t293 * t410 + t333 * t295 + t334 * t297) * t347) / 0.2e1 - t347 * ((-t324 * t403 + t335 * t325 + t336 * t326) * t354 + (-t292 * t403 + t335 * t294 + t336 * t296) * t346 - (-t293 * t403 + t335 * t295 + t336 * t297) * t347) / 0.2e1 + ((-t422 * t372 + ((t303 * t365 - t305 * t366 + t314 * t374 - t316 * t376) * t377 + (-t302 * t365 + t304 * t366 - t313 * t374 + t315 * t376) * t375) * t371) * t397 + (t423 * t372 + (-t330 * t365 + t331 * t366 - t338 * t374 + t339 * t376) * t371) * t361) * t361 / 0.2e1 + (((-t341 * t303 - t342 * t305 - t348 * t314 - t349 * t316) * t377 + (t341 * t302 + t342 * t304 + t348 * t313 + t349 * t315 + t424) * t375) * t397 + (t341 * t330 + t342 * t331 + t348 * t338 + t349 * t339 - t423 * t410) * t361) * t394 / 0.2e1 - (((-t343 * t303 - t344 * t305 - t350 * t314 - t351 * t316 - t424) * t377 + (t343 * t302 + t344 * t304 + t350 * t313 + t351 * t315) * t375) * t397 + (t343 * t330 + t344 * t331 + t350 * t338 + t351 * t339 + t423 * t403) * t361) * t393 / 0.2e1 + (m(2) * (t357 ^ 2 + t359 ^ 2) + Icges(2,3) + Icges(3,2) * t372 ^ 2 + (Icges(3,1) * t371 + 0.2e1 * Icges(3,4) * t372) * t371) * qJD(1) ^ 2 / 0.2e1;
T = t1;
