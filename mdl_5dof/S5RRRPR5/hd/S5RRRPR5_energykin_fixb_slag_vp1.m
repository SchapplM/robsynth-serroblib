% Calculate kinetic energy for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:12:59
% EndTime: 2019-12-31 21:13:01
% DurationCPUTime: 1.74s
% Computational Cost: add. (1217->231), mult. (1297->362), div. (0->0), fcn. (1194->10), ass. (0->132)
t424 = Icges(4,3) + Icges(5,3);
t360 = qJ(2) + qJ(3);
t352 = pkin(9) + t360;
t349 = sin(t352);
t350 = cos(t352);
t356 = sin(t360);
t357 = cos(t360);
t423 = Icges(4,5) * t357 + Icges(5,5) * t350 - Icges(4,6) * t356 - Icges(5,6) * t349;
t363 = sin(qJ(1));
t366 = cos(qJ(1));
t408 = Icges(5,4) * t350;
t382 = -Icges(5,2) * t349 + t408;
t294 = -Icges(5,6) * t366 + t363 * t382;
t295 = Icges(5,6) * t363 + t366 * t382;
t409 = Icges(5,4) * t349;
t385 = Icges(5,1) * t350 - t409;
t296 = -Icges(5,5) * t366 + t363 * t385;
t297 = Icges(5,5) * t363 + t366 * t385;
t410 = Icges(4,4) * t357;
t383 = -Icges(4,2) * t356 + t410;
t305 = -Icges(4,6) * t366 + t363 * t383;
t306 = Icges(4,6) * t363 + t366 * t383;
t411 = Icges(4,4) * t356;
t386 = Icges(4,1) * t357 - t411;
t307 = -Icges(4,5) * t366 + t363 * t386;
t308 = Icges(4,5) * t363 + t366 * t386;
t328 = Icges(5,2) * t350 + t409;
t329 = Icges(5,1) * t349 + t408;
t334 = Icges(4,2) * t357 + t411;
t335 = Icges(4,1) * t356 + t410;
t355 = qJD(2) * t363;
t340 = qJD(3) * t363 + t355;
t341 = (-qJD(2) - qJD(3)) * t366;
t422 = (-t294 * t349 + t296 * t350 - t305 * t356 + t307 * t357) * t341 + (-t295 * t349 + t297 * t350 - t306 * t356 + t308 * t357) * t340 + (-t328 * t349 + t329 * t350 - t334 * t356 + t335 * t357) * qJD(1);
t421 = (t423 * t363 - t424 * t366) * t341 + (t424 * t363 + t423 * t366) * t340 + (Icges(4,5) * t356 + Icges(5,5) * t349 + Icges(4,6) * t357 + Icges(5,6) * t350) * qJD(1);
t417 = pkin(3) * t356;
t365 = cos(qJ(2));
t415 = t365 * pkin(2);
t362 = sin(qJ(2));
t413 = Icges(3,4) * t362;
t412 = Icges(3,4) * t365;
t407 = t349 * t363;
t406 = t349 * t366;
t361 = sin(qJ(5));
t405 = t361 * t363;
t404 = t361 * t366;
t364 = cos(qJ(5));
t403 = t363 * t364;
t402 = t364 * t366;
t301 = -pkin(7) * t366 + t363 * t415;
t302 = pkin(7) * t363 + t366 * t415;
t397 = qJD(2) * t366;
t401 = t301 * t355 + t302 * t397;
t348 = pkin(1) * t363 - pkin(6) * t366;
t400 = -t301 - t348;
t399 = pkin(3) * t357;
t396 = qJD(5) * t349;
t395 = pkin(2) * qJD(2) * t362;
t284 = -qJ(4) * t366 + t363 * t399;
t394 = t340 * t284 + t401;
t393 = -t284 + t400;
t392 = t366 * t395;
t391 = pkin(4) * t350 + pkin(8) * t349;
t390 = rSges(3,1) * t365 - rSges(3,2) * t362;
t389 = rSges(4,1) * t357 - rSges(4,2) * t356;
t388 = rSges(5,1) * t350 - rSges(5,2) * t349;
t387 = Icges(3,1) * t365 - t413;
t384 = -Icges(3,2) * t362 + t412;
t381 = Icges(3,5) * t365 - Icges(3,6) * t362;
t317 = -Icges(3,6) * t366 + t363 * t384;
t319 = -Icges(3,5) * t366 + t363 * t387;
t378 = t317 * t362 - t319 * t365;
t318 = Icges(3,6) * t363 + t366 * t384;
t320 = Icges(3,5) * t363 + t366 * t387;
t377 = -t318 * t362 + t320 * t365;
t343 = Icges(3,2) * t365 + t413;
t344 = Icges(3,1) * t362 + t412;
t376 = -t343 * t362 + t344 * t365;
t338 = qJD(1) * (pkin(1) * t366 + pkin(6) * t363);
t375 = qJD(1) * t302 - t363 * t395 + t338;
t374 = qJD(4) * t363 + t341 * t417 - t392;
t285 = qJ(4) * t363 + t366 * t399;
t371 = qJD(1) * t285 - qJD(4) * t366 + t375;
t347 = rSges(2,1) * t366 - rSges(2,2) * t363;
t346 = rSges(2,1) * t363 + rSges(2,2) * t366;
t345 = rSges(3,1) * t362 + rSges(3,2) * t365;
t342 = Icges(3,5) * t362 + Icges(3,6) * t365;
t339 = -qJD(5) * t350 + qJD(1);
t336 = rSges(4,1) * t356 + rSges(4,2) * t357;
t332 = pkin(4) * t349 - pkin(8) * t350;
t330 = rSges(5,1) * t349 + rSges(5,2) * t350;
t326 = t350 * t402 + t405;
t325 = -t350 * t404 + t403;
t324 = t350 * t403 - t404;
t323 = -t350 * t405 - t402;
t322 = rSges(3,3) * t363 + t366 * t390;
t321 = -rSges(3,3) * t366 + t363 * t390;
t316 = Icges(3,3) * t363 + t366 * t381;
t315 = -Icges(3,3) * t366 + t363 * t381;
t314 = t363 * t396 + t341;
t313 = t366 * t396 + t340;
t312 = t391 * t366;
t311 = t391 * t363;
t310 = rSges(4,3) * t363 + t366 * t389;
t309 = -rSges(4,3) * t366 + t363 * t389;
t299 = rSges(5,3) * t363 + t366 * t388;
t298 = -rSges(5,3) * t366 + t363 * t388;
t289 = -rSges(6,3) * t350 + (rSges(6,1) * t364 - rSges(6,2) * t361) * t349;
t288 = -Icges(6,5) * t350 + (Icges(6,1) * t364 - Icges(6,4) * t361) * t349;
t287 = -Icges(6,6) * t350 + (Icges(6,4) * t364 - Icges(6,2) * t361) * t349;
t286 = -Icges(6,3) * t350 + (Icges(6,5) * t364 - Icges(6,6) * t361) * t349;
t282 = qJD(1) * t322 - t345 * t355 + t338;
t281 = -t345 * t397 + (-t321 - t348) * qJD(1);
t280 = (t321 * t363 + t322 * t366) * qJD(2);
t278 = rSges(6,1) * t326 + rSges(6,2) * t325 + rSges(6,3) * t406;
t277 = rSges(6,1) * t324 + rSges(6,2) * t323 + rSges(6,3) * t407;
t276 = Icges(6,1) * t326 + Icges(6,4) * t325 + Icges(6,5) * t406;
t275 = Icges(6,1) * t324 + Icges(6,4) * t323 + Icges(6,5) * t407;
t274 = Icges(6,4) * t326 + Icges(6,2) * t325 + Icges(6,6) * t406;
t273 = Icges(6,4) * t324 + Icges(6,2) * t323 + Icges(6,6) * t407;
t272 = Icges(6,5) * t326 + Icges(6,6) * t325 + Icges(6,3) * t406;
t271 = Icges(6,5) * t324 + Icges(6,6) * t323 + Icges(6,3) * t407;
t270 = qJD(1) * t310 - t336 * t340 + t375;
t269 = -t392 + t336 * t341 + (-t309 + t400) * qJD(1);
t268 = t309 * t340 - t310 * t341 + t401;
t267 = qJD(1) * t299 + (-t330 - t417) * t340 + t371;
t266 = t330 * t341 + (-t298 + t393) * qJD(1) + t374;
t265 = t298 * t340 + (-t285 - t299) * t341 + t394;
t264 = qJD(1) * t312 + t278 * t339 - t289 * t313 + (-t332 - t417) * t340 + t371;
t263 = -t277 * t339 + t289 * t314 + t332 * t341 + (-t311 + t393) * qJD(1) + t374;
t262 = t277 * t313 - t278 * t314 + t311 * t340 + (-t285 - t312) * t341 + t394;
t1 = m(3) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + ((t363 * t342 + t366 * t376) * qJD(1) + (t363 ^ 2 * t316 + (t378 * t366 + (-t315 + t377) * t363) * t366) * qJD(2)) * t355 / 0.2e1 - ((-t366 * t342 + t376 * t363) * qJD(1) + (t366 ^ 2 * t315 + (t377 * t363 + (-t316 + t378) * t366) * t363) * qJD(2)) * t397 / 0.2e1 + m(4) * (t268 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + m(5) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + m(6) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + t313 * ((t272 * t406 + t325 * t274 + t326 * t276) * t313 + (t271 * t406 + t273 * t325 + t275 * t326) * t314 + (t286 * t406 + t287 * t325 + t288 * t326) * t339) / 0.2e1 + t314 * ((t272 * t407 + t274 * t323 + t276 * t324) * t313 + (t271 * t407 + t323 * t273 + t324 * t275) * t314 + (t286 * t407 + t287 * t323 + t288 * t324) * t339) / 0.2e1 + t339 * ((-t271 * t314 - t272 * t313 - t286 * t339) * t350 + ((-t274 * t361 + t276 * t364) * t313 + (-t273 * t361 + t275 * t364) * t314 + (-t287 * t361 + t288 * t364) * t339) * t349) / 0.2e1 + (t421 * t363 + t422 * t366) * t340 / 0.2e1 + (t422 * t363 - t421 * t366) * t341 / 0.2e1 + (m(2) * (t346 ^ 2 + t347 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t318 * t365 + t320 * t362) * t363 - (t317 * t365 + t319 * t362) * t366) * qJD(2) + (t294 * t350 + t296 * t349 + t305 * t357 + t307 * t356) * t341 + (t295 * t350 + t297 * t349 + t306 * t357 + t308 * t356) * t340 + (t350 * t328 + t349 * t329 + t357 * t334 + t356 * t335 + t365 * t343 + t362 * t344) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
