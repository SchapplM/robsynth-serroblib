% Calculate kinetic energy for
% S5RRPRR8
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:51
% EndTime: 2019-12-31 20:16:53
% DurationCPUTime: 1.65s
% Computational Cost: add. (1199->229), mult. (1279->365), div. (0->0), fcn. (1176->10), ass. (0->132)
t431 = Icges(3,3) + Icges(4,3);
t360 = qJ(2) + pkin(9);
t352 = sin(t360);
t353 = cos(t360);
t363 = sin(qJ(2));
t366 = cos(qJ(2));
t430 = Icges(3,5) * t366 + Icges(4,5) * t353 - Icges(3,6) * t363 - Icges(4,6) * t352;
t364 = sin(qJ(1));
t367 = cos(qJ(1));
t429 = t430 * t364 - t431 * t367;
t428 = t431 * t364 + t430 * t367;
t427 = Icges(3,5) * t363 + Icges(4,5) * t352 + Icges(3,6) * t366 + Icges(4,6) * t353;
t414 = Icges(4,4) * t352;
t334 = Icges(4,2) * t353 + t414;
t413 = Icges(4,4) * t353;
t335 = Icges(4,1) * t352 + t413;
t416 = Icges(3,4) * t363;
t343 = Icges(3,2) * t366 + t416;
t415 = Icges(3,4) * t366;
t344 = Icges(3,1) * t363 + t415;
t426 = -t334 * t352 + t335 * t353 - t343 * t363 + t344 * t366;
t384 = -Icges(4,2) * t352 + t413;
t307 = Icges(4,6) * t364 + t384 * t367;
t387 = Icges(4,1) * t353 - t414;
t309 = Icges(4,5) * t364 + t387 * t367;
t385 = -Icges(3,2) * t363 + t415;
t319 = Icges(3,6) * t364 + t385 * t367;
t388 = Icges(3,1) * t366 - t416;
t321 = Icges(3,5) * t364 + t388 * t367;
t425 = -t307 * t352 + t309 * t353 - t319 * t363 + t321 * t366;
t306 = -Icges(4,6) * t367 + t384 * t364;
t308 = -Icges(4,5) * t367 + t387 * t364;
t318 = -Icges(3,6) * t367 + t385 * t364;
t320 = -Icges(3,5) * t367 + t388 * t364;
t424 = t306 * t352 - t308 * t353 + t318 * t363 - t320 * t366;
t420 = pkin(2) * t363;
t418 = t366 * pkin(2);
t354 = qJ(4) + t360;
t349 = sin(t354);
t412 = Icges(5,4) * t349;
t350 = cos(t354);
t411 = Icges(5,4) * t350;
t410 = t349 * t364;
t409 = t349 * t367;
t362 = sin(qJ(5));
t408 = t364 * t362;
t365 = cos(qJ(5));
t407 = t364 * t365;
t406 = t367 * t362;
t405 = t367 * t365;
t302 = -qJ(3) * t367 + t418 * t364;
t303 = qJ(3) * t364 + t418 * t367;
t357 = qJD(2) * t364;
t400 = qJD(2) * t367;
t404 = t302 * t357 + t303 * t400;
t348 = t364 * pkin(1) - t367 * pkin(6);
t403 = -t302 - t348;
t402 = pkin(3) * t353;
t340 = qJD(4) * t364 + t357;
t399 = qJD(5) * t349;
t285 = -pkin(7) * t367 + t402 * t364;
t398 = -t285 + t403;
t341 = (-qJD(2) - qJD(4)) * t367;
t286 = pkin(7) * t364 + t402 * t367;
t395 = t285 * t357 + t286 * t400 + t404;
t394 = pkin(4) * t350 + pkin(8) * t349;
t338 = qJD(1) * (t367 * pkin(1) + t364 * pkin(6));
t393 = qJD(1) * t303 - qJD(3) * t367 + t338;
t392 = rSges(3,1) * t366 - rSges(3,2) * t363;
t391 = rSges(4,1) * t353 - rSges(4,2) * t352;
t390 = rSges(5,1) * t350 - rSges(5,2) * t349;
t389 = qJD(2) * (-t352 * rSges(4,1) - t353 * rSges(4,2) - t420);
t386 = Icges(5,1) * t350 - t412;
t383 = -Icges(5,2) * t349 + t411;
t380 = Icges(5,5) * t350 - Icges(5,6) * t349;
t373 = qJD(2) * (-pkin(3) * t352 - t420);
t372 = (Icges(5,5) * t349 + Icges(5,6) * t350) * qJD(1) + (-Icges(5,3) * t367 + t380 * t364) * t341 + (Icges(5,3) * t364 + t380 * t367) * t340;
t356 = qJD(3) * t364;
t371 = t367 * t373 + t356;
t370 = qJD(1) * t286 + t364 * t373 + t393;
t295 = -Icges(5,6) * t367 + t383 * t364;
t296 = Icges(5,6) * t364 + t383 * t367;
t297 = -Icges(5,5) * t367 + t386 * t364;
t298 = Icges(5,5) * t364 + t386 * t367;
t329 = Icges(5,2) * t350 + t412;
t330 = Icges(5,1) * t349 + t411;
t369 = (-t296 * t349 + t298 * t350) * t340 + (-t295 * t349 + t297 * t350) * t341 + (-t329 * t349 + t330 * t350) * qJD(1);
t347 = t367 * rSges(2,1) - t364 * rSges(2,2);
t346 = t364 * rSges(2,1) + t367 * rSges(2,2);
t345 = t363 * rSges(3,1) + t366 * rSges(3,2);
t339 = -qJD(5) * t350 + qJD(1);
t332 = t349 * pkin(4) - t350 * pkin(8);
t331 = t349 * rSges(5,1) + t350 * rSges(5,2);
t327 = t350 * t405 + t408;
t326 = -t350 * t406 + t407;
t325 = t350 * t407 - t406;
t324 = -t350 * t408 - t405;
t323 = t364 * rSges(3,3) + t392 * t367;
t322 = -t367 * rSges(3,3) + t392 * t364;
t315 = t364 * t399 + t341;
t314 = t367 * t399 + t340;
t313 = t394 * t367;
t312 = t394 * t364;
t311 = t364 * rSges(4,3) + t391 * t367;
t310 = -t367 * rSges(4,3) + t391 * t364;
t301 = t364 * rSges(5,3) + t390 * t367;
t300 = -t367 * rSges(5,3) + t390 * t364;
t290 = -t350 * rSges(6,3) + (rSges(6,1) * t365 - rSges(6,2) * t362) * t349;
t289 = -Icges(6,5) * t350 + (Icges(6,1) * t365 - Icges(6,4) * t362) * t349;
t288 = -Icges(6,6) * t350 + (Icges(6,4) * t365 - Icges(6,2) * t362) * t349;
t287 = -Icges(6,3) * t350 + (Icges(6,5) * t365 - Icges(6,6) * t362) * t349;
t281 = qJD(1) * t323 - t345 * t357 + t338;
t280 = -t345 * t400 + (-t322 - t348) * qJD(1);
t279 = (t322 * t364 + t323 * t367) * qJD(2);
t278 = t327 * rSges(6,1) + t326 * rSges(6,2) + rSges(6,3) * t409;
t277 = t325 * rSges(6,1) + t324 * rSges(6,2) + rSges(6,3) * t410;
t276 = Icges(6,1) * t327 + Icges(6,4) * t326 + Icges(6,5) * t409;
t275 = Icges(6,1) * t325 + Icges(6,4) * t324 + Icges(6,5) * t410;
t274 = Icges(6,4) * t327 + Icges(6,2) * t326 + Icges(6,6) * t409;
t273 = Icges(6,4) * t325 + Icges(6,2) * t324 + Icges(6,6) * t410;
t272 = Icges(6,5) * t327 + Icges(6,6) * t326 + Icges(6,3) * t409;
t271 = Icges(6,5) * t325 + Icges(6,6) * t324 + Icges(6,3) * t410;
t270 = qJD(1) * t311 + t364 * t389 + t393;
t269 = t356 + t367 * t389 + (-t310 + t403) * qJD(1);
t268 = (t310 * t364 + t311 * t367) * qJD(2) + t404;
t267 = qJD(1) * t301 - t340 * t331 + t370;
t266 = t341 * t331 + (-t300 + t398) * qJD(1) + t371;
t265 = t340 * t300 - t341 * t301 + t395;
t264 = qJD(1) * t313 + t339 * t278 - t314 * t290 - t340 * t332 + t370;
t263 = -t339 * t277 + t315 * t290 + t341 * t332 + (-t312 + t398) * qJD(1) + t371;
t262 = t314 * t277 - t315 * t278 + t340 * t312 - t341 * t313 + t395;
t1 = m(3) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(4) * (t268 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + m(5) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + t340 * (t372 * t364 + t369 * t367) / 0.2e1 + t341 * (t369 * t364 - t372 * t367) / 0.2e1 + m(6) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + t314 * ((t272 * t409 + t326 * t274 + t327 * t276) * t314 + (t271 * t409 + t326 * t273 + t327 * t275) * t315 + (t287 * t409 + t326 * t288 + t327 * t289) * t339) / 0.2e1 + t315 * ((t272 * t410 + t324 * t274 + t325 * t276) * t314 + (t271 * t410 + t324 * t273 + t325 * t275) * t315 + (t287 * t410 + t324 * t288 + t325 * t289) * t339) / 0.2e1 + t339 * ((-t271 * t315 - t272 * t314 - t287 * t339) * t350 + ((-t274 * t362 + t276 * t365) * t314 + (-t273 * t362 + t275 * t365) * t315 + (-t288 * t362 + t289 * t365) * t339) * t349) / 0.2e1 + (m(2) * (t346 ^ 2 + t347 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t428 * t364 ^ 2 + (t424 * t367 + (t425 - t429) * t364) * t367) * qJD(2) + (t364 * t427 + t367 * t426) * qJD(1)) * t357 / 0.2e1 - ((t429 * t367 ^ 2 + (t425 * t364 + (t424 - t428) * t367) * t364) * qJD(2) + (t364 * t426 - t427 * t367) * qJD(1)) * t400 / 0.2e1 + ((t350 * t296 + t349 * t298) * t340 + (t350 * t295 + t349 * t297) * t341 + ((-t353 * t306 - t352 * t308 - t366 * t318 - t363 * t320) * t367 + (t353 * t307 + t352 * t309 + t366 * t319 + t363 * t321) * t364) * qJD(2) + (t350 * t329 + t349 * t330 + t353 * t334 + t352 * t335 + t366 * t343 + t363 * t344) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
