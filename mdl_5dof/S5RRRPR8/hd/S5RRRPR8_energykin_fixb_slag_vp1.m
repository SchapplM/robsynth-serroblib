% Calculate kinetic energy for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:05
% EndTime: 2019-12-31 21:19:07
% DurationCPUTime: 2.40s
% Computational Cost: add. (935->212), mult. (1277->333), div. (0->0), fcn. (1174->8), ass. (0->124)
t439 = Icges(4,4) + Icges(5,6);
t438 = Icges(4,1) + Icges(5,2);
t437 = -Icges(4,2) - Icges(5,3);
t361 = qJ(2) + qJ(3);
t360 = cos(t361);
t436 = t439 * t360;
t359 = sin(t361);
t435 = t439 * t359;
t434 = Icges(5,4) - Icges(4,5);
t433 = Icges(5,5) - Icges(4,6);
t432 = t437 * t359 + t436;
t431 = t438 * t360 - t435;
t430 = Icges(5,1) + Icges(4,3);
t364 = sin(qJ(1));
t367 = cos(qJ(1));
t429 = t432 * t364 + t433 * t367;
t428 = -t433 * t364 + t432 * t367;
t427 = t431 * t364 + t434 * t367;
t426 = -t434 * t364 + t431 * t367;
t425 = t437 * t360 - t435;
t424 = t438 * t359 + t436;
t423 = t433 * t359 - t434 * t360;
t358 = qJD(2) * t364;
t344 = qJD(3) * t364 + t358;
t345 = (-qJD(2) - qJD(3)) * t367;
t422 = (-t429 * t359 + t427 * t360) * t345 + (-t428 * t359 + t426 * t360) * t344 + (t425 * t359 + t424 * t360) * qJD(1);
t421 = (t423 * t364 - t430 * t367) * t345 + (t430 * t364 + t423 * t367) * t344 + (-t434 * t359 - t433 * t360) * qJD(1);
t417 = pkin(8) * t359;
t366 = cos(qJ(2));
t415 = pkin(2) * t366;
t363 = sin(qJ(2));
t413 = Icges(3,4) * t363;
t412 = Icges(3,4) * t366;
t407 = t360 * t364;
t406 = t360 * t367;
t362 = sin(qJ(5));
t405 = t362 * t364;
t404 = t362 * t367;
t365 = cos(qJ(5));
t403 = t364 * t365;
t402 = t365 * t367;
t297 = -pkin(7) * t367 + t364 * t415;
t298 = pkin(7) * t364 + t367 * t415;
t399 = qJD(2) * t367;
t401 = t297 * t358 + t298 * t399;
t354 = pkin(1) * t364 - pkin(6) * t367;
t400 = -t297 - t354;
t398 = qJD(4) * t359;
t397 = qJD(5) * t360;
t396 = pkin(2) * qJD(2) * t363;
t390 = pkin(3) * t360 + qJ(4) * t359;
t326 = t390 * t364;
t395 = -t326 + t400;
t394 = t367 * t396;
t393 = rSges(3,1) * t366 - rSges(3,2) * t363;
t392 = rSges(4,1) * t360 - rSges(4,2) * t359;
t391 = -rSges(5,2) * t360 + rSges(5,3) * t359;
t389 = Icges(3,1) * t366 - t413;
t387 = -Icges(3,2) * t363 + t412;
t384 = Icges(3,5) * t366 - Icges(3,6) * t363;
t318 = -Icges(3,6) * t367 + t364 * t387;
t320 = -Icges(3,5) * t367 + t364 * t389;
t380 = t318 * t363 - t320 * t366;
t319 = Icges(3,6) * t364 + t367 * t387;
t321 = Icges(3,5) * t364 + t367 * t389;
t379 = -t319 * t363 + t321 * t366;
t347 = Icges(3,2) * t366 + t413;
t348 = Icges(3,1) * t363 + t412;
t378 = -t347 * t363 + t348 * t366;
t377 = -qJD(4) * t360 + t344 * t326 + t401;
t343 = qJD(1) * (pkin(1) * t367 + pkin(6) * t364);
t376 = qJD(1) * t298 - t364 * t396 + t343;
t340 = pkin(3) * t359 - qJ(4) * t360;
t375 = t345 * t340 + t367 * t398 - t394;
t327 = t390 * t367;
t372 = qJD(1) * t327 + t364 * t398 + t376;
t355 = qJD(5) * t359 + qJD(1);
t351 = rSges(2,1) * t367 - rSges(2,2) * t364;
t350 = rSges(2,1) * t364 + rSges(2,2) * t367;
t349 = rSges(3,1) * t363 + rSges(3,2) * t366;
t346 = Icges(3,5) * t363 + Icges(3,6) * t366;
t342 = rSges(4,1) * t359 + rSges(4,2) * t360;
t341 = -rSges(5,2) * t359 - rSges(5,3) * t360;
t333 = -pkin(4) * t367 + pkin(8) * t407;
t332 = pkin(4) * t364 + pkin(8) * t406;
t331 = t359 * t405 - t402;
t330 = t359 * t403 + t404;
t329 = t359 * t404 + t403;
t328 = t359 * t402 - t405;
t325 = rSges(3,3) * t364 + t367 * t393;
t324 = -rSges(3,3) * t367 + t364 * t393;
t323 = t364 * t397 + t345;
t322 = t367 * t397 + t344;
t317 = Icges(3,3) * t364 + t367 * t384;
t316 = -Icges(3,3) * t367 + t364 * t384;
t314 = -rSges(5,1) * t367 + t364 * t391;
t313 = rSges(5,1) * t364 + t367 * t391;
t312 = rSges(4,3) * t364 + t367 * t392;
t311 = -rSges(4,3) * t367 + t364 * t392;
t295 = rSges(6,3) * t359 + (-rSges(6,1) * t362 - rSges(6,2) * t365) * t360;
t294 = Icges(6,5) * t359 + (-Icges(6,1) * t362 - Icges(6,4) * t365) * t360;
t293 = Icges(6,6) * t359 + (-Icges(6,4) * t362 - Icges(6,2) * t365) * t360;
t292 = Icges(6,3) * t359 + (-Icges(6,5) * t362 - Icges(6,6) * t365) * t360;
t287 = qJD(1) * t325 - t349 * t358 + t343;
t286 = -t349 * t399 + (-t324 - t354) * qJD(1);
t285 = rSges(6,1) * t331 + rSges(6,2) * t330 + rSges(6,3) * t407;
t284 = rSges(6,1) * t329 + rSges(6,2) * t328 + rSges(6,3) * t406;
t283 = Icges(6,1) * t331 + Icges(6,4) * t330 + Icges(6,5) * t407;
t282 = Icges(6,1) * t329 + Icges(6,4) * t328 + Icges(6,5) * t406;
t281 = Icges(6,4) * t331 + Icges(6,2) * t330 + Icges(6,6) * t407;
t280 = Icges(6,4) * t329 + Icges(6,2) * t328 + Icges(6,6) * t406;
t279 = Icges(6,5) * t331 + Icges(6,6) * t330 + Icges(6,3) * t407;
t278 = Icges(6,5) * t329 + Icges(6,6) * t328 + Icges(6,3) * t406;
t277 = (t324 * t364 + t325 * t367) * qJD(2);
t276 = qJD(1) * t312 - t342 * t344 + t376;
t275 = -t394 + t342 * t345 + (-t311 + t400) * qJD(1);
t274 = t311 * t344 - t312 * t345 + t401;
t273 = qJD(1) * t313 + (-t340 - t341) * t344 + t372;
t272 = t341 * t345 + (-t314 + t395) * qJD(1) + t375;
t271 = t314 * t344 + (-t313 - t327) * t345 + t377;
t270 = qJD(1) * t332 + t284 * t355 - t295 * t322 + (-t340 - t417) * t344 + t372;
t269 = t345 * t417 - t285 * t355 + t295 * t323 + (-t333 + t395) * qJD(1) + t375;
t268 = -t284 * t323 + t285 * t322 + t333 * t344 + (-t327 - t332) * t345 + t377;
t1 = m(3) * (t277 ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1 + ((t364 * t346 + t367 * t378) * qJD(1) + (t364 ^ 2 * t317 + (t380 * t367 + (-t316 + t379) * t364) * t367) * qJD(2)) * t358 / 0.2e1 - ((-t367 * t346 + t364 * t378) * qJD(1) + (t367 ^ 2 * t316 + (t379 * t364 + (-t317 + t380) * t367) * t364) * qJD(2)) * t399 / 0.2e1 + m(4) * (t274 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + m(5) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + m(6) * (t268 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + t322 * ((t278 * t406 + t328 * t280 + t329 * t282) * t322 + (t279 * t406 + t281 * t328 + t283 * t329) * t323 + (t292 * t406 + t293 * t328 + t294 * t329) * t355) / 0.2e1 + t323 * ((t278 * t407 + t280 * t330 + t282 * t331) * t322 + (t279 * t407 + t330 * t281 + t331 * t283) * t323 + (t292 * t407 + t293 * t330 + t294 * t331) * t355) / 0.2e1 + t355 * ((t278 * t322 + t279 * t323 + t292 * t355) * t359 + ((-t280 * t365 - t282 * t362) * t322 + (-t281 * t365 - t283 * t362) * t323 + (-t293 * t365 - t294 * t362) * t355) * t360) / 0.2e1 + (t421 * t364 + t422 * t367) * t344 / 0.2e1 + (t422 * t364 - t421 * t367) * t345 / 0.2e1 + (m(2) * (t350 ^ 2 + t351 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t319 * t366 + t321 * t363) * t364 - (t318 * t366 + t320 * t363) * t367) * qJD(2) + (t427 * t359 + t429 * t360) * t345 + (t426 * t359 + t428 * t360) * t344 + (t366 * t347 + t363 * t348 + t424 * t359 - t425 * t360) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
