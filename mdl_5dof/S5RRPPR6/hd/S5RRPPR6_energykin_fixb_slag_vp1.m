% Calculate kinetic energy for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:33
% EndTime: 2019-12-31 19:31:35
% DurationCPUTime: 1.97s
% Computational Cost: add. (1204->254), mult. (1469->387), div. (0->0), fcn. (1422->10), ass. (0->133)
t435 = Icges(3,3) + Icges(4,3);
t364 = qJ(2) + pkin(8);
t359 = sin(t364);
t361 = cos(t364);
t369 = sin(qJ(2));
t371 = cos(qJ(2));
t434 = Icges(3,5) * t371 + Icges(4,5) * t361 - Icges(3,6) * t369 - Icges(4,6) * t359;
t370 = sin(qJ(1));
t372 = cos(qJ(1));
t433 = t434 * t370 - t435 * t372;
t432 = t435 * t370 + t434 * t372;
t431 = Icges(3,5) * t369 + Icges(4,5) * t359 + Icges(3,6) * t371 + Icges(4,6) * t361;
t417 = Icges(4,4) * t359;
t341 = Icges(4,2) * t361 + t417;
t416 = Icges(4,4) * t361;
t342 = Icges(4,1) * t359 + t416;
t419 = Icges(3,4) * t369;
t347 = Icges(3,2) * t371 + t419;
t418 = Icges(3,4) * t371;
t348 = Icges(3,1) * t369 + t418;
t430 = -t341 * t359 + t342 * t361 - t347 * t369 + t348 * t371;
t387 = -Icges(4,2) * t359 + t416;
t312 = Icges(4,6) * t370 + t387 * t372;
t389 = Icges(4,1) * t361 - t417;
t314 = Icges(4,5) * t370 + t389 * t372;
t388 = -Icges(3,2) * t369 + t418;
t327 = Icges(3,6) * t370 + t388 * t372;
t390 = Icges(3,1) * t371 - t419;
t329 = Icges(3,5) * t370 + t390 * t372;
t429 = -t312 * t359 + t314 * t361 - t327 * t369 + t329 * t371;
t311 = -Icges(4,6) * t372 + t387 * t370;
t313 = -Icges(4,5) * t372 + t389 * t370;
t326 = -Icges(3,6) * t372 + t388 * t370;
t328 = -Icges(3,5) * t372 + t390 * t370;
t428 = t311 * t359 - t313 * t361 + t326 * t369 - t328 * t371;
t424 = pkin(2) * t369;
t422 = pkin(2) * t371;
t366 = cos(pkin(9));
t421 = pkin(4) * t366;
t415 = t359 * t370;
t414 = t359 * t372;
t413 = t361 * t370;
t412 = t361 * t372;
t365 = sin(pkin(9));
t411 = t365 * t370;
t410 = t365 * t372;
t409 = t366 * t370;
t408 = t366 * t372;
t307 = -qJ(3) * t372 + t422 * t370;
t308 = qJ(3) * t370 + t422 * t372;
t402 = qJD(2) * t372;
t403 = qJD(2) * t370;
t406 = t307 * t403 + t308 * t402;
t354 = pkin(1) * t370 - pkin(6) * t372;
t405 = -t307 - t354;
t362 = qJD(3) * t370;
t401 = qJD(4) * t359;
t404 = t372 * t401 + t362;
t400 = qJD(5) * t359;
t392 = pkin(3) * t361 + qJ(4) * t359;
t330 = t392 * t370;
t399 = -t330 + t405;
t396 = -pkin(3) * t359 + qJ(4) * t361 - t424;
t345 = qJD(1) * (pkin(1) * t372 + pkin(6) * t370);
t395 = qJD(1) * t308 - qJD(3) * t372 + t345;
t394 = rSges(3,1) * t371 - rSges(3,2) * t369;
t393 = rSges(4,1) * t361 - rSges(4,2) * t359;
t391 = qJD(2) * (-rSges(4,1) * t359 - rSges(4,2) * t361 - t424);
t378 = qJD(2) * (pkin(7) * t361 - t421 * t359 + t396);
t377 = qJD(2) * (rSges(5,3) * t361 - (rSges(5,1) * t366 - rSges(5,2) * t365) * t359 + t396);
t331 = t392 * t372;
t376 = qJD(1) * t331 + t370 * t401 + t395;
t375 = -qJD(4) * t361 + t330 * t403 + t331 * t402 + t406;
t374 = pkin(7) * t359 + t421 * t361;
t363 = pkin(9) + qJ(5);
t360 = cos(t363);
t358 = sin(t363);
t355 = -qJD(5) * t361 + qJD(1);
t353 = rSges(2,1) * t372 - rSges(2,2) * t370;
t352 = rSges(2,1) * t370 + rSges(2,2) * t372;
t351 = rSges(3,1) * t369 + rSges(3,2) * t371;
t339 = t370 * t400 - t402;
t338 = t372 * t400 + t403;
t337 = t361 * t408 + t411;
t336 = -t361 * t410 + t409;
t335 = t361 * t409 - t410;
t334 = -t361 * t411 - t408;
t333 = rSges(3,3) * t370 + t394 * t372;
t332 = -rSges(3,3) * t372 + t394 * t370;
t322 = t358 * t370 + t360 * t412;
t321 = -t358 * t412 + t360 * t370;
t320 = -t358 * t372 + t360 * t413;
t319 = -t358 * t413 - t360 * t372;
t318 = rSges(4,3) * t370 + t393 * t372;
t317 = -rSges(4,3) * t372 + t393 * t370;
t304 = -Icges(5,5) * t361 + (Icges(5,1) * t366 - Icges(5,4) * t365) * t359;
t303 = -Icges(5,6) * t361 + (Icges(5,4) * t366 - Icges(5,2) * t365) * t359;
t302 = -Icges(5,3) * t361 + (Icges(5,5) * t366 - Icges(5,6) * t365) * t359;
t299 = -rSges(6,3) * t361 + (rSges(6,1) * t360 - rSges(6,2) * t358) * t359;
t298 = -Icges(6,5) * t361 + (Icges(6,1) * t360 - Icges(6,4) * t358) * t359;
t297 = -Icges(6,6) * t361 + (Icges(6,4) * t360 - Icges(6,2) * t358) * t359;
t296 = -Icges(6,3) * t361 + (Icges(6,5) * t360 - Icges(6,6) * t358) * t359;
t294 = qJD(1) * t333 - t351 * t403 + t345;
t293 = -t351 * t402 + (-t332 - t354) * qJD(1);
t292 = (t332 * t370 + t333 * t372) * qJD(2);
t291 = rSges(5,1) * t337 + rSges(5,2) * t336 + rSges(5,3) * t414;
t290 = rSges(5,1) * t335 + rSges(5,2) * t334 + rSges(5,3) * t415;
t289 = Icges(5,1) * t337 + Icges(5,4) * t336 + Icges(5,5) * t414;
t288 = Icges(5,1) * t335 + Icges(5,4) * t334 + Icges(5,5) * t415;
t287 = Icges(5,4) * t337 + Icges(5,2) * t336 + Icges(5,6) * t414;
t286 = Icges(5,4) * t335 + Icges(5,2) * t334 + Icges(5,6) * t415;
t285 = Icges(5,5) * t337 + Icges(5,6) * t336 + Icges(5,3) * t414;
t284 = Icges(5,5) * t335 + Icges(5,6) * t334 + Icges(5,3) * t415;
t283 = pkin(4) * t411 + t374 * t372;
t282 = -pkin(4) * t410 + t374 * t370;
t281 = rSges(6,1) * t322 + rSges(6,2) * t321 + rSges(6,3) * t414;
t280 = rSges(6,1) * t320 + rSges(6,2) * t319 + rSges(6,3) * t415;
t279 = Icges(6,1) * t322 + Icges(6,4) * t321 + Icges(6,5) * t414;
t278 = Icges(6,1) * t320 + Icges(6,4) * t319 + Icges(6,5) * t415;
t277 = Icges(6,4) * t322 + Icges(6,2) * t321 + Icges(6,6) * t414;
t276 = Icges(6,4) * t320 + Icges(6,2) * t319 + Icges(6,6) * t415;
t275 = Icges(6,5) * t322 + Icges(6,6) * t321 + Icges(6,3) * t414;
t274 = Icges(6,5) * t320 + Icges(6,6) * t319 + Icges(6,3) * t415;
t273 = qJD(1) * t318 + t370 * t391 + t395;
t272 = t362 + t372 * t391 + (-t317 + t405) * qJD(1);
t271 = (t317 * t370 + t318 * t372) * qJD(2) + t406;
t270 = qJD(1) * t291 + t370 * t377 + t376;
t269 = t372 * t377 + (-t290 + t399) * qJD(1) + t404;
t268 = (t290 * t370 + t291 * t372) * qJD(2) + t375;
t267 = qJD(1) * t283 + t281 * t355 - t299 * t338 + t370 * t378 + t376;
t266 = -t280 * t355 + t299 * t339 + t372 * t378 + (-t282 + t399) * qJD(1) + t404;
t265 = t280 * t338 - t281 * t339 + (t282 * t370 + t283 * t372) * qJD(2) + t375;
t1 = m(3) * (t292 ^ 2 + t293 ^ 2 + t294 ^ 2) / 0.2e1 + m(4) * (t271 ^ 2 + t272 ^ 2 + t273 ^ 2) / 0.2e1 + m(5) * (t268 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + m(6) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + t338 * ((t275 * t414 + t321 * t277 + t322 * t279) * t338 + (t274 * t414 + t276 * t321 + t278 * t322) * t339 + (t296 * t414 + t297 * t321 + t298 * t322) * t355) / 0.2e1 + t339 * ((t275 * t415 + t277 * t319 + t279 * t320) * t338 + (t274 * t415 + t319 * t276 + t320 * t278) * t339 + (t296 * t415 + t297 * t319 + t298 * t320) * t355) / 0.2e1 + t355 * ((-t274 * t339 - t275 * t338 - t296 * t355) * t361 + ((-t277 * t358 + t279 * t360) * t338 + (-t276 * t358 + t278 * t360) * t339 + (-t297 * t358 + t298 * t360) * t355) * t359) / 0.2e1 + (m(2) * (t352 ^ 2 + t353 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t326 * t371 - t328 * t369 + (-t311 + t284) * t361 + (t286 * t365 - t288 * t366 - t313) * t359) * t372 + (t327 * t371 + t329 * t369 + (t312 - t285) * t361 + (-t287 * t365 + t289 * t366 + t314) * t359) * t370) * qJD(2) + (t371 * t347 + t369 * t348 + (t341 - t302) * t361 + (-t303 * t365 + t304 * t366 + t342) * t359) * qJD(1)) * qJD(1) / 0.2e1 + (((-t284 * t414 - t286 * t336 - t288 * t337 + t428 * t372) * t372 + (t285 * t414 + t336 * t287 + t337 * t289 + (t429 - t433) * t372 + t432 * t370) * t370) * qJD(2) + (t302 * t414 + t303 * t336 + t304 * t337 + t431 * t370 + t430 * t372) * qJD(1)) * t403 / 0.2e1 - (((t285 * t415 + t287 * t334 + t289 * t335 + t429 * t370) * t370 + (-t284 * t415 - t334 * t286 - t335 * t288 + (t428 - t432) * t370 + t433 * t372) * t372) * qJD(2) + (t302 * t415 + t303 * t334 + t304 * t335 + t430 * t370 - t431 * t372) * qJD(1)) * t402 / 0.2e1;
T = t1;
