% Calculate kinetic energy for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR10_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR10_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:00
% EndTime: 2019-12-31 19:43:02
% DurationCPUTime: 1.86s
% Computational Cost: add. (834->238), mult. (2056->360), div. (0->0), fcn. (2229->8), ass. (0->118)
t412 = Icges(4,1) + Icges(5,1);
t411 = Icges(4,4) - Icges(5,5);
t410 = Icges(5,4) + Icges(4,5);
t409 = Icges(4,2) + Icges(5,3);
t408 = -Icges(5,6) + Icges(4,6);
t407 = -Icges(4,3) - Icges(5,2);
t356 = sin(pkin(8));
t357 = cos(pkin(8));
t363 = cos(qJ(1));
t360 = sin(qJ(1));
t362 = cos(qJ(2));
t391 = t362 * t360;
t332 = t356 * t391 + t357 * t363;
t333 = -t356 * t363 + t357 * t391;
t359 = sin(qJ(2));
t393 = t359 * t360;
t406 = -t409 * t332 + t411 * t333 + t408 * t393;
t390 = t362 * t363;
t334 = t356 * t390 - t360 * t357;
t335 = t356 * t360 + t357 * t390;
t392 = t359 * t363;
t405 = t409 * t334 - t411 * t335 - t408 * t392;
t404 = t408 * t332 - t410 * t333 + t407 * t393;
t403 = -t408 * t334 + t410 * t335 - t407 * t392;
t402 = t411 * t332 - t412 * t333 - t410 * t393;
t401 = -t411 * t334 + t412 * t335 + t410 * t392;
t400 = t408 * t362 + (t409 * t356 - t411 * t357) * t359;
t399 = t407 * t362 + (-t408 * t356 + t410 * t357) * t359;
t398 = -t410 * t362 + (-t411 * t356 + t412 * t357) * t359;
t395 = Icges(3,4) * t359;
t394 = Icges(3,4) * t362;
t384 = qJD(3) * t359;
t353 = t363 * t384;
t389 = qJD(4) * t334 + t353;
t346 = pkin(2) * t359 - qJ(3) * t362;
t388 = -(pkin(3) * t357 + qJ(4) * t356) * t359 - t346;
t373 = pkin(2) * t362 + qJ(3) * t359;
t337 = t373 * t360;
t350 = pkin(1) * t360 - pkin(6) * t363;
t387 = -t337 - t350;
t386 = qJD(2) * t360;
t385 = qJD(2) * t363;
t383 = qJD(5) * t359;
t305 = pkin(3) * t333 + qJ(4) * t332;
t382 = -t305 + t387;
t338 = t373 * t363;
t342 = qJD(1) * (pkin(1) * t363 + pkin(6) * t360);
t381 = qJD(1) * t338 + t360 * t384 + t342;
t378 = qJD(2) * (rSges(4,3) * t362 - (rSges(4,1) * t357 - rSges(4,2) * t356) * t359 - t346);
t377 = qJD(2) * (rSges(5,2) * t362 - (rSges(5,1) * t357 + rSges(5,3) * t356) * t359 + t388);
t376 = qJD(2) * (-pkin(4) * t357 * t359 - pkin(7) * t362 + t388);
t375 = -qJD(3) * t362 + t337 * t386 + t338 * t385;
t374 = rSges(3,1) * t362 - rSges(3,2) * t359;
t306 = pkin(3) * t335 + qJ(4) * t334;
t372 = qJD(1) * t306 + qJD(4) * t332 + t381;
t371 = Icges(3,1) * t362 - t395;
t370 = -Icges(3,2) * t359 + t394;
t369 = Icges(3,5) * t362 - Icges(3,6) * t359;
t319 = -Icges(3,6) * t363 + t360 * t370;
t321 = -Icges(3,5) * t363 + t360 * t371;
t368 = t319 * t359 - t321 * t362;
t320 = Icges(3,6) * t360 + t363 * t370;
t322 = Icges(3,5) * t360 + t363 * t371;
t367 = -t320 * t359 + t322 * t362;
t344 = Icges(3,2) * t362 + t395;
t345 = Icges(3,1) * t359 + t394;
t366 = -t344 * t359 + t345 * t362;
t365 = qJD(4) * t359 * t356 + t305 * t386 + t306 * t385 + t375;
t361 = cos(qJ(5));
t358 = sin(qJ(5));
t354 = qJD(5) * t362 + qJD(1);
t349 = rSges(2,1) * t363 - rSges(2,2) * t360;
t348 = rSges(2,1) * t360 + rSges(2,2) * t363;
t347 = rSges(3,1) * t359 + rSges(3,2) * t362;
t343 = Icges(3,5) * t359 + Icges(3,6) * t362;
t340 = -t360 * t383 - t385;
t339 = -t363 * t383 + t386;
t328 = (t356 * t358 + t357 * t361) * t359;
t327 = (t356 * t361 - t357 * t358) * t359;
t326 = rSges(3,3) * t360 + t363 * t374;
t325 = -rSges(3,3) * t363 + t360 * t374;
t318 = Icges(3,3) * t360 + t363 * t369;
t317 = -Icges(3,3) * t363 + t360 * t369;
t308 = pkin(4) * t335 - pkin(7) * t392;
t307 = pkin(4) * t333 - pkin(7) * t393;
t303 = t334 * t358 + t335 * t361;
t302 = t334 * t361 - t335 * t358;
t301 = t332 * t358 + t333 * t361;
t300 = t332 * t361 - t333 * t358;
t297 = rSges(4,1) * t335 - rSges(4,2) * t334 + rSges(4,3) * t392;
t296 = rSges(5,1) * t335 + rSges(5,2) * t392 + rSges(5,3) * t334;
t295 = rSges(4,1) * t333 - rSges(4,2) * t332 + rSges(4,3) * t393;
t294 = rSges(5,1) * t333 + rSges(5,2) * t393 + rSges(5,3) * t332;
t281 = rSges(6,1) * t328 + rSges(6,2) * t327 + rSges(6,3) * t362;
t280 = Icges(6,1) * t328 + Icges(6,4) * t327 + Icges(6,5) * t362;
t279 = Icges(6,4) * t328 + Icges(6,2) * t327 + Icges(6,6) * t362;
t278 = Icges(6,5) * t328 + Icges(6,6) * t327 + Icges(6,3) * t362;
t277 = qJD(1) * t326 - t347 * t386 + t342;
t276 = -t347 * t385 + (-t325 - t350) * qJD(1);
t275 = (t325 * t360 + t326 * t363) * qJD(2);
t274 = rSges(6,1) * t303 + rSges(6,2) * t302 - rSges(6,3) * t392;
t273 = rSges(6,1) * t301 + rSges(6,2) * t300 - rSges(6,3) * t393;
t272 = Icges(6,1) * t303 + Icges(6,4) * t302 - Icges(6,5) * t392;
t271 = Icges(6,1) * t301 + Icges(6,4) * t300 - Icges(6,5) * t393;
t270 = Icges(6,4) * t303 + Icges(6,2) * t302 - Icges(6,6) * t392;
t269 = Icges(6,4) * t301 + Icges(6,2) * t300 - Icges(6,6) * t393;
t268 = Icges(6,5) * t303 + Icges(6,6) * t302 - Icges(6,3) * t392;
t267 = Icges(6,5) * t301 + Icges(6,6) * t300 - Icges(6,3) * t393;
t266 = qJD(1) * t297 + t360 * t378 + t381;
t265 = t353 + t363 * t378 + (-t295 + t387) * qJD(1);
t264 = (t295 * t360 + t297 * t363) * qJD(2) + t375;
t263 = qJD(1) * t296 + t360 * t377 + t372;
t262 = t363 * t377 + (-t294 + t382) * qJD(1) + t389;
t261 = (t294 * t360 + t296 * t363) * qJD(2) + t365;
t260 = qJD(1) * t308 + t274 * t354 - t281 * t339 + t360 * t376 + t372;
t259 = -t273 * t354 + t281 * t340 + t363 * t376 + (-t307 + t382) * qJD(1) + t389;
t258 = t273 * t339 - t274 * t340 + (t307 * t360 + t308 * t363) * qJD(2) + t365;
t1 = m(3) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + m(4) * (t264 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(5) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(6) * (t258 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + t339 * ((-t268 * t392 + t302 * t270 + t303 * t272) * t339 + (-t267 * t392 + t269 * t302 + t271 * t303) * t340 + (-t278 * t392 + t279 * t302 + t280 * t303) * t354) / 0.2e1 + t340 * ((-t268 * t393 + t270 * t300 + t272 * t301) * t339 + (-t267 * t393 + t300 * t269 + t301 * t271) * t340 + (-t278 * t393 + t279 * t300 + t280 * t301) * t354) / 0.2e1 + t354 * ((t268 * t362 + t270 * t327 + t272 * t328) * t339 + (t267 * t362 + t269 * t327 + t271 * t328) * t340 + (t362 * t278 + t327 * t279 + t328 * t280) * t354) / 0.2e1 + (m(2) * (t348 ^ 2 + t349 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((((-t319 - t404) * t363 + (t320 - t403) * t360) * t362 + ((t406 * t356 + t402 * t357 - t321) * t363 + (t405 * t356 + t401 * t357 + t322) * t360) * t359) * qJD(2) + ((t344 - t399) * t362 + (t400 * t356 + t398 * t357 + t345) * t359) * qJD(1)) * qJD(1) / 0.2e1 + (((t406 * t334 + t402 * t335 + t368 * t363 + t404 * t392) * t363 + ((-t317 + t367) * t363 + t318 * t360 + t403 * t392 + t401 * t335 + t405 * t334) * t360) * qJD(2) + (t400 * t334 + t398 * t335 + t360 * t343 + t366 * t363 + t399 * t392) * qJD(1)) * t386 / 0.2e1 - (((t317 * t363 + t406 * t332 + t402 * t333 + t404 * t393) * t363 + (t367 * t360 + (-t318 + t368) * t363 + t403 * t393 + t401 * t333 + t405 * t332) * t360) * qJD(2) + (t400 * t332 + t398 * t333 - t363 * t343 + t366 * t360 + t399 * t393) * qJD(1)) * t385 / 0.2e1;
T = t1;
