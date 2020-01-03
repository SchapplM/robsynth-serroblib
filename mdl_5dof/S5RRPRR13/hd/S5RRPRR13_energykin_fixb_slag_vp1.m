% Calculate kinetic energy for
% S5RRPRR13
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR13_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR13_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:06
% EndTime: 2019-12-31 20:32:08
% DurationCPUTime: 1.99s
% Computational Cost: add. (1252->283), mult. (1818->453), div. (0->0), fcn. (1827->10), ass. (0->137)
t372 = cos(pkin(9));
t415 = t372 * pkin(3);
t374 = sin(qJ(2));
t414 = Icges(3,4) * t374;
t376 = cos(qJ(2));
t413 = Icges(3,4) * t376;
t371 = sin(pkin(9));
t375 = sin(qJ(1));
t412 = t371 * t375;
t377 = cos(qJ(1));
t411 = t371 * t377;
t410 = t374 * t375;
t409 = t374 * t377;
t408 = t375 * t376;
t407 = t376 * t377;
t390 = pkin(2) * t376 + qJ(3) * t374;
t340 = t390 * t375;
t355 = pkin(1) * t375 - pkin(6) * t377;
t405 = -t340 - t355;
t370 = pkin(9) + qJ(4);
t365 = cos(t370);
t404 = pkin(4) * t365;
t367 = qJD(2) * t375;
t400 = qJD(4) * t374;
t342 = t377 * t400 + t367;
t402 = qJD(2) * t377;
t401 = qJD(3) * t374;
t399 = qJD(5) * t374;
t341 = t390 * t377;
t347 = qJD(1) * (pkin(1) * t377 + pkin(6) * t375);
t398 = qJD(1) * t341 + t375 * t401 + t347;
t364 = sin(t370);
t395 = pkin(4) * t364;
t351 = pkin(2) * t374 - qJ(3) * t376;
t394 = qJD(2) * (pkin(7) * t376 - t374 * t415 - t351);
t393 = qJD(2) * (rSges(4,3) * t376 - (rSges(4,1) * t372 - rSges(4,2) * t371) * t374 - t351);
t343 = t375 * t400 - t402;
t392 = -qJD(3) * t376 + t340 * t367 + t341 * t402;
t391 = rSges(3,1) * t376 - rSges(3,2) * t374;
t389 = Icges(3,1) * t376 - t414;
t388 = -Icges(3,2) * t374 + t413;
t387 = Icges(3,5) * t376 - Icges(3,6) * t374;
t327 = -Icges(3,6) * t377 + t375 * t388;
t329 = -Icges(3,5) * t377 + t375 * t389;
t386 = t327 * t374 - t329 * t376;
t328 = Icges(3,6) * t375 + t377 * t388;
t330 = Icges(3,5) * t375 + t377 * t389;
t385 = -t328 * t374 + t330 * t376;
t349 = Icges(3,2) * t376 + t414;
t350 = Icges(3,1) * t374 + t413;
t384 = -t349 * t374 + t350 * t376;
t382 = pkin(7) * t374 + t376 * t415;
t299 = -pkin(3) * t411 + t375 * t382;
t300 = pkin(3) * t412 + t377 * t382;
t383 = t299 * t367 + t300 * t402 + t392;
t381 = pkin(8) * t374 + t376 * t404;
t380 = qJD(1) * t300 + t375 * t394 + t398;
t359 = t377 * t401;
t379 = t359 + (-t299 + t405) * qJD(1) + t377 * t394;
t366 = qJ(5) + t370;
t362 = cos(t366);
t361 = sin(t366);
t360 = -qJD(4) * t376 + qJD(1);
t354 = rSges(2,1) * t377 - rSges(2,2) * t375;
t353 = rSges(2,1) * t375 + rSges(2,2) * t377;
t352 = rSges(3,1) * t374 + rSges(3,2) * t376;
t348 = Icges(3,5) * t374 + Icges(3,6) * t376;
t346 = qJD(1) + (-qJD(4) - qJD(5)) * t376;
t339 = t372 * t407 + t412;
t338 = -t371 * t407 + t372 * t375;
t337 = t372 * t408 - t411;
t336 = -t371 * t408 - t372 * t377;
t334 = rSges(3,3) * t375 + t377 * t391;
t333 = -rSges(3,3) * t377 + t375 * t391;
t326 = Icges(3,3) * t375 + t377 * t387;
t325 = -Icges(3,3) * t377 + t375 * t387;
t324 = t364 * t375 + t365 * t407;
t323 = -t364 * t407 + t365 * t375;
t322 = -t364 * t377 + t365 * t408;
t321 = -t364 * t408 - t365 * t377;
t319 = t375 * t399 + t343;
t318 = t377 * t399 + t342;
t317 = -Icges(4,5) * t376 + (Icges(4,1) * t372 - Icges(4,4) * t371) * t374;
t316 = -Icges(4,6) * t376 + (Icges(4,4) * t372 - Icges(4,2) * t371) * t374;
t315 = -Icges(4,3) * t376 + (Icges(4,5) * t372 - Icges(4,6) * t371) * t374;
t314 = t361 * t375 + t362 * t407;
t313 = -t361 * t407 + t362 * t375;
t312 = -t361 * t377 + t362 * t408;
t311 = -t361 * t408 - t362 * t377;
t310 = -rSges(5,3) * t376 + (rSges(5,1) * t365 - rSges(5,2) * t364) * t374;
t309 = -Icges(5,5) * t376 + (Icges(5,1) * t365 - Icges(5,4) * t364) * t374;
t308 = -Icges(5,6) * t376 + (Icges(5,4) * t365 - Icges(5,2) * t364) * t374;
t307 = -Icges(5,3) * t376 + (Icges(5,5) * t365 - Icges(5,6) * t364) * t374;
t305 = -rSges(6,3) * t376 + (rSges(6,1) * t362 - rSges(6,2) * t361) * t374;
t304 = -Icges(6,5) * t376 + (Icges(6,1) * t362 - Icges(6,4) * t361) * t374;
t303 = -Icges(6,6) * t376 + (Icges(6,4) * t362 - Icges(6,2) * t361) * t374;
t302 = -Icges(6,3) * t376 + (Icges(6,5) * t362 - Icges(6,6) * t361) * t374;
t301 = -pkin(8) * t376 + t374 * t404;
t298 = rSges(4,1) * t339 + rSges(4,2) * t338 + rSges(4,3) * t409;
t297 = rSges(4,1) * t337 + rSges(4,2) * t336 + rSges(4,3) * t410;
t296 = Icges(4,1) * t339 + Icges(4,4) * t338 + Icges(4,5) * t409;
t295 = Icges(4,1) * t337 + Icges(4,4) * t336 + Icges(4,5) * t410;
t294 = Icges(4,4) * t339 + Icges(4,2) * t338 + Icges(4,6) * t409;
t293 = Icges(4,4) * t337 + Icges(4,2) * t336 + Icges(4,6) * t410;
t292 = Icges(4,5) * t339 + Icges(4,6) * t338 + Icges(4,3) * t409;
t291 = Icges(4,5) * t337 + Icges(4,6) * t336 + Icges(4,3) * t410;
t289 = qJD(1) * t334 - t352 * t367 + t347;
t288 = -t352 * t402 + (-t333 - t355) * qJD(1);
t285 = (t333 * t375 + t334 * t377) * qJD(2);
t284 = rSges(5,1) * t324 + rSges(5,2) * t323 + rSges(5,3) * t409;
t283 = rSges(5,1) * t322 + rSges(5,2) * t321 + rSges(5,3) * t410;
t282 = Icges(5,1) * t324 + Icges(5,4) * t323 + Icges(5,5) * t409;
t281 = Icges(5,1) * t322 + Icges(5,4) * t321 + Icges(5,5) * t410;
t280 = Icges(5,4) * t324 + Icges(5,2) * t323 + Icges(5,6) * t409;
t279 = Icges(5,4) * t322 + Icges(5,2) * t321 + Icges(5,6) * t410;
t278 = Icges(5,5) * t324 + Icges(5,6) * t323 + Icges(5,3) * t409;
t277 = Icges(5,5) * t322 + Icges(5,6) * t321 + Icges(5,3) * t410;
t276 = rSges(6,1) * t314 + rSges(6,2) * t313 + rSges(6,3) * t409;
t275 = rSges(6,1) * t312 + rSges(6,2) * t311 + rSges(6,3) * t410;
t274 = Icges(6,1) * t314 + Icges(6,4) * t313 + Icges(6,5) * t409;
t273 = Icges(6,1) * t312 + Icges(6,4) * t311 + Icges(6,5) * t410;
t272 = Icges(6,4) * t314 + Icges(6,2) * t313 + Icges(6,6) * t409;
t271 = Icges(6,4) * t312 + Icges(6,2) * t311 + Icges(6,6) * t410;
t270 = Icges(6,5) * t314 + Icges(6,6) * t313 + Icges(6,3) * t409;
t269 = Icges(6,5) * t312 + Icges(6,6) * t311 + Icges(6,3) * t410;
t268 = t375 * t395 + t377 * t381;
t267 = t375 * t381 - t377 * t395;
t266 = qJD(1) * t298 + t375 * t393 + t398;
t265 = t359 + t377 * t393 + (-t297 + t405) * qJD(1);
t264 = (t297 * t375 + t298 * t377) * qJD(2) + t392;
t263 = t284 * t360 - t310 * t342 + t380;
t262 = -t283 * t360 + t310 * t343 + t379;
t261 = t283 * t342 - t284 * t343 + t383;
t260 = t268 * t360 + t276 * t346 - t301 * t342 - t305 * t318 + t380;
t259 = -t267 * t360 - t275 * t346 + t301 * t343 + t305 * t319 + t379;
t258 = t267 * t342 - t268 * t343 + t275 * t318 - t276 * t319 + t383;
t1 = m(3) * (t285 ^ 2 + t288 ^ 2 + t289 ^ 2) / 0.2e1 + m(4) * (t264 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(5) * (t261 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + t342 * ((t278 * t409 + t280 * t323 + t282 * t324) * t342 + (t277 * t409 + t279 * t323 + t281 * t324) * t343 + (t307 * t409 + t308 * t323 + t309 * t324) * t360) / 0.2e1 + t343 * ((t278 * t410 + t280 * t321 + t282 * t322) * t342 + (t277 * t410 + t279 * t321 + t281 * t322) * t343 + (t307 * t410 + t308 * t321 + t309 * t322) * t360) / 0.2e1 + t360 * ((-t277 * t343 - t278 * t342 - t307 * t360) * t376 + ((-t280 * t364 + t282 * t365) * t342 + (-t279 * t364 + t281 * t365) * t343 + (-t308 * t364 + t309 * t365) * t360) * t374) / 0.2e1 + m(6) * (t258 ^ 2 + t259 ^ 2 + t260 ^ 2) / 0.2e1 + t318 * ((t270 * t409 + t272 * t313 + t274 * t314) * t318 + (t269 * t409 + t271 * t313 + t273 * t314) * t319 + (t302 * t409 + t303 * t313 + t304 * t314) * t346) / 0.2e1 + t319 * ((t270 * t410 + t272 * t311 + t274 * t312) * t318 + (t269 * t410 + t271 * t311 + t273 * t312) * t319 + (t302 * t410 + t303 * t311 + t304 * t312) * t346) / 0.2e1 + t346 * ((-t269 * t319 - t270 * t318 - t302 * t346) * t376 + ((-t272 * t361 + t274 * t362) * t318 + (-t271 * t361 + t273 * t362) * t319 + (-t303 * t361 + t304 * t362) * t346) * t374) / 0.2e1 + (m(2) * (t353 ^ 2 + t354 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t328 * t376 + t330 * t374) * t375 - (t327 * t376 + t374 * t329) * t377 + (t291 * t377 - t292 * t375) * t376 + ((-t294 * t371 + t296 * t372) * t375 - (-t293 * t371 + t295 * t372) * t377) * t374) * qJD(2) + ((t349 - t315) * t376 + (-t316 * t371 + t317 * t372 + t350) * t374) * qJD(1)) * qJD(1) / 0.2e1 + (((-t291 * t409 - t293 * t338 - t295 * t339 + t386 * t377) * t377 + ((-t325 + t385) * t377 + t292 * t409 + t294 * t338 + t296 * t339 + t326 * t375) * t375) * qJD(2) + (t315 * t409 + t316 * t338 + t317 * t339 + t375 * t348 + t377 * t384) * qJD(1)) * t367 / 0.2e1 - (((-t291 * t410 - t293 * t336 - t295 * t337 + t325 * t377) * t377 + ((-t326 + t386) * t377 + t292 * t410 + t294 * t336 + t296 * t337 + t385 * t375) * t375) * qJD(2) + (t315 * t410 + t316 * t336 + t317 * t337 - t377 * t348 + t375 * t384) * qJD(1)) * t402 / 0.2e1;
T = t1;
