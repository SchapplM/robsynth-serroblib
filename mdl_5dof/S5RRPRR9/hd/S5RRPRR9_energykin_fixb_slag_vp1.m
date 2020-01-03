% Calculate kinetic energy for
% S5RRPRR9
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:00
% EndTime: 2019-12-31 20:20:02
% DurationCPUTime: 1.97s
% Computational Cost: add. (1276->254), mult. (1529->405), div. (0->0), fcn. (1482->10), ass. (0->137)
t447 = Icges(3,3) + Icges(4,3);
t377 = qJ(2) + pkin(9);
t371 = sin(t377);
t372 = cos(t377);
t381 = sin(qJ(2));
t384 = cos(qJ(2));
t446 = Icges(3,5) * t384 + Icges(4,5) * t372 - Icges(3,6) * t381 - Icges(4,6) * t371;
t382 = sin(qJ(1));
t385 = cos(qJ(1));
t445 = t446 * t382 - t447 * t385;
t444 = t447 * t382 + t446 * t385;
t443 = Icges(3,5) * t381 + Icges(4,5) * t371 + Icges(3,6) * t384 + Icges(4,6) * t372;
t428 = Icges(4,4) * t371;
t354 = Icges(4,2) * t372 + t428;
t427 = Icges(4,4) * t372;
t355 = Icges(4,1) * t371 + t427;
t430 = Icges(3,4) * t381;
t360 = Icges(3,2) * t384 + t430;
t429 = Icges(3,4) * t384;
t361 = Icges(3,1) * t381 + t429;
t442 = -t354 * t371 + t355 * t372 - t360 * t381 + t361 * t384;
t399 = -Icges(4,2) * t371 + t427;
t324 = Icges(4,6) * t382 + t385 * t399;
t401 = Icges(4,1) * t372 - t428;
t326 = Icges(4,5) * t382 + t385 * t401;
t400 = -Icges(3,2) * t381 + t429;
t339 = Icges(3,6) * t382 + t385 * t400;
t402 = Icges(3,1) * t384 - t430;
t341 = Icges(3,5) * t382 + t385 * t402;
t441 = -t324 * t371 + t326 * t372 - t339 * t381 + t341 * t384;
t323 = -Icges(4,6) * t385 + t382 * t399;
t325 = -Icges(4,5) * t385 + t382 * t401;
t338 = -Icges(3,6) * t385 + t382 * t400;
t340 = -Icges(3,5) * t385 + t382 * t402;
t440 = t323 * t371 - t325 * t372 + t338 * t381 - t340 * t384;
t436 = pkin(2) * t381;
t434 = pkin(2) * t384;
t383 = cos(qJ(4));
t433 = pkin(4) * t383;
t426 = t371 * t382;
t425 = t371 * t385;
t378 = qJ(4) + qJ(5);
t375 = sin(t378);
t424 = t375 * t382;
t423 = t375 * t385;
t376 = cos(t378);
t422 = t376 * t382;
t421 = t376 * t385;
t380 = sin(qJ(4));
t420 = t380 * t382;
t419 = t380 * t385;
t418 = t382 * t383;
t417 = t383 * t385;
t319 = -qJ(3) * t385 + t382 * t434;
t320 = qJ(3) * t382 + t385 * t434;
t374 = qJD(2) * t382;
t414 = qJD(2) * t385;
t416 = t319 * t374 + t320 * t414;
t367 = pkin(1) * t382 - pkin(6) * t385;
t415 = -t319 - t367;
t413 = qJD(4) * t371;
t351 = t385 * t413 + t374;
t412 = qJD(5) * t371;
t352 = t382 * t413 - t414;
t408 = pkin(3) * t372 + pkin(7) * t371;
t344 = t408 * t382;
t345 = t408 * t385;
t409 = t344 * t374 + t345 * t414 + t416;
t358 = qJD(1) * (pkin(1) * t385 + pkin(6) * t382);
t407 = qJD(1) * t320 - qJD(3) * t385 + t358;
t406 = rSges(3,1) * t384 - rSges(3,2) * t381;
t405 = rSges(4,1) * t372 - rSges(4,2) * t371;
t404 = qJD(2) * (-rSges(4,1) * t371 - rSges(4,2) * t372 - t436);
t403 = qJD(2) * (-pkin(3) * t371 + pkin(7) * t372 - t436);
t390 = pkin(8) * t371 + t372 * t433;
t389 = qJD(1) * t345 + t382 * t403 + t407;
t373 = qJD(3) * t382;
t388 = t373 + (-t344 + t415) * qJD(1) + t385 * t403;
t368 = -qJD(4) * t372 + qJD(1);
t366 = rSges(2,1) * t385 - rSges(2,2) * t382;
t365 = rSges(2,1) * t382 + rSges(2,2) * t385;
t364 = rSges(3,1) * t381 + rSges(3,2) * t384;
t350 = qJD(1) + (-qJD(4) - qJD(5)) * t372;
t349 = t372 * t417 + t420;
t348 = -t372 * t419 + t418;
t347 = t372 * t418 - t419;
t346 = -t372 * t420 - t417;
t343 = rSges(3,3) * t382 + t385 * t406;
t342 = -rSges(3,3) * t385 + t382 * t406;
t334 = t372 * t421 + t424;
t333 = -t372 * t423 + t422;
t332 = t372 * t422 - t423;
t331 = -t372 * t424 - t421;
t328 = rSges(4,3) * t382 + t385 * t405;
t327 = -rSges(4,3) * t385 + t382 * t405;
t318 = t382 * t412 + t352;
t317 = t385 * t412 + t351;
t316 = -rSges(5,3) * t372 + (rSges(5,1) * t383 - rSges(5,2) * t380) * t371;
t314 = -Icges(5,5) * t372 + (Icges(5,1) * t383 - Icges(5,4) * t380) * t371;
t313 = -Icges(5,6) * t372 + (Icges(5,4) * t383 - Icges(5,2) * t380) * t371;
t312 = -Icges(5,3) * t372 + (Icges(5,5) * t383 - Icges(5,6) * t380) * t371;
t309 = -rSges(6,3) * t372 + (rSges(6,1) * t376 - rSges(6,2) * t375) * t371;
t308 = -Icges(6,5) * t372 + (Icges(6,1) * t376 - Icges(6,4) * t375) * t371;
t307 = -Icges(6,6) * t372 + (Icges(6,4) * t376 - Icges(6,2) * t375) * t371;
t306 = -Icges(6,3) * t372 + (Icges(6,5) * t376 - Icges(6,6) * t375) * t371;
t305 = -pkin(8) * t372 + t371 * t433;
t304 = qJD(1) * t343 - t364 * t374 + t358;
t303 = -t364 * t414 + (-t342 - t367) * qJD(1);
t302 = (t342 * t382 + t343 * t385) * qJD(2);
t301 = rSges(5,1) * t349 + rSges(5,2) * t348 + rSges(5,3) * t425;
t300 = rSges(5,1) * t347 + rSges(5,2) * t346 + rSges(5,3) * t426;
t299 = Icges(5,1) * t349 + Icges(5,4) * t348 + Icges(5,5) * t425;
t298 = Icges(5,1) * t347 + Icges(5,4) * t346 + Icges(5,5) * t426;
t297 = Icges(5,4) * t349 + Icges(5,2) * t348 + Icges(5,6) * t425;
t296 = Icges(5,4) * t347 + Icges(5,2) * t346 + Icges(5,6) * t426;
t295 = Icges(5,5) * t349 + Icges(5,6) * t348 + Icges(5,3) * t425;
t294 = Icges(5,5) * t347 + Icges(5,6) * t346 + Icges(5,3) * t426;
t293 = pkin(4) * t420 + t385 * t390;
t292 = -pkin(4) * t419 + t382 * t390;
t291 = rSges(6,1) * t334 + rSges(6,2) * t333 + rSges(6,3) * t425;
t290 = rSges(6,1) * t332 + rSges(6,2) * t331 + rSges(6,3) * t426;
t289 = Icges(6,1) * t334 + Icges(6,4) * t333 + Icges(6,5) * t425;
t288 = Icges(6,1) * t332 + Icges(6,4) * t331 + Icges(6,5) * t426;
t287 = Icges(6,4) * t334 + Icges(6,2) * t333 + Icges(6,6) * t425;
t286 = Icges(6,4) * t332 + Icges(6,2) * t331 + Icges(6,6) * t426;
t285 = Icges(6,5) * t334 + Icges(6,6) * t333 + Icges(6,3) * t425;
t284 = Icges(6,5) * t332 + Icges(6,6) * t331 + Icges(6,3) * t426;
t283 = qJD(1) * t328 + t382 * t404 + t407;
t282 = t373 + t385 * t404 + (-t327 + t415) * qJD(1);
t281 = (t327 * t382 + t328 * t385) * qJD(2) + t416;
t280 = t301 * t368 - t316 * t351 + t389;
t279 = -t300 * t368 + t316 * t352 + t388;
t278 = t300 * t351 - t301 * t352 + t409;
t277 = t291 * t350 + t293 * t368 - t305 * t351 - t309 * t317 + t389;
t276 = -t290 * t350 - t292 * t368 + t305 * t352 + t309 * t318 + t388;
t275 = t290 * t317 - t291 * t318 + t292 * t351 - t293 * t352 + t409;
t1 = m(3) * (t302 ^ 2 + t303 ^ 2 + t304 ^ 2) / 0.2e1 + m(4) * (t281 ^ 2 + t282 ^ 2 + t283 ^ 2) / 0.2e1 + m(5) * (t278 ^ 2 + t279 ^ 2 + t280 ^ 2) / 0.2e1 + t351 * ((t295 * t425 + t348 * t297 + t349 * t299) * t351 + (t294 * t425 + t296 * t348 + t298 * t349) * t352 + (t312 * t425 + t313 * t348 + t314 * t349) * t368) / 0.2e1 + t352 * ((t295 * t426 + t297 * t346 + t299 * t347) * t351 + (t294 * t426 + t346 * t296 + t347 * t298) * t352 + (t312 * t426 + t313 * t346 + t314 * t347) * t368) / 0.2e1 + t368 * ((-t294 * t352 - t295 * t351 - t312 * t368) * t372 + ((-t297 * t380 + t299 * t383) * t351 + (-t296 * t380 + t298 * t383) * t352 + (-t313 * t380 + t314 * t383) * t368) * t371) / 0.2e1 + m(6) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + t317 * ((t285 * t425 + t333 * t287 + t334 * t289) * t317 + (t284 * t425 + t286 * t333 + t288 * t334) * t318 + (t306 * t425 + t307 * t333 + t308 * t334) * t350) / 0.2e1 + t318 * ((t285 * t426 + t287 * t331 + t289 * t332) * t317 + (t284 * t426 + t331 * t286 + t332 * t288) * t318 + (t306 * t426 + t307 * t331 + t308 * t332) * t350) / 0.2e1 + t350 * ((-t284 * t318 - t285 * t317 - t306 * t350) * t372 + ((-t287 * t375 + t289 * t376) * t317 + (-t286 * t375 + t288 * t376) * t318 + (-t307 * t375 + t308 * t376) * t350) * t371) / 0.2e1 + (m(2) * (t365 ^ 2 + t366 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t323 * t372 - t325 * t371 - t338 * t384 - t381 * t340) * t385 + (t324 * t372 + t371 * t326 + t339 * t384 + t341 * t381) * t382) * qJD(2) + (t372 * t354 + t371 * t355 + t384 * t360 + t381 * t361) * qJD(1)) * qJD(1) / 0.2e1 + ((t444 * t382 ^ 2 + (t440 * t385 + (t441 - t445) * t382) * t385) * qJD(2) + (t443 * t382 + t442 * t385) * qJD(1)) * t374 / 0.2e1 - ((t445 * t385 ^ 2 + (t441 * t382 + (t440 - t444) * t385) * t382) * qJD(2) + (t442 * t382 - t443 * t385) * qJD(1)) * t414 / 0.2e1;
T = t1;
