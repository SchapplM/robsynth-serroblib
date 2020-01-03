% Calculate kinetic energy for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:37
% EndTime: 2019-12-31 21:59:39
% DurationCPUTime: 2.19s
% Computational Cost: add. (1160->238), mult. (1898->378), div. (0->0), fcn. (1907->8), ass. (0->125)
t440 = Icges(5,1) + Icges(6,1);
t439 = Icges(5,4) + Icges(6,4);
t438 = -Icges(6,5) - Icges(5,5);
t437 = Icges(5,2) + Icges(6,2);
t436 = -Icges(6,6) - Icges(5,6);
t435 = -Icges(6,3) - Icges(5,3);
t372 = qJ(3) + qJ(4);
t368 = sin(t372);
t369 = cos(t372);
t378 = cos(qJ(1));
t375 = sin(qJ(1));
t377 = cos(qJ(2));
t408 = t375 * t377;
t336 = -t368 * t408 - t369 * t378;
t337 = -t368 * t378 + t369 * t408;
t374 = sin(qJ(2));
t410 = t374 * t375;
t434 = -t436 * t336 - t438 * t337 - t435 * t410;
t407 = t377 * t378;
t338 = -t368 * t407 + t369 * t375;
t339 = t368 * t375 + t369 * t407;
t409 = t374 * t378;
t433 = -t436 * t338 - t438 * t339 - t435 * t409;
t432 = t437 * t336 + t439 * t337 - t436 * t410;
t431 = t437 * t338 + t439 * t339 - t436 * t409;
t430 = t439 * t336 + t440 * t337 - t438 * t410;
t429 = t439 * t338 + t440 * t339 - t438 * t409;
t376 = cos(qJ(3));
t418 = t376 * pkin(3);
t383 = pkin(8) * t374 + t377 * t418;
t373 = sin(qJ(3));
t412 = t373 * t375;
t310 = pkin(3) * t412 + t378 * t383;
t313 = -pkin(8) * t377 + t374 * t418;
t367 = qJD(2) * t375;
t398 = qJD(3) * t374;
t349 = t378 * t398 + t367;
t365 = -qJD(3) * t377 + qJD(1);
t428 = t365 * t310 - t313 * t349;
t411 = t373 * t378;
t309 = -pkin(3) * t411 + t375 * t383;
t399 = qJD(2) * t378;
t350 = t375 * t398 - t399;
t427 = -t309 * t365 + t350 * t313;
t426 = t435 * t377 + (t436 * t368 - t438 * t369) * t374;
t425 = t436 * t377 + (-t437 * t368 + t439 * t369) * t374;
t424 = t438 * t377 + (-t439 * t368 + t440 * t369) * t374;
t416 = Icges(3,4) * t374;
t415 = Icges(3,4) * t377;
t401 = pkin(4) * t369;
t381 = qJ(5) * t374 + t377 * t401;
t396 = pkin(4) * t368;
t406 = rSges(6,1) * t337 + rSges(6,2) * t336 + rSges(6,3) * t410 + t375 * t381 - t378 * t396;
t405 = rSges(6,1) * t339 + rSges(6,2) * t338 + rSges(6,3) * t409 + t375 * t396 + t378 * t381;
t404 = (-qJ(5) - rSges(6,3)) * t377 + (rSges(6,1) * t369 - rSges(6,2) * t368 + t401) * t374;
t394 = pkin(2) * t377 + pkin(7) * t374;
t347 = t394 * t375;
t348 = t394 * t378;
t403 = t347 * t367 + t348 * t399;
t353 = qJD(1) * (pkin(1) * t378 + pkin(6) * t375);
t402 = qJD(1) * t348 + t353;
t397 = qJD(4) * t374;
t362 = pkin(1) * t375 - pkin(6) * t378;
t395 = (-t347 - t362) * qJD(1);
t393 = rSges(3,1) * t377 - rSges(3,2) * t374;
t392 = Icges(3,1) * t377 - t416;
t391 = -Icges(3,2) * t374 + t415;
t390 = Icges(3,5) * t377 - Icges(3,6) * t374;
t328 = -Icges(3,6) * t378 + t375 * t391;
t331 = -Icges(3,5) * t378 + t375 * t392;
t389 = t328 * t374 - t331 * t377;
t329 = Icges(3,6) * t375 + t378 * t391;
t332 = Icges(3,5) * t375 + t378 * t392;
t388 = -t329 * t374 + t332 * t377;
t356 = Icges(3,2) * t377 + t416;
t357 = Icges(3,1) * t374 + t415;
t387 = -t356 * t374 + t357 * t377;
t361 = pkin(2) * t374 - pkin(7) * t377;
t386 = -qJD(2) * t361 + qJD(5) * t374;
t385 = t349 * t309 - t310 * t350 + t403;
t384 = -t361 * t367 + t402;
t382 = -t361 * t399 + t395;
t360 = rSges(2,1) * t378 - rSges(2,2) * t375;
t359 = rSges(2,1) * t375 + rSges(2,2) * t378;
t358 = rSges(3,1) * t374 + rSges(3,2) * t377;
t355 = Icges(3,5) * t374 + Icges(3,6) * t377;
t352 = qJD(1) + (-qJD(3) - qJD(4)) * t377;
t346 = t376 * t407 + t412;
t345 = -t373 * t407 + t375 * t376;
t344 = t376 * t408 - t411;
t343 = -t373 * t408 - t376 * t378;
t335 = rSges(3,3) * t375 + t378 * t393;
t334 = -rSges(3,3) * t378 + t375 * t393;
t333 = -rSges(4,3) * t377 + (rSges(4,1) * t376 - rSges(4,2) * t373) * t374;
t330 = -Icges(4,5) * t377 + (Icges(4,1) * t376 - Icges(4,4) * t373) * t374;
t327 = -Icges(4,6) * t377 + (Icges(4,4) * t376 - Icges(4,2) * t373) * t374;
t326 = Icges(3,3) * t375 + t378 * t390;
t325 = -Icges(3,3) * t378 + t375 * t390;
t324 = -Icges(4,3) * t377 + (Icges(4,5) * t376 - Icges(4,6) * t373) * t374;
t323 = t375 * t397 + t350;
t322 = t378 * t397 + t349;
t321 = -rSges(5,3) * t377 + (rSges(5,1) * t369 - rSges(5,2) * t368) * t374;
t308 = rSges(4,1) * t346 + rSges(4,2) * t345 + rSges(4,3) * t409;
t307 = rSges(4,1) * t344 + rSges(4,2) * t343 + rSges(4,3) * t410;
t306 = Icges(4,1) * t346 + Icges(4,4) * t345 + Icges(4,5) * t409;
t305 = Icges(4,1) * t344 + Icges(4,4) * t343 + Icges(4,5) * t410;
t304 = Icges(4,4) * t346 + Icges(4,2) * t345 + Icges(4,6) * t409;
t303 = Icges(4,4) * t344 + Icges(4,2) * t343 + Icges(4,6) * t410;
t302 = Icges(4,5) * t346 + Icges(4,6) * t345 + Icges(4,3) * t409;
t301 = Icges(4,5) * t344 + Icges(4,6) * t343 + Icges(4,3) * t410;
t300 = qJD(1) * t335 - t358 * t367 + t353;
t299 = -t358 * t399 + (-t334 - t362) * qJD(1);
t298 = (t334 * t375 + t335 * t378) * qJD(2);
t296 = rSges(5,1) * t339 + rSges(5,2) * t338 + rSges(5,3) * t409;
t294 = rSges(5,1) * t337 + rSges(5,2) * t336 + rSges(5,3) * t410;
t277 = t308 * t365 - t333 * t349 + t384;
t276 = -t307 * t365 + t333 * t350 + t382;
t275 = t307 * t349 - t308 * t350 + t403;
t274 = t296 * t352 - t321 * t322 + t384 + t428;
t273 = -t294 * t352 + t321 * t323 + t382 + t427;
t272 = t294 * t322 - t296 * t323 + t385;
t271 = -t322 * t404 + t352 * t405 + t375 * t386 + t402 + t428;
t270 = t323 * t404 - t352 * t406 + t378 * t386 + t395 + t427;
t269 = -qJD(5) * t377 + t322 * t406 - t323 * t405 + t385;
t1 = m(3) * (t298 ^ 2 + t299 ^ 2 + t300 ^ 2) / 0.2e1 + ((t375 * t355 + t378 * t387) * qJD(1) + (t375 ^ 2 * t326 + (t389 * t378 + (-t325 + t388) * t375) * t378) * qJD(2)) * t367 / 0.2e1 - ((-t378 * t355 + t375 * t387) * qJD(1) + (t378 ^ 2 * t325 + (t388 * t375 + (-t326 + t389) * t378) * t375) * qJD(2)) * t399 / 0.2e1 + qJD(1) * ((t377 * t356 + t374 * t357) * qJD(1) + ((t329 * t377 + t332 * t374) * t375 - (t328 * t377 + t374 * t331) * t378) * qJD(2)) / 0.2e1 + m(4) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + t349 * ((t302 * t409 + t345 * t304 + t346 * t306) * t349 + (t301 * t409 + t303 * t345 + t305 * t346) * t350 + (t324 * t409 + t327 * t345 + t330 * t346) * t365) / 0.2e1 + t350 * ((t302 * t410 + t304 * t343 + t306 * t344) * t349 + (t301 * t410 + t343 * t303 + t344 * t305) * t350 + (t324 * t410 + t327 * t343 + t330 * t344) * t365) / 0.2e1 + t365 * ((-t301 * t350 - t302 * t349 - t324 * t365) * t377 + ((-t304 * t373 + t306 * t376) * t349 + (-t303 * t373 + t305 * t376) * t350 + (-t327 * t373 + t330 * t376) * t365) * t374) / 0.2e1 + m(5) * (t272 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + m(6) * (t269 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + ((t338 * t425 + t339 * t424 + t409 * t426) * t352 + (t432 * t338 + t430 * t339 + t409 * t434) * t323 + (t431 * t338 + t429 * t339 + t433 * t409) * t322) * t322 / 0.2e1 + ((t336 * t425 + t337 * t424 + t410 * t426) * t352 + (t432 * t336 + t430 * t337 + t434 * t410) * t323 + (t336 * t431 + t337 * t429 + t410 * t433) * t322) * t323 / 0.2e1 + ((-t433 * t322 - t323 * t434 - t426 * t352) * t377 + ((-t368 * t425 + t369 * t424) * t352 + (-t368 * t432 + t369 * t430) * t323 + (-t368 * t431 + t369 * t429) * t322) * t374) * t352 / 0.2e1 + (m(2) * (t359 ^ 2 + t360 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
