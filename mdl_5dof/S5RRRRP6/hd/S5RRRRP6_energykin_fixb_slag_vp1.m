% Calculate kinetic energy for
% S5RRRRP6
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:47
% EndTime: 2019-12-31 21:52:49
% DurationCPUTime: 1.81s
% Computational Cost: add. (1148->212), mult. (1546->334), div. (0->0), fcn. (1499->8), ass. (0->122)
t437 = Icges(5,1) + Icges(6,1);
t436 = Icges(5,4) + Icges(6,4);
t435 = -Icges(6,5) - Icges(5,5);
t434 = Icges(5,2) + Icges(6,2);
t433 = -Icges(6,6) - Icges(5,6);
t432 = -Icges(6,3) - Icges(5,3);
t363 = qJ(2) + qJ(3);
t362 = cos(t363);
t368 = cos(qJ(4));
t370 = cos(qJ(1));
t402 = t368 * t370;
t365 = sin(qJ(4));
t367 = sin(qJ(1));
t405 = t365 * t367;
t337 = -t362 * t405 - t402;
t403 = t367 * t368;
t404 = t365 * t370;
t338 = t362 * t403 - t404;
t361 = sin(t363);
t407 = t361 * t367;
t431 = -t433 * t337 - t435 * t338 - t432 * t407;
t339 = -t362 * t404 + t403;
t340 = t362 * t402 + t405;
t406 = t361 * t370;
t430 = -t433 * t339 - t435 * t340 - t432 * t406;
t429 = t434 * t337 + t436 * t338 - t433 * t407;
t428 = t434 * t339 + t436 * t340 - t433 * t406;
t427 = t436 * t337 + t437 * t338 - t435 * t407;
t426 = t436 * t339 + t437 * t340 - t435 * t406;
t425 = t432 * t362 + (t433 * t365 - t435 * t368) * t361;
t424 = t433 * t362 + (-t434 * t365 + t436 * t368) * t361;
t423 = t435 * t362 + (-t436 * t365 + t437 * t368) * t361;
t391 = pkin(3) * t362 + pkin(8) * t361;
t336 = t391 * t370;
t345 = pkin(3) * t361 - pkin(8) * t362;
t360 = qJD(2) * t367;
t347 = qJD(3) * t367 + t360;
t422 = qJD(1) * t336 - t345 * t347;
t369 = cos(qJ(2));
t416 = pkin(2) * t369;
t415 = pkin(4) * t368;
t366 = sin(qJ(2));
t412 = Icges(3,4) * t366;
t411 = Icges(3,4) * t369;
t410 = Icges(4,4) * t361;
t409 = Icges(4,4) * t362;
t375 = qJ(5) * t361 + t362 * t415;
t401 = rSges(6,1) * t338 + rSges(6,2) * t337 + rSges(6,3) * t407 - pkin(4) * t404 + t367 * t375;
t400 = rSges(6,1) * t340 + rSges(6,2) * t339 + rSges(6,3) * t406 + pkin(4) * t405 + t370 * t375;
t399 = (-qJ(5) - rSges(6,3)) * t362 + (rSges(6,1) * t368 - rSges(6,2) * t365 + t415) * t361;
t313 = -pkin(7) * t370 + t367 * t416;
t314 = pkin(7) * t367 + t370 * t416;
t395 = qJD(2) * t370;
t398 = t313 * t360 + t314 * t395;
t346 = qJD(1) * (pkin(1) * t370 + pkin(6) * t367);
t397 = qJD(1) * t314 + t346;
t355 = pkin(1) * t367 - pkin(6) * t370;
t396 = -t313 - t355;
t394 = qJD(4) * t361;
t393 = pkin(2) * qJD(2) * t366;
t348 = (-qJD(2) - qJD(3)) * t370;
t392 = t370 * t393;
t390 = rSges(3,1) * t369 - rSges(3,2) * t366;
t389 = rSges(4,1) * t362 - rSges(4,2) * t361;
t388 = Icges(3,1) * t369 - t412;
t387 = Icges(4,1) * t362 - t410;
t386 = -Icges(3,2) * t366 + t411;
t385 = -Icges(4,2) * t361 + t409;
t384 = Icges(3,5) * t369 - Icges(3,6) * t366;
t383 = Icges(4,5) * t362 - Icges(4,6) * t361;
t327 = -Icges(3,6) * t370 + t367 * t386;
t329 = -Icges(3,5) * t370 + t367 * t388;
t382 = t327 * t366 - t329 * t369;
t328 = Icges(3,6) * t367 + t370 * t386;
t330 = Icges(3,5) * t367 + t370 * t388;
t381 = -t328 * t366 + t330 * t369;
t350 = Icges(3,2) * t369 + t412;
t351 = Icges(3,1) * t366 + t411;
t380 = -t350 * t366 + t351 * t369;
t335 = t391 * t367;
t379 = t347 * t335 - t336 * t348 + t398;
t378 = qJD(5) * t361 - t393;
t377 = t348 * t345 + (-t335 + t396) * qJD(1);
t376 = -t367 * t393 + t397;
t374 = qJD(1) * (Icges(4,5) * t361 + Icges(4,6) * t362) + (-Icges(4,3) * t370 + t367 * t383) * t348 + (Icges(4,3) * t367 + t370 * t383) * t347;
t318 = -Icges(4,6) * t370 + t367 * t385;
t319 = Icges(4,6) * t367 + t370 * t385;
t320 = -Icges(4,5) * t370 + t367 * t387;
t321 = Icges(4,5) * t367 + t370 * t387;
t342 = Icges(4,2) * t362 + t410;
t343 = Icges(4,1) * t361 + t409;
t373 = (-t319 * t361 + t321 * t362) * t347 + (-t318 * t361 + t320 * t362) * t348 + (-t342 * t361 + t343 * t362) * qJD(1);
t356 = -qJD(4) * t362 + qJD(1);
t354 = rSges(2,1) * t370 - rSges(2,2) * t367;
t353 = rSges(2,1) * t367 + rSges(2,2) * t370;
t352 = rSges(3,1) * t366 + rSges(3,2) * t369;
t349 = Icges(3,5) * t366 + Icges(3,6) * t369;
t344 = rSges(4,1) * t361 + rSges(4,2) * t362;
t334 = rSges(3,3) * t367 + t370 * t390;
t333 = -rSges(3,3) * t370 + t367 * t390;
t332 = t367 * t394 + t348;
t331 = t370 * t394 + t347;
t326 = Icges(3,3) * t367 + t370 * t384;
t325 = -Icges(3,3) * t370 + t367 * t384;
t323 = rSges(4,3) * t367 + t370 * t389;
t322 = -rSges(4,3) * t370 + t367 * t389;
t312 = -rSges(5,3) * t362 + (rSges(5,1) * t368 - rSges(5,2) * t365) * t361;
t299 = qJD(1) * t334 - t352 * t360 + t346;
t298 = -t352 * t395 + (-t333 - t355) * qJD(1);
t297 = rSges(5,1) * t340 + rSges(5,2) * t339 + rSges(5,3) * t406;
t295 = rSges(5,1) * t338 + rSges(5,2) * t337 + rSges(5,3) * t407;
t281 = (t333 * t367 + t334 * t370) * qJD(2);
t278 = qJD(1) * t323 - t344 * t347 + t376;
t277 = -t392 + t344 * t348 + (-t322 + t396) * qJD(1);
t276 = t322 * t347 - t323 * t348 + t398;
t275 = t297 * t356 - t312 * t331 + t376 + t422;
t274 = -t295 * t356 + t312 * t332 + t377 - t392;
t273 = t295 * t331 - t297 * t332 + t379;
t272 = -t331 * t399 + t356 * t400 + t367 * t378 + t397 + t422;
t271 = t332 * t399 - t356 * t401 + t370 * t378 + t377;
t270 = -qJD(5) * t362 + t331 * t401 - t332 * t400 + t379;
t1 = m(3) * (t281 ^ 2 + t298 ^ 2 + t299 ^ 2) / 0.2e1 + ((t367 * t349 + t370 * t380) * qJD(1) + (t367 ^ 2 * t326 + (t382 * t370 + (-t325 + t381) * t367) * t370) * qJD(2)) * t360 / 0.2e1 - ((-t370 * t349 + t367 * t380) * qJD(1) + (t370 ^ 2 * t325 + (t381 * t367 + (-t326 + t382) * t370) * t367) * qJD(2)) * t395 / 0.2e1 + m(4) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + t347 * (t374 * t367 + t373 * t370) / 0.2e1 + t348 * (t373 * t367 - t374 * t370) / 0.2e1 + m(5) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + m(6) * (t270 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + ((t424 * t339 + t423 * t340 + t425 * t406) * t356 + (t429 * t339 + t427 * t340 + t431 * t406) * t332 + (t428 * t339 + t426 * t340 + t430 * t406) * t331) * t331 / 0.2e1 + ((t424 * t337 + t423 * t338 + t425 * t407) * t356 + (t429 * t337 + t427 * t338 + t431 * t407) * t332 + (t428 * t337 + t426 * t338 + t430 * t407) * t331) * t332 / 0.2e1 + ((-t430 * t331 - t431 * t332 - t425 * t356) * t362 + ((-t424 * t365 + t423 * t368) * t356 + (-t429 * t365 + t427 * t368) * t332 + (-t428 * t365 + t426 * t368) * t331) * t361) * t356 / 0.2e1 + (m(2) * (t353 ^ 2 + t354 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t328 * t369 + t330 * t366) * t367 - (t327 * t369 + t329 * t366) * t370) * qJD(2) + (t319 * t362 + t321 * t361) * t347 + (t318 * t362 + t320 * t361) * t348 + (t362 * t342 + t361 * t343 + t369 * t350 + t366 * t351) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
