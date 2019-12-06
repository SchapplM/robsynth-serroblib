% Calculate kinetic energy for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:43
% EndTime: 2019-12-05 16:47:45
% DurationCPUTime: 2.05s
% Computational Cost: add. (1074->192), mult. (1763->313), div. (0->0), fcn. (1799->8), ass. (0->105)
t435 = Icges(5,1) + Icges(6,1);
t434 = Icges(5,4) + Icges(6,4);
t433 = -Icges(6,5) - Icges(5,5);
t432 = Icges(5,2) + Icges(6,2);
t431 = -Icges(6,6) - Icges(5,6);
t430 = -Icges(6,3) - Icges(5,3);
t368 = cos(pkin(8));
t416 = t368 ^ 2;
t367 = sin(pkin(8));
t417 = t367 ^ 2;
t418 = t416 + t417;
t429 = qJD(2) * t418;
t428 = t367 * t368;
t366 = qJ(3) + qJ(4);
t362 = sin(t366);
t363 = cos(t366);
t372 = cos(qJ(2));
t406 = t367 * t372;
t338 = -t362 * t406 - t363 * t368;
t339 = -t362 * t368 + t363 * t406;
t370 = sin(qJ(2));
t407 = t367 * t370;
t427 = -t431 * t338 - t433 * t339 - t430 * t407;
t403 = t368 * t372;
t340 = -t362 * t403 + t363 * t367;
t341 = t362 * t367 + t363 * t403;
t404 = t368 * t370;
t426 = -t431 * t340 - t433 * t341 - t430 * t404;
t425 = t432 * t338 + t434 * t339 - t431 * t407;
t424 = t432 * t340 + t434 * t341 - t431 * t404;
t423 = t434 * t338 + t435 * t339 - t433 * t407;
t422 = t434 * t340 + t435 * t341 - t433 * t404;
t421 = t430 * t372 + (t431 * t362 - t433 * t363) * t370;
t420 = t431 * t372 + (-t432 * t362 + t434 * t363) * t370;
t419 = t433 * t372 + (-t434 * t362 + t435 * t363) * t370;
t415 = qJD(2) ^ 2;
t371 = cos(qJ(3));
t410 = t371 * pkin(3);
t369 = sin(qJ(3));
t408 = t367 * t369;
t405 = t368 * t369;
t402 = t369 * t372;
t401 = t371 * t372;
t396 = pkin(4) * t363;
t375 = qJ(5) * t370 + t372 * t396;
t387 = pkin(4) * t362;
t400 = rSges(6,1) * t339 + rSges(6,2) * t338 + rSges(6,3) * t407 + t367 * t375 - t368 * t387;
t399 = rSges(6,1) * t341 + rSges(6,2) * t340 + rSges(6,3) * t404 + t367 * t387 + t368 * t375;
t376 = pkin(7) * t370 + t372 * t410;
t315 = -pkin(3) * t405 + t367 * t376;
t319 = -pkin(7) * t372 + t370 * t410;
t393 = qJD(3) * t370;
t394 = qJD(2) * t368;
t351 = t367 * t393 - t394;
t392 = qJD(3) * t372;
t398 = t315 * t392 + t351 * t319;
t397 = (-qJ(5) - rSges(6,3)) * t372 + (rSges(6,1) * t363 - rSges(6,2) * t362 + t396) * t370;
t361 = qJD(2) * t367;
t350 = t368 * t393 + t361;
t391 = qJD(4) * t370;
t357 = pkin(2) * t370 - pkin(6) * t372;
t390 = t357 * t361;
t389 = t357 * t394;
t388 = qJD(1) + (pkin(2) * t372 + pkin(6) * t370) * t429;
t384 = Icges(3,5) * t372 - Icges(3,6) * t370;
t381 = -qJD(2) * t357 + qJD(5) * t370;
t316 = pkin(3) * t408 + t368 * t376;
t379 = -t316 * t392 - t319 * t350;
t377 = t350 * t315 - t316 * t351 + t388;
t356 = rSges(3,1) * t370 + rSges(3,2) * t372;
t354 = (-qJD(3) - qJD(4)) * t372;
t349 = t368 * t401 + t408;
t348 = t367 * t371 - t368 * t402;
t347 = t367 * t401 - t405;
t346 = -t367 * t402 - t368 * t371;
t345 = -rSges(4,3) * t372 + (rSges(4,1) * t371 - rSges(4,2) * t369) * t370;
t344 = -Icges(4,5) * t372 + (Icges(4,1) * t371 - Icges(4,4) * t369) * t370;
t343 = -Icges(4,6) * t372 + (Icges(4,4) * t371 - Icges(4,2) * t369) * t370;
t342 = -Icges(4,3) * t372 + (Icges(4,5) * t371 - Icges(4,6) * t369) * t370;
t331 = Icges(3,3) * t367 + t368 * t384;
t330 = -Icges(3,3) * t368 + t367 * t384;
t329 = t367 * t391 + t351;
t328 = t368 * t391 + t350;
t327 = -rSges(5,3) * t372 + (rSges(5,1) * t363 - rSges(5,2) * t362) * t370;
t314 = rSges(4,1) * t349 + rSges(4,2) * t348 + rSges(4,3) * t404;
t313 = rSges(4,1) * t347 + rSges(4,2) * t346 + rSges(4,3) * t407;
t312 = Icges(4,1) * t349 + Icges(4,4) * t348 + Icges(4,5) * t404;
t311 = Icges(4,1) * t347 + Icges(4,4) * t346 + Icges(4,5) * t407;
t310 = Icges(4,4) * t349 + Icges(4,2) * t348 + Icges(4,6) * t404;
t309 = Icges(4,4) * t347 + Icges(4,2) * t346 + Icges(4,6) * t407;
t308 = Icges(4,5) * t349 + Icges(4,6) * t348 + Icges(4,3) * t404;
t307 = Icges(4,5) * t347 + Icges(4,6) * t346 + Icges(4,3) * t407;
t305 = rSges(5,1) * t341 + rSges(5,2) * t340 + rSges(5,3) * t404;
t303 = rSges(5,1) * t339 + rSges(5,2) * t338 + rSges(5,3) * t407;
t289 = qJD(1) + (rSges(3,1) * t372 - rSges(3,2) * t370) * t429;
t285 = t313 * t392 + t345 * t351 - t389;
t284 = -t314 * t392 - t345 * t350 - t390;
t283 = t313 * t350 - t314 * t351 + t388;
t282 = -t303 * t354 + t327 * t329 - t389 + t398;
t281 = t305 * t354 - t327 * t328 + t379 - t390;
t280 = t303 * t328 - t305 * t329 + t377;
t279 = t329 * t397 - t354 * t400 + t368 * t381 + t398;
t278 = -t328 * t397 + t354 * t399 + t367 * t381 + t379;
t277 = -qJD(5) * t372 + t328 * t400 - t329 * t399 + t377;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t418 * t415 * t356 ^ 2 + t289 ^ 2) / 0.2e1 + t415 * t367 * (-t330 * t428 + t417 * t331) / 0.2e1 - t415 * t368 * (t416 * t330 - t331 * t428) / 0.2e1 + m(4) * (t283 ^ 2 + t284 ^ 2 + t285 ^ 2) / 0.2e1 + t350 * ((t308 * t404 + t348 * t310 + t349 * t312) * t350 + (t307 * t404 + t309 * t348 + t311 * t349) * t351 - (t342 * t404 + t343 * t348 + t344 * t349) * t392) / 0.2e1 + t351 * ((t308 * t407 + t310 * t346 + t312 * t347) * t350 + (t307 * t407 + t346 * t309 + t347 * t311) * t351 - (t342 * t407 + t343 * t346 + t344 * t347) * t392) / 0.2e1 - ((-t307 * t351 - t308 * t350 + t342 * t392) * t372 + ((-t310 * t369 + t312 * t371) * t350 + (-t309 * t369 + t311 * t371) * t351 - (-t343 * t369 + t344 * t371) * t392) * t370) * t392 / 0.2e1 + m(5) * (t280 ^ 2 + t281 ^ 2 + t282 ^ 2) / 0.2e1 + m(6) * (t277 ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + ((t420 * t340 + t419 * t341 + t421 * t404) * t354 + (t425 * t340 + t423 * t341 + t427 * t404) * t329 + (t424 * t340 + t422 * t341 + t426 * t404) * t328) * t328 / 0.2e1 + ((t420 * t338 + t419 * t339 + t421 * t407) * t354 + (t425 * t338 + t423 * t339 + t427 * t407) * t329 + (t424 * t338 + t422 * t339 + t426 * t407) * t328) * t329 / 0.2e1 + ((-t426 * t328 - t427 * t329 - t421 * t354) * t372 + ((-t420 * t362 + t419 * t363) * t354 + (-t425 * t362 + t423 * t363) * t329 + (-t424 * t362 + t422 * t363) * t328) * t370) * t354 / 0.2e1;
T = t1;
