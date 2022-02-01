% Calculate kinetic energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:52
% EndTime: 2022-01-23 09:24:54
% DurationCPUTime: 2.11s
% Computational Cost: add. (1197->257), mult. (1600->370), div. (0->0), fcn. (1638->10), ass. (0->127)
t436 = Icges(4,3) + Icges(5,3);
t384 = cos(pkin(8));
t382 = qJ(3) + pkin(9);
t374 = cos(t382);
t389 = cos(qJ(1));
t417 = t389 * t374;
t373 = sin(t382);
t387 = sin(qJ(1));
t425 = t387 * t373;
t346 = -t384 * t425 - t417;
t418 = t389 * t373;
t424 = t387 * t374;
t347 = t384 * t424 - t418;
t348 = -t384 * t418 + t424;
t349 = t384 * t417 + t425;
t388 = cos(qJ(3));
t414 = t389 * t388;
t386 = sin(qJ(3));
t422 = t387 * t386;
t355 = -t384 * t422 - t414;
t415 = t389 * t386;
t421 = t387 * t388;
t356 = t384 * t421 - t415;
t357 = -t384 * t415 + t421;
t358 = t384 * t414 + t422;
t383 = sin(pkin(8));
t416 = t389 * t383;
t423 = t387 * t383;
t435 = (Icges(4,5) * t358 + Icges(5,5) * t349 + Icges(4,6) * t357 + Icges(5,6) * t348 + t436 * t416) * t389 + (Icges(4,5) * t356 + Icges(5,5) * t347 + Icges(4,6) * t355 + Icges(5,6) * t346 + t436 * t423) * t387;
t434 = t436 * t384 + (-Icges(4,5) * t388 - Icges(5,5) * t374 + Icges(4,6) * t386 + Icges(5,6) * t373) * t383;
t433 = t435 * t383;
t430 = t386 * pkin(3);
t429 = t388 * pkin(3);
t428 = pkin(2) * t384 + pkin(1);
t385 = qJ(4) + pkin(6);
t375 = qJ(5) + t382;
t370 = sin(t375);
t427 = t387 * t370;
t371 = cos(t375);
t426 = t387 * t371;
t420 = t389 * t370;
t419 = t389 * t371;
t377 = t387 * qJ(2);
t380 = t389 * pkin(1);
t410 = t380 + t377;
t361 = qJD(1) * t410;
t363 = t383 * pkin(6) + t428;
t413 = qJD(1) * (t363 * t389 - t380) + t361;
t350 = t385 * t383 + t384 * t429 + t428;
t412 = t350 - t363;
t372 = qJ(2) + t430;
t411 = pkin(4) * t373 - t372 + t430;
t409 = qJD(3) * t383;
t408 = qJD(4) * t383;
t407 = qJD(3) + qJD(5);
t406 = -pkin(2) - t429;
t405 = t387 * t409;
t313 = (qJ(2) - t372) * t389 + t412 * t387;
t403 = t389 * t313 * t409 - qJD(4) * t384;
t402 = t383 * t407;
t314 = t372 * t387 + t389 * t412 - t377;
t369 = -qJD(3) * t384 + qJD(1);
t401 = t369 * t314 + t387 * t408 + t413;
t398 = rSges(3,1) * t384 - rSges(3,2) * t383;
t379 = t387 * pkin(1);
t364 = -t389 * qJ(2) + t379;
t376 = qJD(2) * t387;
t395 = t376 + (-t363 * t387 - t364 + t379) * qJD(1);
t359 = pkin(4) * t374 - t406;
t381 = -pkin(7) - t385;
t394 = t359 * t384 - t381 * t383 - t350;
t351 = t383 * t429 + (pkin(6) - t385) * t384;
t391 = t351 * t405 + t389 * t408 + t395;
t366 = t389 * rSges(2,1) - t387 * rSges(2,2);
t365 = t387 * rSges(2,1) + t389 * rSges(2,2);
t360 = -t384 * t407 + qJD(1);
t354 = t389 * t402;
t353 = t387 * t402;
t345 = -t384 * rSges(4,3) + (rSges(4,1) * t388 - rSges(4,2) * t386) * t383;
t344 = -Icges(4,5) * t384 + (Icges(4,1) * t388 - Icges(4,4) * t386) * t383;
t343 = -Icges(4,6) * t384 + (Icges(4,4) * t388 - Icges(4,2) * t386) * t383;
t340 = t384 * t419 + t427;
t339 = -t384 * t420 + t426;
t338 = t384 * t426 - t420;
t337 = -t384 * t427 - t419;
t336 = -t384 * rSges(5,3) + (rSges(5,1) * t374 - rSges(5,2) * t373) * t383;
t335 = -Icges(5,5) * t384 + (Icges(5,1) * t374 - Icges(5,4) * t373) * t383;
t334 = -Icges(5,6) * t384 + (Icges(5,4) * t374 - Icges(5,2) * t373) * t383;
t331 = -t384 * rSges(6,3) + (rSges(6,1) * t371 - rSges(6,2) * t370) * t383;
t330 = -Icges(6,5) * t384 + (Icges(6,1) * t371 - Icges(6,4) * t370) * t383;
t329 = -Icges(6,6) * t384 + (Icges(6,4) * t371 - Icges(6,2) * t370) * t383;
t328 = -Icges(6,3) * t384 + (Icges(6,5) * t371 - Icges(6,6) * t370) * t383;
t327 = qJD(1) * t387 * rSges(3,3) + t361 + (qJD(1) * t398 - qJD(2)) * t389;
t326 = t376 + (t389 * rSges(3,3) - t387 * t398 - t364) * qJD(1);
t325 = (t381 + t385) * t384 + (t359 + t406) * t383;
t324 = t358 * rSges(4,1) + t357 * rSges(4,2) + rSges(4,3) * t416;
t323 = t356 * rSges(4,1) + t355 * rSges(4,2) + rSges(4,3) * t423;
t322 = Icges(4,1) * t358 + Icges(4,4) * t357 + Icges(4,5) * t416;
t321 = Icges(4,1) * t356 + Icges(4,4) * t355 + Icges(4,5) * t423;
t320 = Icges(4,4) * t358 + Icges(4,2) * t357 + Icges(4,6) * t416;
t319 = Icges(4,4) * t356 + Icges(4,2) * t355 + Icges(4,6) * t423;
t316 = t349 * rSges(5,1) + t348 * rSges(5,2) + rSges(5,3) * t416;
t315 = t347 * rSges(5,1) + t346 * rSges(5,2) + rSges(5,3) * t423;
t312 = Icges(5,1) * t349 + Icges(5,4) * t348 + Icges(5,5) * t416;
t311 = Icges(5,1) * t347 + Icges(5,4) * t346 + Icges(5,5) * t423;
t310 = Icges(5,4) * t349 + Icges(5,2) * t348 + Icges(5,6) * t416;
t309 = Icges(5,4) * t347 + Icges(5,2) * t346 + Icges(5,6) * t423;
t306 = t340 * rSges(6,1) + t339 * rSges(6,2) + rSges(6,3) * t416;
t305 = t338 * rSges(6,1) + t337 * rSges(6,2) + rSges(6,3) * t423;
t304 = Icges(6,1) * t340 + Icges(6,4) * t339 + Icges(6,5) * t416;
t303 = Icges(6,1) * t338 + Icges(6,4) * t337 + Icges(6,5) * t423;
t302 = Icges(6,4) * t340 + Icges(6,2) * t339 + Icges(6,6) * t416;
t301 = Icges(6,4) * t338 + Icges(6,2) * t337 + Icges(6,6) * t423;
t300 = Icges(6,5) * t340 + Icges(6,6) * t339 + Icges(6,3) * t416;
t299 = Icges(6,5) * t338 + Icges(6,6) * t337 + Icges(6,3) * t423;
t296 = t387 * t411 + t389 * t394 + t410;
t295 = t379 + (-qJ(2) - t411) * t389 + t394 * t387;
t294 = (t323 * t389 - t324 * t387) * t409;
t293 = t369 * t324 + (-t345 * t409 - qJD(2)) * t389 + t413;
t292 = -t369 * t323 + t345 * t405 + t395;
t291 = t369 * t316 + (-qJD(2) + (-t336 - t351) * t409) * t389 + t401;
t290 = t336 * t405 + (-t313 - t315) * t369 + t391;
t289 = (t315 * t389 + (-t314 - t316) * t387) * t409 + t403;
t288 = t369 * t296 + t360 * t306 - t354 * t331 + (-qJD(2) + (-t325 - t351) * t409) * t389 + t401;
t287 = t325 * t405 - t360 * t305 + t353 * t331 + (-t295 - t313) * t369 + t391;
t286 = t354 * t305 - t353 * t306 + (t295 * t389 + (-t296 - t314) * t387) * t409 + t403;
t1 = m(3) * (t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(4) * (t292 ^ 2 + t293 ^ 2 + t294 ^ 2) / 0.2e1 + m(5) * (t289 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + m(6) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + t354 * ((t300 * t416 + t339 * t302 + t340 * t304) * t354 + (t299 * t416 + t339 * t301 + t340 * t303) * t353 + (t328 * t416 + t339 * t329 + t340 * t330) * t360) / 0.2e1 + t353 * ((t300 * t423 + t337 * t302 + t338 * t304) * t354 + (t299 * t423 + t337 * t301 + t338 * t303) * t353 + (t328 * t423 + t337 * t329 + t338 * t330) * t360) / 0.2e1 + t360 * ((-t299 * t353 - t300 * t354 - t328 * t360) * t384 + ((-t302 * t370 + t304 * t371) * t354 + (-t301 * t370 + t303 * t371) * t353 + (-t329 * t370 + t330 * t371) * t360) * t383) / 0.2e1 + ((-t435 * t384 + ((-t310 * t373 + t312 * t374 - t320 * t386 + t322 * t388) * t389 + (-t309 * t373 + t311 * t374 - t319 * t386 + t321 * t388) * t387) * t383) * t409 + (t434 * t384 + (-t334 * t373 + t335 * t374 - t343 * t386 + t344 * t388) * t383) * t369) * t369 / 0.2e1 + (m(2) * (t365 ^ 2 + t366 ^ 2) + Icges(2,3) + Icges(3,2) * t384 ^ 2 + (Icges(3,1) * t383 + 0.2e1 * Icges(3,4) * t384) * t383) * qJD(1) ^ 2 / 0.2e1 + ((((t346 * t310 + t347 * t312 + t355 * t320 + t356 * t322) * t389 + (t346 * t309 + t347 * t311 + t355 * t319 + t356 * t321 + t433) * t387) * t409 + (t346 * t334 + t347 * t335 + t355 * t343 + t356 * t344 - t434 * t423) * t369) * t387 + (((t348 * t310 + t349 * t312 + t357 * t320 + t358 * t322 + t433) * t389 + (t348 * t309 + t349 * t311 + t357 * t319 + t358 * t321) * t387) * t409 + (t348 * t334 + t349 * t335 + t357 * t343 + t358 * t344 - t434 * t416) * t369) * t389) * t409 / 0.2e1;
T = t1;
