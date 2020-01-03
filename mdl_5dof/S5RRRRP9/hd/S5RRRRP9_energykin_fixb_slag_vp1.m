% Calculate kinetic energy for
% S5RRRRP9
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:45
% EndTime: 2019-12-31 22:03:47
% DurationCPUTime: 2.11s
% Computational Cost: add. (1125->230), mult. (1865->369), div. (0->0), fcn. (1884->8), ass. (0->121)
t432 = Icges(5,1) + Icges(6,1);
t431 = -Icges(5,4) + Icges(6,5);
t430 = Icges(6,4) + Icges(5,5);
t429 = Icges(5,2) + Icges(6,3);
t428 = -Icges(6,6) + Icges(5,6);
t427 = -Icges(5,3) - Icges(6,2);
t426 = rSges(6,1) + pkin(4);
t425 = rSges(6,3) + qJ(5);
t372 = qJ(3) + qJ(4);
t370 = sin(t372);
t371 = cos(t372);
t378 = cos(qJ(1));
t375 = sin(qJ(1));
t377 = cos(qJ(2));
t403 = t375 * t377;
t338 = t370 * t403 + t371 * t378;
t339 = -t370 * t378 + t371 * t403;
t374 = sin(qJ(2));
t405 = t374 * t375;
t424 = t429 * t338 + t431 * t339 - t428 * t405;
t402 = t377 * t378;
t340 = t370 * t402 - t375 * t371;
t341 = t370 * t375 + t371 * t402;
t404 = t374 * t378;
t423 = t429 * t340 + t431 * t341 - t428 * t404;
t422 = -t428 * t338 + t430 * t339 - t427 * t405;
t421 = -t428 * t340 + t430 * t341 - t427 * t404;
t420 = t431 * t338 + t432 * t339 + t430 * t405;
t419 = t431 * t340 + t432 * t341 + t430 * t404;
t418 = t428 * t377 + (t429 * t370 + t431 * t371) * t374;
t417 = t427 * t377 + (-t428 * t370 + t430 * t371) * t374;
t416 = -t430 * t377 + (t431 * t370 + t432 * t371) * t374;
t376 = cos(qJ(3));
t411 = pkin(3) * t376;
t409 = Icges(3,4) * t374;
t408 = Icges(3,4) * t377;
t373 = sin(qJ(3));
t407 = t373 * t375;
t406 = t373 * t378;
t401 = rSges(6,2) * t405 + t425 * t338 + t426 * t339;
t400 = rSges(6,2) * t404 + t425 * t340 + t426 * t341;
t399 = -rSges(6,2) * t377 + (t425 * t370 + t426 * t371) * t374;
t394 = pkin(2) * t377 + pkin(7) * t374;
t350 = t394 * t375;
t351 = t394 * t378;
t369 = qJD(2) * t375;
t397 = qJD(2) * t378;
t398 = t350 * t369 + t351 * t397;
t396 = qJD(3) * t374;
t352 = t378 * t396 + t369;
t395 = qJD(4) * t374;
t353 = t375 * t396 - t397;
t393 = rSges(3,1) * t377 - rSges(3,2) * t374;
t392 = Icges(3,1) * t377 - t409;
t391 = -Icges(3,2) * t374 + t408;
t390 = Icges(3,5) * t377 - Icges(3,6) * t374;
t330 = -Icges(3,6) * t378 + t375 * t391;
t333 = -Icges(3,5) * t378 + t375 * t392;
t389 = t330 * t374 - t333 * t377;
t331 = Icges(3,6) * t375 + t378 * t391;
t334 = Icges(3,5) * t375 + t378 * t392;
t388 = -t331 * t374 + t334 * t377;
t357 = Icges(3,2) * t377 + t409;
t358 = Icges(3,1) * t374 + t408;
t387 = -t357 * t374 + t358 * t377;
t384 = pkin(8) * t374 + t377 * t411;
t312 = -pkin(3) * t406 + t375 * t384;
t313 = pkin(3) * t407 + t378 * t384;
t386 = t352 * t312 - t313 * t353 + t398;
t355 = qJD(1) * (pkin(1) * t378 + pkin(6) * t375);
t362 = pkin(2) * t374 - pkin(7) * t377;
t385 = qJD(1) * t351 - t362 * t369 + t355;
t363 = pkin(1) * t375 - pkin(6) * t378;
t383 = (-t350 - t363) * qJD(1) - t362 * t397;
t315 = -pkin(8) * t377 + t374 * t411;
t367 = -qJD(3) * t377 + qJD(1);
t382 = t367 * t313 - t315 * t352 + t385;
t381 = -t312 * t367 + t353 * t315 + t383;
t361 = rSges(2,1) * t378 - rSges(2,2) * t375;
t360 = rSges(2,1) * t375 + rSges(2,2) * t378;
t359 = rSges(3,1) * t374 + rSges(3,2) * t377;
t356 = Icges(3,5) * t374 + Icges(3,6) * t377;
t354 = qJD(1) + (-qJD(3) - qJD(4)) * t377;
t349 = t376 * t402 + t407;
t348 = -t373 * t402 + t375 * t376;
t347 = t376 * t403 - t406;
t346 = -t373 * t403 - t376 * t378;
t337 = rSges(3,3) * t375 + t378 * t393;
t336 = -rSges(3,3) * t378 + t375 * t393;
t335 = -rSges(4,3) * t377 + (rSges(4,1) * t376 - rSges(4,2) * t373) * t374;
t332 = -Icges(4,5) * t377 + (Icges(4,1) * t376 - Icges(4,4) * t373) * t374;
t329 = -Icges(4,6) * t377 + (Icges(4,4) * t376 - Icges(4,2) * t373) * t374;
t328 = Icges(3,3) * t375 + t378 * t390;
t327 = -Icges(3,3) * t378 + t375 * t390;
t326 = -Icges(4,3) * t377 + (Icges(4,5) * t376 - Icges(4,6) * t373) * t374;
t325 = t375 * t395 + t353;
t324 = t378 * t395 + t352;
t323 = -rSges(5,3) * t377 + (rSges(5,1) * t371 - rSges(5,2) * t370) * t374;
t311 = rSges(4,1) * t349 + rSges(4,2) * t348 + rSges(4,3) * t404;
t310 = rSges(4,1) * t347 + rSges(4,2) * t346 + rSges(4,3) * t405;
t309 = Icges(4,1) * t349 + Icges(4,4) * t348 + Icges(4,5) * t404;
t308 = Icges(4,1) * t347 + Icges(4,4) * t346 + Icges(4,5) * t405;
t307 = Icges(4,4) * t349 + Icges(4,2) * t348 + Icges(4,6) * t404;
t306 = Icges(4,4) * t347 + Icges(4,2) * t346 + Icges(4,6) * t405;
t305 = Icges(4,5) * t349 + Icges(4,6) * t348 + Icges(4,3) * t404;
t304 = Icges(4,5) * t347 + Icges(4,6) * t346 + Icges(4,3) * t405;
t301 = qJD(1) * t337 - t359 * t369 + t355;
t300 = -t359 * t397 + (-t336 - t363) * qJD(1);
t299 = (t336 * t375 + t337 * t378) * qJD(2);
t297 = rSges(5,1) * t341 - rSges(5,2) * t340 + rSges(5,3) * t404;
t295 = rSges(5,1) * t339 - rSges(5,2) * t338 + rSges(5,3) * t405;
t280 = t311 * t367 - t335 * t352 + t385;
t279 = -t310 * t367 + t335 * t353 + t383;
t278 = t310 * t352 - t311 * t353 + t398;
t277 = t297 * t354 - t323 * t324 + t382;
t276 = -t295 * t354 + t323 * t325 + t381;
t275 = t295 * t324 - t297 * t325 + t386;
t274 = qJD(5) * t338 - t324 * t399 + t354 * t400 + t382;
t273 = qJD(5) * t340 + t325 * t399 - t354 * t401 + t381;
t272 = qJD(5) * t370 * t374 + t324 * t401 - t325 * t400 + t386;
t1 = m(3) * (t299 ^ 2 + t300 ^ 2 + t301 ^ 2) / 0.2e1 + ((t375 * t356 + t378 * t387) * qJD(1) + (t375 ^ 2 * t328 + (t389 * t378 + (-t327 + t388) * t375) * t378) * qJD(2)) * t369 / 0.2e1 - ((-t378 * t356 + t375 * t387) * qJD(1) + (t378 ^ 2 * t327 + (t388 * t375 + (-t328 + t389) * t378) * t375) * qJD(2)) * t397 / 0.2e1 + qJD(1) * ((t377 * t357 + t374 * t358) * qJD(1) + ((t331 * t377 + t334 * t374) * t375 - (t330 * t377 + t333 * t374) * t378) * qJD(2)) / 0.2e1 + m(4) * (t278 ^ 2 + t279 ^ 2 + t280 ^ 2) / 0.2e1 + t352 * ((t305 * t404 + t348 * t307 + t349 * t309) * t352 + (t304 * t404 + t306 * t348 + t308 * t349) * t353 + (t326 * t404 + t329 * t348 + t332 * t349) * t367) / 0.2e1 + t353 * ((t305 * t405 + t307 * t346 + t309 * t347) * t352 + (t304 * t405 + t346 * t306 + t347 * t308) * t353 + (t326 * t405 + t329 * t346 + t332 * t347) * t367) / 0.2e1 + t367 * ((-t304 * t353 - t305 * t352 - t326 * t367) * t377 + ((-t307 * t373 + t309 * t376) * t352 + (-t306 * t373 + t308 * t376) * t353 + (-t329 * t373 + t332 * t376) * t367) * t374) / 0.2e1 + m(5) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + m(6) * (t272 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + ((t340 * t418 + t341 * t416 + t404 * t417) * t354 + (t340 * t424 + t420 * t341 + t422 * t404) * t325 + (t423 * t340 + t419 * t341 + t421 * t404) * t324) * t324 / 0.2e1 + ((t338 * t418 + t339 * t416 + t405 * t417) * t354 + (t424 * t338 + t420 * t339 + t422 * t405) * t325 + (t338 * t423 + t339 * t419 + t405 * t421) * t324) * t325 / 0.2e1 + ((-t324 * t421 - t325 * t422 - t354 * t417) * t377 + ((t370 * t418 + t371 * t416) * t354 + (t370 * t424 + t420 * t371) * t325 + (t370 * t423 + t371 * t419) * t324) * t374) * t354 / 0.2e1 + (m(2) * (t360 ^ 2 + t361 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
