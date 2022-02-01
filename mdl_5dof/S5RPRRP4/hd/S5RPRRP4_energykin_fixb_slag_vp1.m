% Calculate kinetic energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:46
% EndTime: 2022-01-23 09:31:48
% DurationCPUTime: 1.48s
% Computational Cost: add. (1078->216), mult. (1635->315), div. (0->0), fcn. (1673->8), ass. (0->110)
t432 = Icges(5,1) + Icges(6,1);
t431 = Icges(5,4) + Icges(6,4);
t430 = -Icges(6,5) - Icges(5,5);
t429 = Icges(5,2) + Icges(6,2);
t428 = -Icges(6,6) - Icges(5,6);
t427 = -Icges(6,3) - Icges(5,3);
t372 = qJ(3) + qJ(4);
t365 = sin(t372);
t366 = cos(t372);
t378 = cos(qJ(1));
t374 = cos(pkin(8));
t376 = sin(qJ(1));
t408 = t374 * t376;
t342 = -t365 * t408 - t366 * t378;
t343 = -t365 * t378 + t366 * t408;
t373 = sin(pkin(8));
t410 = t373 * t376;
t426 = -t428 * t342 - t430 * t343 - t427 * t410;
t407 = t374 * t378;
t344 = -t365 * t407 + t366 * t376;
t345 = t365 * t376 + t366 * t407;
t409 = t373 * t378;
t425 = -t428 * t344 - t430 * t345 - t427 * t409;
t424 = t429 * t342 + t431 * t343 - t428 * t410;
t423 = t429 * t344 + t431 * t345 - t428 * t409;
t422 = t431 * t342 + t432 * t343 - t430 * t410;
t421 = t431 * t344 + t432 * t345 - t430 * t409;
t420 = t427 * t374 + (t428 * t365 - t430 * t366) * t373;
t419 = t428 * t374 + (-t429 * t365 + t431 * t366) * t373;
t418 = t430 * t374 + (-t431 * t365 + t432 * t366) * t373;
t379 = pkin(7) + pkin(6);
t375 = sin(qJ(3));
t413 = pkin(3) * t375;
t377 = cos(qJ(3));
t412 = pkin(3) * t377;
t411 = pkin(2) * t374 + pkin(1);
t406 = t375 * t376;
t405 = t375 * t378;
t404 = t376 * t377;
t403 = t377 * t378;
t369 = t376 * pkin(1);
t341 = t373 * t379 + t374 * t412 + t411;
t391 = -pkin(2) - t412;
t355 = pkin(4) * t366 - t391;
t371 = -qJ(5) - t379;
t383 = t355 * t374 - t371 * t373 - t341;
t363 = qJ(2) + t413;
t397 = pkin(4) * t365 - t363 + t413;
t402 = t369 + (-qJ(2) - t397) * t378 + t383 * t376 + rSges(6,1) * t343 + rSges(6,2) * t342 + rSges(6,3) * t410;
t367 = t376 * qJ(2);
t370 = t378 * pkin(1);
t396 = t370 + t367;
t401 = rSges(6,1) * t345 + rSges(6,2) * t344 + rSges(6,3) * t409 + t376 * t397 + t378 * t383 + t396;
t400 = (t371 + t379 - rSges(6,3)) * t374 + (rSges(6,1) * t366 - rSges(6,2) * t365 + t355 + t391) * t373;
t356 = qJD(1) * t396;
t358 = pkin(6) * t373 + t411;
t399 = qJD(1) * (t358 * t378 - t370) + t356;
t398 = t341 - t358;
t395 = qJD(3) * t373;
t394 = qJD(3) * t378;
t393 = qJD(3) + qJD(4);
t299 = t363 * t376 + t378 * t398 - t367;
t362 = -t374 * qJD(3) + qJD(1);
t392 = t362 * t299 + t399;
t390 = t376 * t395;
t388 = t373 * t393;
t387 = rSges(3,1) * t374 - rSges(3,2) * t373;
t350 = -t374 * t406 - t403;
t351 = t374 * t404 - t405;
t352 = -t374 * t405 + t404;
t353 = t374 * t403 + t406;
t386 = (Icges(4,5) * t351 + Icges(4,6) * t350 + Icges(4,3) * t410) * t376 + (Icges(4,5) * t353 + Icges(4,6) * t352 + Icges(4,3) * t409) * t378;
t359 = -qJ(2) * t378 + t369;
t364 = qJD(2) * t376;
t385 = t364 + (-t358 * t376 - t359 + t369) * qJD(1);
t298 = (qJ(2) - t363) * t378 + t398 * t376;
t384 = t373 * t298 * t394 - t299 * t390;
t382 = t386 * t373;
t347 = t373 * t412 + (pkin(6) - t379) * t374;
t381 = -t298 * t362 + t347 * t390 + t385;
t361 = rSges(2,1) * t378 - rSges(2,2) * t376;
t360 = rSges(2,1) * t376 + rSges(2,2) * t378;
t354 = -t374 * t393 + qJD(1);
t349 = t378 * t388;
t348 = t376 * t388;
t340 = -rSges(4,3) * t374 + (rSges(4,1) * t377 - rSges(4,2) * t375) * t373;
t339 = -Icges(4,5) * t374 + (Icges(4,1) * t377 - Icges(4,4) * t375) * t373;
t338 = -Icges(4,6) * t374 + (Icges(4,4) * t377 - Icges(4,2) * t375) * t373;
t337 = -Icges(4,3) * t374 + (Icges(4,5) * t377 - Icges(4,6) * t375) * t373;
t335 = -rSges(5,3) * t374 + (rSges(5,1) * t366 - rSges(5,2) * t365) * t373;
t326 = qJD(1) * t376 * rSges(3,3) + t356 + (qJD(1) * t387 - qJD(2)) * t378;
t325 = t364 + (t378 * rSges(3,3) - t376 * t387 - t359) * qJD(1);
t323 = rSges(4,1) * t353 + rSges(4,2) * t352 + rSges(4,3) * t409;
t322 = rSges(4,1) * t351 + rSges(4,2) * t350 + rSges(4,3) * t410;
t321 = Icges(4,1) * t353 + Icges(4,4) * t352 + Icges(4,5) * t409;
t320 = Icges(4,1) * t351 + Icges(4,4) * t350 + Icges(4,5) * t410;
t319 = Icges(4,4) * t353 + Icges(4,2) * t352 + Icges(4,6) * t409;
t318 = Icges(4,4) * t351 + Icges(4,2) * t350 + Icges(4,6) * t410;
t315 = rSges(5,1) * t345 + rSges(5,2) * t344 + rSges(5,3) * t409;
t313 = rSges(5,1) * t343 + rSges(5,2) * t342 + rSges(5,3) * t410;
t293 = (t322 * t378 - t323 * t376) * t395;
t292 = t323 * t362 + (-t340 * t395 - qJD(2)) * t378 + t399;
t291 = -t322 * t362 + t340 * t390 + t385;
t290 = t315 * t354 - t335 * t349 + (-t347 * t395 - qJD(2)) * t378 + t392;
t289 = -t313 * t354 + t335 * t348 + t381;
t288 = t313 * t349 - t315 * t348 + t384;
t287 = -qJD(2) * t378 + (qJD(5) * t376 - t347 * t394) * t373 + t401 * t354 - t400 * t349 + t392;
t286 = qJD(5) * t409 + t348 * t400 - t354 * t402 + t381;
t285 = -qJD(5) * t374 - t348 * t401 + t349 * t402 + t384;
t1 = m(3) * (t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(4) * (t291 ^ 2 + t292 ^ 2 + t293 ^ 2) / 0.2e1 + t362 * ((-t374 * t337 + (-t338 * t375 + t339 * t377) * t373) * t362 + (((-t319 * t375 + t321 * t377) * t378 + (-t318 * t375 + t320 * t377) * t376) * t373 - t386 * t374) * t395) / 0.2e1 + m(5) * (t288 ^ 2 + t289 ^ 2 + t290 ^ 2) / 0.2e1 + m(6) * (t285 ^ 2 + t286 ^ 2 + t287 ^ 2) / 0.2e1 + ((t342 * t419 + t343 * t418 + t410 * t420) * t354 + (t342 * t423 + t343 * t421 + t410 * t425) * t349 + (t424 * t342 + t422 * t343 + t426 * t410) * t348) * t348 / 0.2e1 + ((t344 * t419 + t345 * t418 + t409 * t420) * t354 + (t423 * t344 + t421 * t345 + t425 * t409) * t349 + (t424 * t344 + t422 * t345 + t409 * t426) * t348) * t349 / 0.2e1 + ((-t348 * t426 - t425 * t349 - t420 * t354) * t374 + ((-t365 * t419 + t366 * t418) * t354 + (-t365 * t423 + t366 * t421) * t349 + (-t365 * t424 + t366 * t422) * t348) * t373) * t354 / 0.2e1 + (t378 * ((t337 * t409 + t338 * t352 + t339 * t353) * t362 + ((t318 * t352 + t320 * t353) * t376 + (t352 * t319 + t353 * t321 + t382) * t378) * t395) + t376 * ((t337 * t410 + t338 * t350 + t339 * t351) * t362 + ((t319 * t350 + t321 * t351) * t378 + (t350 * t318 + t351 * t320 + t382) * t376) * t395)) * t395 / 0.2e1 + (m(2) * (t360 ^ 2 + t361 ^ 2) + Icges(2,3) + Icges(3,2) * t374 ^ 2 + (Icges(3,1) * t373 + 0.2e1 * Icges(3,4) * t374) * t373) * qJD(1) ^ 2 / 0.2e1;
T = t1;
