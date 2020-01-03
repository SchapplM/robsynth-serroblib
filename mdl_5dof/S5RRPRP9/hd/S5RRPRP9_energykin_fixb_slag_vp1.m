% Calculate kinetic energy for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:40
% EndTime: 2019-12-31 20:05:42
% DurationCPUTime: 2.03s
% Computational Cost: add. (1062->230), mult. (1760->355), div. (0->0), fcn. (1779->8), ass. (0->119)
t434 = Icges(5,1) + Icges(6,1);
t433 = -Icges(5,4) + Icges(6,5);
t432 = Icges(6,4) + Icges(5,5);
t431 = Icges(5,2) + Icges(6,3);
t430 = -Icges(6,6) + Icges(5,6);
t429 = -Icges(5,3) - Icges(6,2);
t428 = rSges(6,1) + pkin(4);
t427 = rSges(6,3) + qJ(5);
t368 = pkin(8) + qJ(4);
t366 = sin(t368);
t367 = cos(t368);
t375 = cos(qJ(1));
t373 = sin(qJ(1));
t374 = cos(qJ(2));
t405 = t374 * t373;
t328 = t366 * t405 + t367 * t375;
t329 = -t366 * t375 + t367 * t405;
t372 = sin(qJ(2));
t407 = t372 * t373;
t426 = t431 * t328 + t433 * t329 - t430 * t407;
t404 = t374 * t375;
t330 = t366 * t404 - t373 * t367;
t331 = t366 * t373 + t367 * t404;
t406 = t372 * t375;
t425 = t431 * t330 + t433 * t331 - t430 * t406;
t424 = -t430 * t328 + t432 * t329 - t429 * t407;
t423 = -t430 * t330 + t432 * t331 - t429 * t406;
t422 = t433 * t328 + t434 * t329 + t432 * t407;
t421 = t433 * t330 + t434 * t331 + t432 * t406;
t420 = t430 * t374 + (t431 * t366 + t433 * t367) * t372;
t419 = t429 * t374 + (-t430 * t366 + t432 * t367) * t372;
t418 = -t432 * t374 + (t433 * t366 + t434 * t367) * t372;
t370 = cos(pkin(8));
t412 = pkin(3) * t370;
t411 = Icges(3,4) * t372;
t410 = Icges(3,4) * t374;
t369 = sin(pkin(8));
t409 = t369 * t373;
t408 = t369 * t375;
t402 = rSges(6,2) * t407 + t427 * t328 + t428 * t329;
t401 = rSges(6,2) * t406 + t427 * t330 + t428 * t331;
t400 = -rSges(6,2) * t374 + (t427 * t366 + t428 * t367) * t372;
t387 = pkin(2) * t374 + qJ(3) * t372;
t348 = t387 * t373;
t360 = pkin(1) * t373 - pkin(6) * t375;
t399 = -t348 - t360;
t398 = qJD(2) * t373;
t397 = qJD(2) * t375;
t396 = qJD(3) * t372;
t395 = qJD(4) * t372;
t349 = t387 * t375;
t352 = qJD(1) * (pkin(1) * t375 + pkin(6) * t373);
t394 = qJD(1) * t349 + t373 * t396 + t352;
t356 = pkin(2) * t372 - qJ(3) * t374;
t391 = qJD(2) * (pkin(7) * t374 - t372 * t412 - t356);
t390 = qJD(2) * (rSges(4,3) * t374 - (rSges(4,1) * t370 - rSges(4,2) * t369) * t372 - t356);
t389 = -qJD(3) * t374 + t348 * t398 + t349 * t397;
t388 = rSges(3,1) * t374 - rSges(3,2) * t372;
t386 = Icges(3,1) * t374 - t411;
t385 = -Icges(3,2) * t372 + t410;
t384 = Icges(3,5) * t374 - Icges(3,6) * t372;
t334 = -Icges(3,6) * t375 + t373 * t385;
t336 = -Icges(3,5) * t375 + t373 * t386;
t383 = t334 * t372 - t336 * t374;
t335 = Icges(3,6) * t373 + t375 * t385;
t337 = Icges(3,5) * t373 + t375 * t386;
t382 = -t335 * t372 + t337 * t374;
t354 = Icges(3,2) * t374 + t411;
t355 = Icges(3,1) * t372 + t410;
t381 = -t354 * t372 + t355 * t374;
t379 = pkin(7) * t372 + t374 * t412;
t313 = -pkin(3) * t408 + t373 * t379;
t314 = pkin(3) * t409 + t375 * t379;
t380 = t313 * t398 + t314 * t397 + t389;
t378 = qJD(1) * t314 + t373 * t391 + t394;
t363 = t375 * t396;
t377 = t363 + (-t313 + t399) * qJD(1) + t375 * t391;
t364 = -qJD(4) * t374 + qJD(1);
t359 = rSges(2,1) * t375 - rSges(2,2) * t373;
t358 = rSges(2,1) * t373 + rSges(2,2) * t375;
t357 = rSges(3,1) * t372 + rSges(3,2) * t374;
t353 = Icges(3,5) * t372 + Icges(3,6) * t374;
t351 = t373 * t395 - t397;
t350 = t375 * t395 + t398;
t347 = t370 * t404 + t409;
t346 = -t369 * t404 + t370 * t373;
t345 = t370 * t405 - t408;
t344 = -t369 * t405 - t370 * t375;
t342 = rSges(3,3) * t373 + t375 * t388;
t341 = -rSges(3,3) * t375 + t373 * t388;
t333 = Icges(3,3) * t373 + t375 * t384;
t332 = -Icges(3,3) * t375 + t373 * t384;
t326 = -Icges(4,5) * t374 + (Icges(4,1) * t370 - Icges(4,4) * t369) * t372;
t325 = -Icges(4,6) * t374 + (Icges(4,4) * t370 - Icges(4,2) * t369) * t372;
t324 = -Icges(4,3) * t374 + (Icges(4,5) * t370 - Icges(4,6) * t369) * t372;
t323 = -rSges(5,3) * t374 + (rSges(5,1) * t367 - rSges(5,2) * t366) * t372;
t312 = rSges(4,1) * t347 + rSges(4,2) * t346 + rSges(4,3) * t406;
t311 = rSges(4,1) * t345 + rSges(4,2) * t344 + rSges(4,3) * t407;
t310 = Icges(4,1) * t347 + Icges(4,4) * t346 + Icges(4,5) * t406;
t309 = Icges(4,1) * t345 + Icges(4,4) * t344 + Icges(4,5) * t407;
t308 = Icges(4,4) * t347 + Icges(4,2) * t346 + Icges(4,6) * t406;
t307 = Icges(4,4) * t345 + Icges(4,2) * t344 + Icges(4,6) * t407;
t306 = Icges(4,5) * t347 + Icges(4,6) * t346 + Icges(4,3) * t406;
t305 = Icges(4,5) * t345 + Icges(4,6) * t344 + Icges(4,3) * t407;
t301 = qJD(1) * t342 - t357 * t398 + t352;
t300 = -t357 * t397 + (-t341 - t360) * qJD(1);
t297 = (t341 * t373 + t342 * t375) * qJD(2);
t296 = rSges(5,1) * t331 - rSges(5,2) * t330 + rSges(5,3) * t406;
t294 = rSges(5,1) * t329 - rSges(5,2) * t328 + rSges(5,3) * t407;
t280 = qJD(1) * t312 + t373 * t390 + t394;
t279 = t363 + t375 * t390 + (-t311 + t399) * qJD(1);
t278 = (t311 * t373 + t312 * t375) * qJD(2) + t389;
t277 = t296 * t364 - t323 * t350 + t378;
t276 = -t294 * t364 + t323 * t351 + t377;
t275 = t294 * t350 - t296 * t351 + t380;
t274 = qJD(5) * t328 - t350 * t400 + t364 * t401 + t378;
t273 = qJD(5) * t330 + t351 * t400 - t364 * t402 + t377;
t272 = qJD(5) * t366 * t372 + t350 * t402 - t351 * t401 + t380;
t1 = m(3) * (t297 ^ 2 + t300 ^ 2 + t301 ^ 2) / 0.2e1 + m(4) * (t278 ^ 2 + t279 ^ 2 + t280 ^ 2) / 0.2e1 + m(5) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + m(6) * (t272 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + ((t420 * t330 + t418 * t331 + t419 * t406) * t364 + (t426 * t330 + t422 * t331 + t424 * t406) * t351 + (t425 * t330 + t421 * t331 + t423 * t406) * t350) * t350 / 0.2e1 + ((t420 * t328 + t418 * t329 + t419 * t407) * t364 + (t426 * t328 + t422 * t329 + t424 * t407) * t351 + (t425 * t328 + t421 * t329 + t423 * t407) * t350) * t351 / 0.2e1 + ((-t423 * t350 - t424 * t351 - t419 * t364) * t374 + ((t420 * t366 + t418 * t367) * t364 + (t426 * t366 + t422 * t367) * t351 + (t425 * t366 + t421 * t367) * t350) * t372) * t364 / 0.2e1 + (m(2) * (t358 ^ 2 + t359 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t335 * t374 + t337 * t372) * t373 - (t334 * t374 + t336 * t372) * t375 + (t305 * t375 - t306 * t373) * t374 + ((-t308 * t369 + t310 * t370) * t373 - (-t307 * t369 + t309 * t370) * t375) * t372) * qJD(2) + ((t354 - t324) * t374 + (-t325 * t369 + t326 * t370 + t355) * t372) * qJD(1)) * qJD(1) / 0.2e1 + (((-t305 * t406 - t307 * t346 - t309 * t347 + t383 * t375) * t375 + ((-t332 + t382) * t375 + t306 * t406 + t346 * t308 + t347 * t310 + t373 * t333) * t373) * qJD(2) + (t324 * t406 + t325 * t346 + t326 * t347 + t373 * t353 + t375 * t381) * qJD(1)) * t398 / 0.2e1 - (((-t305 * t407 - t344 * t307 - t345 * t309 + t375 * t332) * t375 + ((-t333 + t383) * t375 + t306 * t407 + t308 * t344 + t310 * t345 + t382 * t373) * t373) * qJD(2) + (t324 * t407 + t325 * t344 + t326 * t345 - t375 * t353 + t373 * t381) * qJD(1)) * t397 / 0.2e1;
T = t1;
