% Calculate kinetic energy for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:47
% EndTime: 2019-03-09 02:22:49
% DurationCPUTime: 1.26s
% Computational Cost: add. (1233->233), mult. (1301->375), div. (0->0), fcn. (1262->10), ass. (0->118)
t388 = sin(qJ(4));
t391 = cos(qJ(4));
t390 = cos(qJ(5));
t423 = t390 * pkin(5);
t426 = -pkin(9) * t391 + t423 * t388;
t389 = sin(qJ(1));
t424 = t389 * pkin(1);
t421 = Icges(5,4) * t388;
t420 = Icges(5,4) * t391;
t385 = qJ(1) + pkin(10);
t380 = sin(t385);
t387 = sin(qJ(5));
t419 = t380 * t387;
t418 = t380 * t391;
t381 = cos(t385);
t417 = t381 * t387;
t416 = t381 * t391;
t386 = qJ(5) + qJ(6);
t383 = sin(t386);
t415 = t383 * t388;
t384 = cos(t386);
t414 = t384 * t388;
t413 = t387 * t388;
t412 = t388 * t390;
t392 = cos(qJ(1));
t379 = qJD(1) * t392 * pkin(1);
t411 = qJD(1) * (t381 * pkin(2) + t380 * qJ(3)) + t379;
t374 = qJD(4) * t380;
t410 = qJD(5) * t391;
t359 = t381 * t410 + t374;
t375 = qJD(4) * t381;
t377 = qJD(5) * t388 + qJD(1);
t409 = qJD(1) * t381 * pkin(7) + t411;
t408 = -t380 * pkin(2) + t381 * qJ(3) - t424;
t407 = pkin(4) * t388 - pkin(8) * t391;
t406 = rSges(5,1) * t388 + rSges(5,2) * t391;
t405 = Icges(5,1) * t388 + t420;
t404 = Icges(5,2) * t391 + t421;
t403 = Icges(5,5) * t388 + Icges(5,6) * t391;
t330 = Icges(5,6) * t381 + t404 * t380;
t332 = Icges(5,5) * t381 + t405 * t380;
t402 = -t330 * t391 - t332 * t388;
t331 = Icges(5,6) * t380 - t404 * t381;
t333 = Icges(5,5) * t380 - t405 * t381;
t401 = t331 * t391 + t333 * t388;
t366 = -Icges(5,2) * t388 + t420;
t367 = Icges(5,1) * t391 - t421;
t400 = t366 * t391 + t367 * t388;
t399 = -pkin(7) * t380 + t408;
t357 = t407 * t380;
t358 = t407 * t381;
t398 = -t357 * t374 - t358 * t375 + qJD(2);
t372 = t391 * pkin(4) + t388 * pkin(8);
t397 = qJD(1) * t357 + (-qJD(4) * t372 - qJD(3)) * t381 + t409;
t376 = qJD(3) * t380;
t396 = t372 * t374 + t376 + (t358 + t399) * qJD(1);
t394 = qJD(2) ^ 2;
t371 = t392 * rSges(2,1) - t389 * rSges(2,2);
t370 = t391 * rSges(5,1) - t388 * rSges(5,2);
t369 = t389 * rSges(2,1) + t392 * rSges(2,2);
t365 = Icges(5,5) * t391 - Icges(5,6) * t388;
t364 = qJD(6) * t388 + t377;
t360 = -t380 * t410 + t375;
t356 = t388 * rSges(6,3) + (rSges(6,1) * t390 - rSges(6,2) * t387) * t391;
t355 = Icges(6,5) * t388 + (Icges(6,1) * t390 - Icges(6,4) * t387) * t391;
t354 = Icges(6,6) * t388 + (Icges(6,4) * t390 - Icges(6,2) * t387) * t391;
t353 = Icges(6,3) * t388 + (Icges(6,5) * t390 - Icges(6,6) * t387) * t391;
t352 = -t381 * t412 + t419;
t351 = t380 * t390 + t381 * t413;
t350 = t380 * t412 + t417;
t349 = -t380 * t413 + t381 * t390;
t347 = t379 + qJD(1) * (t381 * rSges(3,1) - t380 * rSges(3,2));
t346 = (-t380 * rSges(3,1) - t381 * rSges(3,2) - t424) * qJD(1);
t345 = t388 * rSges(7,3) + (rSges(7,1) * t384 - rSges(7,2) * t383) * t391;
t344 = Icges(7,5) * t388 + (Icges(7,1) * t384 - Icges(7,4) * t383) * t391;
t343 = Icges(7,6) * t388 + (Icges(7,4) * t384 - Icges(7,2) * t383) * t391;
t342 = Icges(7,3) * t388 + (Icges(7,5) * t384 - Icges(7,6) * t383) * t391;
t341 = t380 * t383 - t381 * t414;
t340 = t380 * t384 + t381 * t415;
t339 = t380 * t414 + t381 * t383;
t338 = -t380 * t415 + t381 * t384;
t337 = pkin(9) * t388 + t423 * t391;
t335 = t380 * rSges(5,3) - t406 * t381;
t334 = t381 * rSges(5,3) + t406 * t380;
t329 = Icges(5,3) * t380 - t403 * t381;
t328 = Icges(5,3) * t381 + t403 * t380;
t327 = t375 + (-qJD(5) - qJD(6)) * t418;
t326 = qJD(6) * t416 + t359;
t325 = -qJD(3) * t381 + qJD(1) * (-t381 * rSges(4,2) + t380 * rSges(4,3)) + t411;
t324 = t376 + (t380 * rSges(4,2) + t381 * rSges(4,3) + t408) * qJD(1);
t323 = pkin(5) * t419 - t381 * t426;
t322 = pkin(5) * t417 + t380 * t426;
t321 = t352 * rSges(6,1) + t351 * rSges(6,2) + rSges(6,3) * t416;
t320 = t350 * rSges(6,1) + t349 * rSges(6,2) - rSges(6,3) * t418;
t319 = Icges(6,1) * t352 + Icges(6,4) * t351 + Icges(6,5) * t416;
t318 = Icges(6,1) * t350 + Icges(6,4) * t349 - Icges(6,5) * t418;
t317 = Icges(6,4) * t352 + Icges(6,2) * t351 + Icges(6,6) * t416;
t316 = Icges(6,4) * t350 + Icges(6,2) * t349 - Icges(6,6) * t418;
t315 = Icges(6,5) * t352 + Icges(6,6) * t351 + Icges(6,3) * t416;
t314 = Icges(6,5) * t350 + Icges(6,6) * t349 - Icges(6,3) * t418;
t313 = t341 * rSges(7,1) + t340 * rSges(7,2) + rSges(7,3) * t416;
t312 = t339 * rSges(7,1) + t338 * rSges(7,2) - rSges(7,3) * t418;
t311 = Icges(7,1) * t341 + Icges(7,4) * t340 + Icges(7,5) * t416;
t310 = Icges(7,1) * t339 + Icges(7,4) * t338 - Icges(7,5) * t418;
t309 = Icges(7,4) * t341 + Icges(7,2) * t340 + Icges(7,6) * t416;
t308 = Icges(7,4) * t339 + Icges(7,2) * t338 - Icges(7,6) * t418;
t307 = Icges(7,5) * t341 + Icges(7,6) * t340 + Icges(7,3) * t416;
t306 = Icges(7,5) * t339 + Icges(7,6) * t338 - Icges(7,3) * t418;
t305 = qJD(2) + (-t334 * t380 + t335 * t381) * qJD(4);
t304 = qJD(1) * t334 + (-qJD(4) * t370 - qJD(3)) * t381 + t409;
t303 = t370 * t374 + t376 + (-t335 + t399) * qJD(1);
t302 = t377 * t320 - t360 * t356 + t397;
t301 = -t377 * t321 + t359 * t356 + t396;
t300 = -t359 * t320 + t360 * t321 + t398;
t299 = t364 * t312 + t377 * t322 - t327 * t345 - t360 * t337 + t397;
t298 = -t364 * t313 - t377 * t323 + t326 * t345 + t359 * t337 + t396;
t297 = -t326 * t312 + t327 * t313 - t359 * t322 + t360 * t323 + t398;
t1 = ((t381 * t365 + t400 * t380) * qJD(1) + (t381 ^ 2 * t328 + (t401 * t380 + (t329 - t402) * t381) * t380) * qJD(4)) * t375 / 0.2e1 + m(7) * (t297 ^ 2 + t298 ^ 2 + t299 ^ 2) / 0.2e1 + m(6) * (t300 ^ 2 + t301 ^ 2 + t302 ^ 2) / 0.2e1 + m(5) * (t303 ^ 2 + t304 ^ 2 + t305 ^ 2) / 0.2e1 + m(3) * (t346 ^ 2 + t347 ^ 2 + t394) / 0.2e1 + m(4) * (t324 ^ 2 + t325 ^ 2 + t394) / 0.2e1 + t360 * ((-t314 * t418 + t349 * t316 + t350 * t318) * t360 + (-t315 * t418 + t349 * t317 + t350 * t319) * t359 + (t349 * t354 + t350 * t355 - t353 * t418) * t377) / 0.2e1 + t359 * ((t314 * t416 + t351 * t316 + t352 * t318) * t360 + (t315 * t416 + t351 * t317 + t352 * t319) * t359 + (t351 * t354 + t352 * t355 + t353 * t416) * t377) / 0.2e1 + t364 * ((t306 * t327 + t307 * t326 + t342 * t364) * t388 + ((-t308 * t383 + t310 * t384) * t327 + (-t309 * t383 + t311 * t384) * t326 + (-t343 * t383 + t344 * t384) * t364) * t391) / 0.2e1 + t327 * ((-t306 * t418 + t338 * t308 + t339 * t310) * t327 + (-t307 * t418 + t338 * t309 + t339 * t311) * t326 + (t338 * t343 + t339 * t344 - t342 * t418) * t364) / 0.2e1 + t326 * ((t306 * t416 + t340 * t308 + t341 * t310) * t327 + (t307 * t416 + t340 * t309 + t341 * t311) * t326 + (t340 * t343 + t341 * t344 + t342 * t416) * t364) / 0.2e1 + t377 * ((t314 * t360 + t315 * t359 + t353 * t377) * t388 + ((-t316 * t387 + t318 * t390) * t360 + (-t317 * t387 + t319 * t390) * t359 + (-t354 * t387 + t355 * t390) * t377) * t391) / 0.2e1 + ((t380 * t365 - t400 * t381) * qJD(1) + (t380 ^ 2 * t329 + (t402 * t381 + (t328 - t401) * t380) * t381) * qJD(4)) * t374 / 0.2e1 + qJD(1) * ((-t388 * t366 + t391 * t367) * qJD(1) + ((-t388 * t330 + t391 * t332) * t381 + (-t388 * t331 + t391 * t333) * t380) * qJD(4)) / 0.2e1 + (Icges(3,3) + m(2) * (t369 ^ 2 + t371 ^ 2) + Icges(4,1) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
