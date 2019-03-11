% Calculate kinetic energy for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:43:51
% EndTime: 2019-03-09 01:43:52
% DurationCPUTime: 1.38s
% Computational Cost: add. (1108->207), mult. (1003->313), div. (0->0), fcn. (908->10), ass. (0->112)
t433 = Icges(5,3) + Icges(6,3);
t369 = qJ(4) + pkin(10);
t365 = sin(t369);
t367 = cos(t369);
t373 = sin(qJ(4));
t376 = cos(qJ(4));
t432 = Icges(5,5) * t373 + Icges(6,5) * t365 + Icges(5,6) * t376 + Icges(6,6) * t367;
t370 = qJ(1) + pkin(9);
t366 = sin(t370);
t368 = cos(t370);
t431 = t432 * t366 + t433 * t368;
t430 = t433 * t366 - t432 * t368;
t429 = Icges(5,5) * t376 + Icges(6,5) * t367 - Icges(5,6) * t373 - Icges(6,6) * t365;
t414 = Icges(6,4) * t367;
t347 = -Icges(6,2) * t365 + t414;
t415 = Icges(6,4) * t365;
t348 = Icges(6,1) * t367 - t415;
t416 = Icges(5,4) * t376;
t354 = -Icges(5,2) * t373 + t416;
t417 = Icges(5,4) * t373;
t355 = Icges(5,1) * t376 - t417;
t428 = t347 * t367 + t348 * t365 + t354 * t376 + t355 * t373;
t390 = Icges(6,2) * t367 + t415;
t315 = Icges(6,6) * t366 - t368 * t390;
t392 = Icges(6,1) * t365 + t414;
t317 = Icges(6,5) * t366 - t368 * t392;
t391 = Icges(5,2) * t376 + t417;
t326 = Icges(5,6) * t366 - t368 * t391;
t393 = Icges(5,1) * t373 + t416;
t329 = Icges(5,5) * t366 - t368 * t393;
t427 = t315 * t367 + t317 * t365 + t326 * t376 + t329 * t373;
t314 = Icges(6,6) * t368 + t366 * t390;
t316 = Icges(6,5) * t368 + t366 * t392;
t325 = Icges(5,6) * t368 + t366 * t391;
t328 = Icges(5,5) * t368 + t366 * t393;
t426 = -t314 * t367 - t316 * t365 - t325 * t376 - t328 * t373;
t374 = sin(qJ(1));
t422 = pkin(1) * t374;
t421 = pkin(4) * t373;
t420 = pkin(4) * t376;
t413 = t366 * t367;
t372 = sin(qJ(6));
t412 = t366 * t372;
t375 = cos(qJ(6));
t411 = t366 * t375;
t410 = t367 * t368;
t409 = t368 * t372;
t408 = t368 * t375;
t377 = cos(qJ(1));
t364 = qJD(1) * t377 * pkin(1);
t407 = qJD(1) * (pkin(2) * t368 + qJ(3) * t366) + t364;
t406 = qJD(4) * t366;
t405 = qJD(4) * t368;
t404 = qJD(6) * t367;
t333 = qJ(5) * t366 - t368 * t421;
t403 = t333 * t405 + qJD(2);
t402 = qJD(1) * t368 * pkin(7) + t407;
t363 = qJD(3) * t366;
t401 = qJD(5) * t368 + t406 * t420 + t363;
t398 = -pkin(2) * t366 + qJ(3) * t368 - t422;
t397 = pkin(5) * t365 - pkin(8) * t367;
t396 = rSges(5,1) * t373 + rSges(5,2) * t376;
t395 = rSges(6,1) * t365 + rSges(6,2) * t367;
t334 = qJ(5) * t368 + t366 * t421;
t394 = qJD(1) * t334 + qJD(5) * t366 + t402;
t381 = -pkin(7) * t366 + t398;
t380 = -t333 + t381;
t378 = qJD(2) ^ 2;
t359 = qJD(6) * t365 + qJD(1);
t358 = rSges(2,1) * t377 - rSges(2,2) * t374;
t357 = rSges(5,1) * t376 - rSges(5,2) * t373;
t356 = rSges(2,1) * t374 + rSges(2,2) * t377;
t351 = pkin(5) * t367 + pkin(8) * t365;
t350 = rSges(6,1) * t367 - rSges(6,2) * t365;
t344 = -t366 * t404 + t405;
t343 = t368 * t404 + t406;
t342 = -t365 * t408 + t412;
t341 = t365 * t409 + t411;
t340 = t365 * t411 + t409;
t339 = -t365 * t412 + t408;
t338 = t364 + qJD(1) * (rSges(3,1) * t368 - rSges(3,2) * t366);
t337 = (-rSges(3,1) * t366 - rSges(3,2) * t368 - t422) * qJD(1);
t336 = t397 * t368;
t335 = t397 * t366;
t332 = rSges(5,3) * t366 - t368 * t396;
t331 = rSges(7,3) * t365 + (rSges(7,1) * t375 - rSges(7,2) * t372) * t367;
t330 = rSges(5,3) * t368 + t366 * t396;
t327 = Icges(7,5) * t365 + (Icges(7,1) * t375 - Icges(7,4) * t372) * t367;
t324 = Icges(7,6) * t365 + (Icges(7,4) * t375 - Icges(7,2) * t372) * t367;
t321 = Icges(7,3) * t365 + (Icges(7,5) * t375 - Icges(7,6) * t372) * t367;
t319 = rSges(6,3) * t366 - t368 * t395;
t318 = rSges(6,3) * t368 + t366 * t395;
t310 = -qJD(3) * t368 + qJD(1) * (-rSges(4,2) * t368 + rSges(4,3) * t366) + t407;
t309 = t363 + (rSges(4,2) * t366 + rSges(4,3) * t368 + t398) * qJD(1);
t308 = rSges(7,1) * t342 + rSges(7,2) * t341 + rSges(7,3) * t410;
t307 = rSges(7,1) * t340 + rSges(7,2) * t339 - rSges(7,3) * t413;
t306 = Icges(7,1) * t342 + Icges(7,4) * t341 + Icges(7,5) * t410;
t305 = Icges(7,1) * t340 + Icges(7,4) * t339 - Icges(7,5) * t413;
t304 = Icges(7,4) * t342 + Icges(7,2) * t341 + Icges(7,6) * t410;
t303 = Icges(7,4) * t340 + Icges(7,2) * t339 - Icges(7,6) * t413;
t302 = Icges(7,5) * t342 + Icges(7,6) * t341 + Icges(7,3) * t410;
t301 = Icges(7,5) * t340 + Icges(7,6) * t339 - Icges(7,3) * t413;
t300 = qJD(2) + (-t330 * t366 + t332 * t368) * qJD(4);
t299 = qJD(1) * t330 + (-qJD(4) * t357 - qJD(3)) * t368 + t402;
t298 = t357 * t406 + t363 + (-t332 + t381) * qJD(1);
t297 = qJD(1) * t318 + (-qJD(3) + (-t350 - t420) * qJD(4)) * t368 + t394;
t296 = t350 * t406 + (-t319 + t380) * qJD(1) + t401;
t295 = (t319 * t368 + (-t318 - t334) * t366) * qJD(4) + t403;
t294 = qJD(1) * t335 + t307 * t359 - t331 * t344 + (-qJD(3) + (-t351 - t420) * qJD(4)) * t368 + t394;
t293 = t351 * t406 - t308 * t359 + t331 * t343 + (t336 + t380) * qJD(1) + t401;
t292 = -t307 * t343 + t308 * t344 + (-t336 * t368 + (-t334 - t335) * t366) * qJD(4) + t403;
t1 = m(5) * (t298 ^ 2 + t299 ^ 2 + t300 ^ 2) / 0.2e1 + m(6) * (t295 ^ 2 + t296 ^ 2 + t297 ^ 2) / 0.2e1 + m(7) * (t292 ^ 2 + t293 ^ 2 + t294 ^ 2) / 0.2e1 + t344 * ((-t301 * t413 + t339 * t303 + t340 * t305) * t344 + (-t302 * t413 + t304 * t339 + t306 * t340) * t343 + (-t321 * t413 + t324 * t339 + t327 * t340) * t359) / 0.2e1 + t343 * ((t301 * t410 + t303 * t341 + t305 * t342) * t344 + (t302 * t410 + t341 * t304 + t342 * t306) * t343 + (t321 * t410 + t324 * t341 + t327 * t342) * t359) / 0.2e1 + t359 * ((t301 * t344 + t302 * t343 + t321 * t359) * t365 + ((-t303 * t372 + t305 * t375) * t344 + (-t304 * t372 + t306 * t375) * t343 + (-t324 * t372 + t327 * t375) * t359) * t367) / 0.2e1 + m(3) * (t337 ^ 2 + t338 ^ 2 + t378) / 0.2e1 + m(4) * (t309 ^ 2 + t310 ^ 2 + t378) / 0.2e1 + (((-t314 * t365 + t316 * t367 - t325 * t373 + t328 * t376) * t368 + (-t365 * t315 + t367 * t317 - t326 * t373 + t329 * t376) * t366) * qJD(4) + (-t365 * t347 + t367 * t348 - t373 * t354 + t376 * t355) * qJD(1)) * qJD(1) / 0.2e1 + ((t430 * t366 ^ 2 + (t426 * t368 + (-t427 + t431) * t366) * t368) * qJD(4) + (t366 * t429 - t368 * t428) * qJD(1)) * t406 / 0.2e1 + ((t431 * t368 ^ 2 + (t427 * t366 + (-t426 + t430) * t368) * t366) * qJD(4) + (t366 * t428 + t368 * t429) * qJD(1)) * t405 / 0.2e1 + (Icges(2,3) + Icges(3,3) + Icges(4,1) + m(2) * (t356 ^ 2 + t358 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
