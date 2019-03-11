% Calculate kinetic energy for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:46:14
% EndTime: 2019-03-09 01:46:15
% DurationCPUTime: 1.48s
% Computational Cost: add. (1119->207), mult. (1852->320), div. (0->0), fcn. (2174->10), ass. (0->110)
t438 = Icges(5,3) + Icges(6,3);
t375 = qJ(4) + pkin(10);
t372 = sin(t375);
t373 = cos(t375);
t378 = sin(qJ(4));
t380 = cos(qJ(4));
t437 = -Icges(5,5) * t380 - Icges(6,5) * t373 + Icges(5,6) * t378 + Icges(6,6) * t372;
t420 = sin(pkin(9));
t421 = cos(pkin(9));
t426 = sin(qJ(1));
t427 = cos(qJ(1));
t356 = -t420 * t426 - t421 * t427;
t357 = t420 * t427 - t421 * t426;
t436 = t437 * t356 + t438 * t357;
t435 = -t438 * t356 + t437 * t357;
t434 = -Icges(5,5) * t378 - Icges(6,5) * t372 - Icges(5,6) * t380 - Icges(6,6) * t373;
t417 = Icges(6,4) * t372;
t352 = -Icges(6,2) * t373 - t417;
t416 = Icges(6,4) * t373;
t353 = -Icges(6,1) * t372 - t416;
t419 = Icges(5,4) * t378;
t360 = -Icges(5,2) * t380 - t419;
t418 = Icges(5,4) * t380;
t361 = -Icges(5,1) * t378 - t418;
t433 = t352 * t372 - t353 * t373 + t360 * t378 - t361 * t380;
t395 = Icges(6,2) * t372 - t416;
t317 = Icges(6,6) * t357 + t356 * t395;
t397 = -Icges(6,1) * t373 + t417;
t319 = Icges(6,5) * t357 + t356 * t397;
t396 = Icges(5,2) * t378 - t418;
t327 = Icges(5,6) * t357 + t356 * t396;
t398 = -Icges(5,1) * t380 + t419;
t329 = Icges(5,5) * t357 + t356 * t398;
t432 = t317 * t372 - t319 * t373 + t327 * t378 - t329 * t380;
t316 = -Icges(6,6) * t356 + t357 * t395;
t318 = -Icges(6,5) * t356 + t357 * t397;
t326 = -Icges(5,6) * t356 + t357 * t396;
t328 = -Icges(5,5) * t356 + t357 * t398;
t431 = -t316 * t372 + t318 * t373 - t326 * t378 + t328 * t380;
t425 = pkin(4) * t378;
t423 = pkin(4) * t380;
t415 = t356 * t372;
t414 = t357 * t372;
t377 = sin(qJ(6));
t413 = t373 * t377;
t379 = cos(qJ(6));
t412 = t373 * t379;
t374 = qJD(2) * t426;
t411 = qJD(5) * t357 + t374;
t410 = qJD(4) * t356;
t409 = qJD(4) * t357;
t408 = qJD(4) * (-rSges(5,1) * t378 - rSges(5,2) * t380);
t407 = qJD(6) * t372;
t312 = -qJ(5) * t356 - t357 * t423;
t406 = t312 * t409 - qJD(3);
t363 = pkin(1) * t426 - qJ(2) * t427;
t403 = -pkin(2) * t426 - t363;
t402 = -pkin(5) * t373 - pkin(8) * t372;
t401 = -qJD(2) * t427 + qJD(1) * (pkin(1) * t427 + qJ(2) * t426);
t400 = -rSges(5,1) * t380 + rSges(5,2) * t378;
t399 = -rSges(6,1) * t373 + rSges(6,2) * t372;
t386 = pkin(3) * t357 + pkin(7) * t356 + t403;
t385 = qJD(1) * t427 * pkin(2) + t401;
t384 = -t312 + t386;
t383 = qJD(1) * (-pkin(3) * t356 + pkin(7) * t357) + t385;
t313 = qJ(5) * t357 - t356 * t423;
t382 = qJD(1) * t313 - qJD(5) * t356 + t409 * t425 + t383;
t366 = qJD(6) * t373 + qJD(1);
t365 = rSges(2,1) * t427 - rSges(2,2) * t426;
t364 = rSges(2,1) * t426 + rSges(2,2) * t427;
t355 = -pkin(5) * t372 + pkin(8) * t373;
t354 = -rSges(6,1) * t372 - rSges(6,2) * t373;
t347 = t373 * rSges(7,3) + (-rSges(7,1) * t379 + rSges(7,2) * t377) * t372;
t346 = Icges(7,5) * t373 + (-Icges(7,1) * t379 + Icges(7,4) * t377) * t372;
t345 = Icges(7,6) * t373 + (-Icges(7,4) * t379 + Icges(7,2) * t377) * t372;
t344 = Icges(7,3) * t373 + (-Icges(7,5) * t379 + Icges(7,6) * t377) * t372;
t343 = qJD(1) * (rSges(3,1) * t427 + rSges(3,3) * t426) + t401;
t342 = t374 + (-rSges(3,1) * t426 + rSges(3,3) * t427 - t363) * qJD(1);
t339 = -t357 * t407 - t410;
t338 = -t356 * t407 + t409;
t337 = -t356 * t412 + t357 * t377;
t336 = t356 * t413 + t357 * t379;
t335 = -t356 * t377 - t357 * t412;
t334 = -t356 * t379 + t357 * t413;
t333 = t402 * t356;
t332 = t402 * t357;
t331 = t357 * rSges(5,3) + t356 * t400;
t330 = -t356 * rSges(5,3) + t357 * t400;
t323 = qJD(1) * (-rSges(4,1) * t356 - rSges(4,2) * t357) + t385;
t322 = t374 + (t357 * rSges(4,1) - t356 * rSges(4,2) + t403) * qJD(1);
t321 = t357 * rSges(6,3) + t356 * t399;
t320 = -t356 * rSges(6,3) + t357 * t399;
t309 = rSges(7,1) * t337 + rSges(7,2) * t336 - rSges(7,3) * t415;
t308 = rSges(7,1) * t335 + rSges(7,2) * t334 - rSges(7,3) * t414;
t307 = Icges(7,1) * t337 + Icges(7,4) * t336 - Icges(7,5) * t415;
t306 = Icges(7,1) * t335 + Icges(7,4) * t334 - Icges(7,5) * t414;
t305 = Icges(7,4) * t337 + Icges(7,2) * t336 - Icges(7,6) * t415;
t304 = Icges(7,4) * t335 + Icges(7,2) * t334 - Icges(7,6) * t414;
t303 = Icges(7,5) * t337 + Icges(7,6) * t336 - Icges(7,3) * t415;
t302 = Icges(7,5) * t335 + Icges(7,6) * t334 - Icges(7,3) * t414;
t301 = qJD(1) * t331 - t357 * t408 + t383;
t300 = -t356 * t408 + t374 + (-t330 + t386) * qJD(1);
t299 = -qJD(3) + (t330 * t357 + t331 * t356) * qJD(4);
t298 = qJD(1) * t321 - t354 * t409 + t382;
t297 = (-t354 + t425) * t410 + (-t320 + t384) * qJD(1) + t411;
t296 = (t320 * t357 + (t313 + t321) * t356) * qJD(4) + t406;
t295 = qJD(1) * t333 + t366 * t309 - t338 * t347 - t355 * t409 + t382;
t294 = -t366 * t308 + t339 * t347 + (-t355 + t425) * t410 + (-t332 + t384) * qJD(1) + t411;
t293 = t308 * t338 - t309 * t339 + (t332 * t357 + (t313 + t333) * t356) * qJD(4) + t406;
t1 = t338 * ((-t303 * t415 + t336 * t305 + t337 * t307) * t338 + (-t302 * t415 + t304 * t336 + t306 * t337) * t339 + (t336 * t345 + t337 * t346 - t344 * t415) * t366) / 0.2e1 + t339 * ((-t303 * t414 + t305 * t334 + t307 * t335) * t338 + (-t302 * t414 + t334 * t304 + t335 * t306) * t339 + (t334 * t345 + t335 * t346 - t344 * t414) * t366) / 0.2e1 + t366 * ((t302 * t339 + t303 * t338 + t344 * t366) * t373 + ((t305 * t377 - t307 * t379) * t338 + (t304 * t377 - t306 * t379) * t339 + (t345 * t377 - t346 * t379) * t366) * t372) / 0.2e1 + m(7) * (t293 ^ 2 + t294 ^ 2 + t295 ^ 2) / 0.2e1 + m(6) * (t296 ^ 2 + t297 ^ 2 + t298 ^ 2) / 0.2e1 + m(5) * (t299 ^ 2 + t300 ^ 2 + t301 ^ 2) / 0.2e1 + m(3) * (t342 ^ 2 + t343 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + (((-t317 * t373 - t319 * t372 - t327 * t380 - t329 * t378) * t357 + (t316 * t373 + t318 * t372 + t326 * t380 + t328 * t378) * t356) * qJD(4) + (-t373 * t352 - t372 * t353 - t380 * t360 - t378 * t361) * qJD(1)) * qJD(1) / 0.2e1 - ((t435 * t356 ^ 2 + (t432 * t357 + (t431 - t436) * t356) * t357) * qJD(4) + (-t434 * t356 + t433 * t357) * qJD(1)) * t410 / 0.2e1 + ((t436 * t357 ^ 2 + (t431 * t356 + (t432 - t435) * t357) * t356) * qJD(4) + (t433 * t356 + t434 * t357) * qJD(1)) * t409 / 0.2e1 + (Icges(2,3) + Icges(4,3) + Icges(3,2) + m(2) * (t364 ^ 2 + t365 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
