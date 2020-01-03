% Calculate kinetic energy for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR12_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR12_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:00
% EndTime: 2019-12-31 20:29:03
% DurationCPUTime: 2.55s
% Computational Cost: add. (772->224), mult. (1874->354), div. (0->0), fcn. (2005->8), ass. (0->118)
t442 = Icges(3,4) - Icges(4,5);
t441 = Icges(3,1) + Icges(4,1);
t440 = Icges(3,2) + Icges(4,3);
t374 = sin(qJ(2));
t439 = t442 * t374;
t377 = cos(qJ(2));
t438 = t442 * t377;
t437 = Icges(4,4) + Icges(3,5);
t436 = Icges(3,6) - Icges(4,6);
t435 = t440 * t374 - t438;
t434 = t441 * t377 - t439;
t433 = Icges(4,2) + Icges(3,3);
t375 = sin(qJ(1));
t378 = cos(qJ(1));
t432 = t435 * t375 + t436 * t378;
t431 = -t436 * t375 + t435 * t378;
t430 = -t434 * t375 + t437 * t378;
t429 = t437 * t375 + t434 * t378;
t428 = -t440 * t377 - t439;
t427 = t441 * t374 + t438;
t426 = -t436 * t374 + t437 * t377;
t425 = t426 * t375 - t433 * t378;
t424 = t433 * t375 + t426 * t378;
t423 = t437 * t374 + t436 * t377;
t422 = t428 * t374 + t427 * t377;
t421 = t431 * t374 + t429 * t377;
t420 = -t432 * t374 + t430 * t377;
t373 = sin(qJ(4));
t416 = cos(qJ(4));
t403 = t374 * t416;
t348 = -t377 * t373 + t403;
t409 = t377 * t378;
t396 = pkin(2) * t377 + qJ(3) * t374;
t344 = t396 * t375;
t365 = pkin(1) * t375 - pkin(6) * t378;
t408 = -t344 - t365;
t371 = qJD(2) * t375;
t407 = qJD(2) * t378;
t406 = qJD(3) * t374;
t345 = t396 * t378;
t351 = qJD(1) * (pkin(1) * t378 + pkin(6) * t375);
t405 = qJD(1) * t345 + t375 * t406 + t351;
t349 = pkin(3) * t375 * t377 + pkin(7) * t378;
t404 = -t349 + t408;
t360 = pkin(2) * t374 - qJ(3) * t377;
t400 = qJD(2) * (-rSges(4,1) * t374 + rSges(4,3) * t377 - t360);
t353 = qJD(4) * t378 - t407;
t352 = -qJD(4) * t375 + t371;
t399 = -qJD(3) * t377 + t344 * t371 + t345 * t407;
t398 = rSges(3,1) * t377 - rSges(3,2) * t374;
t397 = rSges(4,1) * t377 + rSges(4,3) * t374;
t395 = qJD(2) * (-pkin(3) * t374 - t360);
t347 = t374 * t373 + t377 * t416;
t350 = pkin(3) * t409 - pkin(7) * t375;
t382 = t349 * t371 + t350 * t407 + t399;
t369 = t378 * t406;
t381 = t378 * t395 + t369;
t380 = qJD(1) * t350 + t375 * t395 + t405;
t376 = cos(qJ(5));
t372 = sin(qJ(5));
t364 = rSges(2,1) * t378 - rSges(2,2) * t375;
t363 = rSges(2,1) * t375 + rSges(2,2) * t378;
t362 = rSges(3,1) * t374 + rSges(3,2) * t377;
t340 = qJD(5) * t347 + qJD(1);
t339 = t347 * t378;
t338 = t373 * t409 - t378 * t403;
t337 = t347 * t375;
t336 = t348 * t375;
t335 = rSges(3,3) * t375 + t378 * t398;
t334 = rSges(4,2) * t375 + t378 * t397;
t333 = -rSges(3,3) * t378 + t375 * t398;
t332 = -rSges(4,2) * t378 + t375 * t397;
t317 = t339 * t376 - t372 * t375;
t316 = -t339 * t372 - t375 * t376;
t315 = t337 * t376 + t372 * t378;
t314 = -t337 * t372 + t376 * t378;
t313 = pkin(4) * t348 + pkin(8) * t347;
t312 = rSges(5,1) * t348 - rSges(5,2) * t347;
t311 = Icges(5,1) * t348 - Icges(5,4) * t347;
t310 = Icges(5,4) * t348 - Icges(5,2) * t347;
t309 = Icges(5,5) * t348 - Icges(5,6) * t347;
t308 = -qJD(5) * t336 + t353;
t307 = qJD(5) * t338 + t352;
t306 = pkin(4) * t339 + pkin(8) * t338;
t305 = pkin(4) * t337 - pkin(8) * t336;
t304 = rSges(5,1) * t339 - rSges(5,2) * t338 - rSges(5,3) * t375;
t303 = rSges(5,1) * t337 + rSges(5,2) * t336 + rSges(5,3) * t378;
t302 = Icges(5,1) * t339 - Icges(5,4) * t338 - Icges(5,5) * t375;
t301 = Icges(5,1) * t337 + Icges(5,4) * t336 + Icges(5,5) * t378;
t300 = Icges(5,4) * t339 - Icges(5,2) * t338 - Icges(5,6) * t375;
t299 = Icges(5,4) * t337 + Icges(5,2) * t336 + Icges(5,6) * t378;
t298 = Icges(5,5) * t339 - Icges(5,6) * t338 - Icges(5,3) * t375;
t297 = Icges(5,5) * t337 + Icges(5,6) * t336 + Icges(5,3) * t378;
t296 = qJD(1) * t335 - t362 * t371 + t351;
t295 = -t362 * t407 + (-t333 - t365) * qJD(1);
t294 = (t333 * t375 + t335 * t378) * qJD(2);
t293 = rSges(6,3) * t347 + (rSges(6,1) * t376 - rSges(6,2) * t372) * t348;
t292 = Icges(6,5) * t347 + (Icges(6,1) * t376 - Icges(6,4) * t372) * t348;
t291 = Icges(6,6) * t347 + (Icges(6,4) * t376 - Icges(6,2) * t372) * t348;
t290 = Icges(6,3) * t347 + (Icges(6,5) * t376 - Icges(6,6) * t372) * t348;
t289 = rSges(6,1) * t317 + rSges(6,2) * t316 + rSges(6,3) * t338;
t288 = rSges(6,1) * t315 + rSges(6,2) * t314 - rSges(6,3) * t336;
t287 = Icges(6,1) * t317 + Icges(6,4) * t316 + Icges(6,5) * t338;
t286 = Icges(6,1) * t315 + Icges(6,4) * t314 - Icges(6,5) * t336;
t285 = Icges(6,4) * t317 + Icges(6,2) * t316 + Icges(6,6) * t338;
t284 = Icges(6,4) * t315 + Icges(6,2) * t314 - Icges(6,6) * t336;
t283 = Icges(6,5) * t317 + Icges(6,6) * t316 + Icges(6,3) * t338;
t282 = Icges(6,5) * t315 + Icges(6,6) * t314 - Icges(6,3) * t336;
t281 = qJD(1) * t334 + t375 * t400 + t405;
t280 = t369 + t378 * t400 + (-t332 + t408) * qJD(1);
t279 = (t332 * t375 + t334 * t378) * qJD(2) + t399;
t278 = qJD(1) * t304 - t312 * t352 + t380;
t277 = t312 * t353 + (-t303 + t404) * qJD(1) + t381;
t276 = t303 * t352 - t304 * t353 + t382;
t275 = qJD(1) * t306 + t289 * t340 - t293 * t307 - t313 * t352 + t380;
t274 = -t288 * t340 + t293 * t308 + t313 * t353 + (-t305 + t404) * qJD(1) + t381;
t273 = t288 * t307 - t289 * t308 + t305 * t352 - t306 * t353 + t382;
t1 = m(3) * (t294 ^ 2 + t295 ^ 2 + t296 ^ 2) / 0.2e1 + m(4) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(5) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + t352 * ((-t298 * t375 - t300 * t338 + t302 * t339) * t352 + (-t297 * t375 - t299 * t338 + t301 * t339) * t353 + (-t309 * t375 - t310 * t338 + t311 * t339) * qJD(1)) / 0.2e1 + t353 * ((t298 * t378 + t300 * t336 + t302 * t337) * t352 + (t297 * t378 + t299 * t336 + t301 * t337) * t353 + (t309 * t378 + t310 * t336 + t311 * t337) * qJD(1)) / 0.2e1 + m(6) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + t307 * ((t283 * t338 + t285 * t316 + t287 * t317) * t307 + (t282 * t338 + t284 * t316 + t286 * t317) * t308 + (t290 * t338 + t291 * t316 + t292 * t317) * t340) / 0.2e1 + t308 * ((-t283 * t336 + t285 * t314 + t287 * t315) * t307 + (-t282 * t336 + t284 * t314 + t286 * t315) * t308 + (-t290 * t336 + t291 * t314 + t292 * t315) * t340) / 0.2e1 + t340 * ((t282 * t308 + t283 * t307 + t290 * t340) * t347 + ((-t285 * t372 + t287 * t376) * t307 + (-t284 * t372 + t286 * t376) * t308 + (-t291 * t372 + t292 * t376) * t340) * t348) / 0.2e1 + (m(2) * (t363 ^ 2 + t364 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t424 * t375 ^ 2 + (t420 * t378 + (t421 - t425) * t375) * t378) * qJD(2) + (t423 * t375 + t422 * t378) * qJD(1)) * t371 / 0.2e1 - ((t425 * t378 ^ 2 + (t421 * t375 + (t420 - t424) * t378) * t375) * qJD(2) + (t422 * t375 - t423 * t378) * qJD(1)) * t407 / 0.2e1 + ((-t300 * t347 + t302 * t348) * t352 + (-t299 * t347 + t301 * t348) * t353 + ((t430 * t374 + t432 * t377) * t378 + (t429 * t374 - t431 * t377) * t375) * qJD(2) + (-t310 * t347 + t311 * t348 + t427 * t374 - t428 * t377) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
