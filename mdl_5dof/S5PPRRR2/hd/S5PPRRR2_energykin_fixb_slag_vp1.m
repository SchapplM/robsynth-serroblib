% Calculate kinetic energy for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:09
% EndTime: 2019-12-05 15:14:10
% DurationCPUTime: 1.35s
% Computational Cost: add. (1002->168), mult. (1119->281), div. (0->0), fcn. (1138->8), ass. (0->93)
t343 = cos(pkin(8));
t378 = t343 ^ 2;
t342 = sin(pkin(8));
t379 = t342 ^ 2;
t380 = t378 + t379;
t382 = qJD(3) * t380;
t381 = t342 * t343;
t345 = cos(qJ(4));
t376 = pkin(4) * t345;
t340 = pkin(9) + qJ(3);
t336 = sin(t340);
t374 = t336 * t342;
t373 = t336 * t343;
t341 = qJ(4) + qJ(5);
t338 = sin(t341);
t372 = t338 * t342;
t371 = t338 * t343;
t339 = cos(t341);
t370 = t339 * t342;
t369 = t339 * t343;
t344 = sin(qJ(4));
t368 = t342 * t344;
t367 = t342 * t345;
t366 = t343 * t344;
t365 = t343 * t345;
t334 = qJD(3) * t342;
t362 = qJD(4) * t336;
t326 = t343 * t362 + t334;
t364 = qJD(2) * t343;
t363 = qJD(3) * t343;
t337 = cos(t340);
t361 = qJD(4) * t337;
t360 = qJD(5) * t336;
t359 = qJD(1) + (pkin(3) * t337 + pkin(6) * t336) * t382;
t327 = t342 * t362 - t363;
t330 = pkin(3) * t336 - pkin(6) * t337;
t335 = qJD(2) * t342;
t358 = -t330 * t363 + t335;
t355 = Icges(4,5) * t337 - Icges(4,6) * t336;
t351 = -t330 * t334 - t364;
t349 = pkin(7) * t336 + t337 * t376;
t348 = qJD(1) ^ 2;
t329 = rSges(4,1) * t336 + rSges(4,2) * t337;
t328 = (-qJD(4) - qJD(5)) * t337;
t325 = t337 * t365 + t368;
t324 = -t337 * t366 + t367;
t323 = t337 * t367 - t366;
t322 = -t337 * t368 - t365;
t321 = t337 * t369 + t372;
t320 = -t337 * t371 + t370;
t319 = t337 * t370 - t371;
t318 = -t337 * t372 - t369;
t315 = -t329 * t334 - t364;
t314 = -t329 * t363 + t335;
t309 = Icges(4,3) * t342 + t343 * t355;
t308 = -Icges(4,3) * t343 + t342 * t355;
t307 = t342 * t360 + t327;
t306 = t343 * t360 + t326;
t305 = -rSges(5,3) * t337 + (rSges(5,1) * t345 - rSges(5,2) * t344) * t336;
t304 = -Icges(5,5) * t337 + (Icges(5,1) * t345 - Icges(5,4) * t344) * t336;
t303 = -Icges(5,6) * t337 + (Icges(5,4) * t345 - Icges(5,2) * t344) * t336;
t302 = -Icges(5,3) * t337 + (Icges(5,5) * t345 - Icges(5,6) * t344) * t336;
t301 = -rSges(6,3) * t337 + (rSges(6,1) * t339 - rSges(6,2) * t338) * t336;
t300 = -Icges(6,5) * t337 + (Icges(6,1) * t339 - Icges(6,4) * t338) * t336;
t299 = -Icges(6,6) * t337 + (Icges(6,4) * t339 - Icges(6,2) * t338) * t336;
t298 = -Icges(6,3) * t337 + (Icges(6,5) * t339 - Icges(6,6) * t338) * t336;
t297 = -pkin(7) * t337 + t336 * t376;
t296 = rSges(5,1) * t325 + rSges(5,2) * t324 + rSges(5,3) * t373;
t295 = rSges(5,1) * t323 + rSges(5,2) * t322 + rSges(5,3) * t374;
t294 = Icges(5,1) * t325 + Icges(5,4) * t324 + Icges(5,5) * t373;
t293 = Icges(5,1) * t323 + Icges(5,4) * t322 + Icges(5,5) * t374;
t292 = Icges(5,4) * t325 + Icges(5,2) * t324 + Icges(5,6) * t373;
t291 = Icges(5,4) * t323 + Icges(5,2) * t322 + Icges(5,6) * t374;
t290 = Icges(5,5) * t325 + Icges(5,6) * t324 + Icges(5,3) * t373;
t289 = Icges(5,5) * t323 + Icges(5,6) * t322 + Icges(5,3) * t374;
t288 = pkin(4) * t368 + t343 * t349;
t287 = -pkin(4) * t366 + t342 * t349;
t286 = rSges(6,1) * t321 + rSges(6,2) * t320 + rSges(6,3) * t373;
t285 = rSges(6,1) * t319 + rSges(6,2) * t318 + rSges(6,3) * t374;
t284 = Icges(6,1) * t321 + Icges(6,4) * t320 + Icges(6,5) * t373;
t283 = Icges(6,1) * t319 + Icges(6,4) * t318 + Icges(6,5) * t374;
t282 = Icges(6,4) * t321 + Icges(6,2) * t320 + Icges(6,6) * t373;
t281 = Icges(6,4) * t319 + Icges(6,2) * t318 + Icges(6,6) * t374;
t280 = Icges(6,5) * t321 + Icges(6,6) * t320 + Icges(6,3) * t373;
t279 = Icges(6,5) * t319 + Icges(6,6) * t318 + Icges(6,3) * t374;
t278 = qJD(1) + (rSges(4,1) * t337 - rSges(4,2) * t336) * t382;
t277 = -t296 * t361 - t305 * t326 + t351;
t276 = t295 * t361 + t305 * t327 + t358;
t275 = t295 * t326 - t296 * t327 + t359;
t274 = t286 * t328 - t288 * t361 - t297 * t326 - t301 * t306 + t351;
t273 = -t285 * t328 + t287 * t361 + t297 * t327 + t301 * t307 + t358;
t272 = t285 * t306 - t286 * t307 + t287 * t326 - t288 * t327 + t359;
t1 = m(2) * t348 / 0.2e1 + m(3) * (t380 * qJD(2) ^ 2 + t348) / 0.2e1 + m(4) * (t278 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + m(5) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + t326 * ((t290 * t373 + t292 * t324 + t294 * t325) * t326 + (t289 * t373 + t291 * t324 + t293 * t325) * t327 - (t302 * t373 + t303 * t324 + t304 * t325) * t361) / 0.2e1 + t327 * ((t290 * t374 + t292 * t322 + t294 * t323) * t326 + (t289 * t374 + t291 * t322 + t293 * t323) * t327 - (t302 * t374 + t303 * t322 + t304 * t323) * t361) / 0.2e1 - ((-t289 * t327 - t290 * t326 + t302 * t361) * t337 + ((-t292 * t344 + t294 * t345) * t326 + (-t291 * t344 + t293 * t345) * t327 - (-t303 * t344 + t304 * t345) * t361) * t336) * t361 / 0.2e1 + m(6) * (t272 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + t306 * ((t280 * t373 + t320 * t282 + t321 * t284) * t306 + (t279 * t373 + t281 * t320 + t283 * t321) * t307 + (t298 * t373 + t299 * t320 + t300 * t321) * t328) / 0.2e1 + t307 * ((t280 * t374 + t282 * t318 + t284 * t319) * t306 + (t279 * t374 + t318 * t281 + t319 * t283) * t307 + (t298 * t374 + t299 * t318 + t300 * t319) * t328) / 0.2e1 + t328 * ((-t279 * t307 - t280 * t306 - t298 * t328) * t337 + ((-t282 * t338 + t284 * t339) * t306 + (-t281 * t338 + t283 * t339) * t307 + (-t299 * t338 + t300 * t339) * t328) * t336) / 0.2e1 + (t342 * (-t308 * t381 + t379 * t309) / 0.2e1 - t343 * (t378 * t308 - t309 * t381) / 0.2e1) * qJD(3) ^ 2;
T = t1;
