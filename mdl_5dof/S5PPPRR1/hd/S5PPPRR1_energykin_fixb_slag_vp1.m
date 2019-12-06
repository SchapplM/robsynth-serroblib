% Calculate kinetic energy for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:48
% EndTime: 2019-12-05 14:57:48
% DurationCPUTime: 0.72s
% Computational Cost: add. (893->152), mult. (1260->260), div. (0->0), fcn. (1431->8), ass. (0->69)
t331 = pkin(9) + qJ(4);
t329 = sin(t331);
t332 = sin(pkin(8));
t349 = t329 * t332;
t333 = sin(pkin(7));
t348 = t332 * t333;
t335 = cos(pkin(7));
t347 = t332 * t335;
t336 = sin(qJ(5));
t346 = t332 * t336;
t337 = cos(qJ(5));
t345 = t332 * t337;
t334 = cos(pkin(8));
t344 = t333 * t334;
t343 = t334 * t335;
t342 = qJD(3) * t332;
t322 = qJD(2) * t333 + t335 * t342;
t341 = qJD(4) * t332;
t323 = -qJD(2) * t335 + t333 * t342;
t327 = -qJD(3) * t334 + qJD(1);
t339 = qJD(1) ^ 2;
t330 = cos(t331);
t321 = -qJD(4) * t334 + qJD(5) * t349;
t320 = t330 * t345 - t334 * t336;
t319 = -t330 * t346 - t334 * t337;
t318 = (pkin(4) * t330 + pkin(6) * t329) * t332;
t317 = t329 * t333 + t330 * t343;
t316 = t329 * t343 - t333 * t330;
t315 = -t329 * t335 + t330 * t344;
t314 = t329 * t344 + t330 * t335;
t313 = -rSges(5,3) * t334 + (rSges(5,1) * t330 - rSges(5,2) * t329) * t332;
t312 = -Icges(5,5) * t334 + (Icges(5,1) * t330 - Icges(5,4) * t329) * t332;
t311 = -Icges(5,6) * t334 + (Icges(5,4) * t330 - Icges(5,2) * t329) * t332;
t310 = -Icges(5,3) * t334 + (Icges(5,5) * t330 - Icges(5,6) * t329) * t332;
t309 = qJD(5) * t316 + t335 * t341;
t308 = qJD(5) * t314 + t333 * t341;
t307 = t317 * t337 + t335 * t346;
t306 = -t317 * t336 + t335 * t345;
t305 = t315 * t337 + t333 * t346;
t304 = -t315 * t336 + t333 * t345;
t303 = pkin(4) * t317 + pkin(6) * t316;
t302 = pkin(4) * t315 + pkin(6) * t314;
t301 = rSges(6,1) * t320 + rSges(6,2) * t319 + rSges(6,3) * t349;
t300 = Icges(6,1) * t320 + Icges(6,4) * t319 + Icges(6,5) * t349;
t299 = Icges(6,4) * t320 + Icges(6,2) * t319 + Icges(6,6) * t349;
t298 = Icges(6,5) * t320 + Icges(6,6) * t319 + Icges(6,3) * t349;
t297 = rSges(5,1) * t317 - rSges(5,2) * t316 + rSges(5,3) * t347;
t296 = rSges(5,1) * t315 - rSges(5,2) * t314 + rSges(5,3) * t348;
t295 = Icges(5,1) * t317 - Icges(5,4) * t316 + Icges(5,5) * t347;
t294 = Icges(5,1) * t315 - Icges(5,4) * t314 + Icges(5,5) * t348;
t293 = Icges(5,4) * t317 - Icges(5,2) * t316 + Icges(5,6) * t347;
t292 = Icges(5,4) * t315 - Icges(5,2) * t314 + Icges(5,6) * t348;
t291 = Icges(5,5) * t317 - Icges(5,6) * t316 + Icges(5,3) * t347;
t290 = Icges(5,5) * t315 - Icges(5,6) * t314 + Icges(5,3) * t348;
t289 = rSges(6,1) * t307 + rSges(6,2) * t306 + rSges(6,3) * t316;
t288 = rSges(6,1) * t305 + rSges(6,2) * t304 + rSges(6,3) * t314;
t287 = Icges(6,1) * t307 + Icges(6,4) * t306 + Icges(6,5) * t316;
t286 = Icges(6,1) * t305 + Icges(6,4) * t304 + Icges(6,5) * t314;
t285 = Icges(6,4) * t307 + Icges(6,2) * t306 + Icges(6,6) * t316;
t284 = Icges(6,4) * t305 + Icges(6,2) * t304 + Icges(6,6) * t314;
t283 = Icges(6,5) * t307 + Icges(6,6) * t306 + Icges(6,3) * t316;
t282 = Icges(6,5) * t305 + Icges(6,6) * t304 + Icges(6,3) * t314;
t281 = (-t297 * t334 - t313 * t347) * qJD(4) + t323;
t280 = (t296 * t334 + t313 * t348) * qJD(4) + t322;
t279 = (t296 * t335 - t297 * t333) * t341 + t327;
t278 = t289 * t321 - t301 * t309 + (-t303 * t334 - t318 * t347) * qJD(4) + t323;
t277 = -t288 * t321 + t301 * t308 + (t302 * t334 + t318 * t348) * qJD(4) + t322;
t276 = t288 * t309 - t289 * t308 + (t302 * t335 - t303 * t333) * t341 + t327;
t1 = m(2) * t339 / 0.2e1 + m(3) * (t339 + (t333 ^ 2 + t335 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t322 ^ 2 + t323 ^ 2 + t327 ^ 2) / 0.2e1 + m(5) * (t279 ^ 2 + t280 ^ 2 + t281 ^ 2) / 0.2e1 + m(6) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + t309 * ((t316 * t283 + t306 * t285 + t307 * t287) * t309 + (t282 * t316 + t284 * t306 + t286 * t307) * t308 + (t298 * t316 + t299 * t306 + t300 * t307) * t321) / 0.2e1 + t308 * ((t283 * t314 + t285 * t304 + t287 * t305) * t309 + (t314 * t282 + t304 * t284 + t305 * t286) * t308 + (t298 * t314 + t299 * t304 + t300 * t305) * t321) / 0.2e1 + t321 * ((t283 * t349 + t285 * t319 + t287 * t320) * t309 + (t282 * t349 + t284 * t319 + t286 * t320) * t308 + (t298 * t349 + t319 * t299 + t320 * t300) * t321) / 0.2e1 + (-t334 * (t334 ^ 2 * t310 + (((-t293 * t329 + t295 * t330) * t335 + (-t292 * t329 + t294 * t330) * t333) * t332 + (-t290 * t333 - t291 * t335 + t311 * t329 - t312 * t330) * t334) * t332) / 0.2e1 + (t335 * ((t291 * t347 - t293 * t316 + t295 * t317) * t347 + (t290 * t347 - t292 * t316 + t294 * t317) * t348 - (t310 * t347 - t311 * t316 + t312 * t317) * t334) + t333 * ((t291 * t348 - t293 * t314 + t295 * t315) * t347 + (t290 * t348 - t292 * t314 + t294 * t315) * t348 - (t310 * t348 - t311 * t314 + t312 * t315) * t334)) * t332 / 0.2e1) * qJD(4) ^ 2;
T = t1;
