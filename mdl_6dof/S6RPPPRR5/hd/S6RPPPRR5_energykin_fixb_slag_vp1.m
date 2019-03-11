% Calculate kinetic energy for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:12
% EndTime: 2019-03-09 01:37:13
% DurationCPUTime: 0.67s
% Computational Cost: add. (625->167), mult. (1341->275), div. (0->0), fcn. (1550->8), ass. (0->87)
t354 = sin(qJ(1));
t353 = sin(pkin(9));
t320 = sin(qJ(5));
t352 = Icges(6,4) * t320;
t322 = cos(qJ(5));
t351 = Icges(6,4) * t322;
t318 = cos(pkin(9));
t323 = cos(qJ(1));
t299 = t318 * t323 - t353 * t354;
t350 = t299 * t320;
t300 = t318 * t354 + t323 * t353;
t349 = t300 * t320;
t348 = t322 * t299;
t347 = t322 * t300;
t317 = qJD(2) * t354;
t346 = qJD(3) * t323 + t317;
t345 = qJD(1) * t323;
t344 = qJD(5) * t299;
t343 = qJD(5) * t300;
t342 = qJD(5) * (rSges(6,1) * t320 + rSges(6,2) * t322);
t341 = qJD(5) * (pkin(5) * t320 - pkin(8) * t322);
t340 = qJD(6) * t320;
t339 = pkin(3) * t345 + t346;
t338 = -qJD(2) * t323 + qJD(1) * (t323 * pkin(1) + qJ(2) * t354);
t337 = pkin(5) * t322 + pkin(8) * t320;
t306 = pkin(1) * t354 - t323 * qJ(2);
t336 = -qJ(3) * t354 - t306;
t335 = rSges(6,1) * t322 - rSges(6,2) * t320;
t334 = Icges(6,1) * t322 - t352;
t333 = -Icges(6,2) * t320 + t351;
t332 = Icges(6,5) * t322 - Icges(6,6) * t320;
t275 = -Icges(6,6) * t299 + t300 * t333;
t277 = -Icges(6,5) * t299 + t300 * t334;
t331 = -t275 * t320 + t277 * t322;
t276 = -Icges(6,6) * t300 - t299 * t333;
t278 = -Icges(6,5) * t300 - t299 * t334;
t330 = t276 * t320 - t278 * t322;
t303 = Icges(6,2) * t322 + t352;
t304 = Icges(6,1) * t320 + t351;
t329 = t303 * t320 - t304 * t322;
t328 = pkin(4) * t299 + pkin(7) * t300 + t336;
t327 = qJ(3) * t345 + qJD(3) * t354 + t338;
t326 = qJD(1) * t354 * pkin(3) + t327;
t325 = qJD(1) * (pkin(4) * t300 - pkin(7) * t299) + t326;
t321 = cos(qJ(6));
t319 = sin(qJ(6));
t310 = -qJD(6) * t322 + qJD(1);
t308 = t323 * rSges(2,1) - rSges(2,2) * t354;
t307 = rSges(2,1) * t354 + t323 * rSges(2,2);
t302 = Icges(6,5) * t320 + Icges(6,6) * t322;
t298 = -rSges(7,3) * t322 + (rSges(7,1) * t321 - rSges(7,2) * t319) * t320;
t297 = -Icges(7,5) * t322 + (Icges(7,1) * t321 - Icges(7,4) * t319) * t320;
t296 = -Icges(7,6) * t322 + (Icges(7,4) * t321 - Icges(7,2) * t319) * t320;
t295 = -Icges(7,3) * t322 + (Icges(7,5) * t321 - Icges(7,6) * t319) * t320;
t294 = qJD(1) * (-t323 * rSges(3,2) + rSges(3,3) * t354) + t338;
t293 = t317 + (rSges(3,2) * t354 + t323 * rSges(3,3) - t306) * qJD(1);
t290 = t300 * t340 - t344;
t289 = -t299 * t340 - t343;
t288 = t337 * t299;
t287 = t337 * t300;
t286 = -t300 * t319 - t321 * t348;
t285 = -t300 * t321 + t319 * t348;
t284 = -t299 * t319 + t321 * t347;
t283 = -t299 * t321 - t319 * t347;
t282 = qJD(1) * (rSges(4,1) * t354 + t323 * rSges(4,3)) + t327;
t281 = (t323 * rSges(4,1) - rSges(4,3) * t354 + t336) * qJD(1) + t346;
t280 = -rSges(6,3) * t300 - t299 * t335;
t279 = -rSges(6,3) * t299 + t300 * t335;
t274 = -Icges(6,3) * t300 - t299 * t332;
t273 = -Icges(6,3) * t299 + t300 * t332;
t272 = qJD(1) * (rSges(5,1) * t300 + rSges(5,2) * t299) + t326;
t271 = (t299 * rSges(5,1) - t300 * rSges(5,2) + t336) * qJD(1) + t339;
t270 = rSges(7,1) * t286 + rSges(7,2) * t285 - rSges(7,3) * t350;
t269 = rSges(7,1) * t284 + rSges(7,2) * t283 + rSges(7,3) * t349;
t268 = Icges(7,1) * t286 + Icges(7,4) * t285 - Icges(7,5) * t350;
t267 = Icges(7,1) * t284 + Icges(7,4) * t283 + Icges(7,5) * t349;
t266 = Icges(7,4) * t286 + Icges(7,2) * t285 - Icges(7,6) * t350;
t265 = Icges(7,4) * t284 + Icges(7,2) * t283 + Icges(7,6) * t349;
t264 = Icges(7,5) * t286 + Icges(7,6) * t285 - Icges(7,3) * t350;
t263 = Icges(7,5) * t284 + Icges(7,6) * t283 + Icges(7,3) * t349;
t262 = qJD(4) + (t279 * t300 - t280 * t299) * qJD(5);
t261 = qJD(1) * t279 + t299 * t342 + t325;
t260 = -t300 * t342 + (-t280 + t328) * qJD(1) + t339;
t259 = qJD(1) * t287 + t269 * t310 - t290 * t298 + t299 * t341 + t325;
t258 = -t300 * t341 - t310 * t270 + t289 * t298 + (t288 + t328) * qJD(1) + t339;
t257 = -t269 * t289 + t270 * t290 + qJD(4) + (t287 * t300 + t288 * t299) * qJD(5);
t1 = m(3) * (t293 ^ 2 + t294 ^ 2) / 0.2e1 + m(4) * (t281 ^ 2 + t282 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + m(6) * (t260 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + m(7) * (t257 ^ 2 + t258 ^ 2 + t259 ^ 2) / 0.2e1 - ((-t299 * t302 - t300 * t329) * qJD(1) + (t299 ^ 2 * t273 + (t330 * t300 + (t274 - t331) * t299) * t300) * qJD(5)) * t344 / 0.2e1 - ((t299 * t329 - t300 * t302) * qJD(1) + (t300 ^ 2 * t274 + (t331 * t299 + (t273 - t330) * t300) * t299) * qJD(5)) * t343 / 0.2e1 + qJD(1) * ((t322 * t303 + t320 * t304) * qJD(1) + (-(t275 * t322 + t320 * t277) * t299 - (t276 * t322 + t278 * t320) * t300) * qJD(5)) / 0.2e1 + t290 * ((t263 * t349 + t283 * t265 + t284 * t267) * t290 + (t264 * t349 + t266 * t283 + t268 * t284) * t289 + (t283 * t296 + t284 * t297 + t295 * t349) * t310) / 0.2e1 + t289 * ((-t263 * t350 + t265 * t285 + t267 * t286) * t290 + (-t264 * t350 + t285 * t266 + t286 * t268) * t289 + (t285 * t296 + t286 * t297 - t295 * t350) * t310) / 0.2e1 + t310 * ((-t263 * t290 - t264 * t289 - t295 * t310) * t322 + ((-t265 * t319 + t267 * t321) * t290 + (-t266 * t319 + t268 * t321) * t289 + (-t296 * t319 + t297 * t321) * t310) * t320) / 0.2e1 + (m(2) * (t307 ^ 2 + t308 ^ 2) + Icges(2,3) + Icges(3,1) + Icges(4,2) + Icges(5,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
