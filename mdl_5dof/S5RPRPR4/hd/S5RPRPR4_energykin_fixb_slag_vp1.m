% Calculate kinetic energy for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:12
% EndTime: 2020-01-03 11:38:14
% DurationCPUTime: 1.48s
% Computational Cost: add. (985->168), mult. (765->256), div. (0->0), fcn. (634->10), ass. (0->102)
t377 = -Icges(4,3) - Icges(5,3);
t312 = qJ(3) + pkin(9);
t305 = sin(t312);
t307 = cos(t312);
t315 = sin(qJ(3));
t317 = cos(qJ(3));
t376 = Icges(4,5) * t317 + Icges(5,5) * t307 - Icges(4,6) * t315 - Icges(5,6) * t305;
t313 = qJ(1) + pkin(8);
t306 = sin(t313);
t308 = cos(t313);
t375 = t376 * t306 + t377 * t308;
t374 = t377 * t306 - t376 * t308;
t373 = -Icges(4,5) * t315 - Icges(5,5) * t305 - Icges(4,6) * t317 - Icges(5,6) * t307;
t357 = Icges(5,4) * t305;
t288 = Icges(5,2) * t307 + t357;
t356 = Icges(5,4) * t307;
t289 = Icges(5,1) * t305 + t356;
t359 = Icges(4,4) * t315;
t295 = Icges(4,2) * t317 + t359;
t358 = Icges(4,4) * t317;
t296 = Icges(4,1) * t315 + t358;
t372 = t288 * t305 - t289 * t307 + t295 * t315 - t296 * t317;
t335 = -Icges(5,2) * t305 + t356;
t265 = -Icges(5,6) * t306 - t308 * t335;
t338 = Icges(5,1) * t307 - t357;
t267 = -Icges(5,5) * t306 - t308 * t338;
t336 = -Icges(4,2) * t315 + t358;
t273 = -Icges(4,6) * t306 - t308 * t336;
t339 = Icges(4,1) * t317 - t359;
t275 = -Icges(4,5) * t306 - t308 * t339;
t371 = t265 * t305 - t267 * t307 + t273 * t315 - t275 * t317;
t264 = -Icges(5,6) * t308 + t306 * t335;
t266 = -Icges(5,5) * t308 + t306 * t338;
t272 = -Icges(4,6) * t308 + t306 * t336;
t274 = -Icges(4,5) * t308 + t306 * t339;
t370 = -t264 * t305 + t266 * t307 - t272 * t315 + t274 * t317;
t309 = qJ(5) + t312;
t300 = sin(t309);
t301 = cos(t309);
t354 = Icges(6,4) * t301;
t334 = -Icges(6,2) * t300 + t354;
t256 = -Icges(6,6) * t308 + t306 * t334;
t257 = -Icges(6,6) * t306 - t308 * t334;
t355 = Icges(6,4) * t300;
t337 = Icges(6,1) * t301 - t355;
t258 = -Icges(6,5) * t308 + t306 * t337;
t259 = -Icges(6,5) * t306 - t308 * t337;
t281 = Icges(6,2) * t301 + t355;
t282 = Icges(6,1) * t300 + t354;
t346 = -qJD(3) - qJD(5);
t285 = t346 * t306;
t286 = t346 * t308;
t369 = (t281 * t300 - t282 * t301) * qJD(1) + (t256 * t300 - t258 * t301) * t286 + (t257 * t300 - t259 * t301) * t285;
t365 = pkin(3) * t315;
t364 = pkin(4) * t305;
t362 = t317 * pkin(3);
t360 = pkin(1) * qJD(1);
t253 = -qJ(4) * t306 - t308 * t362;
t291 = -pkin(2) * t308 - pkin(6) * t306;
t353 = -t253 - t291;
t316 = sin(qJ(1));
t303 = t316 * t360;
t352 = qJD(1) * (pkin(2) * t306 - pkin(6) * t308) + t303;
t351 = pkin(4) * t307;
t349 = qJD(3) * t306;
t348 = qJD(3) * t308;
t252 = -qJ(4) * t308 + t306 * t362;
t347 = t252 * t349 + qJD(2);
t318 = cos(qJ(1));
t304 = t318 * t360;
t343 = -qJD(4) * t308 + t304;
t342 = rSges(4,1) * t317 - rSges(4,2) * t315;
t341 = rSges(5,1) * t307 - rSges(5,2) * t305;
t340 = rSges(6,1) * t301 - rSges(6,2) * t300;
t331 = Icges(6,5) * t301 - Icges(6,6) * t300;
t321 = qJD(1) * t252 - qJD(4) * t306 + t348 * t365 + t352;
t320 = -qJD(1) * (Icges(6,5) * t300 + Icges(6,6) * t301) - (-Icges(6,3) * t308 + t306 * t331) * t286 - (-Icges(6,3) * t306 - t308 * t331) * t285;
t299 = -rSges(2,1) * t318 + rSges(2,2) * t316;
t298 = rSges(2,1) * t316 + rSges(2,2) * t318;
t297 = rSges(4,1) * t315 + rSges(4,2) * t317;
t290 = rSges(5,1) * t305 + rSges(5,2) * t307;
t283 = rSges(6,1) * t300 + rSges(6,2) * t301;
t279 = t304 - qJD(1) * (-rSges(3,1) * t308 + rSges(3,2) * t306);
t278 = t303 + qJD(1) * (rSges(3,1) * t306 + rSges(3,2) * t308);
t277 = -rSges(4,3) * t306 - t308 * t342;
t276 = -rSges(4,3) * t308 + t306 * t342;
t269 = -rSges(5,3) * t306 - t308 * t341;
t268 = -rSges(5,3) * t308 + t306 * t341;
t261 = -rSges(6,3) * t306 - t308 * t340;
t260 = -rSges(6,3) * t308 + t306 * t340;
t249 = -pkin(7) * t306 - t308 * t351;
t248 = -pkin(7) * t308 + t306 * t351;
t247 = -t297 * t349 + t304 + (-t277 - t291) * qJD(1);
t246 = qJD(1) * t276 + t297 * t348 + t352;
t245 = qJD(2) + (t276 * t306 - t277 * t308) * qJD(3);
t244 = (-t290 - t365) * t349 + (-t269 + t353) * qJD(1) + t343;
t243 = qJD(1) * t268 + t290 * t348 + t321;
t242 = (t268 * t306 + (-t253 - t269) * t308) * qJD(3) + t347;
t241 = t283 * t285 + (-t364 - t365) * t349 + (-t249 - t261 + t353) * qJD(1) + t343;
t240 = t348 * t364 - t283 * t286 + (t248 + t260) * qJD(1) + t321;
t239 = -t260 * t285 + t261 * t286 + (t248 * t306 + (-t249 - t253) * t308) * qJD(3) + t347;
t1 = m(3) * (qJD(2) ^ 2 + t278 ^ 2 + t279 ^ 2) / 0.2e1 + m(4) * (t245 ^ 2 + t246 ^ 2 + t247 ^ 2) / 0.2e1 + m(5) * (t242 ^ 2 + t243 ^ 2 + t244 ^ 2) / 0.2e1 + m(6) * (t239 ^ 2 + t240 ^ 2 + t241 ^ 2) / 0.2e1 + t286 * (-t369 * t306 + t320 * t308) / 0.2e1 + t285 * (t320 * t306 + t369 * t308) / 0.2e1 - ((t374 * t306 ^ 2 + (t370 * t308 + (-t371 + t375) * t306) * t308) * qJD(3) + (t373 * t306 + t372 * t308) * qJD(1)) * t349 / 0.2e1 - ((t375 * t308 ^ 2 + (t371 * t306 + (-t370 + t374) * t308) * t306) * qJD(3) + (-t372 * t306 + t373 * t308) * qJD(1)) * t348 / 0.2e1 + (m(2) * (t298 ^ 2 + t299 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1 + ((t256 * t301 + t258 * t300) * t286 + (t257 * t301 + t259 * t300) * t285 + ((-t264 * t307 - t305 * t266 - t272 * t317 - t274 * t315) * t308 + (-t265 * t307 - t267 * t305 - t273 * t317 - t275 * t315) * t306) * qJD(3) + (t301 * t281 + t300 * t282 + t307 * t288 + t305 * t289 + t317 * t295 + t315 * t296) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
