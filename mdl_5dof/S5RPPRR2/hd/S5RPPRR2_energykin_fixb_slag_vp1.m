% Calculate kinetic energy for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:32
% EndTime: 2019-12-05 17:39:33
% DurationCPUTime: 0.64s
% Computational Cost: add. (517->137), mult. (563->222), div. (0->0), fcn. (458->8), ass. (0->79)
t300 = sin(qJ(1));
t301 = cos(qJ(1));
t337 = qJD(1) * t301 * qJ(3) + qJD(3) * t300;
t296 = pkin(8) + qJ(4);
t291 = qJ(5) + t296;
t287 = cos(t291);
t286 = sin(t291);
t330 = Icges(6,4) * t286;
t312 = Icges(6,2) * t287 + t330;
t255 = Icges(6,6) * t301 + t312 * t300;
t256 = Icges(6,6) * t300 - t312 * t301;
t329 = Icges(6,4) * t287;
t314 = Icges(6,1) * t286 + t329;
t257 = Icges(6,5) * t301 + t314 * t300;
t258 = Icges(6,5) * t300 - t314 * t301;
t272 = -Icges(6,2) * t286 + t329;
t273 = Icges(6,1) * t287 - t330;
t324 = qJD(4) + qJD(5);
t281 = t324 * t300;
t282 = t324 * t301;
t336 = (t255 * t287 + t257 * t286) * t282 + (t256 * t287 + t258 * t286) * t281 + (t272 * t287 + t273 * t286) * qJD(1);
t297 = sin(pkin(8));
t334 = pkin(3) * t297;
t289 = sin(t296);
t332 = Icges(5,4) * t289;
t290 = cos(t296);
t331 = Icges(5,4) * t290;
t294 = qJD(2) * t300;
t327 = qJD(3) * t301 + t294;
t325 = qJD(4) * t300;
t323 = pkin(4) * qJD(4) * t290;
t322 = pkin(4) * t289;
t280 = qJD(1) * (pkin(1) * t301 + qJ(2) * t300);
t321 = -qJD(2) * t301 + t280;
t320 = qJD(1) * (pkin(6) * t301 + t300 * t334) + t280 + t337;
t283 = pkin(1) * t300 - qJ(2) * t301;
t319 = t301 * t334 - t283 + (-pkin(6) - qJ(3)) * t300;
t298 = cos(pkin(8));
t318 = rSges(4,1) * t297 + rSges(4,2) * t298;
t317 = rSges(5,1) * t289 + rSges(5,2) * t290;
t316 = rSges(6,1) * t286 + rSges(6,2) * t287;
t315 = Icges(5,1) * t289 + t331;
t313 = Icges(5,2) * t290 + t332;
t311 = Icges(5,5) * t289 + Icges(5,6) * t290;
t310 = Icges(6,5) * t286 + Icges(6,6) * t287;
t263 = Icges(5,6) * t301 + t313 * t300;
t265 = Icges(5,5) * t301 + t315 * t300;
t307 = -t263 * t290 - t265 * t289;
t264 = Icges(5,6) * t300 - t313 * t301;
t266 = Icges(5,5) * t300 - t315 * t301;
t306 = t264 * t290 + t266 * t289;
t276 = -Icges(5,2) * t289 + t331;
t277 = Icges(5,1) * t290 - t332;
t304 = t276 * t290 + t277 * t289;
t303 = (Icges(6,5) * t287 - Icges(6,6) * t286) * qJD(1) + (Icges(6,3) * t301 + t310 * t300) * t282 + (Icges(6,3) * t300 - t310 * t301) * t281;
t285 = rSges(2,1) * t301 - rSges(2,2) * t300;
t284 = rSges(2,1) * t300 + rSges(2,2) * t301;
t278 = rSges(5,1) * t290 - rSges(5,2) * t289;
t275 = Icges(5,5) * t290 - Icges(5,6) * t289;
t274 = rSges(6,1) * t287 - rSges(6,2) * t286;
t268 = rSges(5,3) * t300 - t317 * t301;
t267 = rSges(5,3) * t301 + t317 * t300;
t262 = Icges(5,3) * t300 - t311 * t301;
t261 = Icges(5,3) * t301 + t311 * t300;
t260 = rSges(6,3) * t300 - t316 * t301;
t259 = rSges(6,3) * t301 + t316 * t300;
t252 = qJD(1) * (-rSges(3,2) * t301 + rSges(3,3) * t300) + t321;
t251 = t294 + (rSges(3,2) * t300 + rSges(3,3) * t301 - t283) * qJD(1);
t250 = pkin(7) * t301 + t322 * t300;
t249 = pkin(7) * t300 - t322 * t301;
t248 = qJD(1) * (rSges(4,3) * t301 + t318 * t300) + t321 + t337;
t247 = (-t283 + t318 * t301 + (-rSges(4,3) - qJ(3)) * t300) * qJD(1) + t327;
t246 = (-t267 * t300 + t268 * t301) * qJD(4);
t245 = qJD(1) * t267 + (-qJD(4) * t278 - qJD(2)) * t301 + t320;
t244 = t278 * t325 + (-t268 + t319) * qJD(1) + t327;
t243 = -t274 * t282 + (-qJD(2) - t323) * t301 + (t250 + t259) * qJD(1) + t320;
t242 = t300 * t323 + t274 * t281 + (-t249 - t260 + t319) * qJD(1) + t327;
t241 = -t259 * t281 + t260 * t282 + (t249 * t301 - t250 * t300) * qJD(4);
t1 = m(3) * (t251 ^ 2 + t252 ^ 2) / 0.2e1 + m(4) * (t247 ^ 2 + t248 ^ 2) / 0.2e1 + m(5) * (t244 ^ 2 + t245 ^ 2 + t246 ^ 2) / 0.2e1 + qJD(4) * t301 * ((t301 * t275 + t304 * t300) * qJD(1) + (t301 ^ 2 * t261 + (t306 * t300 + (t262 - t307) * t301) * t300) * qJD(4)) / 0.2e1 + ((t300 * t275 - t304 * t301) * qJD(1) + (t300 ^ 2 * t262 + (t307 * t301 + (t261 - t306) * t300) * t301) * qJD(4)) * t325 / 0.2e1 + m(6) * (t241 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + t282 * (t336 * t300 + t303 * t301) / 0.2e1 + t281 * (t303 * t300 - t336 * t301) / 0.2e1 + (((-t289 * t263 + t265 * t290) * t301 + (-t264 * t289 + t266 * t290) * t300) * qJD(4) + (-t255 * t286 + t257 * t287) * t282 + (-t256 * t286 + t258 * t287) * t281 + (-t286 * t272 + t287 * t273 - t289 * t276 + t290 * t277) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t284 ^ 2 + t285 ^ 2) + Icges(2,3) + Icges(3,1) + Icges(4,1) * t298 ^ 2 + (-0.2e1 * Icges(4,4) * t298 + Icges(4,2) * t297) * t297) * qJD(1) ^ 2 / 0.2e1;
T = t1;
