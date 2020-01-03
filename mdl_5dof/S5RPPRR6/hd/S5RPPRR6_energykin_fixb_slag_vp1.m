% Calculate kinetic energy for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:39
% EndTime: 2019-12-31 17:57:40
% DurationCPUTime: 0.70s
% Computational Cost: add. (882->159), mult. (761->261), div. (0->0), fcn. (710->10), ass. (0->88)
t301 = sin(qJ(1));
t332 = t301 * pkin(1);
t298 = cos(pkin(9));
t331 = t298 * pkin(3);
t295 = pkin(9) + qJ(4);
t291 = sin(t295);
t330 = Icges(5,4) * t291;
t293 = cos(t295);
t329 = Icges(5,4) * t293;
t296 = qJ(1) + pkin(8);
t292 = sin(t296);
t328 = t291 * t292;
t294 = cos(t296);
t327 = t291 * t294;
t300 = sin(qJ(5));
t326 = t292 * t300;
t302 = cos(qJ(5));
t325 = t292 * t302;
t324 = t294 * t300;
t323 = t294 * t302;
t303 = cos(qJ(1));
t290 = qJD(1) * t303 * pkin(1);
t321 = qJD(1) * (t294 * pkin(2) + t292 * qJ(3)) + t290;
t320 = qJD(4) * t292;
t319 = qJD(4) * t294;
t318 = qJD(5) * t291;
t317 = -t292 * pkin(2) + t294 * qJ(3) - t332;
t316 = pkin(6) * t294 - t292 * t331 + t317;
t315 = pkin(4) * t293 + pkin(7) * t291;
t297 = sin(pkin(9));
t314 = rSges(4,1) * t298 - rSges(4,2) * t297;
t313 = rSges(5,1) * t293 - rSges(5,2) * t291;
t312 = Icges(5,1) * t293 - t330;
t311 = -Icges(5,2) * t291 + t329;
t310 = Icges(5,5) * t293 - Icges(5,6) * t291;
t258 = -Icges(5,6) * t294 + t292 * t311;
t260 = -Icges(5,5) * t294 + t292 * t312;
t309 = t258 * t291 - t260 * t293;
t259 = Icges(5,6) * t292 + t294 * t311;
t261 = Icges(5,5) * t292 + t294 * t312;
t308 = -t259 * t291 + t261 * t293;
t280 = Icges(5,2) * t293 + t330;
t281 = Icges(5,1) * t291 + t329;
t307 = -t280 * t291 + t281 * t293;
t306 = -qJD(3) * t294 + qJD(1) * (pkin(6) * t292 + t294 * t331) + t321;
t304 = qJD(2) ^ 2;
t288 = qJD(3) * t292;
t287 = -qJD(5) * t293 + qJD(1);
t286 = t303 * rSges(2,1) - t301 * rSges(2,2);
t285 = t301 * rSges(2,1) + t303 * rSges(2,2);
t284 = t291 * pkin(4) - t293 * pkin(7);
t282 = t291 * rSges(5,1) + t293 * rSges(5,2);
t279 = Icges(5,5) * t291 + Icges(5,6) * t293;
t277 = t292 * t318 - t319;
t276 = t294 * t318 + t320;
t275 = t293 * t323 + t326;
t274 = -t293 * t324 + t325;
t273 = t293 * t325 - t324;
t272 = -t293 * t326 - t323;
t271 = t290 + qJD(1) * (t294 * rSges(3,1) - t292 * rSges(3,2));
t270 = (-t292 * rSges(3,1) - t294 * rSges(3,2) - t332) * qJD(1);
t269 = t315 * t294;
t268 = t315 * t292;
t267 = -t293 * rSges(6,3) + (rSges(6,1) * t302 - rSges(6,2) * t300) * t291;
t266 = -Icges(6,5) * t293 + (Icges(6,1) * t302 - Icges(6,4) * t300) * t291;
t265 = -Icges(6,6) * t293 + (Icges(6,4) * t302 - Icges(6,2) * t300) * t291;
t264 = -Icges(6,3) * t293 + (Icges(6,5) * t302 - Icges(6,6) * t300) * t291;
t263 = t292 * rSges(5,3) + t294 * t313;
t262 = -t294 * rSges(5,3) + t292 * t313;
t257 = Icges(5,3) * t292 + t294 * t310;
t256 = -Icges(5,3) * t294 + t292 * t310;
t253 = qJD(1) * t292 * rSges(4,3) + (qJD(1) * t314 - qJD(3)) * t294 + t321;
t252 = t288 + (t294 * rSges(4,3) - t292 * t314 + t317) * qJD(1);
t251 = t275 * rSges(6,1) + t274 * rSges(6,2) + rSges(6,3) * t327;
t250 = t273 * rSges(6,1) + t272 * rSges(6,2) + rSges(6,3) * t328;
t249 = Icges(6,1) * t275 + Icges(6,4) * t274 + Icges(6,5) * t327;
t248 = Icges(6,1) * t273 + Icges(6,4) * t272 + Icges(6,5) * t328;
t247 = Icges(6,4) * t275 + Icges(6,2) * t274 + Icges(6,6) * t327;
t246 = Icges(6,4) * t273 + Icges(6,2) * t272 + Icges(6,6) * t328;
t245 = Icges(6,5) * t275 + Icges(6,6) * t274 + Icges(6,3) * t327;
t244 = Icges(6,5) * t273 + Icges(6,6) * t272 + Icges(6,3) * t328;
t243 = qJD(2) + (t262 * t292 + t263 * t294) * qJD(4);
t242 = qJD(1) * t263 - t282 * t320 + t306;
t241 = -t282 * t319 + t288 + (-t262 + t316) * qJD(1);
t240 = t276 * t250 - t277 * t251 + qJD(2) + (t268 * t292 + t269 * t294) * qJD(4);
t239 = qJD(1) * t269 + t287 * t251 - t276 * t267 - t284 * t320 + t306;
t238 = -t284 * t319 - t287 * t250 + t277 * t267 + t288 + (-t268 + t316) * qJD(1);
t1 = m(3) * (t270 ^ 2 + t271 ^ 2 + t304) / 0.2e1 + m(4) * (t252 ^ 2 + t253 ^ 2 + t304) / 0.2e1 + m(5) * (t241 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + ((t292 * t279 + t294 * t307) * qJD(1) + (t292 ^ 2 * t257 + (t309 * t294 + (-t256 + t308) * t292) * t294) * qJD(4)) * t320 / 0.2e1 - ((-t294 * t279 + t292 * t307) * qJD(1) + (t294 ^ 2 * t256 + (t308 * t292 + (-t257 + t309) * t294) * t292) * qJD(4)) * t319 / 0.2e1 + qJD(1) * ((t293 * t280 + t291 * t281) * qJD(1) + ((t293 * t259 + t291 * t261) * t292 - (t293 * t258 + t291 * t260) * t294) * qJD(4)) / 0.2e1 + m(6) * (t238 ^ 2 + t239 ^ 2 + t240 ^ 2) / 0.2e1 + t276 * ((t245 * t327 + t274 * t247 + t275 * t249) * t276 + (t244 * t327 + t274 * t246 + t275 * t248) * t277 + (t264 * t327 + t274 * t265 + t275 * t266) * t287) / 0.2e1 + t277 * ((t245 * t328 + t272 * t247 + t273 * t249) * t276 + (t244 * t328 + t272 * t246 + t273 * t248) * t277 + (t264 * t328 + t272 * t265 + t273 * t266) * t287) / 0.2e1 + t287 * ((-t244 * t277 - t245 * t276 - t264 * t287) * t293 + ((-t247 * t300 + t249 * t302) * t276 + (-t246 * t300 + t248 * t302) * t277 + (-t265 * t300 + t266 * t302) * t287) * t291) / 0.2e1 + (m(2) * (t285 ^ 2 + t286 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t298 ^ 2 + (Icges(4,1) * t297 + 0.2e1 * Icges(4,4) * t298) * t297) * qJD(1) ^ 2 / 0.2e1;
T = t1;
