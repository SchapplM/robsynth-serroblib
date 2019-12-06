% Calculate kinetic energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:44
% EndTime: 2019-12-05 17:28:45
% DurationCPUTime: 0.54s
% Computational Cost: add. (609->130), mult. (605->214), div. (0->0), fcn. (578->10), ass. (0->63)
t292 = cos(pkin(9));
t293 = cos(pkin(8));
t311 = t292 * t293;
t291 = sin(pkin(8));
t288 = pkin(9) + qJ(5);
t284 = sin(t288);
t286 = cos(t288);
t289 = qJ(1) + pkin(7);
t287 = cos(t289);
t285 = sin(t289);
t316 = t285 * t293;
t266 = t284 * t316 + t286 * t287;
t267 = t284 * t287 - t286 * t316;
t313 = t287 * t293;
t268 = -t284 * t313 + t285 * t286;
t269 = t284 * t285 + t286 * t313;
t314 = t287 * t291;
t317 = t285 * t291;
t299 = (Icges(6,5) * t267 + Icges(6,6) * t266 - Icges(6,3) * t317) * t285 - (Icges(6,5) * t269 + Icges(6,6) * t268 + Icges(6,3) * t314) * t287;
t323 = t299 * t291;
t295 = sin(qJ(1));
t320 = t295 * pkin(1);
t296 = cos(qJ(1));
t319 = t296 * pkin(1);
t290 = sin(pkin(9));
t318 = t285 * t290;
t315 = t287 * t290;
t312 = t290 * t293;
t309 = qJD(1) * (-pkin(2) * t285 + qJ(3) * t287) + qJD(3) * t285;
t308 = qJD(4) * t291;
t307 = qJD(5) * (-t293 * rSges(6,3) + (rSges(6,1) * t286 - rSges(6,2) * t284) * t291);
t306 = qJD(5) * t291;
t305 = pkin(4) * t311;
t304 = -pkin(2) * t287 - qJ(3) * t285 - t319;
t279 = -qJD(4) * t293 + qJD(2);
t300 = pkin(3) * t293 + qJ(4) * t291;
t303 = -qJD(1) * t300 * t285 + t287 * t308 + t309;
t302 = -t300 * t287 + t304;
t301 = -rSges(4,1) * t293 + rSges(4,2) * t291;
t297 = qJD(2) ^ 2;
t282 = qJD(3) * t287;
t280 = -qJD(5) * t293 + qJD(1);
t278 = rSges(2,1) * t296 - rSges(2,2) * t295;
t277 = -rSges(2,1) * t295 - rSges(2,2) * t296;
t271 = (-rSges(3,1) * t287 + rSges(3,2) * t285 - t319) * qJD(1);
t270 = (-rSges(3,1) * t285 - rSges(3,2) * t287 - t320) * qJD(1);
t264 = -Icges(6,5) * t293 + (Icges(6,1) * t286 - Icges(6,4) * t284) * t291;
t263 = -Icges(6,6) * t293 + (Icges(6,4) * t286 - Icges(6,2) * t284) * t291;
t262 = -Icges(6,3) * t293 + (Icges(6,5) * t286 - Icges(6,6) * t284) * t291;
t261 = t282 + (-t285 * rSges(4,3) + t301 * t287 + t304) * qJD(1);
t260 = (t287 * rSges(4,3) + t301 * t285 - t320) * qJD(1) + t309;
t259 = rSges(6,1) * t269 + rSges(6,2) * t268 + rSges(6,3) * t314;
t258 = rSges(6,1) * t267 + rSges(6,2) * t266 - rSges(6,3) * t317;
t257 = Icges(6,1) * t269 + Icges(6,4) * t268 + Icges(6,5) * t314;
t256 = Icges(6,1) * t267 + Icges(6,4) * t266 - Icges(6,5) * t317;
t255 = Icges(6,4) * t269 + Icges(6,2) * t268 + Icges(6,6) * t314;
t254 = Icges(6,4) * t267 + Icges(6,2) * t266 - Icges(6,6) * t317;
t251 = -t285 * t308 + t282 + (-(t287 * t311 + t318) * rSges(5,1) - (t285 * t292 - t287 * t312) * rSges(5,2) - rSges(5,3) * t314 + t302) * qJD(1);
t250 = (-t320 + (-t285 * t311 + t315) * rSges(5,1) + (t285 * t312 + t287 * t292) * rSges(5,2) - rSges(5,3) * t317) * qJD(1) + t303;
t249 = (-t258 * t287 - t259 * t285) * t306 + t279;
t248 = -t280 * t259 + t282 + (-qJD(4) * t285 + t287 * t307) * t291 + (-pkin(4) * t318 + (-pkin(6) * t291 - t305) * t287 + t302) * qJD(1);
t247 = t280 * t258 + (pkin(4) * t315 - t320) * qJD(1) + (-qJD(1) * t305 + (-pkin(6) * qJD(1) + t307) * t291) * t285 + t303;
t1 = m(3) * (t270 ^ 2 + t271 ^ 2 + t297) / 0.2e1 + m(4) * (t260 ^ 2 + t261 ^ 2 + t297) / 0.2e1 + m(5) * (t250 ^ 2 + t251 ^ 2 + t279 ^ 2) / 0.2e1 + m(6) * (t247 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + t280 * ((-t293 * t262 + (-t263 * t284 + t264 * t286) * t291) * t280 + ((-(-t254 * t284 + t256 * t286) * t285 + (-t255 * t284 + t257 * t286) * t287) * t291 + t299 * t293) * t306) / 0.2e1 - t285 * ((-t262 * t317 + t263 * t266 + t264 * t267) * t280 + ((t255 * t266 + t257 * t267) * t287 + (-t266 * t254 - t267 * t256 + t323) * t285) * t306) * t306 / 0.2e1 + t287 * ((t262 * t314 + t263 * t268 + t264 * t269) * t280 + (-(t254 * t268 + t256 * t269) * t285 + (t268 * t255 + t269 * t257 - t323) * t287) * t306) * t306 / 0.2e1 + (m(2) * (t277 ^ 2 + t278 ^ 2) + Icges(2,3) + Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t293 ^ 2 + ((Icges(4,1) + Icges(5,1) * t292 ^ 2 + (-0.2e1 * Icges(5,4) * t292 + Icges(5,2) * t290) * t290) * t291 + 0.2e1 * (-Icges(5,5) * t292 + Icges(5,6) * t290 + Icges(4,4)) * t293) * t291) * qJD(1) ^ 2 / 0.2e1;
T = t1;
