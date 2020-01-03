% Calculate kinetic energy for
% S5RPPRR9
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:18
% EndTime: 2019-12-31 18:02:18
% DurationCPUTime: 0.65s
% Computational Cost: add. (598->155), mult. (1296->260), div. (0->0), fcn. (1524->8), ass. (0->81)
t316 = cos(qJ(1));
t315 = sin(qJ(1));
t314 = cos(pkin(8));
t313 = sin(pkin(8));
t285 = sin(qJ(4));
t312 = Icges(5,4) * t285;
t287 = cos(qJ(4));
t311 = Icges(5,4) * t287;
t267 = -t313 * t315 - t314 * t316;
t310 = t267 * t285;
t268 = t313 * t316 - t314 * t315;
t309 = t268 * t285;
t284 = sin(qJ(5));
t308 = t284 * t287;
t286 = cos(qJ(5));
t307 = t286 * t287;
t306 = qJD(4) * t267;
t305 = qJD(4) * t268;
t304 = qJD(4) * (-t285 * rSges(5,1) - t287 * rSges(5,2));
t303 = qJD(4) * (-t285 * pkin(4) + t287 * pkin(7));
t302 = qJD(5) * t285;
t274 = pkin(1) * t315 - qJ(2) * t316;
t301 = -pkin(2) * t315 - t274;
t300 = -pkin(4) * t287 - pkin(7) * t285;
t299 = -qJD(2) * t316 + qJD(1) * (pkin(1) * t316 + qJ(2) * t315);
t298 = -rSges(5,1) * t287 + rSges(5,2) * t285;
t297 = -Icges(5,1) * t287 + t312;
t296 = Icges(5,2) * t285 - t311;
t295 = -Icges(5,5) * t287 + Icges(5,6) * t285;
t245 = -Icges(5,6) * t267 + t268 * t296;
t247 = -Icges(5,5) * t267 + t268 * t297;
t294 = -t245 * t285 + t247 * t287;
t246 = Icges(5,6) * t268 + t267 * t296;
t248 = Icges(5,5) * t268 + t267 * t297;
t293 = t246 * t285 - t248 * t287;
t271 = -Icges(5,2) * t287 - t312;
t272 = -Icges(5,1) * t285 - t311;
t292 = t271 * t285 - t272 * t287;
t291 = t268 * pkin(3) + t267 * pkin(6) + t301;
t290 = qJD(1) * t316 * pkin(2) + t299;
t289 = qJD(1) * (-t267 * pkin(3) + t268 * pkin(6)) + t290;
t283 = qJD(2) * t315;
t278 = qJD(5) * t287 + qJD(1);
t276 = rSges(2,1) * t316 - rSges(2,2) * t315;
t275 = rSges(2,1) * t315 + rSges(2,2) * t316;
t270 = -Icges(5,5) * t285 - Icges(5,6) * t287;
t266 = t287 * rSges(6,3) + (-rSges(6,1) * t286 + rSges(6,2) * t284) * t285;
t265 = Icges(6,5) * t287 + (-Icges(6,1) * t286 + Icges(6,4) * t284) * t285;
t264 = Icges(6,6) * t287 + (-Icges(6,4) * t286 + Icges(6,2) * t284) * t285;
t263 = Icges(6,3) * t287 + (-Icges(6,5) * t286 + Icges(6,6) * t284) * t285;
t262 = qJD(1) * (rSges(3,1) * t316 + rSges(3,3) * t315) + t299;
t261 = t283 + (-rSges(3,1) * t315 + rSges(3,3) * t316 - t274) * qJD(1);
t258 = -t268 * t302 - t306;
t257 = -t267 * t302 + t305;
t256 = t300 * t267;
t255 = t300 * t268;
t254 = -t267 * t307 + t268 * t284;
t253 = t267 * t308 + t268 * t286;
t252 = -t267 * t284 - t268 * t307;
t251 = -t267 * t286 + t268 * t308;
t250 = t268 * rSges(5,3) + t267 * t298;
t249 = -t267 * rSges(5,3) + t268 * t298;
t244 = Icges(5,3) * t268 + t267 * t295;
t243 = -Icges(5,3) * t267 + t268 * t295;
t242 = qJD(1) * (-t267 * rSges(4,1) - t268 * rSges(4,2)) + t290;
t241 = t283 + (t268 * rSges(4,1) - t267 * rSges(4,2) + t301) * qJD(1);
t240 = t254 * rSges(6,1) + t253 * rSges(6,2) - rSges(6,3) * t310;
t239 = t252 * rSges(6,1) + t251 * rSges(6,2) - rSges(6,3) * t309;
t238 = Icges(6,1) * t254 + Icges(6,4) * t253 - Icges(6,5) * t310;
t237 = Icges(6,1) * t252 + Icges(6,4) * t251 - Icges(6,5) * t309;
t236 = Icges(6,4) * t254 + Icges(6,2) * t253 - Icges(6,6) * t310;
t235 = Icges(6,4) * t252 + Icges(6,2) * t251 - Icges(6,6) * t309;
t234 = Icges(6,5) * t254 + Icges(6,6) * t253 - Icges(6,3) * t310;
t233 = Icges(6,5) * t252 + Icges(6,6) * t251 - Icges(6,3) * t309;
t232 = qJD(1) * t250 - t268 * t304 + t289;
t231 = -t267 * t304 + t283 + (-t249 + t291) * qJD(1);
t230 = -qJD(3) + (t249 * t268 + t250 * t267) * qJD(4);
t229 = qJD(1) * t256 + t278 * t240 - t257 * t266 - t268 * t303 + t289;
t228 = -t267 * t303 - t278 * t239 + t258 * t266 + t283 + (-t255 + t291) * qJD(1);
t227 = t257 * t239 - t258 * t240 - qJD(3) + (t255 * t268 + t256 * t267) * qJD(4);
t1 = m(3) * (t261 ^ 2 + t262 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + m(5) * (t230 ^ 2 + t231 ^ 2 + t232 ^ 2) / 0.2e1 + ((t292 * t267 + t268 * t270) * qJD(1) + (t268 ^ 2 * t244 + (t294 * t267 + (-t243 + t293) * t268) * t267) * qJD(4)) * t305 / 0.2e1 - ((-t267 * t270 + t268 * t292) * qJD(1) + (t267 ^ 2 * t243 + (t293 * t268 + (-t244 + t294) * t267) * t268) * qJD(4)) * t306 / 0.2e1 + qJD(1) * ((-t287 * t271 - t285 * t272) * qJD(1) + ((-t287 * t246 - t285 * t248) * t268 - (-t287 * t245 - t285 * t247) * t267) * qJD(4)) / 0.2e1 + m(6) * (t227 ^ 2 + t228 ^ 2 + t229 ^ 2) / 0.2e1 + t257 * ((-t234 * t310 + t253 * t236 + t254 * t238) * t257 + (-t233 * t310 + t253 * t235 + t254 * t237) * t258 + (t253 * t264 + t254 * t265 - t263 * t310) * t278) / 0.2e1 + t258 * ((-t234 * t309 + t251 * t236 + t252 * t238) * t257 + (-t233 * t309 + t251 * t235 + t252 * t237) * t258 + (t251 * t264 + t252 * t265 - t263 * t309) * t278) / 0.2e1 + t278 * ((t233 * t258 + t234 * t257 + t263 * t278) * t287 + ((t236 * t284 - t238 * t286) * t257 + (t235 * t284 - t237 * t286) * t258 + (t264 * t284 - t265 * t286) * t278) * t285) / 0.2e1 + (m(2) * (t275 ^ 2 + t276 ^ 2) + Icges(2,3) + Icges(3,2) + Icges(4,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
