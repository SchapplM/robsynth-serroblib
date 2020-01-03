% Calculate kinetic energy for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:30
% EndTime: 2019-12-31 17:59:30
% DurationCPUTime: 0.60s
% Computational Cost: add. (609->154), mult. (733->253), div. (0->0), fcn. (682->8), ass. (0->81)
t277 = sin(qJ(1));
t304 = pkin(1) * t277;
t276 = sin(qJ(4));
t303 = Icges(5,4) * t276;
t279 = cos(qJ(4));
t302 = Icges(5,4) * t279;
t274 = qJ(1) + pkin(8);
t272 = sin(t274);
t301 = t272 * t279;
t273 = cos(t274);
t300 = t273 * t279;
t275 = sin(qJ(5));
t299 = t275 * t276;
t278 = cos(qJ(5));
t298 = t276 * t278;
t280 = cos(qJ(1));
t271 = qJD(1) * t280 * pkin(1);
t297 = qJD(1) * (pkin(2) * t273 + qJ(3) * t272) + t271;
t296 = qJD(4) * t272;
t295 = qJD(4) * t273;
t294 = qJD(5) * t279;
t293 = qJD(1) * t273 * pkin(6) + t297;
t292 = -pkin(2) * t272 + qJ(3) * t273 - t304;
t291 = pkin(4) * t276 - pkin(7) * t279;
t290 = rSges(5,1) * t276 + rSges(5,2) * t279;
t289 = Icges(5,1) * t276 + t302;
t288 = Icges(5,2) * t279 + t303;
t287 = Icges(5,5) * t276 + Icges(5,6) * t279;
t239 = Icges(5,6) * t273 + t288 * t272;
t241 = Icges(5,5) * t273 + t289 * t272;
t286 = -t239 * t279 - t241 * t276;
t240 = Icges(5,6) * t272 - t288 * t273;
t242 = Icges(5,5) * t272 - t289 * t273;
t285 = t240 * t279 + t242 * t276;
t262 = -Icges(5,2) * t276 + t302;
t263 = Icges(5,1) * t279 - t303;
t284 = t262 * t279 + t263 * t276;
t283 = -pkin(6) * t272 + t292;
t281 = qJD(2) ^ 2;
t270 = qJD(5) * t276 + qJD(1);
t269 = qJD(3) * t272;
t267 = pkin(4) * t279 + pkin(7) * t276;
t266 = rSges(2,1) * t280 - rSges(2,2) * t277;
t265 = rSges(5,1) * t279 - rSges(5,2) * t276;
t264 = rSges(2,1) * t277 + rSges(2,2) * t280;
t261 = Icges(5,5) * t279 - Icges(5,6) * t276;
t258 = -t272 * t294 + t295;
t257 = t273 * t294 + t296;
t256 = t291 * t273;
t255 = t291 * t272;
t254 = rSges(6,3) * t276 + (rSges(6,1) * t278 - rSges(6,2) * t275) * t279;
t253 = Icges(6,5) * t276 + (Icges(6,1) * t278 - Icges(6,4) * t275) * t279;
t252 = Icges(6,6) * t276 + (Icges(6,4) * t278 - Icges(6,2) * t275) * t279;
t251 = Icges(6,3) * t276 + (Icges(6,5) * t278 - Icges(6,6) * t275) * t279;
t250 = t272 * t275 - t273 * t298;
t249 = t272 * t278 + t273 * t299;
t248 = t272 * t298 + t273 * t275;
t247 = -t272 * t299 + t273 * t278;
t246 = t271 + qJD(1) * (rSges(3,1) * t273 - rSges(3,2) * t272);
t245 = (-rSges(3,1) * t272 - rSges(3,2) * t273 - t304) * qJD(1);
t244 = rSges(5,3) * t272 - t290 * t273;
t243 = rSges(5,3) * t273 + t290 * t272;
t238 = Icges(5,3) * t272 - t287 * t273;
t237 = Icges(5,3) * t273 + t287 * t272;
t236 = -qJD(3) * t273 + qJD(1) * (-rSges(4,2) * t273 + rSges(4,3) * t272) + t297;
t235 = t269 + (rSges(4,2) * t272 + rSges(4,3) * t273 + t292) * qJD(1);
t234 = rSges(6,1) * t250 + rSges(6,2) * t249 + rSges(6,3) * t300;
t233 = rSges(6,1) * t248 + rSges(6,2) * t247 - rSges(6,3) * t301;
t232 = Icges(6,1) * t250 + Icges(6,4) * t249 + Icges(6,5) * t300;
t231 = Icges(6,1) * t248 + Icges(6,4) * t247 - Icges(6,5) * t301;
t230 = Icges(6,4) * t250 + Icges(6,2) * t249 + Icges(6,6) * t300;
t229 = Icges(6,4) * t248 + Icges(6,2) * t247 - Icges(6,6) * t301;
t228 = Icges(6,5) * t250 + Icges(6,6) * t249 + Icges(6,3) * t300;
t227 = Icges(6,5) * t248 + Icges(6,6) * t247 - Icges(6,3) * t301;
t226 = qJD(2) + (-t243 * t272 + t244 * t273) * qJD(4);
t225 = qJD(1) * t243 + (-qJD(4) * t265 - qJD(3)) * t273 + t293;
t224 = t265 * t296 + t269 + (-t244 + t283) * qJD(1);
t223 = qJD(1) * t255 + t233 * t270 - t254 * t258 + (-qJD(4) * t267 - qJD(3)) * t273 + t293;
t222 = t267 * t296 - t234 * t270 + t254 * t257 + t269 + (t256 + t283) * qJD(1);
t221 = -t233 * t257 + t234 * t258 + qJD(2) + (-t255 * t272 - t256 * t273) * qJD(4);
t1 = m(3) * (t245 ^ 2 + t246 ^ 2 + t281) / 0.2e1 + m(4) * (t235 ^ 2 + t236 ^ 2 + t281) / 0.2e1 + m(5) * (t224 ^ 2 + t225 ^ 2 + t226 ^ 2) / 0.2e1 + ((t273 * t261 + t284 * t272) * qJD(1) + (t273 ^ 2 * t237 + (t285 * t272 + (t238 - t286) * t273) * t272) * qJD(4)) * t295 / 0.2e1 + ((t272 * t261 - t284 * t273) * qJD(1) + (t272 ^ 2 * t238 + (t286 * t273 + (t237 - t285) * t272) * t273) * qJD(4)) * t296 / 0.2e1 + qJD(1) * ((-t276 * t262 + t279 * t263) * qJD(1) + ((-t239 * t276 + t241 * t279) * t273 + (-t240 * t276 + t242 * t279) * t272) * qJD(4)) / 0.2e1 + m(6) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + t258 * ((-t227 * t301 + t247 * t229 + t248 * t231) * t258 + (-t228 * t301 + t230 * t247 + t232 * t248) * t257 + (t247 * t252 + t248 * t253 - t251 * t301) * t270) / 0.2e1 + t257 * ((t227 * t300 + t229 * t249 + t231 * t250) * t258 + (t228 * t300 + t249 * t230 + t250 * t232) * t257 + (t249 * t252 + t250 * t253 + t251 * t300) * t270) / 0.2e1 + t270 * ((t227 * t258 + t228 * t257 + t251 * t270) * t276 + ((-t229 * t275 + t231 * t278) * t258 + (-t230 * t275 + t232 * t278) * t257 + (-t252 * t275 + t253 * t278) * t270) * t279) / 0.2e1 + (m(2) * (t264 ^ 2 + t266 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,1)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
