% Calculate kinetic energy for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:37
% EndTime: 2019-12-05 15:42:37
% DurationCPUTime: 0.55s
% Computational Cost: add. (771->128), mult. (533->213), div. (0->0), fcn. (440->8), ass. (0->78)
t287 = cos(pkin(9));
t318 = t287 * pkin(3);
t284 = pkin(9) + qJ(4);
t277 = sin(t284);
t317 = Icges(5,4) * t277;
t279 = cos(t284);
t316 = Icges(5,4) * t279;
t281 = qJ(5) + t284;
t274 = sin(t281);
t315 = Icges(6,4) * t274;
t275 = cos(t281);
t314 = Icges(6,4) * t275;
t285 = pkin(8) + qJ(2);
t278 = sin(t285);
t280 = cos(t285);
t269 = pkin(2) * t278 - qJ(3) * t280;
t312 = pkin(6) * t280 - t318 * t278 - t269;
t311 = pkin(4) * t279;
t309 = qJD(4) * t278;
t308 = qJD(4) * t280;
t307 = qJD(4) + qJD(5);
t306 = pkin(4) * qJD(4) * t277;
t262 = qJD(2) * (pkin(2) * t280 + qJ(3) * t278);
t305 = -qJD(3) * t280 + qJD(2) * (pkin(6) * t278 + t318 * t280) + t262;
t286 = sin(pkin(9));
t304 = rSges(4,1) * t287 - rSges(4,2) * t286;
t303 = rSges(5,1) * t279 - rSges(5,2) * t277;
t302 = rSges(6,1) * t275 - rSges(6,2) * t274;
t301 = Icges(5,1) * t279 - t317;
t300 = Icges(6,1) * t275 - t315;
t299 = -Icges(5,2) * t277 + t316;
t298 = -Icges(6,2) * t274 + t314;
t297 = Icges(5,5) * t279 - Icges(5,6) * t277;
t296 = Icges(6,5) * t275 - Icges(6,6) * t274;
t252 = -Icges(5,6) * t280 + t299 * t278;
t254 = -Icges(5,5) * t280 + t301 * t278;
t295 = t252 * t277 - t254 * t279;
t253 = Icges(5,6) * t278 + t299 * t280;
t255 = Icges(5,5) * t278 + t301 * t280;
t294 = -t253 * t277 + t255 * t279;
t266 = Icges(5,2) * t279 + t317;
t267 = Icges(5,1) * t277 + t316;
t293 = -t266 * t277 + t267 * t279;
t263 = t307 * t278;
t264 = t307 * t280;
t292 = (Icges(6,5) * t274 + Icges(6,6) * t275) * qJD(2) - (-Icges(6,3) * t280 + t296 * t278) * t264 + (Icges(6,3) * t278 + t296 * t280) * t263;
t244 = -Icges(6,6) * t280 + t298 * t278;
t245 = Icges(6,6) * t278 + t298 * t280;
t246 = -Icges(6,5) * t280 + t300 * t278;
t247 = Icges(6,5) * t278 + t300 * t280;
t259 = Icges(6,2) * t275 + t315;
t260 = Icges(6,1) * t274 + t314;
t291 = (-t245 * t274 + t247 * t275) * t263 - (-t244 * t274 + t246 * t275) * t264 + (-t259 * t274 + t260 * t275) * qJD(2);
t290 = qJD(1) ^ 2;
t289 = qJD(2) ^ 2;
t273 = qJD(3) * t278;
t271 = rSges(3,1) * t280 - rSges(3,2) * t278;
t270 = rSges(3,1) * t278 + rSges(3,2) * t280;
t268 = rSges(5,1) * t277 + rSges(5,2) * t279;
t265 = Icges(5,5) * t277 + Icges(5,6) * t279;
t261 = rSges(6,1) * t274 + rSges(6,2) * t275;
t257 = rSges(5,3) * t278 + t303 * t280;
t256 = -rSges(5,3) * t280 + t303 * t278;
t251 = Icges(5,3) * t278 + t297 * t280;
t250 = -Icges(5,3) * t280 + t297 * t278;
t249 = rSges(6,3) * t278 + t302 * t280;
t248 = -rSges(6,3) * t280 + t302 * t278;
t239 = pkin(7) * t278 + t311 * t280;
t238 = -pkin(7) * t280 + t311 * t278;
t237 = qJD(2) * t278 * rSges(4,3) + t262 + (qJD(2) * t304 - qJD(3)) * t280;
t236 = t273 + (t280 * rSges(4,3) - t304 * t278 - t269) * qJD(2);
t235 = qJD(1) + (t256 * t278 + t257 * t280) * qJD(4);
t234 = qJD(2) * t257 - t268 * t309 + t305;
t233 = -t268 * t308 + t273 + (-t256 + t312) * qJD(2);
t232 = -t278 * t306 - t261 * t263 + (t239 + t249) * qJD(2) + t305;
t231 = -t280 * t306 - t261 * t264 + t273 + (-t238 - t248 + t312) * qJD(2);
t230 = t248 * t263 + t249 * t264 + qJD(1) + (t238 * t278 + t239 * t280) * qJD(4);
t1 = m(2) * t290 / 0.2e1 + m(3) * (t290 + (t270 ^ 2 + t271 ^ 2) * t289) / 0.2e1 + m(4) * (t236 ^ 2 + t237 ^ 2 + t290) / 0.2e1 + m(5) * (t233 ^ 2 + t234 ^ 2 + t235 ^ 2) / 0.2e1 + ((t278 * t265 + t293 * t280) * qJD(2) + (t278 ^ 2 * t251 + (t295 * t280 + (-t250 + t294) * t278) * t280) * qJD(4)) * t309 / 0.2e1 - ((-t280 * t265 + t293 * t278) * qJD(2) + (t280 ^ 2 * t250 + (t294 * t278 + (-t251 + t295) * t280) * t278) * qJD(4)) * t308 / 0.2e1 + m(6) * (t230 ^ 2 + t231 ^ 2 + t232 ^ 2) / 0.2e1 + t263 * (t292 * t278 + t291 * t280) / 0.2e1 - t264 * (t291 * t278 - t292 * t280) / 0.2e1 + (Icges(3,3) + Icges(4,2) * t287 ^ 2 + (Icges(4,1) * t286 + 0.2e1 * Icges(4,4) * t287) * t286) * t289 / 0.2e1 + (((t253 * t279 + t255 * t277) * t278 - (t252 * t279 + t254 * t277) * t280) * qJD(4) + (t245 * t275 + t247 * t274) * t263 - (t244 * t275 + t246 * t274) * t264 + (t275 * t259 + t274 * t260 + t279 * t266 + t277 * t267) * qJD(2)) * qJD(2) / 0.2e1;
T = t1;
