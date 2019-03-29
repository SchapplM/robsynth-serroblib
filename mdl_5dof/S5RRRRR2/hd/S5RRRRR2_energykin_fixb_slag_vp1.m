% Calculate kinetic energy for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_energykin_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:25:42
% EndTime: 2019-03-29 15:25:43
% DurationCPUTime: 0.98s
% Computational Cost: add. (1017->175), mult. (914->301), div. (0->0), fcn. (848->10), ass. (0->105)
t273 = qJ(1) + qJ(2);
t268 = sin(t273);
t317 = t268 ^ 2;
t270 = cos(t273);
t316 = t270 ^ 2;
t313 = pkin(1) * qJD(1);
t275 = sin(qJ(3));
t312 = Icges(4,4) * t275;
t278 = cos(qJ(3));
t311 = Icges(4,4) * t278;
t272 = qJ(3) + qJ(4);
t267 = sin(t272);
t310 = Icges(5,4) * t267;
t269 = cos(t272);
t309 = Icges(5,4) * t269;
t308 = t267 * t268;
t307 = t267 * t270;
t274 = sin(qJ(5));
t306 = t268 * t274;
t277 = cos(qJ(5));
t305 = t268 * t277;
t304 = t270 * t274;
t303 = t270 * t277;
t271 = qJD(1) + qJD(2);
t302 = t271 * t278;
t301 = (t316 + t317) * pkin(2) * qJD(3) * t278;
t265 = qJD(3) * t268;
t248 = qJD(4) * t268 + t265;
t300 = qJD(3) * t270;
t299 = qJD(5) * t267;
t298 = t275 * qJD(3);
t276 = sin(qJ(1));
t296 = t276 * t313;
t249 = (-qJD(3) - qJD(4)) * t270;
t295 = rSges(4,1) * t278 - rSges(4,2) * t275;
t294 = rSges(5,1) * t269 - rSges(5,2) * t267;
t293 = Icges(4,1) * t278 - t312;
t292 = Icges(5,1) * t269 - t310;
t291 = -Icges(4,2) * t275 + t311;
t290 = -Icges(5,2) * t267 + t309;
t289 = Icges(4,5) * t278 - Icges(4,6) * t275;
t288 = Icges(5,5) * t269 - Icges(5,6) * t267;
t232 = -Icges(4,6) * t270 + t268 * t291;
t235 = -Icges(4,5) * t270 + t268 * t293;
t287 = t232 * t275 - t235 * t278;
t233 = Icges(4,6) * t268 + t270 * t291;
t236 = Icges(4,5) * t268 + t270 * t293;
t286 = -t233 * t275 + t236 * t278;
t259 = Icges(4,2) * t278 + t312;
t260 = Icges(4,1) * t275 + t311;
t285 = -t259 * t275 + t260 * t278;
t279 = cos(qJ(1));
t266 = t279 * t313;
t284 = t266 + (-t268 * t298 + t270 * t302) * pkin(2);
t283 = (-Icges(5,3) * t270 + t268 * t288) * t249 + (Icges(5,3) * t268 + t270 * t288) * t248 + (Icges(5,5) * t267 + Icges(5,6) * t269) * t271;
t282 = (-t268 * t302 - t270 * t298) * pkin(2) - t296;
t222 = -Icges(5,6) * t270 + t268 * t290;
t223 = Icges(5,6) * t268 + t270 * t290;
t224 = -Icges(5,5) * t270 + t268 * t292;
t225 = Icges(5,5) * t268 + t270 * t292;
t251 = Icges(5,2) * t269 + t310;
t252 = Icges(5,1) * t267 + t309;
t281 = (-t223 * t267 + t225 * t269) * t248 + (-t222 * t267 + t224 * t269) * t249 + (-t251 * t267 + t252 * t269) * t271;
t263 = rSges(2,1) * t279 - rSges(2,2) * t276;
t262 = rSges(2,1) * t276 + rSges(2,2) * t279;
t261 = rSges(4,1) * t275 + rSges(4,2) * t278;
t258 = Icges(4,5) * t275 + Icges(4,6) * t278;
t255 = -qJD(5) * t269 + t271;
t253 = rSges(5,1) * t267 + rSges(5,2) * t269;
t247 = t269 * t303 + t306;
t246 = -t269 * t304 + t305;
t245 = t269 * t305 - t304;
t244 = -t269 * t306 - t303;
t243 = t266 + t271 * (rSges(3,1) * t270 - rSges(3,2) * t268);
t242 = -t296 - t271 * (rSges(3,1) * t268 + rSges(3,2) * t270);
t241 = rSges(4,3) * t268 + t270 * t295;
t240 = -rSges(4,3) * t270 + t268 * t295;
t239 = -rSges(6,3) * t269 + (rSges(6,1) * t277 - rSges(6,2) * t274) * t267;
t238 = t268 * t299 + t249;
t237 = t270 * t299 + t248;
t234 = -Icges(6,5) * t269 + (Icges(6,1) * t277 - Icges(6,4) * t274) * t267;
t231 = -Icges(6,6) * t269 + (Icges(6,4) * t277 - Icges(6,2) * t274) * t267;
t230 = Icges(4,3) * t268 + t270 * t289;
t229 = -Icges(4,3) * t270 + t268 * t289;
t228 = -Icges(6,3) * t269 + (Icges(6,5) * t277 - Icges(6,6) * t274) * t267;
t227 = rSges(5,3) * t268 + t270 * t294;
t226 = -rSges(5,3) * t270 + t268 * t294;
t219 = t241 * t271 - t261 * t265 + t266;
t218 = -t240 * t271 - t261 * t300 - t296;
t217 = rSges(6,1) * t247 + rSges(6,2) * t246 + rSges(6,3) * t307;
t216 = rSges(6,1) * t245 + rSges(6,2) * t244 + rSges(6,3) * t308;
t215 = Icges(6,1) * t247 + Icges(6,4) * t246 + Icges(6,5) * t307;
t214 = Icges(6,1) * t245 + Icges(6,4) * t244 + Icges(6,5) * t308;
t213 = Icges(6,4) * t247 + Icges(6,2) * t246 + Icges(6,6) * t307;
t212 = Icges(6,4) * t245 + Icges(6,2) * t244 + Icges(6,6) * t308;
t211 = Icges(6,5) * t247 + Icges(6,6) * t246 + Icges(6,3) * t307;
t210 = Icges(6,5) * t245 + Icges(6,6) * t244 + Icges(6,3) * t308;
t209 = (t240 * t268 + t241 * t270) * qJD(3);
t208 = t227 * t271 - t248 * t253 + t284;
t207 = -t226 * t271 + t249 * t253 + t282;
t206 = t226 * t248 - t227 * t249 + t301;
t205 = t217 * t255 - t237 * t239 + t284;
t204 = -t216 * t255 + t238 * t239 + t282;
t203 = t216 * t237 - t217 * t238 + t301;
t1 = m(3) * (t242 ^ 2 + t243 ^ 2) / 0.2e1 + t271 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t209 ^ 2 + t218 ^ 2 + t219 ^ 2) / 0.2e1 + ((t268 * t258 + t270 * t285) * t271 + (t317 * t230 + (t287 * t270 + (-t229 + t286) * t268) * t270) * qJD(3)) * t265 / 0.2e1 - ((-t270 * t258 + t268 * t285) * t271 + (t316 * t229 + (t286 * t268 + (-t230 + t287) * t270) * t268) * qJD(3)) * t300 / 0.2e1 + m(5) * (t206 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + t248 * (t283 * t268 + t281 * t270) / 0.2e1 + t249 * (t281 * t268 - t283 * t270) / 0.2e1 + m(6) * (t203 ^ 2 + t204 ^ 2 + t205 ^ 2) / 0.2e1 + t237 * ((t211 * t307 + t246 * t213 + t247 * t215) * t237 + (t210 * t307 + t212 * t246 + t214 * t247) * t238 + (t228 * t307 + t231 * t246 + t234 * t247) * t255) / 0.2e1 + t238 * ((t211 * t308 + t213 * t244 + t215 * t245) * t237 + (t210 * t308 + t244 * t212 + t245 * t214) * t238 + (t228 * t308 + t231 * t244 + t234 * t245) * t255) / 0.2e1 + t255 * ((-t210 * t238 - t211 * t237 - t228 * t255) * t269 + ((-t213 * t274 + t215 * t277) * t237 + (-t212 * t274 + t214 * t277) * t238 + (-t231 * t274 + t234 * t277) * t255) * t267) / 0.2e1 + (((t233 * t278 + t236 * t275) * t268 - (t232 * t278 + t235 * t275) * t270) * qJD(3) + (t223 * t269 + t225 * t267) * t248 + (t222 * t269 + t224 * t267) * t249 + (t269 * t251 + t267 * t252 + t278 * t259 + t275 * t260) * t271) * t271 / 0.2e1 + (m(2) * (t262 ^ 2 + t263 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
