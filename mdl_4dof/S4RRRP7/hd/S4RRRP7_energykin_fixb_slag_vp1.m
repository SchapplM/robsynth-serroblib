% Calculate kinetic energy for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP7_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP7_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:08
% EndTime: 2019-12-31 17:20:09
% DurationCPUTime: 1.13s
% Computational Cost: add. (497->151), mult. (1193->246), div. (0->0), fcn. (1202->6), ass. (0->86)
t317 = Icges(4,1) + Icges(5,1);
t316 = -Icges(4,4) + Icges(5,5);
t315 = Icges(5,4) + Icges(4,5);
t314 = Icges(4,2) + Icges(5,3);
t313 = -Icges(5,6) + Icges(4,6);
t312 = -Icges(4,3) - Icges(5,2);
t311 = rSges(5,1) + pkin(3);
t310 = rSges(5,3) + qJ(4);
t266 = sin(qJ(3));
t269 = cos(qJ(3));
t271 = cos(qJ(1));
t290 = t271 * t269;
t268 = sin(qJ(1));
t270 = cos(qJ(2));
t292 = t268 * t270;
t246 = t266 * t292 + t290;
t291 = t271 * t266;
t247 = t269 * t292 - t291;
t267 = sin(qJ(2));
t294 = t267 * t268;
t309 = t246 * t314 + t247 * t316 - t294 * t313;
t248 = -t268 * t269 + t270 * t291;
t249 = t266 * t268 + t270 * t290;
t293 = t267 * t271;
t308 = t248 * t314 + t249 * t316 - t293 * t313;
t307 = -t246 * t313 + t247 * t315 - t294 * t312;
t306 = -t248 * t313 + t249 * t315 - t293 * t312;
t305 = t246 * t316 + t247 * t317 + t294 * t315;
t304 = t248 * t316 + t249 * t317 + t293 * t315;
t303 = t313 * t270 + (t266 * t314 + t269 * t316) * t267;
t302 = t312 * t270 + (-t266 * t313 + t269 * t315) * t267;
t301 = -t315 * t270 + (t266 * t316 + t269 * t317) * t267;
t296 = Icges(3,4) * t267;
t295 = Icges(3,4) * t270;
t289 = rSges(5,2) * t294 + t246 * t310 + t247 * t311;
t288 = rSges(5,2) * t293 + t248 * t310 + t249 * t311;
t287 = -t270 * rSges(5,2) + (t266 * t310 + t269 * t311) * t267;
t282 = pkin(2) * t270 + pkin(6) * t267;
t251 = t282 * t268;
t252 = t282 * t271;
t284 = qJD(2) * t271;
t285 = qJD(2) * t268;
t286 = t251 * t285 + t252 * t284;
t283 = qJD(3) * t267;
t281 = rSges(3,1) * t270 - rSges(3,2) * t267;
t280 = Icges(3,1) * t270 - t296;
t279 = -Icges(3,2) * t267 + t295;
t278 = Icges(3,5) * t270 - Icges(3,6) * t267;
t233 = -Icges(3,6) * t271 + t268 * t279;
t237 = -Icges(3,5) * t271 + t268 * t280;
t277 = t233 * t267 - t237 * t270;
t234 = Icges(3,6) * t268 + t271 * t279;
t238 = Icges(3,5) * t268 + t271 * t280;
t276 = -t234 * t267 + t238 * t270;
t257 = Icges(3,2) * t270 + t296;
t258 = Icges(3,1) * t267 + t295;
t275 = -t257 * t267 + t258 * t270;
t255 = qJD(1) * (pkin(1) * t271 + pkin(5) * t268);
t262 = pkin(2) * t267 - pkin(6) * t270;
t274 = qJD(1) * t252 - t262 * t285 + t255;
t263 = pkin(1) * t268 - pkin(5) * t271;
t273 = (-t251 - t263) * qJD(1) - t262 * t284;
t264 = -qJD(3) * t270 + qJD(1);
t261 = rSges(2,1) * t271 - rSges(2,2) * t268;
t260 = rSges(2,1) * t268 + rSges(2,2) * t271;
t259 = rSges(3,1) * t267 + rSges(3,2) * t270;
t256 = Icges(3,5) * t267 + Icges(3,6) * t270;
t254 = t268 * t283 - t284;
t253 = t271 * t283 + t285;
t242 = t268 * rSges(3,3) + t271 * t281;
t241 = -t271 * rSges(3,3) + t268 * t281;
t240 = -t270 * rSges(4,3) + (rSges(4,1) * t269 - rSges(4,2) * t266) * t267;
t230 = Icges(3,3) * t268 + t271 * t278;
t229 = -Icges(3,3) * t271 + t268 * t278;
t224 = rSges(4,1) * t249 - rSges(4,2) * t248 + rSges(4,3) * t293;
t222 = rSges(4,1) * t247 - rSges(4,2) * t246 + rSges(4,3) * t294;
t208 = qJD(1) * t242 - t259 * t285 + t255;
t207 = -t259 * t284 + (-t241 - t263) * qJD(1);
t206 = (t241 * t268 + t242 * t271) * qJD(2);
t205 = t224 * t264 - t240 * t253 + t274;
t204 = -t222 * t264 + t240 * t254 + t273;
t203 = t222 * t253 - t224 * t254 + t286;
t202 = qJD(4) * t246 - t253 * t287 + t264 * t288 + t274;
t201 = qJD(4) * t248 + t254 * t287 - t264 * t289 + t273;
t200 = qJD(4) * t267 * t266 + t253 * t289 - t254 * t288 + t286;
t1 = m(3) * (t206 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + ((t268 * t256 + t271 * t275) * qJD(1) + (t268 ^ 2 * t230 + (t277 * t271 + (-t229 + t276) * t268) * t271) * qJD(2)) * t285 / 0.2e1 - ((-t256 * t271 + t268 * t275) * qJD(1) + (t271 ^ 2 * t229 + (t276 * t268 + (-t230 + t277) * t271) * t268) * qJD(2)) * t284 / 0.2e1 + qJD(1) * ((t270 * t257 + t267 * t258) * qJD(1) + ((t270 * t234 + t267 * t238) * t268 - (t233 * t270 + t237 * t267) * t271) * qJD(2)) / 0.2e1 + m(4) * (t203 ^ 2 + t204 ^ 2 + t205 ^ 2) / 0.2e1 + m(5) * (t200 ^ 2 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + ((t248 * t303 + t249 * t301 + t293 * t302) * t264 + (t248 * t309 + t249 * t305 + t293 * t307) * t254 + (t308 * t248 + t304 * t249 + t306 * t293) * t253) * t253 / 0.2e1 + ((t246 * t303 + t247 * t301 + t294 * t302) * t264 + (t309 * t246 + t305 * t247 + t307 * t294) * t254 + (t246 * t308 + t247 * t304 + t294 * t306) * t253) * t254 / 0.2e1 + ((-t253 * t306 - t254 * t307 - t264 * t302) * t270 + ((t266 * t303 + t269 * t301) * t264 + (t266 * t309 + t269 * t305) * t254 + (t266 * t308 + t269 * t304) * t253) * t267) * t264 / 0.2e1 + (m(2) * (t260 ^ 2 + t261 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
