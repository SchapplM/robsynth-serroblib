% Calculate kinetic energy for
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP6_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP6_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:08
% EndTime: 2019-12-31 17:18:09
% DurationCPUTime: 1.15s
% Computational Cost: add. (515->157), mult. (1202->255), div. (0->0), fcn. (1201->6), ass. (0->88)
t327 = Icges(4,1) + Icges(5,1);
t326 = Icges(4,4) + Icges(5,4);
t325 = -Icges(5,5) - Icges(4,5);
t324 = Icges(4,2) + Icges(5,2);
t323 = Icges(4,6) + Icges(5,6);
t322 = -Icges(5,3) - Icges(4,3);
t273 = sin(qJ(3));
t276 = cos(qJ(3));
t278 = cos(qJ(1));
t275 = sin(qJ(1));
t277 = cos(qJ(2));
t300 = t277 * t275;
t253 = -t273 * t300 - t276 * t278;
t303 = t273 * t278;
t254 = t276 * t300 - t303;
t274 = sin(qJ(2));
t302 = t274 * t275;
t321 = t323 * t253 - t325 * t254 - t322 * t302;
t299 = t277 * t278;
t255 = -t273 * t299 + t275 * t276;
t304 = t273 * t275;
t256 = t276 * t299 + t304;
t301 = t274 * t278;
t320 = t323 * t255 - t325 * t256 - t322 * t301;
t319 = t324 * t253 + t326 * t254 + t323 * t302;
t318 = t324 * t255 + t326 * t256 + t323 * t301;
t317 = t326 * t253 + t327 * t254 - t325 * t302;
t316 = t326 * t255 + t327 * t256 - t325 * t301;
t315 = t322 * t277 + (-t323 * t273 - t325 * t276) * t274;
t314 = -t323 * t277 + (-t324 * t273 + t326 * t276) * t274;
t313 = t325 * t277 + (-t326 * t273 + t327 * t276) * t274;
t308 = pkin(3) * t276;
t306 = Icges(3,4) * t274;
t305 = Icges(3,4) * t277;
t280 = qJ(4) * t274 + t277 * t308;
t298 = rSges(5,1) * t254 + rSges(5,2) * t253 + rSges(5,3) * t302 - pkin(3) * t303 + t275 * t280;
t297 = rSges(5,1) * t256 + rSges(5,2) * t255 + rSges(5,3) * t301 + pkin(3) * t304 + t278 * t280;
t296 = (-qJ(4) - rSges(5,3)) * t277 + (rSges(5,1) * t276 - rSges(5,2) * t273 + t308) * t274;
t289 = pkin(2) * t277 + pkin(6) * t274;
t257 = t289 * t275;
t258 = t289 * t278;
t292 = qJD(2) * t278;
t293 = qJD(2) * t275;
t295 = t257 * t293 + t258 * t292;
t261 = qJD(1) * (pkin(1) * t278 + pkin(5) * t275);
t294 = qJD(1) * t258 + t261;
t291 = qJD(3) * t274;
t269 = pkin(1) * t275 - pkin(5) * t278;
t290 = (-t257 - t269) * qJD(1);
t288 = rSges(3,1) * t277 - rSges(3,2) * t274;
t287 = Icges(3,1) * t277 - t306;
t286 = -Icges(3,2) * t274 + t305;
t285 = Icges(3,5) * t277 - Icges(3,6) * t274;
t240 = -Icges(3,6) * t278 + t275 * t286;
t244 = -Icges(3,5) * t278 + t275 * t287;
t284 = t240 * t274 - t244 * t277;
t241 = Icges(3,6) * t275 + t278 * t286;
t245 = Icges(3,5) * t275 + t278 * t287;
t283 = -t241 * t274 + t245 * t277;
t263 = Icges(3,2) * t277 + t306;
t264 = Icges(3,1) * t274 + t305;
t282 = -t263 * t274 + t264 * t277;
t268 = pkin(2) * t274 - pkin(6) * t277;
t281 = -qJD(2) * t268 + qJD(4) * t274;
t270 = -qJD(3) * t277 + qJD(1);
t267 = rSges(2,1) * t278 - rSges(2,2) * t275;
t266 = rSges(2,1) * t275 + rSges(2,2) * t278;
t265 = rSges(3,1) * t274 + rSges(3,2) * t277;
t262 = Icges(3,5) * t274 + Icges(3,6) * t277;
t260 = t275 * t291 - t292;
t259 = t278 * t291 + t293;
t249 = rSges(3,3) * t275 + t278 * t288;
t248 = -rSges(3,3) * t278 + t275 * t288;
t247 = -rSges(4,3) * t277 + (rSges(4,1) * t276 - rSges(4,2) * t273) * t274;
t237 = Icges(3,3) * t275 + t278 * t285;
t236 = -Icges(3,3) * t278 + t275 * t285;
t232 = rSges(4,1) * t256 + rSges(4,2) * t255 + rSges(4,3) * t301;
t230 = rSges(4,1) * t254 + rSges(4,2) * t253 + rSges(4,3) * t302;
t214 = qJD(1) * t249 - t265 * t293 + t261;
t213 = -t265 * t292 + (-t248 - t269) * qJD(1);
t212 = (t248 * t275 + t249 * t278) * qJD(2);
t211 = t232 * t270 - t247 * t259 - t268 * t293 + t294;
t210 = -t230 * t270 + t247 * t260 - t268 * t292 + t290;
t209 = t230 * t259 - t232 * t260 + t295;
t208 = -t259 * t296 + t270 * t297 + t275 * t281 + t294;
t207 = t260 * t296 - t270 * t298 + t278 * t281 + t290;
t206 = -qJD(4) * t277 + t259 * t298 - t260 * t297 + t295;
t1 = m(3) * (t212 ^ 2 + t213 ^ 2 + t214 ^ 2) / 0.2e1 + ((t275 * t262 + t278 * t282) * qJD(1) + (t275 ^ 2 * t237 + (t284 * t278 + (-t236 + t283) * t275) * t278) * qJD(2)) * t293 / 0.2e1 - ((-t278 * t262 + t275 * t282) * qJD(1) + (t278 ^ 2 * t236 + (t283 * t275 + (-t237 + t284) * t278) * t275) * qJD(2)) * t292 / 0.2e1 + qJD(1) * ((t277 * t263 + t274 * t264) * qJD(1) + ((t241 * t277 + t245 * t274) * t275 - (t240 * t277 + t244 * t274) * t278) * qJD(2)) / 0.2e1 + m(4) * (t209 ^ 2 + t210 ^ 2 + t211 ^ 2) / 0.2e1 + m(5) * (t206 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + ((t314 * t255 + t313 * t256 + t315 * t301) * t270 + (t319 * t255 + t317 * t256 + t321 * t301) * t260 + (t318 * t255 + t316 * t256 + t320 * t301) * t259) * t259 / 0.2e1 + ((t314 * t253 + t313 * t254 + t315 * t302) * t270 + (t319 * t253 + t317 * t254 + t321 * t302) * t260 + (t318 * t253 + t316 * t254 + t320 * t302) * t259) * t260 / 0.2e1 + ((-t320 * t259 - t321 * t260 - t315 * t270) * t277 + ((-t314 * t273 + t313 * t276) * t270 + (-t319 * t273 + t317 * t276) * t260 + (-t318 * t273 + t316 * t276) * t259) * t274) * t270 / 0.2e1 + (m(2) * (t266 ^ 2 + t267 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
