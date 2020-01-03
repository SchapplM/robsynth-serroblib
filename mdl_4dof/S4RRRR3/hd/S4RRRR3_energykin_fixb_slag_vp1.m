% Calculate kinetic energy for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:15
% EndTime: 2019-12-31 17:24:16
% DurationCPUTime: 0.80s
% Computational Cost: add. (630->157), mult. (759->262), div. (0->0), fcn. (646->8), ass. (0->98)
t258 = qJ(2) + qJ(3);
t253 = sin(t258);
t302 = pkin(3) * t253;
t261 = cos(qJ(2));
t300 = t261 * pkin(2);
t259 = sin(qJ(2));
t298 = Icges(3,4) * t259;
t297 = Icges(3,4) * t261;
t296 = Icges(4,4) * t253;
t254 = cos(t258);
t295 = Icges(4,4) * t254;
t255 = qJ(4) + t258;
t248 = sin(t255);
t294 = Icges(5,4) * t248;
t249 = cos(t255);
t293 = Icges(5,4) * t249;
t260 = sin(qJ(1));
t262 = cos(qJ(1));
t209 = -pkin(6) * t262 + t260 * t300;
t210 = pkin(6) * t260 + t262 * t300;
t252 = qJD(2) * t260;
t288 = qJD(2) * t262;
t292 = t209 * t252 + t210 * t288;
t247 = pkin(1) * t260 - pkin(5) * t262;
t291 = -t209 - t247;
t290 = pkin(3) * t254;
t239 = qJD(3) * t260 + t252;
t287 = -qJD(2) - qJD(3);
t286 = pkin(2) * qJD(2) * t259;
t285 = t262 * t286;
t284 = rSges(3,1) * t261 - rSges(3,2) * t259;
t283 = rSges(4,1) * t254 - rSges(4,2) * t253;
t282 = rSges(5,1) * t249 - rSges(5,2) * t248;
t281 = Icges(3,1) * t261 - t298;
t280 = Icges(4,1) * t254 - t296;
t279 = Icges(5,1) * t249 - t294;
t278 = -Icges(3,2) * t259 + t297;
t277 = -Icges(4,2) * t253 + t295;
t276 = -Icges(5,2) * t248 + t293;
t275 = Icges(3,5) * t261 - Icges(3,6) * t259;
t274 = Icges(4,5) * t254 - Icges(4,6) * t253;
t273 = Icges(5,5) * t249 - Icges(5,6) * t248;
t221 = -Icges(3,6) * t262 + t260 * t278;
t223 = -Icges(3,5) * t262 + t260 * t281;
t272 = t221 * t259 - t223 * t261;
t222 = Icges(3,6) * t260 + t262 * t278;
t224 = Icges(3,5) * t260 + t262 * t281;
t271 = -t222 * t259 + t224 * t261;
t242 = Icges(3,2) * t261 + t298;
t243 = Icges(3,1) * t259 + t297;
t270 = -t242 * t259 + t243 * t261;
t238 = qJD(1) * (pkin(1) * t262 + pkin(5) * t260);
t269 = qJD(1) * t210 - t260 * t286 + t238;
t231 = qJD(4) * t260 + t239;
t232 = (-qJD(4) + t287) * t262;
t268 = (Icges(5,5) * t248 + Icges(5,6) * t249) * qJD(1) + (-Icges(5,3) * t262 + t260 * t273) * t232 + (Icges(5,3) * t260 + t262 * t273) * t231;
t240 = t287 * t262;
t267 = (Icges(4,5) * t253 + Icges(4,6) * t254) * qJD(1) + (-Icges(4,3) * t262 + t260 * t274) * t240 + (Icges(4,3) * t260 + t262 * t274) * t239;
t203 = -Icges(5,6) * t262 + t260 * t276;
t204 = Icges(5,6) * t260 + t262 * t276;
t205 = -Icges(5,5) * t262 + t260 * t279;
t206 = Icges(5,5) * t260 + t262 * t279;
t228 = Icges(5,2) * t249 + t294;
t229 = Icges(5,1) * t248 + t293;
t266 = (-t204 * t248 + t206 * t249) * t231 + (-t203 * t248 + t205 * t249) * t232 + (-t228 * t248 + t229 * t249) * qJD(1);
t213 = -Icges(4,6) * t262 + t260 * t277;
t214 = Icges(4,6) * t260 + t262 * t277;
t215 = -Icges(4,5) * t262 + t260 * t280;
t216 = Icges(4,5) * t260 + t262 * t280;
t234 = Icges(4,2) * t254 + t296;
t235 = Icges(4,1) * t253 + t295;
t265 = (-t214 * t253 + t216 * t254) * t239 + (-t213 * t253 + t215 * t254) * t240 + (-t234 * t253 + t235 * t254) * qJD(1);
t246 = rSges(2,1) * t262 - rSges(2,2) * t260;
t245 = rSges(2,1) * t260 + rSges(2,2) * t262;
t244 = rSges(3,1) * t259 + rSges(3,2) * t261;
t241 = Icges(3,5) * t259 + Icges(3,6) * t261;
t236 = rSges(4,1) * t253 + rSges(4,2) * t254;
t230 = rSges(5,1) * t248 + rSges(5,2) * t249;
t226 = rSges(3,3) * t260 + t262 * t284;
t225 = -rSges(3,3) * t262 + t260 * t284;
t220 = Icges(3,3) * t260 + t262 * t275;
t219 = -Icges(3,3) * t262 + t260 * t275;
t218 = rSges(4,3) * t260 + t262 * t283;
t217 = -rSges(4,3) * t262 + t260 * t283;
t208 = rSges(5,3) * t260 + t262 * t282;
t207 = -rSges(5,3) * t262 + t260 * t282;
t197 = pkin(7) * t260 + t262 * t290;
t196 = -pkin(7) * t262 + t260 * t290;
t195 = qJD(1) * t226 - t244 * t252 + t238;
t194 = -t244 * t288 + (-t225 - t247) * qJD(1);
t193 = (t225 * t260 + t226 * t262) * qJD(2);
t192 = qJD(1) * t218 - t236 * t239 + t269;
t191 = -t285 + t236 * t240 + (-t217 + t291) * qJD(1);
t190 = t217 * t239 - t218 * t240 + t292;
t189 = -t239 * t302 - t230 * t231 + (t197 + t208) * qJD(1) + t269;
t188 = -t285 + t240 * t302 + t230 * t232 + (-t196 - t207 + t291) * qJD(1);
t187 = t196 * t239 - t197 * t240 + t207 * t231 - t208 * t232 + t292;
t1 = m(3) * (t193 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + ((t260 * t241 + t262 * t270) * qJD(1) + (t260 ^ 2 * t220 + (t272 * t262 + (-t219 + t271) * t260) * t262) * qJD(2)) * t252 / 0.2e1 - ((-t262 * t241 + t260 * t270) * qJD(1) + (t262 ^ 2 * t219 + (t271 * t260 + (-t220 + t272) * t262) * t260) * qJD(2)) * t288 / 0.2e1 + m(4) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + t239 * (t267 * t260 + t265 * t262) / 0.2e1 + t240 * (t265 * t260 - t267 * t262) / 0.2e1 + m(5) * (t187 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + t231 * (t268 * t260 + t266 * t262) / 0.2e1 + t232 * (t266 * t260 - t268 * t262) / 0.2e1 + (m(2) * (t245 ^ 2 + t246 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t222 * t261 + t224 * t259) * t260 - (t221 * t261 + t223 * t259) * t262) * qJD(2) + (t214 * t254 + t216 * t253) * t239 + (t213 * t254 + t215 * t253) * t240 + (t204 * t249 + t206 * t248) * t231 + (t203 * t249 + t205 * t248) * t232 + (t249 * t228 + t248 * t229 + t254 * t234 + t253 * t235 + t261 * t242 + t259 * t243) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
