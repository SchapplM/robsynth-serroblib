% Calculate kinetic energy for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:04
% EndTime: 2019-03-09 09:17:08
% DurationCPUTime: 3.62s
% Computational Cost: add. (1876->344), mult. (4163->487), div. (0->0), fcn. (4711->10), ass. (0->156)
t338 = Icges(3,1) + Icges(4,1) + Icges(5,2);
t337 = Icges(5,1) + Icges(3,2) + Icges(4,3);
t336 = Icges(3,4) + Icges(5,4) - Icges(4,5);
t335 = Icges(4,6) - Icges(3,6) - Icges(5,5);
t334 = Icges(5,6) + Icges(4,4) + Icges(3,5);
t333 = -Icges(5,3) - Icges(3,3) - Icges(4,2);
t281 = cos(pkin(6));
t287 = cos(qJ(2));
t288 = cos(qJ(1));
t310 = t287 * t288;
t284 = sin(qJ(2));
t285 = sin(qJ(1));
t313 = t284 * t285;
t244 = -t281 * t310 + t313;
t311 = t285 * t287;
t312 = t284 * t288;
t245 = t281 * t312 + t311;
t280 = sin(pkin(6));
t314 = t280 * t288;
t332 = -t335 * t244 - t334 * t245 - t333 * t314;
t246 = t281 * t311 + t312;
t247 = -t281 * t313 + t310;
t316 = t280 * t285;
t331 = -t335 * t246 - t334 * t247 + t333 * t316;
t330 = t337 * t244 - t336 * t245 - t335 * t314;
t329 = t337 * t246 - t336 * t247 + t335 * t316;
t328 = -t336 * t244 + t338 * t245 - t334 * t314;
t327 = -t336 * t246 + t338 * t247 + t334 * t316;
t326 = t333 * t281 + (-t334 * t284 + t335 * t287) * t280;
t325 = t335 * t281 + (-t336 * t284 - t337 * t287) * t280;
t324 = t334 * t281 + (t338 * t284 + t336 * t287) * t280;
t320 = cos(qJ(5));
t319 = pkin(8) * t281;
t318 = Icges(2,4) * t285;
t317 = t280 * t284;
t315 = t280 * t287;
t204 = pkin(2) * t245 + qJ(3) * t244;
t221 = t245 * pkin(3) + qJ(4) * t314;
t309 = -t204 - t221;
t205 = pkin(2) * t247 + qJ(3) * t246;
t222 = pkin(3) * t247 - qJ(4) * t316;
t308 = -t205 - t222;
t248 = (pkin(2) * t284 - qJ(3) * t287) * t280;
t251 = pkin(3) * t317 - qJ(4) * t281;
t307 = -t248 - t251;
t306 = qJD(2) * t280;
t305 = qJD(4) * t280;
t304 = V_base(5) * pkin(7) + V_base(1);
t301 = t280 * t320;
t258 = t285 * t306 + V_base(4);
t277 = V_base(6) + qJD(1);
t211 = qJD(5) * t247 + t258;
t259 = qJD(2) * t281 + t277;
t240 = qJD(5) * t317 + t259;
t257 = -t288 * t306 + V_base(5);
t252 = t285 * pkin(1) - pkin(8) * t314;
t300 = -t252 * t277 + V_base(5) * t319 + t304;
t253 = pkin(1) * t288 + pkin(8) * t316;
t299 = V_base(4) * t252 - t253 * V_base(5) + V_base(3);
t210 = qJD(5) * t245 + t257;
t298 = qJD(3) * t246 + t257 * t248 + t300;
t297 = t277 * t253 + V_base(2) + (-pkin(7) - t319) * V_base(4);
t296 = qJD(3) * t244 + t259 * t205 + t297;
t295 = -qJD(3) * t315 + t258 * t204 + t299;
t294 = t259 * t222 + t288 * t305 + t296;
t293 = t257 * t251 - t285 * t305 + t298;
t292 = -qJD(4) * t281 + t258 * t221 + t295;
t207 = pkin(4) * t246 + pkin(9) * t247;
t250 = (-pkin(4) * t287 + pkin(9) * t284) * t280;
t291 = t259 * t207 + (-t250 + t307) * t258 + t294;
t206 = pkin(4) * t244 + pkin(9) * t245;
t290 = t257 * t250 + (-t206 + t309) * t259 + t293;
t289 = t258 * t206 + (-t207 + t308) * t257 + t292;
t286 = cos(qJ(6));
t283 = sin(qJ(5));
t282 = sin(qJ(6));
t278 = Icges(2,4) * t288;
t267 = rSges(2,1) * t288 - t285 * rSges(2,2);
t266 = t285 * rSges(2,1) + rSges(2,2) * t288;
t265 = Icges(2,1) * t288 - t318;
t264 = Icges(2,1) * t285 + t278;
t263 = -Icges(2,2) * t285 + t278;
t262 = Icges(2,2) * t288 + t318;
t256 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t255 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t254 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t243 = -t281 * t283 - t287 * t301;
t242 = -t281 * t320 + t283 * t315;
t235 = -rSges(5,3) * t281 + (-rSges(5,1) * t287 - rSges(5,2) * t284) * t280;
t234 = rSges(3,3) * t281 + (rSges(3,1) * t284 + rSges(3,2) * t287) * t280;
t233 = rSges(4,2) * t281 + (rSges(4,1) * t284 - rSges(4,3) * t287) * t280;
t220 = V_base(5) * rSges(2,3) - t266 * t277 + t304;
t219 = t267 * t277 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t216 = t266 * V_base(4) - t267 * V_base(5) + V_base(3);
t215 = t246 * t320 - t283 * t316;
t214 = t246 * t283 + t285 * t301;
t213 = t244 * t320 + t283 * t314;
t212 = t244 * t283 - t288 * t301;
t209 = t243 * t286 + t282 * t317;
t208 = -t243 * t282 + t286 * t317;
t203 = -qJD(6) * t242 + t240;
t202 = pkin(5) * t243 - pkin(10) * t242;
t199 = rSges(3,1) * t247 - rSges(3,2) * t246 + rSges(3,3) * t316;
t198 = rSges(4,1) * t247 + rSges(4,2) * t316 + rSges(4,3) * t246;
t197 = rSges(5,1) * t246 - rSges(5,2) * t247 - rSges(5,3) * t316;
t196 = t245 * rSges(3,1) - t244 * rSges(3,2) - rSges(3,3) * t314;
t195 = t245 * rSges(4,1) - rSges(4,2) * t314 + t244 * rSges(4,3);
t194 = t244 * rSges(5,1) - t245 * rSges(5,2) + rSges(5,3) * t314;
t175 = rSges(6,1) * t243 + rSges(6,2) * t242 + rSges(6,3) * t317;
t174 = Icges(6,1) * t243 + Icges(6,4) * t242 + Icges(6,5) * t317;
t173 = Icges(6,4) * t243 + Icges(6,2) * t242 + Icges(6,6) * t317;
t172 = Icges(6,5) * t243 + Icges(6,6) * t242 + Icges(6,3) * t317;
t167 = t215 * t286 + t247 * t282;
t166 = -t215 * t282 + t247 * t286;
t165 = t213 * t286 + t245 * t282;
t164 = -t213 * t282 + t245 * t286;
t163 = qJD(6) * t214 + t211;
t162 = qJD(6) * t212 + t210;
t161 = pkin(5) * t215 + pkin(10) * t214;
t160 = pkin(5) * t213 + pkin(10) * t212;
t159 = rSges(6,1) * t215 - rSges(6,2) * t214 + rSges(6,3) * t247;
t158 = rSges(6,1) * t213 - rSges(6,2) * t212 + rSges(6,3) * t245;
t157 = Icges(6,1) * t215 - Icges(6,4) * t214 + Icges(6,5) * t247;
t156 = Icges(6,1) * t213 - Icges(6,4) * t212 + Icges(6,5) * t245;
t155 = Icges(6,4) * t215 - Icges(6,2) * t214 + Icges(6,6) * t247;
t154 = Icges(6,4) * t213 - Icges(6,2) * t212 + Icges(6,6) * t245;
t153 = Icges(6,5) * t215 - Icges(6,6) * t214 + Icges(6,3) * t247;
t152 = Icges(6,5) * t213 - Icges(6,6) * t212 + Icges(6,3) * t245;
t151 = rSges(7,1) * t209 + rSges(7,2) * t208 - rSges(7,3) * t242;
t150 = Icges(7,1) * t209 + Icges(7,4) * t208 - Icges(7,5) * t242;
t149 = Icges(7,4) * t209 + Icges(7,2) * t208 - Icges(7,6) * t242;
t148 = Icges(7,5) * t209 + Icges(7,6) * t208 - Icges(7,3) * t242;
t147 = -t196 * t259 + t234 * t257 + t300;
t146 = t199 * t259 - t234 * t258 + t297;
t145 = rSges(7,1) * t167 + rSges(7,2) * t166 + rSges(7,3) * t214;
t144 = rSges(7,1) * t165 + rSges(7,2) * t164 + rSges(7,3) * t212;
t143 = Icges(7,1) * t167 + Icges(7,4) * t166 + Icges(7,5) * t214;
t142 = Icges(7,1) * t165 + Icges(7,4) * t164 + Icges(7,5) * t212;
t141 = Icges(7,4) * t167 + Icges(7,2) * t166 + Icges(7,6) * t214;
t140 = Icges(7,4) * t165 + Icges(7,2) * t164 + Icges(7,6) * t212;
t139 = Icges(7,5) * t167 + Icges(7,6) * t166 + Icges(7,3) * t214;
t138 = Icges(7,5) * t165 + Icges(7,6) * t164 + Icges(7,3) * t212;
t137 = t196 * t258 - t199 * t257 + t299;
t136 = t233 * t257 + (-t195 - t204) * t259 + t298;
t135 = t198 * t259 + (-t233 - t248) * t258 + t296;
t134 = t195 * t258 + (-t198 - t205) * t257 + t295;
t133 = t235 * t257 + (-t194 + t309) * t259 + t293;
t132 = t197 * t259 + (-t235 + t307) * t258 + t294;
t131 = t194 * t258 + (-t197 + t308) * t257 + t292;
t130 = -t158 * t240 + t175 * t210 + t290;
t129 = t159 * t240 - t175 * t211 + t291;
t128 = t158 * t211 - t159 * t210 + t289;
t127 = -t144 * t203 + t151 * t162 - t160 * t240 + t202 * t210 + t290;
t126 = t145 * t203 - t151 * t163 + t161 * t240 - t202 * t211 + t291;
t125 = t144 * t163 - t145 * t162 + t160 * t211 - t161 * t210 + t289;
t1 = t211 * ((t247 * t153 - t214 * t155 + t215 * t157) * t211 + (t152 * t247 - t154 * t214 + t156 * t215) * t210 + (t172 * t247 - t173 * t214 + t174 * t215) * t240) / 0.2e1 + m(1) * (t254 ^ 2 + t255 ^ 2 + t256 ^ 2) / 0.2e1 + t210 * ((t153 * t245 - t155 * t212 + t157 * t213) * t211 + (t245 * t152 - t212 * t154 + t213 * t156) * t210 + (t172 * t245 - t173 * t212 + t174 * t213) * t240) / 0.2e1 + t203 * ((-t139 * t242 + t141 * t208 + t143 * t209) * t163 + (-t138 * t242 + t140 * t208 + t142 * t209) * t162 + (-t242 * t148 + t208 * t149 + t209 * t150) * t203) / 0.2e1 + m(2) * (t216 ^ 2 + t219 ^ 2 + t220 ^ 2) / 0.2e1 + t163 * ((t214 * t139 + t166 * t141 + t167 * t143) * t163 + (t138 * t214 + t140 * t166 + t142 * t167) * t162 + (t148 * t214 + t149 * t166 + t150 * t167) * t203) / 0.2e1 + t162 * ((t139 * t212 + t141 * t164 + t143 * t165) * t163 + (t212 * t138 + t164 * t140 + t165 * t142) * t162 + (t148 * t212 + t149 * t164 + t150 * t165) * t203) / 0.2e1 + m(3) * (t137 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(4) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(6) * (t128 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(5) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(7) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + t240 * ((t153 * t317 + t155 * t242 + t157 * t243) * t211 + (t152 * t317 + t154 * t242 + t156 * t243) * t210 + (t172 * t317 + t242 * t173 + t243 * t174) * t240) / 0.2e1 + ((-t285 * t262 + t264 * t288 + Icges(1,4)) * V_base(5) + (-t285 * t263 + t288 * t265 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t288 * t262 + t285 * t264 + Icges(1,2)) * V_base(5) + (t263 * t288 + t285 * t265 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t244 * t325 + t245 * t324 + t314 * t326) * t259 + (t244 * t329 + t245 * t327 + t314 * t331) * t258 + (t330 * t244 + t328 * t245 + t332 * t314) * t257) * t257 / 0.2e1 + ((t246 * t325 + t247 * t324 - t316 * t326) * t259 + (t329 * t246 + t327 * t247 - t331 * t316) * t258 + (t330 * t246 + t328 * t247 - t316 * t332) * t257) * t258 / 0.2e1 + ((-t257 * t332 - t331 * t258 - t326 * t259) * t281 + ((t284 * t324 - t287 * t325) * t259 + (t284 * t327 - t287 * t329) * t258 + (t284 * t328 - t287 * t330) * t257) * t280) * t259 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t285 + Icges(2,6) * t288) * V_base(5) + (Icges(2,5) * t288 - Icges(2,6) * t285) * V_base(4) + Icges(2,3) * t277 / 0.2e1) * t277;
T  = t1;
