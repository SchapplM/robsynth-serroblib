% Calculate kinetic energy for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:25:00
% EndTime: 2019-03-09 18:25:05
% DurationCPUTime: 4.51s
% Computational Cost: add. (2481->404), mult. (2801->602), div. (0->0), fcn. (2751->12), ass. (0->182)
t337 = -Icges(5,3) - Icges(4,3);
t276 = qJ(3) + pkin(11);
t265 = sin(t276);
t266 = cos(t276);
t283 = cos(qJ(1));
t280 = sin(qJ(1));
t282 = cos(qJ(2));
t320 = t280 * t282;
t197 = -t265 * t320 - t266 * t283;
t198 = -t265 * t283 + t266 * t320;
t278 = sin(qJ(3));
t281 = cos(qJ(3));
t219 = -t278 * t320 - t281 * t283;
t324 = t278 * t283;
t220 = t281 * t320 - t324;
t279 = sin(qJ(2));
t323 = t279 * t280;
t336 = Icges(4,5) * t220 + Icges(5,5) * t198 + Icges(4,6) * t219 + Icges(5,6) * t197 - t337 * t323;
t319 = t282 * t283;
t199 = -t265 * t319 + t280 * t266;
t200 = t280 * t265 + t266 * t319;
t221 = -t278 * t319 + t280 * t281;
t321 = t280 * t278;
t222 = t281 * t319 + t321;
t322 = t279 * t283;
t335 = Icges(4,5) * t222 + Icges(5,5) * t200 + Icges(4,6) * t221 + Icges(5,6) * t199 - t337 * t322;
t334 = t337 * t282 + (Icges(4,5) * t281 + Icges(5,5) * t266 - Icges(4,6) * t278 - Icges(5,6) * t265) * t279;
t329 = t281 * pkin(3);
t327 = Icges(2,4) * t280;
t326 = Icges(3,4) * t279;
t325 = Icges(3,4) * t282;
t267 = qJ(5) + t276;
t262 = cos(t267);
t318 = pkin(5) * t262;
t261 = sin(t267);
t317 = pkin(5) * t261;
t316 = pkin(4) * t266;
t313 = qJD(3) * t279;
t312 = qJD(4) * t279;
t311 = qJD(5) * t279;
t310 = qJD(6) * t279;
t309 = -qJD(3) - qJD(5);
t308 = V_base(5) * pkin(6) + V_base(1);
t250 = qJD(2) * t280 + V_base(4);
t268 = V_base(6) + qJD(1);
t305 = pkin(4) * t265;
t218 = t283 * t313 + t250;
t304 = pkin(2) * t282 + pkin(8) * t279;
t249 = -qJD(2) * t283 + V_base(5);
t303 = rSges(3,1) * t282 - rSges(3,2) * t279;
t191 = t283 * t311 + t218;
t302 = Icges(3,1) * t282 - t326;
t301 = -Icges(3,2) * t279 + t325;
t300 = Icges(3,5) * t282 - Icges(3,6) * t279;
t217 = t280 * t313 + t249;
t248 = pkin(1) * t283 + t280 * pkin(7);
t299 = -V_base(4) * pkin(6) + t268 * t248 + V_base(2);
t247 = t280 * pkin(1) - pkin(7) * t283;
t298 = V_base(4) * t247 - t248 * V_base(5) + V_base(3);
t190 = t280 * t311 + t217;
t297 = qJ(4) * t279 + t282 * t329;
t224 = t304 * t280;
t246 = pkin(2) * t279 - pkin(8) * t282;
t296 = t249 * t246 + (-t224 - t247) * t268 + t308;
t295 = (-Icges(3,3) * t283 + t280 * t300) * t249 + (Icges(3,3) * t280 + t283 * t300) * t250 + (Icges(3,5) * t279 + Icges(3,6) * t282) * t268;
t294 = pkin(10) * t279 + t282 * t318;
t293 = pkin(9) * t279 + t282 * t316;
t225 = t304 * t283;
t292 = t268 * t225 - t246 * t250 + t299;
t184 = -qJ(4) * t282 + t279 * t329;
t291 = t217 * t184 + t283 * t312 + t296;
t290 = t250 * t224 - t225 * t249 + t298;
t162 = pkin(3) * t321 + t283 * t297;
t242 = -qJD(3) * t282 + t268;
t289 = t242 * t162 + t280 * t312 + t292;
t161 = -pkin(3) * t324 + t280 * t297;
t288 = -qJD(4) * t282 + t218 * t161 + t290;
t123 = t280 * t293 - t283 * t305;
t166 = -pkin(9) * t282 + t279 * t316;
t287 = t217 * t166 + (-t123 - t161) * t242 + t291;
t124 = t280 * t305 + t283 * t293;
t286 = t242 * t124 + (-t166 - t184) * t218 + t289;
t285 = t218 * t123 + (-t124 - t162) * t217 + t288;
t206 = -Icges(3,6) * t283 + t280 * t301;
t207 = Icges(3,6) * t280 + t283 * t301;
t209 = -Icges(3,5) * t283 + t280 * t302;
t210 = Icges(3,5) * t280 + t283 * t302;
t236 = Icges(3,2) * t282 + t326;
t239 = Icges(3,1) * t279 + t325;
t284 = (-t207 * t279 + t210 * t282) * t250 + (-t206 * t279 + t209 * t282) * t249 + (-t236 * t279 + t239 * t282) * t268;
t271 = Icges(2,4) * t283;
t263 = qJ(6) + t267;
t252 = cos(t263);
t251 = sin(t263);
t245 = rSges(2,1) * t283 - t280 * rSges(2,2);
t244 = t280 * rSges(2,1) + rSges(2,2) * t283;
t243 = rSges(3,1) * t279 + rSges(3,2) * t282;
t241 = Icges(2,1) * t283 - t327;
t240 = Icges(2,1) * t280 + t271;
t238 = -Icges(2,2) * t280 + t271;
t237 = Icges(2,2) * t283 + t327;
t231 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t230 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t229 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t226 = t282 * t309 + t268;
t215 = (-qJD(6) + t309) * t282 + t268;
t213 = t280 * rSges(3,3) + t283 * t303;
t212 = -rSges(3,3) * t283 + t280 * t303;
t211 = -rSges(4,3) * t282 + (rSges(4,1) * t281 - rSges(4,2) * t278) * t279;
t208 = -Icges(4,5) * t282 + (Icges(4,1) * t281 - Icges(4,4) * t278) * t279;
t205 = -Icges(4,6) * t282 + (Icges(4,4) * t281 - Icges(4,2) * t278) * t279;
t195 = t280 * t261 + t262 * t319;
t194 = -t261 * t319 + t280 * t262;
t193 = -t261 * t283 + t262 * t320;
t192 = -t261 * t320 - t262 * t283;
t188 = -rSges(5,3) * t282 + (rSges(5,1) * t266 - rSges(5,2) * t265) * t279;
t187 = -Icges(5,5) * t282 + (Icges(5,1) * t266 - Icges(5,4) * t265) * t279;
t186 = -Icges(5,6) * t282 + (Icges(5,4) * t266 - Icges(5,2) * t265) * t279;
t183 = t280 * t251 + t252 * t319;
t182 = -t251 * t319 + t280 * t252;
t181 = -t251 * t283 + t252 * t320;
t180 = -t251 * t320 - t252 * t283;
t179 = -rSges(6,3) * t282 + (rSges(6,1) * t262 - rSges(6,2) * t261) * t279;
t178 = V_base(5) * rSges(2,3) - t244 * t268 + t308;
t177 = t245 * t268 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t176 = -Icges(6,5) * t282 + (Icges(6,1) * t262 - Icges(6,4) * t261) * t279;
t175 = -Icges(6,6) * t282 + (Icges(6,4) * t262 - Icges(6,2) * t261) * t279;
t174 = -Icges(6,3) * t282 + (Icges(6,5) * t262 - Icges(6,6) * t261) * t279;
t173 = t244 * V_base(4) - t245 * V_base(5) + V_base(3);
t172 = -rSges(7,3) * t282 + (rSges(7,1) * t252 - rSges(7,2) * t251) * t279;
t171 = -Icges(7,5) * t282 + (Icges(7,1) * t252 - Icges(7,4) * t251) * t279;
t170 = -Icges(7,6) * t282 + (Icges(7,4) * t252 - Icges(7,2) * t251) * t279;
t169 = -Icges(7,3) * t282 + (Icges(7,5) * t252 - Icges(7,6) * t251) * t279;
t168 = t283 * t310 + t191;
t167 = t280 * t310 + t190;
t164 = t222 * rSges(4,1) + t221 * rSges(4,2) + rSges(4,3) * t322;
t163 = rSges(4,1) * t220 + rSges(4,2) * t219 + rSges(4,3) * t323;
t160 = Icges(4,1) * t222 + Icges(4,4) * t221 + Icges(4,5) * t322;
t159 = Icges(4,1) * t220 + Icges(4,4) * t219 + Icges(4,5) * t323;
t158 = Icges(4,4) * t222 + Icges(4,2) * t221 + Icges(4,6) * t322;
t157 = Icges(4,4) * t220 + Icges(4,2) * t219 + Icges(4,6) * t323;
t154 = t200 * rSges(5,1) + t199 * rSges(5,2) + rSges(5,3) * t322;
t153 = rSges(5,1) * t198 + rSges(5,2) * t197 + rSges(5,3) * t323;
t152 = Icges(5,1) * t200 + Icges(5,4) * t199 + Icges(5,5) * t322;
t151 = Icges(5,1) * t198 + Icges(5,4) * t197 + Icges(5,5) * t323;
t150 = Icges(5,4) * t200 + Icges(5,2) * t199 + Icges(5,6) * t322;
t149 = Icges(5,4) * t198 + Icges(5,2) * t197 + Icges(5,6) * t323;
t145 = -pkin(10) * t282 + t279 * t318;
t143 = t195 * rSges(6,1) + t194 * rSges(6,2) + rSges(6,3) * t322;
t142 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t323;
t141 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t322;
t140 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t323;
t139 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t322;
t138 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t323;
t137 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t322;
t136 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t323;
t134 = t183 * rSges(7,1) + t182 * rSges(7,2) + rSges(7,3) * t322;
t133 = rSges(7,1) * t181 + rSges(7,2) * t180 + rSges(7,3) * t323;
t132 = Icges(7,1) * t183 + Icges(7,4) * t182 + Icges(7,5) * t322;
t131 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t323;
t130 = Icges(7,4) * t183 + Icges(7,2) * t182 + Icges(7,6) * t322;
t129 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t323;
t128 = Icges(7,5) * t183 + Icges(7,6) * t182 + Icges(7,3) * t322;
t127 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t323;
t126 = t243 * t249 + (-t212 - t247) * t268 + t308;
t125 = t213 * t268 - t243 * t250 + t299;
t121 = t212 * t250 - t213 * t249 + t298;
t119 = t280 * t317 + t283 * t294;
t118 = t280 * t294 - t283 * t317;
t117 = -t163 * t242 + t211 * t217 + t296;
t116 = t164 * t242 - t211 * t218 + t292;
t115 = t163 * t218 - t164 * t217 + t290;
t114 = t188 * t217 + (-t153 - t161) * t242 + t291;
t113 = t154 * t242 + (-t184 - t188) * t218 + t289;
t112 = t153 * t218 + (-t154 - t162) * t217 + t288;
t111 = -t142 * t226 + t179 * t190 + t287;
t110 = t143 * t226 - t179 * t191 + t286;
t109 = t142 * t191 - t143 * t190 + t285;
t108 = -t118 * t226 - t133 * t215 + t145 * t190 + t167 * t172 + t287;
t107 = t119 * t226 + t134 * t215 - t145 * t191 - t168 * t172 + t286;
t106 = t118 * t191 - t119 * t190 + t133 * t168 - t134 * t167 + t285;
t1 = t190 * ((t137 * t323 + t139 * t192 + t141 * t193) * t191 + (t136 * t323 + t138 * t192 + t140 * t193) * t190 + (t174 * t323 + t175 * t192 + t176 * t193) * t226) / 0.2e1 + t167 * ((t128 * t323 + t130 * t180 + t132 * t181) * t168 + (t127 * t323 + t180 * t129 + t181 * t131) * t167 + (t169 * t323 + t170 * t180 + t171 * t181) * t215) / 0.2e1 + t191 * ((t137 * t322 + t194 * t139 + t195 * t141) * t191 + (t136 * t322 + t194 * t138 + t195 * t140) * t190 + (t174 * t322 + t194 * t175 + t195 * t176) * t226) / 0.2e1 + t168 * ((t128 * t322 + t182 * t130 + t183 * t132) * t168 + (t127 * t322 + t182 * t129 + t183 * t131) * t167 + (t169 * t322 + t182 * t170 + t183 * t171) * t215) / 0.2e1 + t250 * (t295 * t280 + t284 * t283) / 0.2e1 + t249 * (t284 * t280 - t295 * t283) / 0.2e1 + m(1) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + m(2) * (t173 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(3) * (t121 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(7) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + t215 * ((-t127 * t167 - t128 * t168 - t169 * t215) * t282 + ((-t130 * t251 + t132 * t252) * t168 + (-t129 * t251 + t131 * t252) * t167 + (-t170 * t251 + t171 * t252) * t215) * t279) / 0.2e1 + t226 * ((-t136 * t190 - t137 * t191 - t174 * t226) * t282 + ((-t139 * t261 + t141 * t262) * t191 + (-t138 * t261 + t140 * t262) * t190 + (-t175 * t261 + t176 * t262) * t226) * t279) / 0.2e1 + ((t186 * t197 + t187 * t198 + t205 * t219 + t208 * t220 + t323 * t334) * t242 + (t150 * t197 + t152 * t198 + t158 * t219 + t160 * t220 + t323 * t335) * t218 + (t149 * t197 + t151 * t198 + t157 * t219 + t159 * t220 + t323 * t336) * t217) * t217 / 0.2e1 + ((t199 * t186 + t200 * t187 + t221 * t205 + t222 * t208 + t322 * t334) * t242 + (t199 * t150 + t200 * t152 + t221 * t158 + t222 * t160 + t322 * t335) * t218 + (t199 * t149 + t200 * t151 + t221 * t157 + t222 * t159 + t322 * t336) * t217) * t218 / 0.2e1 + ((-t217 * t336 - t335 * t218 - t334 * t242) * t282 + ((-t186 * t265 + t187 * t266 - t205 * t278 + t208 * t281) * t242 + (-t150 * t265 + t152 * t266 - t158 * t278 + t160 * t281) * t218 + (-t149 * t265 + t151 * t266 - t157 * t278 + t159 * t281) * t217) * t279) * t242 / 0.2e1 + ((t207 * t282 + t210 * t279) * t250 + (t206 * t282 + t209 * t279) * t249 + (t236 * t282 + t239 * t279 + Icges(2,3)) * t268) * t268 / 0.2e1 + ((-t280 * t237 + t240 * t283 + Icges(1,4)) * V_base(5) + (-t280 * t238 + t241 * t283 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t237 * t283 + t280 * t240 + Icges(1,2)) * V_base(5) + (t238 * t283 + t280 * t241 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t268 * (Icges(2,5) * t283 - Icges(2,6) * t280) + V_base(5) * t268 * (Icges(2,5) * t280 + Icges(2,6) * t283) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
