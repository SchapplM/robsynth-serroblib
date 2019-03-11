% Calculate kinetic energy for
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:18:33
% EndTime: 2019-03-09 22:18:38
% DurationCPUTime: 4.52s
% Computational Cost: add. (2508->404), mult. (2846->602), div. (0->0), fcn. (2796->12), ass. (0->182)
t337 = -Icges(6,3) - Icges(5,3);
t276 = qJ(3) + qJ(4);
t265 = pkin(11) + t276;
t259 = sin(t265);
t260 = cos(t265);
t282 = cos(qJ(1));
t279 = sin(qJ(1));
t281 = cos(qJ(2));
t320 = t279 * t281;
t192 = -t259 * t320 - t260 * t282;
t193 = -t259 * t282 + t260 * t320;
t268 = sin(t276);
t269 = cos(t276);
t210 = -t268 * t320 - t269 * t282;
t211 = -t268 * t282 + t269 * t320;
t278 = sin(qJ(2));
t322 = t278 * t279;
t336 = Icges(5,5) * t211 + Icges(6,5) * t193 + Icges(5,6) * t210 + Icges(6,6) * t192 - t322 * t337;
t319 = t281 * t282;
t194 = -t259 * t319 + t260 * t279;
t195 = t259 * t279 + t260 * t319;
t212 = -t268 * t319 + t269 * t279;
t213 = t268 * t279 + t269 * t319;
t321 = t278 * t282;
t335 = Icges(5,5) * t213 + Icges(6,5) * t195 + Icges(5,6) * t212 + Icges(6,6) * t194 - t321 * t337;
t334 = t337 * t281 + (Icges(5,5) * t269 + Icges(6,5) * t260 - Icges(5,6) * t268 - Icges(6,6) * t259) * t278;
t280 = cos(qJ(3));
t329 = t280 * pkin(3);
t327 = Icges(2,4) * t279;
t326 = Icges(3,4) * t278;
t325 = Icges(3,4) * t281;
t277 = sin(qJ(3));
t324 = t277 * t279;
t323 = t277 * t282;
t318 = pkin(5) * t260;
t317 = pkin(5) * t259;
t316 = pkin(4) * t269;
t313 = qJD(3) * t278;
t312 = qJD(4) * t278;
t311 = qJD(5) * t278;
t310 = qJD(6) * t278;
t309 = -qJD(3) - qJD(4);
t308 = V_base(5) * pkin(6) + V_base(1);
t250 = qJD(2) * t279 + V_base(4);
t266 = V_base(6) + qJD(1);
t305 = pkin(4) * t268;
t218 = t282 * t313 + t250;
t304 = pkin(2) * t281 + pkin(8) * t278;
t249 = -qJD(2) * t282 + V_base(5);
t303 = rSges(3,1) * t281 - rSges(3,2) * t278;
t191 = t282 * t312 + t218;
t302 = Icges(3,1) * t281 - t326;
t301 = -Icges(3,2) * t278 + t325;
t300 = Icges(3,5) * t281 - Icges(3,6) * t278;
t217 = t279 * t313 + t249;
t248 = pkin(1) * t282 + pkin(7) * t279;
t299 = -V_base(4) * pkin(6) + t266 * t248 + V_base(2);
t247 = pkin(1) * t279 - pkin(7) * t282;
t298 = V_base(4) * t247 - t248 * V_base(5) + V_base(3);
t190 = t279 * t312 + t217;
t297 = pkin(9) * t278 + t281 * t329;
t224 = t304 * t279;
t246 = t278 * pkin(2) - t281 * pkin(8);
t296 = t249 * t246 + (-t224 - t247) * t266 + t308;
t295 = (-Icges(3,3) * t282 + t279 * t300) * t249 + (Icges(3,3) * t279 + t282 * t300) * t250 + (Icges(3,5) * t278 + Icges(3,6) * t281) * t266;
t294 = pkin(10) * t278 + t281 * t318;
t293 = qJ(5) * t278 + t281 * t316;
t225 = t304 * t282;
t292 = t266 * t225 - t246 * t250 + t299;
t291 = t250 * t224 - t225 * t249 + t298;
t163 = -pkin(3) * t323 + t279 * t297;
t184 = -pkin(9) * t281 + t278 * t329;
t242 = -qJD(3) * t281 + t266;
t290 = -t163 * t242 + t217 * t184 + t296;
t164 = pkin(3) * t324 + t282 * t297;
t289 = t242 * t164 - t184 * t218 + t292;
t166 = -qJ(5) * t281 + t278 * t316;
t288 = t190 * t166 + t282 * t311 + t290;
t287 = t218 * t163 - t164 * t217 + t291;
t124 = t279 * t305 + t282 * t293;
t226 = t281 * t309 + t266;
t286 = t226 * t124 + t279 * t311 + t289;
t123 = t279 * t293 - t282 * t305;
t285 = -qJD(5) * t281 + t191 * t123 + t287;
t202 = -Icges(3,6) * t282 + t279 * t301;
t203 = Icges(3,6) * t279 + t282 * t301;
t205 = -Icges(3,5) * t282 + t279 * t302;
t206 = Icges(3,5) * t279 + t282 * t302;
t236 = Icges(3,2) * t281 + t326;
t239 = Icges(3,1) * t278 + t325;
t284 = (-t203 * t278 + t206 * t281) * t250 + (-t202 * t278 + t205 * t281) * t249 + (-t236 * t278 + t239 * t281) * t266;
t271 = Icges(2,4) * t282;
t263 = qJ(6) + t265;
t252 = cos(t263);
t251 = sin(t263);
t245 = rSges(2,1) * t282 - rSges(2,2) * t279;
t244 = rSges(2,1) * t279 + rSges(2,2) * t282;
t243 = rSges(3,1) * t278 + rSges(3,2) * t281;
t241 = Icges(2,1) * t282 - t327;
t240 = Icges(2,1) * t279 + t271;
t238 = -Icges(2,2) * t279 + t271;
t237 = Icges(2,2) * t282 + t327;
t231 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t230 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t229 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t222 = t280 * t319 + t324;
t221 = -t277 * t319 + t279 * t280;
t220 = t280 * t320 - t323;
t219 = -t277 * t320 - t280 * t282;
t214 = (-qJD(6) + t309) * t281 + t266;
t209 = rSges(3,3) * t279 + t282 * t303;
t208 = -rSges(3,3) * t282 + t279 * t303;
t207 = -rSges(4,3) * t281 + (rSges(4,1) * t280 - rSges(4,2) * t277) * t278;
t204 = -Icges(4,5) * t281 + (Icges(4,1) * t280 - Icges(4,4) * t277) * t278;
t201 = -Icges(4,6) * t281 + (Icges(4,4) * t280 - Icges(4,2) * t277) * t278;
t198 = -Icges(4,3) * t281 + (Icges(4,5) * t280 - Icges(4,6) * t277) * t278;
t189 = -rSges(5,3) * t281 + (rSges(5,1) * t269 - rSges(5,2) * t268) * t278;
t187 = -Icges(5,5) * t281 + (Icges(5,1) * t269 - Icges(5,4) * t268) * t278;
t186 = -Icges(5,6) * t281 + (Icges(5,4) * t269 - Icges(5,2) * t268) * t278;
t183 = t251 * t279 + t252 * t319;
t182 = -t251 * t319 + t252 * t279;
t181 = -t251 * t282 + t252 * t320;
t180 = -t251 * t320 - t252 * t282;
t179 = -rSges(6,3) * t281 + (rSges(6,1) * t260 - rSges(6,2) * t259) * t278;
t178 = V_base(5) * rSges(2,3) - t244 * t266 + t308;
t177 = t245 * t266 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t176 = -Icges(6,5) * t281 + (Icges(6,1) * t260 - Icges(6,4) * t259) * t278;
t175 = -Icges(6,6) * t281 + (Icges(6,4) * t260 - Icges(6,2) * t259) * t278;
t173 = t244 * V_base(4) - t245 * V_base(5) + V_base(3);
t172 = -rSges(7,3) * t281 + (rSges(7,1) * t252 - rSges(7,2) * t251) * t278;
t171 = -Icges(7,5) * t281 + (Icges(7,1) * t252 - Icges(7,4) * t251) * t278;
t170 = -Icges(7,6) * t281 + (Icges(7,4) * t252 - Icges(7,2) * t251) * t278;
t169 = -Icges(7,3) * t281 + (Icges(7,5) * t252 - Icges(7,6) * t251) * t278;
t168 = t282 * t310 + t191;
t167 = t279 * t310 + t190;
t162 = rSges(4,1) * t222 + rSges(4,2) * t221 + rSges(4,3) * t321;
t161 = rSges(4,1) * t220 + rSges(4,2) * t219 + rSges(4,3) * t322;
t160 = Icges(4,1) * t222 + Icges(4,4) * t221 + Icges(4,5) * t321;
t159 = Icges(4,1) * t220 + Icges(4,4) * t219 + Icges(4,5) * t322;
t158 = Icges(4,4) * t222 + Icges(4,2) * t221 + Icges(4,6) * t321;
t157 = Icges(4,4) * t220 + Icges(4,2) * t219 + Icges(4,6) * t322;
t156 = Icges(4,5) * t222 + Icges(4,6) * t221 + Icges(4,3) * t321;
t155 = Icges(4,5) * t220 + Icges(4,6) * t219 + Icges(4,3) * t322;
t154 = rSges(5,1) * t213 + rSges(5,2) * t212 + rSges(5,3) * t321;
t153 = rSges(5,1) * t211 + rSges(5,2) * t210 + rSges(5,3) * t322;
t152 = Icges(5,1) * t213 + Icges(5,4) * t212 + Icges(5,5) * t321;
t151 = Icges(5,1) * t211 + Icges(5,4) * t210 + Icges(5,5) * t322;
t150 = Icges(5,4) * t213 + Icges(5,2) * t212 + Icges(5,6) * t321;
t149 = Icges(5,4) * t211 + Icges(5,2) * t210 + Icges(5,6) * t322;
t145 = -pkin(10) * t281 + t278 * t318;
t144 = rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t321;
t143 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t322;
t142 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t321;
t141 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t322;
t140 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t321;
t139 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t322;
t135 = rSges(7,1) * t183 + rSges(7,2) * t182 + rSges(7,3) * t321;
t134 = rSges(7,1) * t181 + rSges(7,2) * t180 + rSges(7,3) * t322;
t133 = Icges(7,1) * t183 + Icges(7,4) * t182 + Icges(7,5) * t321;
t132 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t322;
t131 = Icges(7,4) * t183 + Icges(7,2) * t182 + Icges(7,6) * t321;
t130 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t322;
t129 = Icges(7,5) * t183 + Icges(7,6) * t182 + Icges(7,3) * t321;
t128 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t322;
t126 = t243 * t249 + (-t208 - t247) * t266 + t308;
t125 = t209 * t266 - t243 * t250 + t299;
t122 = t208 * t250 - t209 * t249 + t298;
t119 = t279 * t317 + t282 * t294;
t118 = t279 * t294 - t282 * t317;
t117 = -t161 * t242 + t207 * t217 + t296;
t116 = t162 * t242 - t207 * t218 + t292;
t115 = t161 * t218 - t162 * t217 + t291;
t114 = -t153 * t226 + t189 * t190 + t290;
t113 = t154 * t226 - t189 * t191 + t289;
t112 = t153 * t191 - t154 * t190 + t287;
t111 = t179 * t190 + (-t123 - t143) * t226 + t288;
t110 = t144 * t226 + (-t166 - t179) * t191 + t286;
t109 = t143 * t191 + (-t124 - t144) * t190 + t285;
t108 = -t134 * t214 + t145 * t190 + t167 * t172 + (-t118 - t123) * t226 + t288;
t107 = t119 * t226 + t135 * t214 - t168 * t172 + (-t145 - t166) * t191 + t286;
t106 = t118 * t191 + t134 * t168 - t135 * t167 + (-t119 - t124) * t190 + t285;
t1 = m(2) * (t173 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(3) * (t122 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + t217 * ((t156 * t322 + t158 * t219 + t160 * t220) * t218 + (t155 * t322 + t157 * t219 + t159 * t220) * t217 + (t198 * t322 + t201 * t219 + t204 * t220) * t242) / 0.2e1 + t167 * ((t129 * t322 + t131 * t180 + t133 * t181) * t168 + (t128 * t322 + t180 * t130 + t181 * t132) * t167 + (t169 * t322 + t170 * t180 + t171 * t181) * t214) / 0.2e1 + t250 * (t295 * t279 + t284 * t282) / 0.2e1 + t242 * ((-t155 * t217 - t156 * t218 - t198 * t242) * t281 + ((-t158 * t277 + t160 * t280) * t218 + (-t157 * t277 + t159 * t280) * t217 + (-t201 * t277 + t204 * t280) * t242) * t278) / 0.2e1 + t249 * (t284 * t279 - t295 * t282) / 0.2e1 + t218 * ((t156 * t321 + t158 * t221 + t160 * t222) * t218 + (t155 * t321 + t157 * t221 + t159 * t222) * t217 + (t198 * t321 + t201 * t221 + t204 * t222) * t242) / 0.2e1 + t168 * ((t129 * t321 + t182 * t131 + t183 * t133) * t168 + (t128 * t321 + t130 * t182 + t132 * t183) * t167 + (t169 * t321 + t170 * t182 + t171 * t183) * t214) / 0.2e1 + m(1) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + t214 * ((-t128 * t167 - t129 * t168 - t169 * t214) * t281 + ((-t131 * t251 + t133 * t252) * t168 + (-t130 * t251 + t132 * t252) * t167 + (-t170 * t251 + t171 * t252) * t214) * t278) / 0.2e1 + ((t175 * t192 + t176 * t193 + t186 * t210 + t187 * t211 + t322 * t334) * t226 + (t140 * t192 + t142 * t193 + t150 * t210 + t152 * t211 + t322 * t335) * t191 + (t139 * t192 + t141 * t193 + t149 * t210 + t151 * t211 + t336 * t322) * t190) * t190 / 0.2e1 + ((t175 * t194 + t176 * t195 + t186 * t212 + t187 * t213 + t321 * t334) * t226 + (t140 * t194 + t142 * t195 + t150 * t212 + t152 * t213 + t335 * t321) * t191 + (t139 * t194 + t141 * t195 + t149 * t212 + t151 * t213 + t321 * t336) * t190) * t191 / 0.2e1 + ((-t190 * t336 - t335 * t191 - t334 * t226) * t281 + ((-t175 * t259 + t176 * t260 - t186 * t268 + t187 * t269) * t226 + (-t140 * t259 + t142 * t260 - t150 * t268 + t152 * t269) * t191 + (-t139 * t259 + t141 * t260 - t149 * t268 + t151 * t269) * t190) * t278) * t226 / 0.2e1 + ((t203 * t281 + t206 * t278) * t250 + (t202 * t281 + t205 * t278) * t249 + (t236 * t281 + t239 * t278 + Icges(2,3)) * t266) * t266 / 0.2e1 + ((-t237 * t279 + t240 * t282 + Icges(1,4)) * V_base(5) + (-t238 * t279 + t241 * t282 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t237 * t282 + t240 * t279 + Icges(1,2)) * V_base(5) + (t238 * t282 + t241 * t279 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t266 * (Icges(2,5) * t282 - Icges(2,6) * t279) + V_base(5) * t266 * (Icges(2,5) * t279 + Icges(2,6) * t282) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
