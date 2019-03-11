% Calculate kinetic energy for
% S6RRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:46:37
% EndTime: 2019-03-10 03:46:41
% DurationCPUTime: 3.98s
% Computational Cost: add. (2544->410), mult. (2906->635), div. (0->0), fcn. (2856->12), ass. (0->188)
t283 = cos(qJ(3));
t333 = t283 * pkin(3);
t282 = sin(qJ(1));
t331 = Icges(2,4) * t282;
t281 = sin(qJ(2));
t330 = Icges(3,4) * t281;
t284 = cos(qJ(2));
t329 = Icges(3,4) * t284;
t280 = sin(qJ(3));
t328 = t280 * t282;
t285 = cos(qJ(1));
t327 = t280 * t285;
t326 = t281 * t282;
t325 = t281 * t285;
t324 = t282 * t284;
t323 = t284 * t285;
t279 = qJ(3) + qJ(4);
t273 = qJ(5) + t279;
t265 = cos(t273);
t322 = pkin(5) * t265;
t264 = sin(t273);
t321 = pkin(5) * t264;
t271 = cos(t279);
t320 = pkin(4) * t271;
t317 = qJD(3) * t281;
t316 = qJD(4) * t281;
t315 = qJD(5) * t281;
t314 = qJD(6) * t281;
t313 = -qJD(3) - qJD(4);
t312 = V_base(5) * pkin(6) + V_base(1);
t253 = qJD(2) * t282 + V_base(4);
t268 = V_base(6) + qJD(1);
t270 = sin(t279);
t309 = pkin(4) * t270;
t308 = -qJD(5) + t313;
t221 = t285 * t317 + t253;
t307 = pkin(2) * t284 + pkin(8) * t281;
t252 = -qJD(2) * t285 + V_base(5);
t306 = rSges(3,1) * t284 - rSges(3,2) * t281;
t194 = t285 * t316 + t221;
t305 = Icges(3,1) * t284 - t330;
t304 = -Icges(3,2) * t281 + t329;
t303 = Icges(3,5) * t284 - Icges(3,6) * t281;
t220 = t282 * t317 + t252;
t251 = pkin(1) * t285 + pkin(7) * t282;
t302 = -V_base(4) * pkin(6) + t268 * t251 + V_base(2);
t170 = t285 * t315 + t194;
t250 = pkin(1) * t282 - pkin(7) * t285;
t301 = V_base(4) * t250 - t251 * V_base(5) + V_base(3);
t193 = t282 * t316 + t220;
t300 = pkin(9) * t281 + t284 * t333;
t227 = t307 * t282;
t249 = t281 * pkin(2) - t284 * pkin(8);
t299 = t252 * t249 + (-t227 - t250) * t268 + t312;
t298 = (-Icges(3,3) * t285 + t282 * t303) * t252 + (Icges(3,3) * t282 + t285 * t303) * t253 + (Icges(3,5) * t281 + Icges(3,6) * t284) * t268;
t169 = t282 * t315 + t193;
t297 = pkin(11) * t281 + t284 * t322;
t296 = pkin(10) * t281 + t284 * t320;
t228 = t307 * t285;
t295 = t268 * t228 - t249 * t253 + t302;
t294 = t253 * t227 - t228 * t252 + t301;
t163 = -pkin(3) * t327 + t282 * t300;
t182 = -pkin(9) * t284 + t281 * t333;
t245 = -qJD(3) * t284 + t268;
t293 = -t163 * t245 + t220 * t182 + t299;
t164 = pkin(3) * t328 + t285 * t300;
t292 = t245 * t164 - t182 * t221 + t295;
t291 = t221 * t163 - t164 * t220 + t294;
t123 = t282 * t296 - t285 * t309;
t168 = -pkin(10) * t284 + t281 * t320;
t229 = t284 * t313 + t268;
t290 = -t123 * t229 + t193 * t168 + t293;
t124 = t282 * t309 + t285 * t296;
t289 = t229 * t124 - t168 * t194 + t292;
t288 = t194 * t123 - t124 * t193 + t291;
t204 = -Icges(3,6) * t285 + t282 * t304;
t205 = Icges(3,6) * t282 + t285 * t304;
t207 = -Icges(3,5) * t285 + t282 * t305;
t208 = Icges(3,5) * t282 + t285 * t305;
t239 = Icges(3,2) * t284 + t330;
t242 = Icges(3,1) * t281 + t329;
t287 = (-t205 * t281 + t208 * t284) * t253 + (-t204 * t281 + t207 * t284) * t252 + (-t239 * t281 + t242 * t284) * t268;
t272 = Icges(2,4) * t285;
t267 = qJ(6) + t273;
t261 = cos(t267);
t260 = sin(t267);
t248 = rSges(2,1) * t285 - rSges(2,2) * t282;
t247 = rSges(2,1) * t282 + rSges(2,2) * t285;
t246 = rSges(3,1) * t281 + rSges(3,2) * t284;
t244 = Icges(2,1) * t285 - t331;
t243 = Icges(2,1) * t282 + t272;
t241 = -Icges(2,2) * t282 + t272;
t240 = Icges(2,2) * t285 + t331;
t234 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t233 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t232 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t225 = t283 * t323 + t328;
t224 = -t280 * t323 + t282 * t283;
t223 = t283 * t324 - t327;
t222 = -t280 * t324 - t283 * t285;
t217 = t284 * t308 + t268;
t216 = t270 * t282 + t271 * t323;
t215 = -t270 * t323 + t271 * t282;
t214 = -t270 * t285 + t271 * t324;
t213 = -t270 * t324 - t271 * t285;
t212 = rSges(3,3) * t282 + t285 * t306;
t211 = -rSges(3,3) * t285 + t282 * t306;
t210 = -rSges(4,3) * t284 + (rSges(4,1) * t283 - rSges(4,2) * t280) * t281;
t206 = -Icges(4,5) * t284 + (Icges(4,1) * t283 - Icges(4,4) * t280) * t281;
t203 = -Icges(4,6) * t284 + (Icges(4,4) * t283 - Icges(4,2) * t280) * t281;
t200 = -Icges(4,3) * t284 + (Icges(4,5) * t283 - Icges(4,6) * t280) * t281;
t198 = t264 * t282 + t265 * t323;
t197 = -t264 * t323 + t265 * t282;
t196 = -t264 * t285 + t265 * t324;
t195 = -t264 * t324 - t265 * t285;
t192 = -rSges(5,3) * t284 + (rSges(5,1) * t271 - rSges(5,2) * t270) * t281;
t190 = -Icges(5,5) * t284 + (Icges(5,1) * t271 - Icges(5,4) * t270) * t281;
t189 = -Icges(5,6) * t284 + (Icges(5,4) * t271 - Icges(5,2) * t270) * t281;
t188 = -Icges(5,3) * t284 + (Icges(5,5) * t271 - Icges(5,6) * t270) * t281;
t187 = (-qJD(6) + t308) * t284 + t268;
t186 = t260 * t282 + t261 * t323;
t185 = -t260 * t323 + t261 * t282;
t184 = -t260 * t285 + t261 * t324;
t183 = -t260 * t324 - t261 * t285;
t181 = -rSges(6,3) * t284 + (rSges(6,1) * t265 - rSges(6,2) * t264) * t281;
t180 = -Icges(6,5) * t284 + (Icges(6,1) * t265 - Icges(6,4) * t264) * t281;
t179 = -Icges(6,6) * t284 + (Icges(6,4) * t265 - Icges(6,2) * t264) * t281;
t178 = -Icges(6,3) * t284 + (Icges(6,5) * t265 - Icges(6,6) * t264) * t281;
t177 = V_base(5) * rSges(2,3) - t247 * t268 + t312;
t176 = t248 * t268 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t175 = -rSges(7,3) * t284 + (rSges(7,1) * t261 - rSges(7,2) * t260) * t281;
t174 = t247 * V_base(4) - t248 * V_base(5) + V_base(3);
t173 = -Icges(7,5) * t284 + (Icges(7,1) * t261 - Icges(7,4) * t260) * t281;
t172 = -Icges(7,6) * t284 + (Icges(7,4) * t261 - Icges(7,2) * t260) * t281;
t171 = -Icges(7,3) * t284 + (Icges(7,5) * t261 - Icges(7,6) * t260) * t281;
t167 = t285 * t314 + t170;
t166 = t282 * t314 + t169;
t162 = rSges(4,1) * t225 + rSges(4,2) * t224 + rSges(4,3) * t325;
t161 = rSges(4,1) * t223 + rSges(4,2) * t222 + rSges(4,3) * t326;
t160 = Icges(4,1) * t225 + Icges(4,4) * t224 + Icges(4,5) * t325;
t159 = Icges(4,1) * t223 + Icges(4,4) * t222 + Icges(4,5) * t326;
t158 = Icges(4,4) * t225 + Icges(4,2) * t224 + Icges(4,6) * t325;
t157 = Icges(4,4) * t223 + Icges(4,2) * t222 + Icges(4,6) * t326;
t156 = Icges(4,5) * t225 + Icges(4,6) * t224 + Icges(4,3) * t325;
t155 = Icges(4,5) * t223 + Icges(4,6) * t222 + Icges(4,3) * t326;
t154 = rSges(5,1) * t216 + rSges(5,2) * t215 + rSges(5,3) * t325;
t153 = rSges(5,1) * t214 + rSges(5,2) * t213 + rSges(5,3) * t326;
t152 = Icges(5,1) * t216 + Icges(5,4) * t215 + Icges(5,5) * t325;
t151 = Icges(5,1) * t214 + Icges(5,4) * t213 + Icges(5,5) * t326;
t150 = Icges(5,4) * t216 + Icges(5,2) * t215 + Icges(5,6) * t325;
t149 = Icges(5,4) * t214 + Icges(5,2) * t213 + Icges(5,6) * t326;
t148 = Icges(5,5) * t216 + Icges(5,6) * t215 + Icges(5,3) * t325;
t147 = Icges(5,5) * t214 + Icges(5,6) * t213 + Icges(5,3) * t326;
t145 = -pkin(11) * t284 + t281 * t322;
t144 = rSges(6,1) * t198 + rSges(6,2) * t197 + rSges(6,3) * t325;
t143 = rSges(6,1) * t196 + rSges(6,2) * t195 + rSges(6,3) * t326;
t142 = Icges(6,1) * t198 + Icges(6,4) * t197 + Icges(6,5) * t325;
t141 = Icges(6,1) * t196 + Icges(6,4) * t195 + Icges(6,5) * t326;
t140 = Icges(6,4) * t198 + Icges(6,2) * t197 + Icges(6,6) * t325;
t139 = Icges(6,4) * t196 + Icges(6,2) * t195 + Icges(6,6) * t326;
t138 = Icges(6,5) * t198 + Icges(6,6) * t197 + Icges(6,3) * t325;
t137 = Icges(6,5) * t196 + Icges(6,6) * t195 + Icges(6,3) * t326;
t135 = rSges(7,1) * t186 + rSges(7,2) * t185 + rSges(7,3) * t325;
t134 = rSges(7,1) * t184 + rSges(7,2) * t183 + rSges(7,3) * t326;
t133 = Icges(7,1) * t186 + Icges(7,4) * t185 + Icges(7,5) * t325;
t132 = Icges(7,1) * t184 + Icges(7,4) * t183 + Icges(7,5) * t326;
t131 = Icges(7,4) * t186 + Icges(7,2) * t185 + Icges(7,6) * t325;
t130 = Icges(7,4) * t184 + Icges(7,2) * t183 + Icges(7,6) * t326;
t129 = Icges(7,5) * t186 + Icges(7,6) * t185 + Icges(7,3) * t325;
t128 = Icges(7,5) * t184 + Icges(7,6) * t183 + Icges(7,3) * t326;
t126 = t246 * t252 + (-t211 - t250) * t268 + t312;
t125 = t212 * t268 - t246 * t253 + t302;
t122 = t211 * t253 - t212 * t252 + t301;
t119 = t282 * t321 + t285 * t297;
t118 = t282 * t297 - t285 * t321;
t117 = -t161 * t245 + t210 * t220 + t299;
t116 = t162 * t245 - t210 * t221 + t295;
t115 = t161 * t221 - t162 * t220 + t294;
t114 = -t153 * t229 + t192 * t193 + t293;
t113 = t154 * t229 - t192 * t194 + t292;
t112 = t153 * t194 - t154 * t193 + t291;
t111 = -t143 * t217 + t169 * t181 + t290;
t110 = t144 * t217 - t170 * t181 + t289;
t109 = t143 * t170 - t144 * t169 + t288;
t108 = -t118 * t217 - t134 * t187 + t145 * t169 + t166 * t175 + t290;
t107 = t119 * t217 + t135 * t187 - t145 * t170 - t167 * t175 + t289;
t106 = t118 * t170 - t119 * t169 + t134 * t167 - t135 * t166 + t288;
t1 = (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((t240 * t285 + t243 * t282 + Icges(1,2)) * V_base(5) + (t241 * t285 + t244 * t282 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + m(2) * (t174 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + m(3) * (t122 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + ((t205 * t284 + t208 * t281) * t253 + (t204 * t284 + t207 * t281) * t252 + (t239 * t284 + t242 * t281 + Icges(2,3)) * t268) * t268 / 0.2e1 + t253 * (t298 * t282 + t287 * t285) / 0.2e1 + t252 * (t287 * t282 - t298 * t285) / 0.2e1 + t268 * V_base(4) * (Icges(2,5) * t285 - Icges(2,6) * t282) + t187 * ((-t128 * t166 - t129 * t167 - t171 * t187) * t284 + ((-t131 * t260 + t133 * t261) * t167 + (-t130 * t260 + t132 * t261) * t166 + (-t172 * t260 + t173 * t261) * t187) * t281) / 0.2e1 + t217 * ((-t137 * t169 - t138 * t170 - t178 * t217) * t284 + ((-t140 * t264 + t142 * t265) * t170 + (-t139 * t264 + t141 * t265) * t169 + (-t179 * t264 + t180 * t265) * t217) * t281) / 0.2e1 + t245 * ((-t155 * t220 - t156 * t221 - t200 * t245) * t284 + ((-t158 * t280 + t160 * t283) * t221 + (-t157 * t280 + t159 * t283) * t220 + (-t203 * t280 + t206 * t283) * t245) * t281) / 0.2e1 + t229 * ((-t147 * t193 - t148 * t194 - t188 * t229) * t284 + ((-t150 * t270 + t152 * t271) * t194 + (-t149 * t270 + t151 * t271) * t193 + (-t189 * t270 + t190 * t271) * t229) * t281) / 0.2e1 + m(1) * (t232 ^ 2 + t233 ^ 2 + t234 ^ 2) / 0.2e1 + V_base(5) * t268 * (Icges(2,5) * t282 + Icges(2,6) * t285) + ((-t240 * t282 + t243 * t285 + Icges(1,4)) * V_base(5) + (-t241 * t282 + t244 * t285 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + t221 * ((t156 * t325 + t158 * t224 + t160 * t225) * t221 + (t155 * t325 + t157 * t224 + t159 * t225) * t220 + (t200 * t325 + t203 * t224 + t206 * t225) * t245) / 0.2e1 + t194 * ((t148 * t325 + t150 * t215 + t152 * t216) * t194 + (t147 * t325 + t149 * t215 + t151 * t216) * t193 + (t188 * t325 + t189 * t215 + t190 * t216) * t229) / 0.2e1 + t170 * ((t138 * t325 + t197 * t140 + t198 * t142) * t170 + (t137 * t325 + t139 * t197 + t141 * t198) * t169 + (t178 * t325 + t179 * t197 + t180 * t198) * t217) / 0.2e1 + t167 * ((t129 * t325 + t185 * t131 + t186 * t133) * t167 + (t128 * t325 + t130 * t185 + t132 * t186) * t166 + (t171 * t325 + t172 * t185 + t173 * t186) * t187) / 0.2e1 + t220 * ((t156 * t326 + t158 * t222 + t160 * t223) * t221 + (t155 * t326 + t157 * t222 + t159 * t223) * t220 + (t200 * t326 + t203 * t222 + t206 * t223) * t245) / 0.2e1 + t193 * ((t148 * t326 + t150 * t213 + t152 * t214) * t194 + (t147 * t326 + t149 * t213 + t151 * t214) * t193 + (t188 * t326 + t189 * t213 + t190 * t214) * t229) / 0.2e1 + t169 * ((t138 * t326 + t140 * t195 + t142 * t196) * t170 + (t137 * t326 + t195 * t139 + t196 * t141) * t169 + (t178 * t326 + t179 * t195 + t180 * t196) * t217) / 0.2e1 + t166 * ((t129 * t326 + t131 * t183 + t133 * t184) * t167 + (t128 * t326 + t183 * t130 + t184 * t132) * t166 + (t171 * t326 + t172 * t183 + t173 * t184) * t187) / 0.2e1;
T  = t1;
