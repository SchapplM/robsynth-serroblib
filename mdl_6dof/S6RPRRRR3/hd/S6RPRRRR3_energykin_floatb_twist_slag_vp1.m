% Calculate kinetic energy for
% S6RPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:12
% EndTime: 2019-03-09 07:00:16
% DurationCPUTime: 3.35s
% Computational Cost: add. (2435->370), mult. (2152->554), div. (0->0), fcn. (2044->12), ass. (0->177)
t259 = sin(qJ(1));
t314 = pkin(1) * t259;
t262 = cos(qJ(1));
t313 = pkin(1) * t262;
t260 = cos(qJ(4));
t311 = t260 * pkin(4);
t310 = -pkin(6) - qJ(2);
t308 = Icges(2,4) * t259;
t254 = qJ(1) + pkin(11);
t244 = sin(t254);
t307 = Icges(3,4) * t244;
t258 = sin(qJ(3));
t306 = Icges(4,4) * t258;
t261 = cos(qJ(3));
t305 = Icges(4,4) * t261;
t257 = sin(qJ(4));
t304 = t244 * t257;
t303 = t244 * t258;
t302 = t244 * t261;
t245 = cos(t254);
t301 = t245 * t257;
t300 = t245 * t258;
t299 = t245 * t261;
t256 = qJ(4) + qJ(5);
t247 = sin(t256);
t298 = t247 * t261;
t248 = cos(t256);
t297 = t248 * t261;
t296 = t257 * t261;
t295 = t260 * t261;
t294 = pkin(5) * t248;
t292 = qJD(4) * t258;
t291 = qJD(5) * t258;
t290 = qJD(6) * t258;
t289 = -qJD(4) - qJD(5);
t246 = V_base(6) + qJD(1);
t288 = t246 * t313 + V_base(2);
t287 = V_base(5) * pkin(6) + V_base(1);
t218 = qJD(3) * t244 + V_base(4);
t210 = pkin(2) * t244 - pkin(7) * t245;
t284 = -t210 - t314;
t283 = pkin(5) * t247;
t282 = V_base(5) * qJ(2) + t287;
t281 = t314 * V_base(4) + qJD(2) + V_base(3);
t187 = t245 * t292 + t218;
t280 = pkin(3) * t261 + pkin(8) * t258;
t217 = -qJD(3) * t245 + V_base(5);
t279 = rSges(4,1) * t261 - rSges(4,2) * t258;
t155 = t245 * t291 + t187;
t278 = Icges(4,1) * t261 - t306;
t277 = -Icges(4,2) * t258 + t305;
t276 = Icges(4,5) * t261 - Icges(4,6) * t258;
t186 = t244 * t292 + t217;
t154 = t244 * t291 + t186;
t275 = pkin(9) * t258 + t261 * t311;
t274 = (-Icges(4,3) * t245 + t244 * t276) * t217 + (Icges(4,3) * t244 + t245 * t276) * t218 + (Icges(4,5) * t258 + Icges(4,6) * t261) * t246;
t273 = pkin(10) * t258 + t261 * t294;
t211 = pkin(2) * t245 + pkin(7) * t244;
t272 = t211 * t246 + t310 * V_base(4) + t288;
t198 = t280 * t244;
t236 = t258 * pkin(3) - pkin(8) * t261;
t271 = t217 * t236 + (-t198 + t284) * t246 + t282;
t270 = V_base(4) * t210 + (-t211 - t313) * V_base(5) + t281;
t199 = t280 * t245;
t269 = t199 * t246 - t218 * t236 + t272;
t144 = -pkin(4) * t301 + t244 * t275;
t176 = -pkin(9) * t261 + t258 * t311;
t232 = -qJD(4) * t261 + t246;
t268 = -t144 * t232 + t176 * t186 + t271;
t267 = t198 * t218 - t199 * t217 + t270;
t145 = pkin(4) * t304 + t245 * t275;
t266 = t145 * t232 - t176 * t187 + t269;
t265 = t144 * t187 - t145 * t186 + t267;
t159 = -Icges(4,6) * t245 + t244 * t277;
t160 = Icges(4,6) * t244 + t245 * t277;
t161 = -Icges(4,5) * t245 + t244 * t278;
t162 = Icges(4,5) * t244 + t245 * t278;
t222 = Icges(4,2) * t261 + t306;
t225 = Icges(4,1) * t258 + t305;
t264 = (-t160 * t258 + t162 * t261) * t218 + (-t159 * t258 + t161 * t261) * t217 + (-t222 * t258 + t225 * t261) * t246;
t251 = qJ(6) + t256;
t250 = Icges(2,4) * t262;
t242 = cos(t251);
t241 = sin(t251);
t240 = Icges(3,4) * t245;
t235 = rSges(2,1) * t262 - rSges(2,2) * t259;
t234 = rSges(2,1) * t259 + rSges(2,2) * t262;
t233 = rSges(4,1) * t258 + rSges(4,2) * t261;
t227 = Icges(2,1) * t262 - t308;
t226 = Icges(2,1) * t259 + t250;
t224 = -Icges(2,2) * t259 + t250;
t223 = Icges(2,2) * t262 + t308;
t215 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t214 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t213 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t209 = rSges(3,1) * t245 - rSges(3,2) * t244;
t208 = rSges(3,1) * t244 + rSges(3,2) * t245;
t207 = t261 * t289 + t246;
t206 = Icges(3,1) * t245 - t307;
t205 = Icges(3,1) * t244 + t240;
t204 = -Icges(3,2) * t244 + t240;
t203 = Icges(3,2) * t245 + t307;
t197 = (-qJD(6) + t289) * t261 + t246;
t195 = -rSges(5,3) * t261 + (rSges(5,1) * t260 - rSges(5,2) * t257) * t258;
t194 = -Icges(5,5) * t261 + (Icges(5,1) * t260 - Icges(5,4) * t257) * t258;
t193 = -Icges(5,6) * t261 + (Icges(5,4) * t260 - Icges(5,2) * t257) * t258;
t192 = -Icges(5,3) * t261 + (Icges(5,5) * t260 - Icges(5,6) * t257) * t258;
t191 = t245 * t295 + t304;
t190 = t244 * t260 - t245 * t296;
t189 = t244 * t295 - t301;
t188 = -t244 * t296 - t245 * t260;
t184 = -rSges(6,3) * t261 + (rSges(6,1) * t248 - rSges(6,2) * t247) * t258;
t183 = -Icges(6,5) * t261 + (Icges(6,1) * t248 - Icges(6,4) * t247) * t258;
t182 = -Icges(6,6) * t261 + (Icges(6,4) * t248 - Icges(6,2) * t247) * t258;
t181 = -Icges(6,3) * t261 + (Icges(6,5) * t248 - Icges(6,6) * t247) * t258;
t180 = t244 * t247 + t245 * t297;
t179 = t244 * t248 - t245 * t298;
t178 = t244 * t297 - t245 * t247;
t177 = -t244 * t298 - t245 * t248;
t174 = -rSges(7,3) * t261 + (rSges(7,1) * t242 - rSges(7,2) * t241) * t258;
t173 = -Icges(7,5) * t261 + (Icges(7,1) * t242 - Icges(7,4) * t241) * t258;
t172 = -Icges(7,6) * t261 + (Icges(7,4) * t242 - Icges(7,2) * t241) * t258;
t171 = -Icges(7,3) * t261 + (Icges(7,5) * t242 - Icges(7,6) * t241) * t258;
t170 = rSges(4,3) * t244 + t245 * t279;
t169 = -rSges(4,3) * t245 + t244 * t279;
t168 = t241 * t244 + t242 * t299;
t167 = -t241 * t299 + t242 * t244;
t166 = -t241 * t245 + t242 * t302;
t165 = -t241 * t302 - t242 * t245;
t164 = V_base(5) * rSges(2,3) - t234 * t246 + t287;
t163 = t235 * t246 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t156 = t234 * V_base(4) - t235 * V_base(5) + V_base(3);
t152 = -pkin(10) * t261 + t258 * t294;
t151 = t245 * t290 + t155;
t150 = t244 * t290 + t154;
t149 = V_base(5) * rSges(3,3) + (-t208 - t314) * t246 + t282;
t148 = t209 * t246 + (-rSges(3,3) + t310) * V_base(4) + t288;
t146 = t208 * V_base(4) + (-t209 - t313) * V_base(5) + t281;
t143 = rSges(5,1) * t191 + rSges(5,2) * t190 + rSges(5,3) * t300;
t142 = rSges(5,1) * t189 + rSges(5,2) * t188 + rSges(5,3) * t303;
t141 = Icges(5,1) * t191 + Icges(5,4) * t190 + Icges(5,5) * t300;
t140 = Icges(5,1) * t189 + Icges(5,4) * t188 + Icges(5,5) * t303;
t139 = Icges(5,4) * t191 + Icges(5,2) * t190 + Icges(5,6) * t300;
t138 = Icges(5,4) * t189 + Icges(5,2) * t188 + Icges(5,6) * t303;
t137 = Icges(5,5) * t191 + Icges(5,6) * t190 + Icges(5,3) * t300;
t136 = Icges(5,5) * t189 + Icges(5,6) * t188 + Icges(5,3) * t303;
t135 = rSges(6,1) * t180 + rSges(6,2) * t179 + rSges(6,3) * t300;
t134 = rSges(6,1) * t178 + rSges(6,2) * t177 + rSges(6,3) * t303;
t133 = Icges(6,1) * t180 + Icges(6,4) * t179 + Icges(6,5) * t300;
t132 = Icges(6,1) * t178 + Icges(6,4) * t177 + Icges(6,5) * t303;
t131 = Icges(6,4) * t180 + Icges(6,2) * t179 + Icges(6,6) * t300;
t130 = Icges(6,4) * t178 + Icges(6,2) * t177 + Icges(6,6) * t303;
t129 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t300;
t128 = Icges(6,5) * t178 + Icges(6,6) * t177 + Icges(6,3) * t303;
t126 = rSges(7,1) * t168 + rSges(7,2) * t167 + rSges(7,3) * t300;
t125 = rSges(7,1) * t166 + rSges(7,2) * t165 + rSges(7,3) * t303;
t124 = Icges(7,1) * t168 + Icges(7,4) * t167 + Icges(7,5) * t300;
t123 = Icges(7,1) * t166 + Icges(7,4) * t165 + Icges(7,5) * t303;
t122 = Icges(7,4) * t168 + Icges(7,2) * t167 + Icges(7,6) * t300;
t121 = Icges(7,4) * t166 + Icges(7,2) * t165 + Icges(7,6) * t303;
t120 = Icges(7,5) * t168 + Icges(7,6) * t167 + Icges(7,3) * t300;
t119 = Icges(7,5) * t166 + Icges(7,6) * t165 + Icges(7,3) * t303;
t117 = t244 * t283 + t245 * t273;
t116 = t244 * t273 - t245 * t283;
t115 = t217 * t233 + (-t169 + t284) * t246 + t282;
t114 = t170 * t246 - t218 * t233 + t272;
t113 = t169 * t218 - t170 * t217 + t270;
t112 = -t142 * t232 + t186 * t195 + t271;
t111 = t143 * t232 - t187 * t195 + t269;
t110 = t142 * t187 - t143 * t186 + t267;
t109 = -t134 * t207 + t154 * t184 + t268;
t108 = t135 * t207 - t155 * t184 + t266;
t107 = t134 * t155 - t135 * t154 + t265;
t106 = -t116 * t207 - t125 * t197 + t150 * t174 + t152 * t154 + t268;
t105 = t117 * t207 + t126 * t197 - t151 * t174 - t152 * t155 + t266;
t104 = t116 * t155 - t117 * t154 + t125 * t151 - t126 * t150 + t265;
t1 = t197 * ((-t119 * t150 - t120 * t151 - t171 * t197) * t261 + ((-t122 * t241 + t124 * t242) * t151 + (-t121 * t241 + t123 * t242) * t150 + (-t172 * t241 + t173 * t242) * t197) * t258) / 0.2e1 + m(2) * (t156 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + m(3) * (t146 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(4) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(7) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(6) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + t207 * ((-t128 * t154 - t129 * t155 - t181 * t207) * t261 + ((-t131 * t247 + t133 * t248) * t155 + (-t130 * t247 + t132 * t248) * t154 + (-t182 * t247 + t183 * t248) * t207) * t258) / 0.2e1 + t218 * (t244 * t274 + t245 * t264) / 0.2e1 + t217 * (t244 * t264 - t274 * t245) / 0.2e1 + m(1) * (t213 ^ 2 + t214 ^ 2 + t215 ^ 2) / 0.2e1 + t187 * ((t137 * t300 + t190 * t139 + t191 * t141) * t187 + (t136 * t300 + t138 * t190 + t140 * t191) * t186 + (t190 * t193 + t191 * t194 + t192 * t300) * t232) / 0.2e1 + t155 * ((t129 * t300 + t179 * t131 + t180 * t133) * t155 + (t128 * t300 + t130 * t179 + t132 * t180) * t154 + (t179 * t182 + t180 * t183 + t181 * t300) * t207) / 0.2e1 + t151 * ((t120 * t300 + t167 * t122 + t168 * t124) * t151 + (t119 * t300 + t121 * t167 + t123 * t168) * t150 + (t167 * t172 + t168 * t173 + t171 * t300) * t197) / 0.2e1 + t186 * ((t137 * t303 + t139 * t188 + t141 * t189) * t187 + (t136 * t303 + t188 * t138 + t189 * t140) * t186 + (t188 * t193 + t189 * t194 + t192 * t303) * t232) / 0.2e1 + t154 * ((t129 * t303 + t131 * t177 + t133 * t178) * t155 + (t128 * t303 + t177 * t130 + t178 * t132) * t154 + (t177 * t182 + t178 * t183 + t181 * t303) * t207) / 0.2e1 + t150 * ((t120 * t303 + t122 * t165 + t124 * t166) * t151 + (t119 * t303 + t165 * t121 + t166 * t123) * t150 + (t165 * t172 + t166 * t173 + t171 * t303) * t197) / 0.2e1 + t232 * ((-t136 * t186 - t137 * t187 - t192 * t232) * t261 + ((-t139 * t257 + t141 * t260) * t187 + (-t138 * t257 + t140 * t260) * t186 + (-t193 * t257 + t194 * t260) * t232) * t258) / 0.2e1 + ((t160 * t261 + t162 * t258) * t218 + (t159 * t261 + t161 * t258) * t217 + (t222 * t261 + t225 * t258 + Icges(2,3) + Icges(3,3)) * t246) * t246 / 0.2e1 + ((-t203 * t244 + t205 * t245 - t223 * t259 + t226 * t262 + Icges(1,4)) * V_base(5) + (-t204 * t244 + t206 * t245 - t224 * t259 + t227 * t262 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t203 * t245 + t205 * t244 + t223 * t262 + t226 * t259 + Icges(1,2)) * V_base(5) + (t204 * t245 + t206 * t244 + t224 * t262 + t227 * t259 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t246 * (Icges(2,5) * t259 + Icges(3,5) * t244 + Icges(2,6) * t262 + Icges(3,6) * t245) + V_base(4) * t246 * (Icges(2,5) * t262 + Icges(3,5) * t245 - Icges(2,6) * t259 - Icges(3,6) * t244) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
