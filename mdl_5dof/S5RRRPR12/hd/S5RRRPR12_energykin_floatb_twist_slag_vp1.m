% Calculate kinetic energy for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR12_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:36:59
% EndTime: 2019-12-31 21:37:02
% DurationCPUTime: 3.41s
% Computational Cost: add. (2157->355), mult. (4560->509), div. (0->0), fcn. (5454->12), ass. (0->153)
t298 = Icges(4,2) + Icges(5,3);
t256 = cos(pkin(5));
t260 = sin(qJ(1));
t261 = cos(qJ(2));
t280 = t260 * t261;
t259 = sin(qJ(2));
t262 = cos(qJ(1));
t281 = t259 * t262;
t219 = t256 * t281 + t280;
t258 = sin(qJ(3));
t254 = sin(pkin(5));
t283 = t254 * t262;
t291 = cos(qJ(3));
t200 = t219 * t291 - t258 * t283;
t279 = t261 * t262;
t282 = t259 * t260;
t218 = -t256 * t279 + t282;
t253 = sin(pkin(10));
t255 = cos(pkin(10));
t169 = -t200 * t253 + t218 * t255;
t287 = t218 * t253;
t170 = t200 * t255 + t287;
t272 = t254 * t291;
t199 = t219 * t258 + t262 * t272;
t297 = -Icges(4,4) * t200 + Icges(5,5) * t170 - Icges(4,6) * t218 + Icges(5,6) * t169 + t199 * t298;
t221 = -t256 * t282 + t279;
t285 = t254 * t260;
t202 = t221 * t291 + t258 * t285;
t220 = t256 * t280 + t281;
t171 = -t202 * t253 + t220 * t255;
t286 = t220 * t253;
t172 = t202 * t255 + t286;
t201 = t221 * t258 - t260 * t272;
t296 = -Icges(4,4) * t202 + Icges(5,5) * t172 - Icges(4,6) * t220 + Icges(5,6) * t171 + t201 * t298;
t217 = t256 * t258 + t259 * t272;
t284 = t254 * t261;
t195 = -t217 * t253 - t255 * t284;
t273 = t253 * t284;
t196 = t217 * t255 - t273;
t216 = t254 * t258 * t259 - t256 * t291;
t295 = -Icges(4,4) * t217 + Icges(5,5) * t196 + Icges(4,6) * t284 + Icges(5,6) * t195 + t216 * t298;
t290 = pkin(7) * t256;
t289 = pkin(4) * t255;
t288 = Icges(2,4) * t260;
t277 = qJD(2) * t254;
t276 = V_base(5) * pkin(6) + V_base(1);
t230 = t260 * t277 + V_base(4);
t249 = V_base(6) + qJD(1);
t198 = qJD(3) * t220 + t230;
t231 = qJD(2) * t256 + t249;
t229 = -t262 * t277 + V_base(5);
t224 = pkin(1) * t260 - pkin(7) * t283;
t271 = -t224 * t249 + V_base(5) * t290 + t276;
t225 = pkin(1) * t262 + pkin(7) * t285;
t270 = t224 * V_base(4) - t225 * V_base(5) + V_base(3);
t197 = qJD(3) * t218 + t229;
t214 = -qJD(3) * t284 + t231;
t269 = t249 * t225 + V_base(2) + (-pkin(6) - t290) * V_base(4);
t189 = pkin(2) * t219 + pkin(8) * t218;
t223 = (pkin(2) * t259 - pkin(8) * t261) * t254;
t268 = -t189 * t231 + t223 * t229 + t271;
t190 = pkin(2) * t221 + pkin(8) * t220;
t267 = t189 * t230 - t190 * t229 + t270;
t187 = pkin(3) * t217 + qJ(4) * t216;
t266 = qJD(4) * t201 + t187 * t197 + t268;
t161 = pkin(3) * t200 + qJ(4) * t199;
t265 = qJD(4) * t216 + t161 * t198 + t267;
t264 = t190 * t231 - t223 * t230 + t269;
t162 = pkin(3) * t202 + qJ(4) * t201;
t263 = qJD(4) * t199 + t162 * t214 + t264;
t252 = pkin(10) + qJ(5);
t250 = Icges(2,4) * t262;
t248 = cos(t252);
t247 = sin(t252);
t239 = rSges(2,1) * t262 - rSges(2,2) * t260;
t238 = rSges(2,1) * t260 + rSges(2,2) * t262;
t237 = Icges(2,1) * t262 - t288;
t236 = Icges(2,1) * t260 + t250;
t235 = -Icges(2,2) * t260 + t250;
t234 = Icges(2,2) * t262 + t288;
t228 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t227 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t226 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t210 = rSges(3,3) * t256 + (rSges(3,1) * t259 + rSges(3,2) * t261) * t254;
t209 = Icges(3,5) * t256 + (Icges(3,1) * t259 + Icges(3,4) * t261) * t254;
t208 = Icges(3,6) * t256 + (Icges(3,4) * t259 + Icges(3,2) * t261) * t254;
t207 = Icges(3,3) * t256 + (Icges(3,5) * t259 + Icges(3,6) * t261) * t254;
t206 = V_base(5) * rSges(2,3) - t238 * t249 + t276;
t205 = t239 * t249 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t203 = t238 * V_base(4) - t239 * V_base(5) + V_base(3);
t192 = t217 * t248 - t247 * t284;
t191 = -t217 * t247 - t248 * t284;
t188 = qJD(5) * t216 + t214;
t186 = rSges(3,1) * t221 - rSges(3,2) * t220 + rSges(3,3) * t285;
t185 = rSges(3,1) * t219 - rSges(3,2) * t218 - rSges(3,3) * t283;
t184 = Icges(3,1) * t221 - Icges(3,4) * t220 + Icges(3,5) * t285;
t183 = Icges(3,1) * t219 - Icges(3,4) * t218 - Icges(3,5) * t283;
t182 = Icges(3,4) * t221 - Icges(3,2) * t220 + Icges(3,6) * t285;
t181 = Icges(3,4) * t219 - Icges(3,2) * t218 - Icges(3,6) * t283;
t180 = Icges(3,5) * t221 - Icges(3,6) * t220 + Icges(3,3) * t285;
t179 = Icges(3,5) * t219 - Icges(3,6) * t218 - Icges(3,3) * t283;
t178 = rSges(4,1) * t217 - rSges(4,2) * t216 - rSges(4,3) * t284;
t177 = Icges(4,1) * t217 - Icges(4,4) * t216 - Icges(4,5) * t284;
t175 = Icges(4,5) * t217 - Icges(4,6) * t216 - Icges(4,3) * t284;
t168 = t202 * t248 + t220 * t247;
t167 = -t202 * t247 + t220 * t248;
t166 = t200 * t248 + t218 * t247;
t165 = -t200 * t247 + t218 * t248;
t164 = qJD(5) * t201 + t198;
t163 = qJD(5) * t199 + t197;
t159 = rSges(4,1) * t202 - rSges(4,2) * t201 + rSges(4,3) * t220;
t158 = rSges(4,1) * t200 - rSges(4,2) * t199 + rSges(4,3) * t218;
t157 = Icges(4,1) * t202 - Icges(4,4) * t201 + Icges(4,5) * t220;
t156 = Icges(4,1) * t200 - Icges(4,4) * t199 + Icges(4,5) * t218;
t153 = Icges(4,5) * t202 - Icges(4,6) * t201 + Icges(4,3) * t220;
t152 = Icges(4,5) * t200 - Icges(4,6) * t199 + Icges(4,3) * t218;
t150 = rSges(5,1) * t196 + rSges(5,2) * t195 + rSges(5,3) * t216;
t149 = Icges(5,1) * t196 + Icges(5,4) * t195 + Icges(5,5) * t216;
t148 = Icges(5,4) * t196 + Icges(5,2) * t195 + Icges(5,6) * t216;
t146 = -pkin(4) * t273 + pkin(9) * t216 + t217 * t289;
t145 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t216;
t144 = Icges(6,1) * t192 + Icges(6,4) * t191 + Icges(6,5) * t216;
t143 = Icges(6,4) * t192 + Icges(6,2) * t191 + Icges(6,6) * t216;
t142 = Icges(6,5) * t192 + Icges(6,6) * t191 + Icges(6,3) * t216;
t140 = -t185 * t231 + t210 * t229 + t271;
t139 = t186 * t231 - t210 * t230 + t269;
t138 = rSges(5,1) * t172 + rSges(5,2) * t171 + rSges(5,3) * t201;
t137 = rSges(5,1) * t170 + rSges(5,2) * t169 + rSges(5,3) * t199;
t136 = Icges(5,1) * t172 + Icges(5,4) * t171 + Icges(5,5) * t201;
t135 = Icges(5,1) * t170 + Icges(5,4) * t169 + Icges(5,5) * t199;
t134 = Icges(5,4) * t172 + Icges(5,2) * t171 + Icges(5,6) * t201;
t133 = Icges(5,4) * t170 + Icges(5,2) * t169 + Icges(5,6) * t199;
t130 = t185 * t230 - t186 * t229 + t270;
t129 = rSges(6,1) * t168 + rSges(6,2) * t167 + rSges(6,3) * t201;
t128 = rSges(6,1) * t166 + rSges(6,2) * t165 + rSges(6,3) * t199;
t127 = Icges(6,1) * t168 + Icges(6,4) * t167 + Icges(6,5) * t201;
t126 = Icges(6,1) * t166 + Icges(6,4) * t165 + Icges(6,5) * t199;
t125 = Icges(6,4) * t168 + Icges(6,2) * t167 + Icges(6,6) * t201;
t124 = Icges(6,4) * t166 + Icges(6,2) * t165 + Icges(6,6) * t199;
t123 = Icges(6,5) * t168 + Icges(6,6) * t167 + Icges(6,3) * t201;
t122 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t199;
t121 = pkin(4) * t286 + pkin(9) * t201 + t202 * t289;
t120 = pkin(4) * t287 + pkin(9) * t199 + t200 * t289;
t119 = -t158 * t214 + t178 * t197 + t268;
t118 = t159 * t214 - t178 * t198 + t264;
t117 = t158 * t198 - t159 * t197 + t267;
t116 = t150 * t197 + (-t137 - t161) * t214 + t266;
t115 = t138 * t214 + (-t150 - t187) * t198 + t263;
t114 = t137 * t198 + (-t138 - t162) * t197 + t265;
t113 = -t128 * t188 + t145 * t163 + t146 * t197 + (-t120 - t161) * t214 + t266;
t112 = t121 * t214 + t129 * t188 - t145 * t164 + (-t146 - t187) * t198 + t263;
t111 = t120 * t198 + t128 * t164 - t129 * t163 + (-t121 - t162) * t197 + t265;
t1 = m(1) * (t226 ^ 2 + t227 ^ 2 + t228 ^ 2) / 0.2e1 + m(2) * (t203 ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + m(3) * (t130 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + t230 * ((t180 * t285 - t220 * t182 + t221 * t184) * t230 + (t179 * t285 - t181 * t220 + t183 * t221) * t229 + (t207 * t285 - t208 * t220 + t209 * t221) * t231) / 0.2e1 + t229 * ((-t180 * t283 - t182 * t218 + t184 * t219) * t230 + (-t179 * t283 - t218 * t181 + t219 * t183) * t229 + (-t207 * t283 - t208 * t218 + t209 * t219) * t231) / 0.2e1 + t231 * ((t179 * t229 + t180 * t230 + t207 * t231) * t256 + ((t182 * t261 + t184 * t259) * t230 + (t181 * t261 + t183 * t259) * t229 + (t208 * t261 + t209 * t259) * t231) * t254) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(5) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(6) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + t164 * ((t201 * t123 + t167 * t125 + t168 * t127) * t164 + (t122 * t201 + t124 * t167 + t126 * t168) * t163 + (t142 * t201 + t143 * t167 + t144 * t168) * t188) / 0.2e1 + t163 * ((t123 * t199 + t125 * t165 + t127 * t166) * t164 + (t199 * t122 + t165 * t124 + t166 * t126) * t163 + (t142 * t199 + t143 * t165 + t144 * t166) * t188) / 0.2e1 + t188 * ((t123 * t216 + t125 * t191 + t127 * t192) * t164 + (t122 * t216 + t124 * t191 + t126 * t192) * t163 + (t216 * t142 + t191 * t143 + t192 * t144) * t188) / 0.2e1 + ((t148 * t169 + t149 * t170 + t175 * t218 + t177 * t200 + t199 * t295) * t214 + (t134 * t169 + t136 * t170 + t153 * t218 + t157 * t200 + t199 * t296) * t198 + (t169 * t133 + t170 * t135 + t218 * t152 + t200 * t156 + t297 * t199) * t197) * t197 / 0.2e1 + ((t148 * t171 + t149 * t172 + t175 * t220 + t177 * t202 + t201 * t295) * t214 + (t171 * t134 + t172 * t136 + t220 * t153 + t202 * t157 + t296 * t201) * t198 + (t133 * t171 + t135 * t172 + t152 * t220 + t202 * t156 + t201 * t297) * t197) * t198 / 0.2e1 + ((t195 * t148 + t196 * t149 - t175 * t284 + t217 * t177 + t295 * t216) * t214 + (t134 * t195 + t136 * t196 - t153 * t284 + t157 * t217 + t216 * t296) * t198 + (t133 * t195 + t135 * t196 - t152 * t284 + t156 * t217 + t216 * t297) * t197) * t214 / 0.2e1 + ((-t260 * t234 + t236 * t262 + Icges(1,4)) * V_base(5) + (-t260 * t235 + t262 * t237 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t262 * t234 + t260 * t236 + Icges(1,2)) * V_base(5) + (t235 * t262 + t260 * t237 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t260 + Icges(2,6) * t262) * V_base(5) + (Icges(2,5) * t262 - Icges(2,6) * t260) * V_base(4) + Icges(2,3) * t249 / 0.2e1) * t249;
T = t1;
