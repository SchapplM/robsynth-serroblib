% Calculate kinetic energy for
% S6RPRRRR1
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:54:06
% EndTime: 2019-03-09 06:54:08
% DurationCPUTime: 2.08s
% Computational Cost: add. (2237->322), mult. (1550->467), div. (0->0), fcn. (1330->12), ass. (0->170)
t234 = sin(qJ(1));
t299 = pkin(1) * t234;
t237 = cos(qJ(1));
t298 = pkin(1) * t237;
t233 = sin(qJ(3));
t297 = pkin(3) * t233;
t231 = qJ(3) + qJ(4);
t222 = sin(t231);
t296 = pkin(4) * t222;
t236 = cos(qJ(3));
t295 = t236 * pkin(3);
t294 = -pkin(6) - qJ(2);
t292 = Icges(2,4) * t234;
t229 = qJ(1) + pkin(11);
t219 = sin(t229);
t291 = Icges(3,4) * t219;
t290 = Icges(4,4) * t233;
t289 = Icges(4,4) * t236;
t288 = Icges(5,4) * t222;
t223 = cos(t231);
t287 = Icges(5,4) * t223;
t226 = qJ(5) + t231;
t216 = sin(t226);
t286 = Icges(6,4) * t216;
t217 = cos(t226);
t285 = Icges(6,4) * t217;
t284 = t216 * t219;
t220 = cos(t229);
t283 = t216 * t220;
t232 = sin(qJ(6));
t282 = t219 * t232;
t235 = cos(qJ(6));
t281 = t219 * t235;
t280 = t220 * t232;
t279 = t220 * t235;
t278 = pkin(4) * t223;
t276 = qJD(6) * t216;
t275 = -qJD(3) - qJD(4);
t221 = V_base(6) + qJD(1);
t274 = t221 * t298 + V_base(2);
t273 = V_base(5) * pkin(6) + V_base(1);
t197 = qJD(3) * t219 + V_base(4);
t185 = t219 * pkin(2) - t220 * pkin(7);
t270 = -t185 - t299;
t269 = V_base(5) * qJ(2) + t273;
t268 = V_base(4) * t299 + qJD(2) + V_base(3);
t169 = qJD(4) * t219 + t197;
t125 = -pkin(8) * t220 + t219 * t295;
t267 = -t125 + t270;
t196 = -qJD(3) * t220 + V_base(5);
t266 = t196 * t297 + t269;
t265 = pkin(5) * t217 + pkin(10) * t216;
t264 = rSges(4,1) * t236 - rSges(4,2) * t233;
t263 = rSges(5,1) * t223 - rSges(5,2) * t222;
t262 = rSges(6,1) * t217 - rSges(6,2) * t216;
t161 = qJD(5) * t219 + t169;
t261 = Icges(4,1) * t236 - t290;
t260 = Icges(5,1) * t223 - t288;
t259 = Icges(6,1) * t217 - t286;
t258 = -Icges(4,2) * t233 + t289;
t257 = -Icges(5,2) * t222 + t287;
t256 = -Icges(6,2) * t216 + t285;
t255 = Icges(4,5) * t236 - Icges(4,6) * t233;
t254 = Icges(5,5) * t223 - Icges(5,6) * t222;
t253 = Icges(6,5) * t217 - Icges(6,6) * t216;
t117 = -pkin(9) * t220 + t219 * t278;
t252 = -t117 + t267;
t168 = t220 * t275 + V_base(5);
t251 = t168 * t296 + t266;
t160 = V_base(5) + (-qJD(5) + t275) * t220;
t250 = (-Icges(6,3) * t220 + t219 * t253) * t160 + (Icges(6,3) * t219 + t220 * t253) * t161 + (Icges(6,5) * t216 + Icges(6,6) * t217) * t221;
t249 = (-Icges(5,3) * t220 + t219 * t254) * t168 + (Icges(5,3) * t219 + t220 * t254) * t169 + (Icges(5,5) * t222 + Icges(5,6) * t223) * t221;
t248 = (-Icges(4,3) * t220 + t219 * t255) * t196 + (Icges(4,3) * t219 + t220 * t255) * t197 + (Icges(4,5) * t233 + Icges(4,6) * t236) * t221;
t186 = t220 * pkin(2) + t219 * pkin(7);
t247 = t221 * t186 + t294 * V_base(4) + t274;
t246 = V_base(4) * t185 + (-t186 - t298) * V_base(5) + t268;
t126 = pkin(8) * t219 + t220 * t295;
t245 = t221 * t126 - t197 * t297 + t247;
t244 = t197 * t125 - t126 * t196 + t246;
t118 = pkin(9) * t219 + t220 * t278;
t243 = t221 * t118 - t169 * t296 + t245;
t242 = t169 * t117 - t118 * t168 + t244;
t129 = -Icges(6,6) * t220 + t219 * t256;
t130 = Icges(6,6) * t219 + t220 * t256;
t131 = -Icges(6,5) * t220 + t219 * t259;
t132 = Icges(6,5) * t219 + t220 * t259;
t172 = Icges(6,2) * t217 + t286;
t173 = Icges(6,1) * t216 + t285;
t241 = (-t130 * t216 + t132 * t217) * t161 + (-t129 * t216 + t131 * t217) * t160 + (-t172 * t216 + t173 * t217) * t221;
t137 = -Icges(5,6) * t220 + t219 * t257;
t138 = Icges(5,6) * t219 + t220 * t257;
t139 = -Icges(5,5) * t220 + t219 * t260;
t140 = Icges(5,5) * t219 + t220 * t260;
t188 = Icges(5,2) * t223 + t288;
t189 = Icges(5,1) * t222 + t287;
t240 = (-t138 * t222 + t140 * t223) * t169 + (-t137 * t222 + t139 * t223) * t168 + (-t188 * t222 + t189 * t223) * t221;
t150 = -Icges(4,6) * t220 + t219 * t258;
t151 = Icges(4,6) * t219 + t220 * t258;
t152 = -Icges(4,5) * t220 + t219 * t261;
t153 = Icges(4,5) * t219 + t220 * t261;
t201 = Icges(4,2) * t236 + t290;
t204 = Icges(4,1) * t233 + t289;
t239 = (-t151 * t233 + t153 * t236) * t197 + (-t150 * t233 + t152 * t236) * t196 + (-t201 * t233 + t204 * t236) * t221;
t225 = Icges(2,4) * t237;
t215 = Icges(3,4) * t220;
t209 = rSges(2,1) * t237 - rSges(2,2) * t234;
t208 = rSges(2,1) * t234 + rSges(2,2) * t237;
t207 = rSges(4,1) * t233 + rSges(4,2) * t236;
t206 = Icges(2,1) * t237 - t292;
t205 = Icges(2,1) * t234 + t225;
t203 = -Icges(2,2) * t234 + t225;
t202 = Icges(2,2) * t237 + t292;
t195 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t194 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t193 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t191 = -qJD(6) * t217 + t221;
t190 = rSges(5,1) * t222 + rSges(5,2) * t223;
t184 = rSges(3,1) * t220 - rSges(3,2) * t219;
t183 = rSges(3,1) * t219 + rSges(3,2) * t220;
t182 = Icges(3,1) * t220 - t291;
t181 = Icges(3,1) * t219 + t215;
t180 = -Icges(3,2) * t219 + t215;
t179 = Icges(3,2) * t220 + t291;
t175 = pkin(5) * t216 - pkin(10) * t217;
t174 = rSges(6,1) * t216 + rSges(6,2) * t217;
t166 = t217 * t279 + t282;
t165 = -t217 * t280 + t281;
t164 = t217 * t281 - t280;
t163 = -t217 * t282 - t279;
t159 = t265 * t220;
t158 = t265 * t219;
t157 = rSges(4,3) * t219 + t220 * t264;
t156 = -rSges(4,3) * t220 + t219 * t264;
t155 = V_base(5) * rSges(2,3) - t208 * t221 + t273;
t154 = t209 * t221 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t147 = t208 * V_base(4) - t209 * V_base(5) + V_base(3);
t146 = -rSges(7,3) * t217 + (rSges(7,1) * t235 - rSges(7,2) * t232) * t216;
t145 = -Icges(7,5) * t217 + (Icges(7,1) * t235 - Icges(7,4) * t232) * t216;
t144 = -Icges(7,6) * t217 + (Icges(7,4) * t235 - Icges(7,2) * t232) * t216;
t143 = -Icges(7,3) * t217 + (Icges(7,5) * t235 - Icges(7,6) * t232) * t216;
t142 = rSges(5,3) * t219 + t220 * t263;
t141 = -rSges(5,3) * t220 + t219 * t263;
t134 = rSges(6,3) * t219 + t220 * t262;
t133 = -rSges(6,3) * t220 + t219 * t262;
t124 = t220 * t276 + t161;
t123 = t219 * t276 + t160;
t120 = V_base(5) * rSges(3,3) + (-t183 - t299) * t221 + t269;
t119 = t184 * t221 + (-rSges(3,3) + t294) * V_base(4) + t274;
t115 = t183 * V_base(4) + (-t184 - t298) * V_base(5) + t268;
t114 = rSges(7,1) * t166 + rSges(7,2) * t165 + rSges(7,3) * t283;
t113 = rSges(7,1) * t164 + rSges(7,2) * t163 + rSges(7,3) * t284;
t112 = Icges(7,1) * t166 + Icges(7,4) * t165 + Icges(7,5) * t283;
t111 = Icges(7,1) * t164 + Icges(7,4) * t163 + Icges(7,5) * t284;
t110 = Icges(7,4) * t166 + Icges(7,2) * t165 + Icges(7,6) * t283;
t109 = Icges(7,4) * t164 + Icges(7,2) * t163 + Icges(7,6) * t284;
t108 = Icges(7,5) * t166 + Icges(7,6) * t165 + Icges(7,3) * t283;
t107 = Icges(7,5) * t164 + Icges(7,6) * t163 + Icges(7,3) * t284;
t105 = t196 * t207 + (-t156 + t270) * t221 + t269;
t104 = t157 * t221 - t197 * t207 + t247;
t103 = t156 * t197 - t157 * t196 + t246;
t102 = t168 * t190 + (-t141 + t267) * t221 + t266;
t101 = t142 * t221 - t169 * t190 + t245;
t100 = t141 * t169 - t142 * t168 + t244;
t99 = t160 * t174 + (-t133 + t252) * t221 + t251;
t98 = t134 * t221 - t161 * t174 + t243;
t97 = -t113 * t191 + t123 * t146 + t160 * t175 + (-t158 + t252) * t221 + t251;
t96 = t114 * t191 - t124 * t146 + t159 * t221 - t161 * t175 + t243;
t95 = t133 * t161 - t134 * t160 + t242;
t94 = t113 * t124 - t114 * t123 + t158 * t161 - t159 * t160 + t242;
t1 = t191 * ((-t107 * t123 - t108 * t124 - t143 * t191) * t217 + ((-t110 * t232 + t112 * t235) * t124 + (-t109 * t232 + t111 * t235) * t123 + (-t144 * t232 + t145 * t235) * t191) * t216) / 0.2e1 + m(1) * (t193 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + m(2) * (t147 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(7) * (t94 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(6) * (t95 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(4) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + t197 * (t219 * t248 + t220 * t239) / 0.2e1 + t196 * (t219 * t239 - t220 * t248) / 0.2e1 + t169 * (t219 * t249 + t220 * t240) / 0.2e1 + t168 * (t219 * t240 - t220 * t249) / 0.2e1 + t161 * (t219 * t250 + t220 * t241) / 0.2e1 + t160 * (t219 * t241 - t220 * t250) / 0.2e1 + t124 * ((t108 * t283 + t165 * t110 + t166 * t112) * t124 + (t107 * t283 + t109 * t165 + t111 * t166) * t123 + (t143 * t283 + t144 * t165 + t145 * t166) * t191) / 0.2e1 + t123 * ((t108 * t284 + t110 * t163 + t112 * t164) * t124 + (t107 * t284 + t163 * t109 + t164 * t111) * t123 + (t143 * t284 + t144 * t163 + t145 * t164) * t191) / 0.2e1 + ((-t179 * t219 + t181 * t220 - t202 * t234 + t205 * t237 + Icges(1,4)) * V_base(5) + (-t180 * t219 + t182 * t220 - t203 * t234 + t206 * t237 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t179 * t220 + t181 * t219 + t202 * t237 + t205 * t234 + Icges(1,2)) * V_base(5) + (t180 * t220 + t182 * t219 + t203 * t237 + t206 * t234 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t151 * t236 + t153 * t233) * t197 + (t150 * t236 + t152 * t233) * t196 + (t138 * t223 + t140 * t222) * t169 + (t137 * t223 + t139 * t222) * t168 + (t130 * t217 + t132 * t216) * t161 + (t129 * t217 + t131 * t216) * t160 + (t172 * t217 + t173 * t216 + t188 * t223 + t189 * t222 + t201 * t236 + t204 * t233 + Icges(2,3) + Icges(3,3)) * t221) * t221 / 0.2e1 + t221 * V_base(5) * (Icges(2,5) * t234 + Icges(3,5) * t219 + Icges(2,6) * t237 + Icges(3,6) * t220) + t221 * V_base(4) * (Icges(2,5) * t237 + Icges(3,5) * t220 - Icges(2,6) * t234 - Icges(3,6) * t219) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
