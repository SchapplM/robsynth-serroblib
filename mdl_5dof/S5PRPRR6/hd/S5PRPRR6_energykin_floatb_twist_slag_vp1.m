% Calculate kinetic energy for
% S5PRPRR6
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:22
% EndTime: 2019-12-05 15:56:25
% DurationCPUTime: 3.46s
% Computational Cost: add. (2189->360), mult. (3934->517), div. (0->0), fcn. (4590->12), ass. (0->157)
t252 = sin(pkin(9));
t255 = cos(pkin(9));
t261 = cos(qJ(2));
t256 = cos(pkin(5));
t259 = sin(qJ(2));
t284 = t256 * t259;
t218 = t252 * t261 + t255 * t284;
t251 = sin(pkin(10));
t254 = cos(pkin(10));
t253 = sin(pkin(5));
t287 = t253 * t255;
t196 = -t218 * t251 - t254 * t287;
t273 = t251 * t287;
t197 = t218 * t254 - t273;
t283 = t256 * t261;
t217 = t252 * t259 - t255 * t283;
t149 = Icges(4,5) * t197 + Icges(4,6) * t196 + Icges(4,3) * t217;
t177 = Icges(3,4) * t218 - Icges(3,2) * t217 - Icges(3,6) * t287;
t298 = t149 - t177;
t220 = -t252 * t284 + t255 * t261;
t288 = t252 * t253;
t198 = -t220 * t251 + t254 * t288;
t274 = t251 * t288;
t199 = t220 * t254 + t274;
t219 = t252 * t283 + t255 * t259;
t150 = Icges(4,5) * t199 + Icges(4,6) * t198 + Icges(4,3) * t219;
t178 = Icges(3,4) * t220 - Icges(3,2) * t219 + Icges(3,6) * t288;
t297 = t150 - t178;
t286 = t253 * t259;
t215 = -t251 * t286 + t254 * t256;
t289 = t251 * t256;
t216 = t254 * t286 + t289;
t285 = t253 * t261;
t172 = Icges(4,5) * t216 + Icges(4,6) * t215 - Icges(4,3) * t285;
t206 = Icges(3,6) * t256 + (Icges(3,4) * t259 + Icges(3,2) * t261) * t253;
t296 = t172 - t206;
t292 = pkin(6) * t256;
t291 = pkin(3) * t254;
t290 = Icges(2,4) * t252;
t281 = qJD(2) * t253;
t280 = pkin(10) + qJ(4);
t279 = V_base(5) * qJ(1) + V_base(1);
t275 = qJD(1) + V_base(3);
t232 = t252 * t281 + V_base(4);
t243 = qJD(2) * t256 + V_base(6);
t272 = cos(t280);
t201 = qJD(4) * t219 + t232;
t271 = t253 * t272;
t231 = -t255 * t281 + V_base(5);
t200 = qJD(4) * t217 + t231;
t221 = -qJD(4) * t285 + t243;
t225 = pkin(1) * t252 - pkin(6) * t287;
t270 = -t225 * V_base(6) + V_base(5) * t292 + t279;
t226 = pkin(1) * t255 + pkin(6) * t288;
t269 = V_base(4) * t225 - V_base(5) * t226 + t275;
t268 = V_base(6) * t226 + V_base(2) + (-qJ(1) - t292) * V_base(4);
t224 = (pkin(2) * t259 - qJ(3) * t261) * t253;
t267 = qJD(3) * t219 + t231 * t224 + t270;
t188 = pkin(2) * t220 + qJ(3) * t219;
t266 = qJD(3) * t217 + t243 * t188 + t268;
t187 = pkin(2) * t218 + qJ(3) * t217;
t265 = -qJD(3) * t285 + t232 * t187 + t269;
t147 = -pkin(3) * t273 + pkin(7) * t217 + t218 * t291;
t183 = pkin(3) * t289 + (-pkin(7) * t261 + t259 * t291) * t253;
t264 = t231 * t183 + (-t147 - t187) * t243 + t267;
t148 = pkin(3) * t274 + pkin(7) * t219 + t220 * t291;
t263 = t243 * t148 + (-t183 - t224) * t232 + t266;
t262 = t232 * t147 + (-t148 - t188) * t231 + t265;
t260 = cos(qJ(5));
t258 = sin(qJ(5));
t249 = Icges(2,4) * t255;
t248 = sin(t280);
t240 = rSges(2,1) * t255 - rSges(2,2) * t252;
t239 = rSges(2,1) * t252 + rSges(2,2) * t255;
t238 = Icges(2,1) * t255 - t290;
t237 = Icges(2,1) * t252 + t249;
t236 = -Icges(2,2) * t252 + t249;
t235 = Icges(2,2) * t255 + t290;
t230 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t229 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t228 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t210 = t256 * rSges(3,3) + (rSges(3,1) * t259 + rSges(3,2) * t261) * t253;
t209 = t256 * t248 + t259 * t271;
t208 = t248 * t286 - t256 * t272;
t207 = Icges(3,5) * t256 + (Icges(3,1) * t259 + Icges(3,4) * t261) * t253;
t205 = Icges(3,3) * t256 + (Icges(3,5) * t259 + Icges(3,6) * t261) * t253;
t204 = V_base(5) * rSges(2,3) - t239 * V_base(6) + t279;
t203 = t240 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t195 = t239 * V_base(4) - t240 * V_base(5) + t275;
t194 = t209 * t260 - t258 * t285;
t193 = -t209 * t258 - t260 * t285;
t192 = t220 * t272 + t248 * t288;
t191 = t220 * t248 - t252 * t271;
t190 = t218 * t272 - t248 * t287;
t189 = t218 * t248 + t255 * t271;
t186 = qJD(5) * t208 + t221;
t185 = rSges(3,1) * t220 - rSges(3,2) * t219 + rSges(3,3) * t288;
t184 = rSges(3,1) * t218 - rSges(3,2) * t217 - rSges(3,3) * t287;
t182 = pkin(4) * t209 + pkin(8) * t208;
t181 = t216 * rSges(4,1) + t215 * rSges(4,2) - rSges(4,3) * t285;
t180 = Icges(3,1) * t220 - Icges(3,4) * t219 + Icges(3,5) * t288;
t179 = Icges(3,1) * t218 - Icges(3,4) * t217 - Icges(3,5) * t287;
t176 = Icges(3,5) * t220 - Icges(3,6) * t219 + Icges(3,3) * t288;
t175 = Icges(3,5) * t218 - Icges(3,6) * t217 - Icges(3,3) * t287;
t174 = Icges(4,1) * t216 + Icges(4,4) * t215 - Icges(4,5) * t285;
t173 = Icges(4,4) * t216 + Icges(4,2) * t215 - Icges(4,6) * t285;
t169 = t209 * rSges(5,1) - t208 * rSges(5,2) - rSges(5,3) * t285;
t168 = Icges(5,1) * t209 - Icges(5,4) * t208 - Icges(5,5) * t285;
t167 = Icges(5,4) * t209 - Icges(5,2) * t208 - Icges(5,6) * t285;
t166 = Icges(5,5) * t209 - Icges(5,6) * t208 - Icges(5,3) * t285;
t164 = t192 * t260 + t219 * t258;
t163 = -t192 * t258 + t219 * t260;
t162 = t190 * t260 + t217 * t258;
t161 = -t190 * t258 + t217 * t260;
t160 = qJD(5) * t191 + t201;
t159 = qJD(5) * t189 + t200;
t158 = pkin(4) * t192 + pkin(8) * t191;
t157 = pkin(4) * t190 + pkin(8) * t189;
t156 = rSges(4,1) * t199 + rSges(4,2) * t198 + rSges(4,3) * t219;
t155 = rSges(4,1) * t197 + rSges(4,2) * t196 + rSges(4,3) * t217;
t154 = Icges(4,1) * t199 + Icges(4,4) * t198 + Icges(4,5) * t219;
t153 = Icges(4,1) * t197 + Icges(4,4) * t196 + Icges(4,5) * t217;
t152 = Icges(4,4) * t199 + Icges(4,2) * t198 + Icges(4,6) * t219;
t151 = Icges(4,4) * t197 + Icges(4,2) * t196 + Icges(4,6) * t217;
t146 = rSges(5,1) * t192 - rSges(5,2) * t191 + rSges(5,3) * t219;
t145 = rSges(5,1) * t190 - rSges(5,2) * t189 + rSges(5,3) * t217;
t144 = Icges(5,1) * t192 - Icges(5,4) * t191 + Icges(5,5) * t219;
t143 = Icges(5,1) * t190 - Icges(5,4) * t189 + Icges(5,5) * t217;
t142 = Icges(5,4) * t192 - Icges(5,2) * t191 + Icges(5,6) * t219;
t141 = Icges(5,4) * t190 - Icges(5,2) * t189 + Icges(5,6) * t217;
t140 = Icges(5,5) * t192 - Icges(5,6) * t191 + Icges(5,3) * t219;
t139 = Icges(5,5) * t190 - Icges(5,6) * t189 + Icges(5,3) * t217;
t138 = rSges(6,1) * t194 + rSges(6,2) * t193 + rSges(6,3) * t208;
t137 = Icges(6,1) * t194 + Icges(6,4) * t193 + Icges(6,5) * t208;
t136 = Icges(6,4) * t194 + Icges(6,2) * t193 + Icges(6,6) * t208;
t135 = Icges(6,5) * t194 + Icges(6,6) * t193 + Icges(6,3) * t208;
t132 = -t184 * t243 + t210 * t231 + t270;
t131 = t185 * t243 - t210 * t232 + t268;
t130 = t184 * t232 - t185 * t231 + t269;
t129 = rSges(6,1) * t164 + rSges(6,2) * t163 + rSges(6,3) * t191;
t128 = rSges(6,1) * t162 + rSges(6,2) * t161 + rSges(6,3) * t189;
t127 = Icges(6,1) * t164 + Icges(6,4) * t163 + Icges(6,5) * t191;
t126 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t189;
t125 = Icges(6,4) * t164 + Icges(6,2) * t163 + Icges(6,6) * t191;
t124 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t189;
t123 = Icges(6,5) * t164 + Icges(6,6) * t163 + Icges(6,3) * t191;
t122 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t189;
t121 = t181 * t231 + (-t155 - t187) * t243 + t267;
t120 = t156 * t243 + (-t181 - t224) * t232 + t266;
t119 = t232 * t155 + (-t156 - t188) * t231 + t265;
t118 = -t145 * t221 + t169 * t200 + t264;
t117 = t146 * t221 - t169 * t201 + t263;
t116 = t201 * t145 - t200 * t146 + t262;
t115 = -t128 * t186 + t138 * t159 - t157 * t221 + t182 * t200 + t264;
t114 = t129 * t186 - t138 * t160 + t158 * t221 - t182 * t201 + t263;
t113 = t160 * t128 - t159 * t129 + t201 * t157 - t200 * t158 + t262;
t1 = m(1) * (t228 ^ 2 + t229 ^ 2 + t230 ^ 2) / 0.2e1 + m(2) * (t195 ^ 2 + t203 ^ 2 + t204 ^ 2) / 0.2e1 + m(3) * (t130 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(4) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(5) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + t201 * ((t140 * t219 - t142 * t191 + t144 * t192) * t201 + (t139 * t219 - t141 * t191 + t143 * t192) * t200 + (t166 * t219 - t167 * t191 + t168 * t192) * t221) / 0.2e1 + t200 * ((t140 * t217 - t142 * t189 + t144 * t190) * t201 + (t139 * t217 - t141 * t189 + t143 * t190) * t200 + (t166 * t217 - t167 * t189 + t168 * t190) * t221) / 0.2e1 + t221 * ((-t140 * t285 - t208 * t142 + t209 * t144) * t201 + (-t139 * t285 - t208 * t141 + t209 * t143) * t200 + (-t166 * t285 - t208 * t167 + t209 * t168) * t221) / 0.2e1 + m(6) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + t160 * ((t123 * t191 + t125 * t163 + t127 * t164) * t160 + (t122 * t191 + t124 * t163 + t126 * t164) * t159 + (t135 * t191 + t136 * t163 + t137 * t164) * t186) / 0.2e1 + t159 * ((t123 * t189 + t125 * t161 + t127 * t162) * t160 + (t122 * t189 + t124 * t161 + t126 * t162) * t159 + (t135 * t189 + t136 * t161 + t137 * t162) * t186) / 0.2e1 + t186 * ((t123 * t208 + t125 * t193 + t127 * t194) * t160 + (t122 * t208 + t124 * t193 + t126 * t194) * t159 + (t135 * t208 + t136 * t193 + t137 * t194) * t186) / 0.2e1 + ((t173 * t196 + t174 * t197 - t205 * t287 + t207 * t218 + t217 * t296) * t243 + (t152 * t196 + t154 * t197 - t176 * t287 + t180 * t218 + t217 * t297) * t232 + (t151 * t196 + t153 * t197 - t175 * t287 + t179 * t218 + t217 * t298) * t231) * t231 / 0.2e1 + ((t173 * t198 + t174 * t199 + t205 * t288 + t207 * t220 + t219 * t296) * t243 + (t152 * t198 + t154 * t199 + t176 * t288 + t180 * t220 + t219 * t297) * t232 + (t151 * t198 + t153 * t199 + t175 * t288 + t179 * t220 + t219 * t298) * t231) * t232 / 0.2e1 + ((t175 * t231 + t176 * t232 + t205 * t243) * t256 + ((t178 * t261 + t180 * t259) * t232 + (t177 * t261 + t179 * t259) * t231 + (t206 * t261 + t207 * t259) * t243) * t253 + (-t150 * t285 + t215 * t152 + t216 * t154) * t232 + (-t149 * t285 + t215 * t151 + t216 * t153) * t231 + (-t172 * t285 + t215 * t173 + t216 * t174) * t243) * t243 / 0.2e1 + ((-t235 * t252 + t237 * t255 + Icges(1,4)) * V_base(5) + (-t236 * t252 + t238 * t255 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t235 * t255 + t237 * t252 + Icges(1,2)) * V_base(5) + (t236 * t255 + t238 * t252 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t252 + Icges(2,6) * t255 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t255 - Icges(2,6) * t252 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
