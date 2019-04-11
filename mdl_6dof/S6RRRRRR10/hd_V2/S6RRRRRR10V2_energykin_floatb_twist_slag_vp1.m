% Calculate kinetic energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR10V2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR10V2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:32
% EndTime: 2019-04-11 14:41:35
% DurationCPUTime: 3.49s
% Computational Cost: add. (2573->368), mult. (3404->563), div. (0->0), fcn. (3784->12), ass. (0->172)
t306 = cos(qJ(5));
t258 = sin(qJ(1));
t305 = pkin(1) * t258;
t262 = cos(qJ(1));
t304 = pkin(1) * t262;
t257 = sin(qJ(2));
t303 = pkin(2) * t257;
t261 = cos(qJ(2));
t302 = pkin(2) * t261;
t301 = Icges(2,4) * t258;
t300 = Icges(3,4) * t257;
t299 = Icges(3,4) * t261;
t253 = qJ(2) + qJ(3);
t249 = sin(t253);
t298 = Icges(4,4) * t249;
t250 = cos(t253);
t297 = Icges(4,4) * t250;
t256 = sin(qJ(4));
t296 = t249 * t256;
t295 = t249 * t258;
t294 = t249 * t262;
t293 = t256 * t262;
t292 = t258 * t256;
t260 = cos(qJ(4));
t291 = t258 * t260;
t290 = t260 * t262;
t289 = qJD(4) * t249;
t288 = V_base(5) * pkin(4) + V_base(1);
t285 = t249 * t306;
t242 = qJD(2) * t258 + V_base(4);
t246 = V_base(6) + qJD(1);
t209 = t302 * t258;
t284 = -t209 - t305;
t241 = -qJD(2) * t262 + V_base(5);
t283 = t241 * t303 + t288;
t217 = qJD(3) * t258 + t242;
t282 = pkin(3) * t250 + pkin(5) * t249;
t281 = rSges(3,1) * t261 - rSges(3,2) * t257;
t280 = rSges(4,1) * t250 - rSges(4,2) * t249;
t188 = t262 * t289 + t217;
t279 = Icges(3,1) * t261 - t300;
t278 = Icges(4,1) * t250 - t298;
t277 = -Icges(3,2) * t257 + t299;
t276 = -Icges(4,2) * t249 + t297;
t275 = Icges(3,5) * t261 - Icges(3,6) * t257;
t274 = Icges(4,5) * t250 - Icges(4,6) * t249;
t273 = -V_base(4) * pkin(4) + t246 * t304 + V_base(2);
t222 = -qJD(4) * t250 + t246;
t207 = t250 * t293 - t291;
t161 = qJD(5) * t207 + t188;
t216 = V_base(5) + (-qJD(2) - qJD(3)) * t262;
t272 = -t304 * V_base(5) + V_base(4) * t305 + V_base(3);
t189 = qJD(5) * t296 + t222;
t187 = t258 * t289 + t216;
t271 = (-Icges(4,3) * t262 + t258 * t274) * t216 + (Icges(4,3) * t258 + t262 * t274) * t217 + (Icges(4,5) * t249 + Icges(4,6) * t250) * t246;
t270 = (-Icges(3,3) * t262 + t258 * t275) * t241 + (Icges(3,3) * t258 + t262 * t275) * t242 + (Icges(3,5) * t257 + Icges(3,6) * t261) * t246;
t205 = t250 * t292 + t290;
t160 = qJD(5) * t205 + t187;
t210 = t302 * t262;
t269 = t246 * t210 - t242 * t303 + t273;
t268 = t242 * t209 - t241 * t210 + t272;
t203 = t282 * t258;
t215 = pkin(3) * t249 - pkin(5) * t250;
t267 = t216 * t215 + (-t203 + t284) * t246 + t283;
t204 = t282 * t262;
t266 = t246 * t204 - t215 * t217 + t269;
t265 = t217 * t203 - t216 * t204 + t268;
t180 = -Icges(4,6) * t262 + t258 * t276;
t181 = Icges(4,6) * t258 + t262 * t276;
t182 = -Icges(4,5) * t262 + t258 * t278;
t183 = Icges(4,5) * t258 + t262 * t278;
t212 = Icges(4,2) * t250 + t298;
t213 = Icges(4,1) * t249 + t297;
t264 = (-t181 * t249 + t183 * t250) * t217 + (-t180 * t249 + t182 * t250) * t216 + (-t212 * t249 + t213 * t250) * t246;
t192 = -Icges(3,6) * t262 + t258 * t277;
t193 = Icges(3,6) * t258 + t262 * t277;
t194 = -Icges(3,5) * t262 + t258 * t279;
t195 = Icges(3,5) * t258 + t262 * t279;
t228 = Icges(3,2) * t261 + t300;
t231 = Icges(3,1) * t257 + t299;
t263 = (-t193 * t257 + t195 * t261) * t242 + (-t192 * t257 + t194 * t261) * t241 + (-t228 * t257 + t231 * t261) * t246;
t259 = cos(qJ(6));
t255 = sin(qJ(5));
t254 = sin(qJ(6));
t251 = Icges(2,4) * t262;
t236 = rSges(2,1) * t262 - t258 * rSges(2,2);
t235 = t258 * rSges(2,1) + rSges(2,2) * t262;
t234 = rSges(3,1) * t257 + rSges(3,2) * t261;
t233 = Icges(2,1) * t262 - t301;
t232 = Icges(2,1) * t258 + t251;
t230 = -Icges(2,2) * t258 + t251;
t229 = Icges(2,2) * t262 + t301;
t221 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t220 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t219 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t214 = rSges(4,1) * t249 + rSges(4,2) * t250;
t208 = t250 * t290 + t292;
t206 = t250 * t291 - t293;
t202 = -t250 * t255 + t260 * t285;
t201 = t249 * t260 * t255 + t250 * t306;
t200 = t258 * rSges(3,3) + t262 * t281;
t199 = -rSges(3,3) * t262 + t258 * t281;
t185 = t258 * rSges(4,3) + t262 * t280;
t184 = -rSges(4,3) * t262 + t258 * t280;
t176 = -rSges(5,3) * t250 + (rSges(5,1) * t260 - rSges(5,2) * t256) * t249;
t175 = -Icges(5,5) * t250 + (Icges(5,1) * t260 - Icges(5,4) * t256) * t249;
t174 = -Icges(5,6) * t250 + (Icges(5,4) * t260 - Icges(5,2) * t256) * t249;
t173 = -Icges(5,3) * t250 + (Icges(5,5) * t260 - Icges(5,6) * t256) * t249;
t172 = V_base(5) * rSges(2,3) - t235 * t246 + t288;
t171 = t236 * t246 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t169 = t235 * V_base(4) - t236 * V_base(5) + V_base(3);
t167 = t208 * t306 + t255 * t294;
t166 = t208 * t255 - t262 * t285;
t165 = t206 * t306 + t255 * t295;
t164 = t206 * t255 - t258 * t285;
t163 = t202 * t259 + t254 * t296;
t162 = -t202 * t254 + t259 * t296;
t159 = qJD(6) * t201 + t189;
t158 = t208 * rSges(5,1) - t207 * rSges(5,2) + rSges(5,3) * t294;
t157 = rSges(5,1) * t206 - rSges(5,2) * t205 + rSges(5,3) * t295;
t156 = Icges(5,1) * t208 - Icges(5,4) * t207 + Icges(5,5) * t294;
t155 = Icges(5,1) * t206 - Icges(5,4) * t205 + Icges(5,5) * t295;
t154 = Icges(5,4) * t208 - Icges(5,2) * t207 + Icges(5,6) * t294;
t153 = Icges(5,4) * t206 - Icges(5,2) * t205 + Icges(5,6) * t295;
t152 = Icges(5,5) * t208 - Icges(5,6) * t207 + Icges(5,3) * t294;
t151 = Icges(5,5) * t206 - Icges(5,6) * t205 + Icges(5,3) * t295;
t150 = rSges(6,1) * t202 - rSges(6,2) * t201 + rSges(6,3) * t296;
t149 = Icges(6,1) * t202 - Icges(6,4) * t201 + Icges(6,5) * t296;
t148 = Icges(6,4) * t202 - Icges(6,2) * t201 + Icges(6,6) * t296;
t147 = Icges(6,5) * t202 - Icges(6,6) * t201 + Icges(6,3) * t296;
t146 = t167 * t259 + t207 * t254;
t145 = -t167 * t254 + t207 * t259;
t144 = t165 * t259 + t205 * t254;
t143 = -t165 * t254 + t205 * t259;
t142 = t234 * t241 + (-t199 - t305) * t246 + t288;
t141 = t200 * t246 - t234 * t242 + t273;
t140 = t242 * t199 - t241 * t200 + t272;
t139 = qJD(6) * t166 + t161;
t138 = qJD(6) * t164 + t160;
t137 = rSges(6,1) * t167 - rSges(6,2) * t166 + rSges(6,3) * t207;
t136 = rSges(6,1) * t165 - rSges(6,2) * t164 + rSges(6,3) * t205;
t135 = Icges(6,1) * t167 - Icges(6,4) * t166 + Icges(6,5) * t207;
t134 = Icges(6,1) * t165 - Icges(6,4) * t164 + Icges(6,5) * t205;
t133 = Icges(6,4) * t167 - Icges(6,2) * t166 + Icges(6,6) * t207;
t132 = Icges(6,4) * t165 - Icges(6,2) * t164 + Icges(6,6) * t205;
t131 = Icges(6,5) * t167 - Icges(6,6) * t166 + Icges(6,3) * t207;
t130 = Icges(6,5) * t165 - Icges(6,6) * t164 + Icges(6,3) * t205;
t129 = rSges(7,1) * t163 + rSges(7,2) * t162 + rSges(7,3) * t201;
t128 = Icges(7,1) * t163 + Icges(7,4) * t162 + Icges(7,5) * t201;
t127 = Icges(7,4) * t163 + Icges(7,2) * t162 + Icges(7,6) * t201;
t126 = Icges(7,5) * t163 + Icges(7,6) * t162 + Icges(7,3) * t201;
t125 = t214 * t216 + (-t184 + t284) * t246 + t283;
t124 = t185 * t246 - t214 * t217 + t269;
t123 = t217 * t184 - t216 * t185 + t268;
t122 = rSges(7,1) * t146 + rSges(7,2) * t145 + rSges(7,3) * t166;
t121 = rSges(7,1) * t144 + rSges(7,2) * t143 + rSges(7,3) * t164;
t120 = Icges(7,1) * t146 + Icges(7,4) * t145 + Icges(7,5) * t166;
t119 = Icges(7,1) * t144 + Icges(7,4) * t143 + Icges(7,5) * t164;
t118 = Icges(7,4) * t146 + Icges(7,2) * t145 + Icges(7,6) * t166;
t117 = Icges(7,4) * t144 + Icges(7,2) * t143 + Icges(7,6) * t164;
t116 = Icges(7,5) * t146 + Icges(7,6) * t145 + Icges(7,3) * t166;
t115 = Icges(7,5) * t144 + Icges(7,6) * t143 + Icges(7,3) * t164;
t114 = -t157 * t222 + t176 * t187 + t267;
t113 = t158 * t222 - t176 * t188 + t266;
t112 = t188 * t157 - t187 * t158 + t265;
t111 = -t136 * t189 + t150 * t160 + t267;
t110 = t137 * t189 - t150 * t161 + t266;
t109 = t161 * t136 - t160 * t137 + t265;
t108 = -t121 * t159 + t129 * t138 + (t160 * t201 - t164 * t189) * pkin(6) + t267;
t107 = t122 * t159 - t129 * t139 + (-t161 * t201 + t166 * t189) * pkin(6) + t266;
t106 = t139 * t121 - t138 * t122 + (-t160 * t166 + t161 * t164) * pkin(6) + t265;
t1 = t139 * ((t166 * t116 + t145 * t118 + t146 * t120) * t139 + (t115 * t166 + t117 * t145 + t119 * t146) * t138 + (t126 * t166 + t127 * t145 + t128 * t146) * t159) / 0.2e1 + m(2) * (t169 ^ 2 + t171 ^ 2 + t172 ^ 2) / 0.2e1 + t138 * ((t116 * t164 + t118 * t143 + t120 * t144) * t139 + (t164 * t115 + t143 * t117 + t144 * t119) * t138 + (t126 * t164 + t127 * t143 + t128 * t144) * t159) / 0.2e1 + m(3) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(4) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + ((-t258 * t229 + t232 * t262 + Icges(1,4)) * V_base(5) + (-t258 * t230 + t233 * t262 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t229 * t262 + t258 * t232 + Icges(1,2)) * V_base(5) + (t230 * t262 + t258 * t233 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t188 * ((t152 * t294 - t207 * t154 + t208 * t156) * t188 + (t151 * t294 - t207 * t153 + t208 * t155) * t187 + (t173 * t294 - t207 * t174 + t208 * t175) * t222) / 0.2e1 + t242 * (t258 * t270 + t262 * t263) / 0.2e1 + t241 * (t258 * t263 - t262 * t270) / 0.2e1 + t217 * (t258 * t271 + t262 * t264) / 0.2e1 + t216 * (t258 * t264 - t271 * t262) / 0.2e1 + ((t193 * t261 + t195 * t257) * t242 + (t192 * t261 + t194 * t257) * t241 + (t181 * t250 + t183 * t249) * t217 + (t180 * t250 + t182 * t249) * t216 + (t212 * t250 + t213 * t249 + t228 * t261 + t231 * t257 + Icges(2,3)) * t246) * t246 / 0.2e1 + t246 * V_base(5) * (Icges(2,5) * t258 + Icges(2,6) * t262) + t222 * ((-t151 * t187 - t152 * t188 - t173 * t222) * t250 + ((-t154 * t256 + t156 * t260) * t188 + (-t153 * t256 + t155 * t260) * t187 + (-t174 * t256 + t175 * t260) * t222) * t249) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + t187 * ((t152 * t295 - t154 * t205 + t156 * t206) * t188 + (t151 * t295 - t153 * t205 + t155 * t206) * t187 + (t173 * t295 - t174 * t205 + t175 * t206) * t222) / 0.2e1 + t189 * ((t131 * t296 - t133 * t201 + t135 * t202) * t161 + (t130 * t296 - t132 * t201 + t134 * t202) * t160 + (t147 * t296 - t148 * t201 + t149 * t202) * t189) / 0.2e1 + t246 * V_base(4) * (Icges(2,5) * t262 - Icges(2,6) * t258) + t159 * ((t116 * t201 + t118 * t162 + t120 * t163) * t139 + (t115 * t201 + t117 * t162 + t119 * t163) * t138 + (t201 * t126 + t162 * t127 + t163 * t128) * t159) / 0.2e1 + t160 * ((t131 * t205 - t133 * t164 + t135 * t165) * t161 + (t205 * t130 - t164 * t132 + t165 * t134) * t160 + (t147 * t205 - t148 * t164 + t149 * t165) * t189) / 0.2e1 + t161 * ((t207 * t131 - t166 * t133 + t167 * t135) * t161 + (t130 * t207 - t132 * t166 + t134 * t167) * t160 + (t147 * t207 - t148 * t166 + t149 * t167) * t189) / 0.2e1 + m(1) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1;
T  = t1;
