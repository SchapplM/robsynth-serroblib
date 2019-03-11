% Calculate kinetic energy for
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:33:48
% EndTime: 2019-03-09 03:33:50
% DurationCPUTime: 2.42s
% Computational Cost: add. (2153->319), mult. (1508->445), div. (0->0), fcn. (1288->12), ass. (0->165)
t305 = Icges(4,3) + Icges(5,3);
t228 = qJ(3) + pkin(11);
t217 = sin(t228);
t219 = cos(t228);
t232 = sin(qJ(3));
t235 = cos(qJ(3));
t304 = Icges(4,5) * t235 + Icges(5,5) * t219 - Icges(4,6) * t232 - Icges(5,6) * t217;
t229 = qJ(1) + pkin(10);
t218 = sin(t229);
t220 = cos(t229);
t284 = Icges(5,4) * t219;
t256 = -Icges(5,2) * t217 + t284;
t135 = -Icges(5,6) * t220 + t218 * t256;
t136 = Icges(5,6) * t218 + t220 * t256;
t285 = Icges(5,4) * t217;
t259 = Icges(5,1) * t219 - t285;
t137 = -Icges(5,5) * t220 + t218 * t259;
t138 = Icges(5,5) * t218 + t220 * t259;
t286 = Icges(4,4) * t235;
t257 = -Icges(4,2) * t232 + t286;
t150 = -Icges(4,6) * t220 + t218 * t257;
t151 = Icges(4,6) * t218 + t220 * t257;
t287 = Icges(4,4) * t232;
t260 = Icges(4,1) * t235 - t287;
t152 = -Icges(4,5) * t220 + t218 * t260;
t153 = Icges(4,5) * t218 + t220 * t260;
t178 = Icges(5,2) * t219 + t285;
t181 = Icges(5,1) * t217 + t284;
t194 = -qJD(3) * t220 + V_base(5);
t195 = qJD(3) * t218 + V_base(4);
t199 = Icges(4,2) * t235 + t287;
t202 = Icges(4,1) * t232 + t286;
t222 = V_base(6) + qJD(1);
t301 = (-t178 * t217 + t181 * t219 - t199 * t232 + t202 * t235) * t222 + (-t136 * t217 + t138 * t219 - t151 * t232 + t153 * t235) * t195 + (-t135 * t217 + t137 * t219 - t150 * t232 + t152 * t235) * t194;
t300 = (Icges(4,5) * t232 + Icges(5,5) * t217 + Icges(4,6) * t235 + Icges(5,6) * t219) * t222 + (t305 * t218 + t304 * t220) * t195 + (t304 * t218 - t305 * t220) * t194;
t233 = sin(qJ(1));
t296 = pkin(1) * t233;
t236 = cos(qJ(1));
t295 = pkin(1) * t236;
t294 = pkin(3) * t232;
t293 = pkin(4) * t217;
t292 = t235 * pkin(3);
t291 = -pkin(6) - qJ(2);
t289 = Icges(2,4) * t233;
t288 = Icges(3,4) * t218;
t221 = qJ(5) + t228;
t214 = sin(t221);
t283 = Icges(6,4) * t214;
t215 = cos(t221);
t282 = Icges(6,4) * t215;
t281 = t214 * t218;
t280 = t214 * t220;
t231 = sin(qJ(6));
t279 = t218 * t231;
t234 = cos(qJ(6));
t278 = t218 * t234;
t277 = t220 * t231;
t276 = t220 * t234;
t275 = pkin(4) * t219;
t273 = qJD(6) * t214;
t272 = t222 * t295 + V_base(2);
t271 = V_base(5) * pkin(6) + V_base(1);
t187 = pkin(2) * t218 - pkin(7) * t220;
t268 = -t187 - t296;
t267 = V_base(5) * qJ(2) + t271;
t266 = V_base(4) * t296 + qJD(2) + V_base(3);
t166 = qJD(5) * t218 + t195;
t123 = -qJ(4) * t220 + t218 * t292;
t265 = -t123 + t268;
t264 = pkin(5) * t215 + pkin(9) * t214;
t263 = rSges(4,1) * t235 - rSges(4,2) * t232;
t262 = rSges(5,1) * t219 - rSges(5,2) * t217;
t261 = rSges(6,1) * t215 - rSges(6,2) * t214;
t258 = Icges(6,1) * t215 - t283;
t255 = -Icges(6,2) * t214 + t282;
t252 = Icges(6,5) * t215 - Icges(6,6) * t214;
t117 = -pkin(8) * t220 + t218 * t275;
t251 = -t117 + t265;
t250 = qJD(4) * t218 + t194 * t294 + t267;
t165 = V_base(5) + (-qJD(3) - qJD(5)) * t220;
t249 = t194 * t293 + t250;
t248 = (-Icges(6,3) * t220 + t218 * t252) * t165 + (Icges(6,3) * t218 + t220 * t252) * t166 + (Icges(6,5) * t214 + Icges(6,6) * t215) * t222;
t188 = pkin(2) * t220 + pkin(7) * t218;
t245 = t222 * t188 + t291 * V_base(4) + t272;
t244 = V_base(4) * t187 + (-t188 - t295) * V_base(5) + t266;
t243 = t195 * t123 + t244;
t124 = qJ(4) * t218 + t220 * t292;
t242 = -qJD(4) * t220 + t222 * t124 + t245;
t118 = pkin(8) * t218 + t220 * t275;
t241 = t195 * t117 + (-t118 - t124) * t194 + t243;
t240 = t222 * t118 + (-t293 - t294) * t195 + t242;
t127 = -Icges(6,6) * t220 + t218 * t255;
t128 = Icges(6,6) * t218 + t220 * t255;
t129 = -Icges(6,5) * t220 + t218 * t258;
t130 = Icges(6,5) * t218 + t220 * t258;
t169 = Icges(6,2) * t215 + t283;
t170 = Icges(6,1) * t214 + t282;
t239 = (-t128 * t214 + t130 * t215) * t166 + (-t127 * t214 + t129 * t215) * t165 + (-t169 * t214 + t170 * t215) * t222;
t224 = Icges(2,4) * t236;
t213 = Icges(3,4) * t220;
t207 = rSges(2,1) * t236 - t233 * rSges(2,2);
t206 = t233 * rSges(2,1) + rSges(2,2) * t236;
t205 = rSges(4,1) * t232 + rSges(4,2) * t235;
t204 = Icges(2,1) * t236 - t289;
t203 = Icges(2,1) * t233 + t224;
t201 = -Icges(2,2) * t233 + t224;
t200 = Icges(2,2) * t236 + t289;
t193 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t192 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t191 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t189 = -qJD(6) * t215 + t222;
t186 = rSges(3,1) * t220 - rSges(3,2) * t218;
t185 = rSges(3,1) * t218 + rSges(3,2) * t220;
t184 = rSges(5,1) * t217 + rSges(5,2) * t219;
t183 = Icges(3,1) * t220 - t288;
t182 = Icges(3,1) * t218 + t213;
t180 = -Icges(3,2) * t218 + t213;
t179 = Icges(3,2) * t220 + t288;
t173 = pkin(5) * t214 - pkin(9) * t215;
t172 = rSges(6,1) * t214 + rSges(6,2) * t215;
t163 = t215 * t276 + t279;
t162 = -t215 * t277 + t278;
t161 = t215 * t278 - t277;
t160 = -t215 * t279 - t276;
t159 = t264 * t220;
t158 = t264 * t218;
t157 = rSges(4,3) * t218 + t220 * t263;
t156 = -rSges(4,3) * t220 + t218 * t263;
t155 = V_base(5) * rSges(2,3) - t206 * t222 + t271;
t154 = t207 * t222 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t147 = t206 * V_base(4) - t207 * V_base(5) + V_base(3);
t146 = -rSges(7,3) * t215 + (rSges(7,1) * t234 - rSges(7,2) * t231) * t214;
t145 = -Icges(7,5) * t215 + (Icges(7,1) * t234 - Icges(7,4) * t231) * t214;
t144 = -Icges(7,6) * t215 + (Icges(7,4) * t234 - Icges(7,2) * t231) * t214;
t143 = -Icges(7,3) * t215 + (Icges(7,5) * t234 - Icges(7,6) * t231) * t214;
t142 = rSges(5,3) * t218 + t220 * t262;
t141 = -rSges(5,3) * t220 + t218 * t262;
t140 = t220 * t273 + t166;
t139 = t218 * t273 + t165;
t132 = rSges(6,3) * t218 + t220 * t261;
t131 = -rSges(6,3) * t220 + t218 * t261;
t121 = V_base(5) * rSges(3,3) + (-t185 - t296) * t222 + t267;
t120 = t186 * t222 + (-rSges(3,3) + t291) * V_base(4) + t272;
t116 = V_base(4) * t185 + (-t186 - t295) * V_base(5) + t266;
t113 = rSges(7,1) * t163 + rSges(7,2) * t162 + rSges(7,3) * t280;
t112 = rSges(7,1) * t161 + rSges(7,2) * t160 + rSges(7,3) * t281;
t111 = Icges(7,1) * t163 + Icges(7,4) * t162 + Icges(7,5) * t280;
t110 = Icges(7,1) * t161 + Icges(7,4) * t160 + Icges(7,5) * t281;
t109 = Icges(7,4) * t163 + Icges(7,2) * t162 + Icges(7,6) * t280;
t108 = Icges(7,4) * t161 + Icges(7,2) * t160 + Icges(7,6) * t281;
t107 = Icges(7,5) * t163 + Icges(7,6) * t162 + Icges(7,3) * t280;
t106 = Icges(7,5) * t161 + Icges(7,6) * t160 + Icges(7,3) * t281;
t105 = t194 * t205 + (-t156 + t268) * t222 + t267;
t104 = t157 * t222 - t195 * t205 + t245;
t103 = t195 * t156 - t194 * t157 + t244;
t102 = t184 * t194 + (-t141 + t265) * t222 + t250;
t101 = t142 * t222 + (-t184 - t294) * t195 + t242;
t100 = t195 * t141 + (-t124 - t142) * t194 + t243;
t99 = t165 * t172 + (-t131 + t251) * t222 + t249;
t98 = t132 * t222 - t166 * t172 + t240;
t97 = t166 * t131 - t165 * t132 + t241;
t96 = -t112 * t189 + t139 * t146 + t165 * t173 + (-t158 + t251) * t222 + t249;
t95 = t113 * t189 - t140 * t146 + t159 * t222 - t166 * t173 + t240;
t94 = t140 * t112 - t139 * t113 + t166 * t158 - t165 * t159 + t241;
t1 = t189 * ((-t106 * t139 - t107 * t140 - t143 * t189) * t215 + ((-t109 * t231 + t111 * t234) * t140 + (-t108 * t231 + t110 * t234) * t139 + (-t144 * t231 + t145 * t234) * t189) * t214) / 0.2e1 + t140 * ((t107 * t280 + t109 * t162 + t111 * t163) * t140 + (t106 * t280 + t108 * t162 + t110 * t163) * t139 + (t143 * t280 + t144 * t162 + t145 * t163) * t189) / 0.2e1 + t139 * ((t107 * t281 + t109 * t160 + t111 * t161) * t140 + (t106 * t281 + t108 * t160 + t110 * t161) * t139 + (t143 * t281 + t144 * t160 + t145 * t161) * t189) / 0.2e1 + t166 * (t248 * t218 + t239 * t220) / 0.2e1 + t165 * (t239 * t218 - t248 * t220) / 0.2e1 + m(1) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + m(2) * (t147 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(3) * (t116 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(4) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(7) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + (t301 * t218 - t300 * t220) * t194 / 0.2e1 + (t300 * t218 + t301 * t220) * t195 / 0.2e1 + ((-t179 * t218 + t182 * t220 - t233 * t200 + t203 * t236 + Icges(1,4)) * V_base(5) + (-t180 * t218 + t183 * t220 - t233 * t201 + t204 * t236 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t179 * t220 + t182 * t218 + t200 * t236 + t233 * t203 + Icges(1,2)) * V_base(5) + (t180 * t220 + t183 * t218 + t201 * t236 + t233 * t204 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t128 * t215 + t130 * t214) * t166 + (t127 * t215 + t129 * t214) * t165 + (t136 * t219 + t138 * t217 + t151 * t235 + t153 * t232) * t195 + (t135 * t219 + t137 * t217 + t150 * t235 + t152 * t232) * t194 + (t169 * t215 + t170 * t214 + t178 * t219 + t181 * t217 + t199 * t235 + t202 * t232 + Icges(2,3) + Icges(3,3)) * t222) * t222 / 0.2e1 + t222 * V_base(5) * (Icges(2,5) * t233 + Icges(3,5) * t218 + Icges(2,6) * t236 + Icges(3,6) * t220) + t222 * V_base(4) * (Icges(2,5) * t236 + Icges(3,5) * t220 - Icges(2,6) * t233 - Icges(3,6) * t218) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
