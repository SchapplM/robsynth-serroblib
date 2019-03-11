% Calculate kinetic energy for
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:27
% EndTime: 2019-03-09 04:57:30
% DurationCPUTime: 2.40s
% Computational Cost: add. (2189->319), mult. (1526->445), div. (0->0), fcn. (1306->12), ass. (0->165)
t305 = Icges(5,3) + Icges(6,3);
t229 = qJ(3) + qJ(4);
t219 = pkin(11) + t229;
t214 = sin(t219);
t215 = cos(t219);
t221 = sin(t229);
t222 = cos(t229);
t304 = Icges(5,5) * t222 + Icges(6,5) * t215 - Icges(5,6) * t221 - Icges(6,6) * t214;
t228 = qJ(1) + pkin(10);
t217 = sin(t228);
t218 = cos(t228);
t282 = Icges(6,4) * t215;
t254 = -Icges(6,2) * t214 + t282;
t125 = -Icges(6,6) * t218 + t217 * t254;
t126 = Icges(6,6) * t217 + t218 * t254;
t283 = Icges(6,4) * t214;
t257 = Icges(6,1) * t215 - t283;
t127 = -Icges(6,5) * t218 + t217 * t257;
t128 = Icges(6,5) * t217 + t218 * t257;
t284 = Icges(5,4) * t222;
t255 = -Icges(5,2) * t221 + t284;
t137 = -Icges(5,6) * t218 + t217 * t255;
t138 = Icges(5,6) * t217 + t218 * t255;
t285 = Icges(5,4) * t221;
t258 = Icges(5,1) * t222 - t285;
t139 = -Icges(5,5) * t218 + t217 * t258;
t140 = Icges(5,5) * t217 + t218 * t258;
t166 = V_base(5) + (-qJD(3) - qJD(4)) * t218;
t195 = qJD(3) * t217 + V_base(4);
t167 = qJD(4) * t217 + t195;
t170 = Icges(6,2) * t215 + t283;
t171 = Icges(6,1) * t214 + t282;
t186 = Icges(5,2) * t222 + t285;
t187 = Icges(5,1) * t221 + t284;
t220 = V_base(6) + qJD(1);
t301 = (-t170 * t214 + t171 * t215 - t186 * t221 + t187 * t222) * t220 + (-t126 * t214 + t128 * t215 - t138 * t221 + t140 * t222) * t167 + (-t125 * t214 + t127 * t215 - t137 * t221 + t139 * t222) * t166;
t300 = (Icges(5,5) * t221 + Icges(6,5) * t214 + Icges(5,6) * t222 + Icges(6,6) * t215) * t220 + (t305 * t217 + t304 * t218) * t167 + (t304 * t217 - t305 * t218) * t166;
t232 = sin(qJ(1));
t296 = pkin(1) * t232;
t235 = cos(qJ(1));
t295 = pkin(1) * t235;
t231 = sin(qJ(3));
t294 = pkin(3) * t231;
t293 = pkin(4) * t221;
t234 = cos(qJ(3));
t292 = t234 * pkin(3);
t291 = -pkin(6) - qJ(2);
t289 = Icges(2,4) * t232;
t288 = Icges(3,4) * t217;
t287 = Icges(4,4) * t231;
t286 = Icges(4,4) * t234;
t281 = t214 * t217;
t280 = t214 * t218;
t230 = sin(qJ(6));
t279 = t217 * t230;
t233 = cos(qJ(6));
t278 = t217 * t233;
t277 = t218 * t230;
t276 = t218 * t233;
t275 = pkin(4) * t222;
t273 = qJD(6) * t214;
t272 = t220 * t295 + V_base(2);
t271 = V_base(5) * pkin(6) + V_base(1);
t183 = t217 * pkin(2) - t218 * pkin(7);
t268 = -t183 - t296;
t267 = V_base(5) * qJ(2) + t271;
t266 = V_base(4) * t296 + qJD(2) + V_base(3);
t129 = -pkin(8) * t218 + t217 * t292;
t265 = -t129 + t268;
t194 = -qJD(3) * t218 + V_base(5);
t264 = t194 * t294 + t267;
t263 = pkin(5) * t215 + pkin(9) * t214;
t262 = rSges(4,1) * t234 - rSges(4,2) * t231;
t261 = rSges(5,1) * t222 - rSges(5,2) * t221;
t260 = rSges(6,1) * t215 - rSges(6,2) * t214;
t259 = Icges(4,1) * t234 - t287;
t256 = -Icges(4,2) * t231 + t286;
t253 = Icges(4,5) * t234 - Icges(4,6) * t231;
t117 = -qJ(5) * t218 + t217 * t275;
t250 = -t117 + t265;
t249 = qJD(5) * t217 + t166 * t293 + t264;
t246 = (-Icges(4,3) * t218 + t217 * t253) * t194 + (Icges(4,3) * t217 + t218 * t253) * t195 + (Icges(4,5) * t231 + Icges(4,6) * t234) * t220;
t184 = t218 * pkin(2) + t217 * pkin(7);
t245 = t220 * t184 + t291 * V_base(4) + t272;
t244 = V_base(4) * t183 + (-t184 - t295) * V_base(5) + t266;
t130 = pkin(8) * t217 + t218 * t292;
t243 = t220 * t130 - t195 * t294 + t245;
t242 = t195 * t129 - t130 * t194 + t244;
t241 = t167 * t117 + t242;
t118 = qJ(5) * t217 + t218 * t275;
t240 = -qJD(5) * t218 + t220 * t118 + t243;
t150 = -Icges(4,6) * t218 + t217 * t256;
t151 = Icges(4,6) * t217 + t218 * t256;
t152 = -Icges(4,5) * t218 + t217 * t259;
t153 = Icges(4,5) * t217 + t218 * t259;
t199 = Icges(4,2) * t234 + t287;
t202 = Icges(4,1) * t231 + t286;
t237 = (-t151 * t231 + t153 * t234) * t195 + (-t150 * t231 + t152 * t234) * t194 + (-t199 * t231 + t202 * t234) * t220;
t224 = Icges(2,4) * t235;
t213 = Icges(3,4) * t218;
t207 = rSges(2,1) * t235 - rSges(2,2) * t232;
t206 = rSges(2,1) * t232 + rSges(2,2) * t235;
t205 = rSges(4,1) * t231 + rSges(4,2) * t234;
t204 = Icges(2,1) * t235 - t289;
t203 = Icges(2,1) * t232 + t224;
t201 = -Icges(2,2) * t232 + t224;
t200 = Icges(2,2) * t235 + t289;
t193 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t192 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t191 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t189 = -qJD(6) * t215 + t220;
t188 = rSges(5,1) * t221 + rSges(5,2) * t222;
t182 = rSges(3,1) * t218 - rSges(3,2) * t217;
t181 = rSges(3,1) * t217 + rSges(3,2) * t218;
t180 = Icges(3,1) * t218 - t288;
t179 = Icges(3,1) * t217 + t213;
t178 = -Icges(3,2) * t217 + t213;
t177 = Icges(3,2) * t218 + t288;
t173 = pkin(5) * t214 - pkin(9) * t215;
t172 = rSges(6,1) * t214 + rSges(6,2) * t215;
t164 = t215 * t276 + t279;
t163 = -t215 * t277 + t278;
t162 = t215 * t278 - t277;
t161 = -t215 * t279 - t276;
t159 = t263 * t218;
t158 = t263 * t217;
t157 = rSges(4,3) * t217 + t218 * t262;
t156 = -rSges(4,3) * t218 + t217 * t262;
t155 = V_base(5) * rSges(2,3) - t206 * t220 + t271;
t154 = t207 * t220 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t147 = t206 * V_base(4) - t207 * V_base(5) + V_base(3);
t146 = -rSges(7,3) * t215 + (rSges(7,1) * t233 - rSges(7,2) * t230) * t214;
t145 = rSges(5,3) * t217 + t218 * t261;
t144 = -rSges(5,3) * t218 + t217 * t261;
t143 = -Icges(7,5) * t215 + (Icges(7,1) * t233 - Icges(7,4) * t230) * t214;
t142 = -Icges(7,6) * t215 + (Icges(7,4) * t233 - Icges(7,2) * t230) * t214;
t141 = -Icges(7,3) * t215 + (Icges(7,5) * t233 - Icges(7,6) * t230) * t214;
t134 = t218 * t273 + t167;
t133 = t217 * t273 + t166;
t132 = rSges(6,3) * t217 + t218 * t260;
t131 = -rSges(6,3) * t218 + t217 * t260;
t120 = V_base(5) * rSges(3,3) + (-t181 - t296) * t220 + t267;
t119 = t182 * t220 + (-rSges(3,3) + t291) * V_base(4) + t272;
t115 = t181 * V_base(4) + (-t182 - t295) * V_base(5) + t266;
t114 = rSges(7,1) * t164 + rSges(7,2) * t163 + rSges(7,3) * t280;
t113 = rSges(7,1) * t162 + rSges(7,2) * t161 + rSges(7,3) * t281;
t112 = Icges(7,1) * t164 + Icges(7,4) * t163 + Icges(7,5) * t280;
t111 = Icges(7,1) * t162 + Icges(7,4) * t161 + Icges(7,5) * t281;
t110 = Icges(7,4) * t164 + Icges(7,2) * t163 + Icges(7,6) * t280;
t109 = Icges(7,4) * t162 + Icges(7,2) * t161 + Icges(7,6) * t281;
t108 = Icges(7,5) * t164 + Icges(7,6) * t163 + Icges(7,3) * t280;
t107 = Icges(7,5) * t162 + Icges(7,6) * t161 + Icges(7,3) * t281;
t105 = t194 * t205 + (-t156 + t268) * t220 + t267;
t104 = t157 * t220 - t195 * t205 + t245;
t103 = t156 * t195 - t157 * t194 + t244;
t102 = t166 * t188 + (-t144 + t265) * t220 + t264;
t101 = t145 * t220 - t167 * t188 + t243;
t100 = t144 * t167 - t145 * t166 + t242;
t99 = t166 * t172 + (-t131 + t250) * t220 + t249;
t98 = t132 * t220 + (-t172 - t293) * t167 + t240;
t97 = t131 * t167 + (-t118 - t132) * t166 + t241;
t96 = -t113 * t189 + t133 * t146 + t166 * t173 + (-t158 + t250) * t220 + t249;
t95 = t114 * t189 - t134 * t146 + t159 * t220 + (-t173 - t293) * t167 + t240;
t94 = t113 * t134 - t114 * t133 + t158 * t167 + (-t118 - t159) * t166 + t241;
t1 = t189 * ((-t107 * t133 - t108 * t134 - t141 * t189) * t215 + ((-t110 * t230 + t112 * t233) * t134 + (-t109 * t230 + t111 * t233) * t133 + (-t142 * t230 + t143 * t233) * t189) * t214) / 0.2e1 + m(1) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + m(2) * (t147 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + t195 * (t246 * t217 + t237 * t218) / 0.2e1 + t194 * (t237 * t217 - t246 * t218) / 0.2e1 + m(3) * (t115 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + t134 * ((t108 * t280 + t110 * t163 + t112 * t164) * t134 + (t107 * t280 + t109 * t163 + t111 * t164) * t133 + (t141 * t280 + t142 * t163 + t143 * t164) * t189) / 0.2e1 + t133 * ((t108 * t281 + t110 * t161 + t112 * t162) * t134 + (t107 * t281 + t109 * t161 + t111 * t162) * t133 + (t141 * t281 + t142 * t161 + t143 * t162) * t189) / 0.2e1 + m(4) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(7) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + (t301 * t217 - t300 * t218) * t166 / 0.2e1 + (t217 * t300 + t301 * t218) * t167 / 0.2e1 + ((-t177 * t217 + t179 * t218 - t200 * t232 + t203 * t235 + Icges(1,4)) * V_base(5) + (-t178 * t217 + t180 * t218 - t201 * t232 + t204 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t177 * t218 + t179 * t217 + t200 * t235 + t203 * t232 + Icges(1,2)) * V_base(5) + (t178 * t218 + t180 * t217 + t201 * t235 + t204 * t232 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t151 * t234 + t153 * t231) * t195 + (t150 * t234 + t152 * t231) * t194 + (t126 * t215 + t128 * t214 + t138 * t222 + t140 * t221) * t167 + (t125 * t215 + t127 * t214 + t137 * t222 + t139 * t221) * t166 + (t170 * t215 + t171 * t214 + t186 * t222 + t187 * t221 + t199 * t234 + t202 * t231 + Icges(2,3) + Icges(3,3)) * t220) * t220 / 0.2e1 + t220 * V_base(5) * (Icges(2,5) * t232 + Icges(3,5) * t217 + Icges(2,6) * t235 + Icges(3,6) * t218) + t220 * V_base(4) * (Icges(2,5) * t235 + Icges(3,5) * t218 - Icges(2,6) * t232 - Icges(3,6) * t217) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
