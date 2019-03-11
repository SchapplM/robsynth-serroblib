% Calculate kinetic energy for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:18
% EndTime: 2019-03-09 02:25:22
% DurationCPUTime: 3.41s
% Computational Cost: add. (1568->320), mult. (2668->457), div. (0->0), fcn. (3058->10), ass. (0->156)
t306 = Icges(2,4) - Icges(3,5);
t305 = Icges(2,1) + Icges(3,1);
t304 = Icges(3,4) + Icges(2,5);
t303 = Icges(2,2) + Icges(3,3);
t302 = Icges(2,6) - Icges(3,6);
t291 = sin(qJ(1));
t301 = t306 * t291;
t292 = cos(qJ(1));
t300 = t306 * t292;
t299 = -t303 * t292 - t301;
t298 = t303 * t291 - t300;
t297 = t305 * t291 + t300;
t296 = t305 * t292 - t301;
t242 = cos(qJ(5));
t290 = pkin(5) * t242;
t288 = cos(pkin(10));
t287 = sin(pkin(10));
t194 = -t287 * t291 - t288 * t292;
t286 = Icges(4,4) * t194;
t241 = sin(qJ(4));
t285 = Icges(5,4) * t241;
t243 = cos(qJ(4));
t284 = Icges(5,4) * t243;
t240 = sin(qJ(5));
t283 = t194 * t240;
t282 = t194 * t241;
t195 = t287 * t292 - t288 * t291;
t281 = t195 * t240;
t280 = t195 * t241;
t239 = qJ(5) + qJ(6);
t233 = sin(t239);
t279 = t233 * t243;
t234 = cos(t239);
t278 = t234 * t243;
t277 = t240 * t243;
t276 = t242 * t243;
t275 = qJD(5) * t241;
t274 = qJD(6) * t241;
t217 = pkin(1) * t291 - qJ(2) * t292;
t273 = V_base(4) * t217 + V_base(3);
t272 = V_base(5) * pkin(6) + V_base(1);
t269 = t292 * pkin(2);
t268 = t291 * pkin(2);
t186 = qJD(4) * t195 + V_base(4);
t185 = -qJD(4) * t194 + V_base(5);
t230 = V_base(6) + qJD(1);
t265 = qJD(2) * t291 + t272;
t153 = -t194 * t275 + t186;
t152 = -t195 * t275 + t185;
t215 = qJD(5) * t243 + t230;
t264 = -t217 - t268;
t220 = pkin(1) * t292 + qJ(2) * t291;
t263 = -t220 - t269;
t262 = -pkin(4) * t243 - pkin(8) * t241;
t261 = V_base(4) * t268 - qJD(3) + t273;
t260 = -rSges(5,1) * t243 + rSges(5,2) * t241;
t259 = -Icges(5,1) * t243 + t285;
t258 = Icges(5,2) * t241 - t284;
t257 = -Icges(5,5) * t243 + Icges(5,6) * t241;
t170 = -pkin(3) * t195 - pkin(7) * t194;
t256 = -t170 + t264;
t255 = -qJD(2) * t292 + t230 * t220 + V_base(2);
t254 = -V_base(5) * qJ(3) + t265;
t253 = -pkin(9) * t241 - t243 * t290;
t252 = (-Icges(5,3) * t194 + t195 * t257) * t185 + (Icges(5,3) * t195 + t194 * t257) * t186 + (-Icges(5,5) * t241 - Icges(5,6) * t243) * t230;
t251 = V_base(4) * qJ(3) + t230 * t269 + t255;
t171 = -pkin(3) * t194 + pkin(7) * t195;
t250 = -V_base(4) * pkin(6) + t230 * t171 + t251;
t249 = V_base(4) * t170 + (-t171 + t263) * V_base(5) + t261;
t160 = t262 * t194;
t223 = -t241 * pkin(4) + t243 * pkin(8);
t248 = t230 * t160 - t186 * t223 + t250;
t159 = t262 * t195;
t247 = t185 * t223 + (-t159 + t256) * t230 + t254;
t246 = t186 * t159 - t185 * t160 + t249;
t139 = -Icges(5,6) * t194 + t195 * t258;
t140 = Icges(5,6) * t195 + t194 * t258;
t141 = -Icges(5,5) * t194 + t195 * t259;
t142 = Icges(5,5) * t195 + t194 * t259;
t205 = -Icges(5,2) * t243 - t285;
t210 = -Icges(5,1) * t241 - t284;
t245 = (t140 * t241 - t142 * t243) * t186 + (t139 * t241 - t141 * t243) * t185 + (t205 * t241 - t210 * t243) * t230;
t222 = rSges(2,1) * t292 - rSges(2,2) * t291;
t221 = rSges(3,1) * t292 + rSges(3,3) * t291;
t219 = rSges(2,1) * t291 + rSges(2,2) * t292;
t218 = rSges(3,1) * t291 - rSges(3,3) * t292;
t216 = -rSges(5,1) * t241 - rSges(5,2) * t243;
t199 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t198 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t197 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t193 = qJD(6) * t243 + t215;
t191 = Icges(4,4) * t195;
t184 = rSges(6,3) * t243 + (-rSges(6,1) * t242 + rSges(6,2) * t240) * t241;
t183 = Icges(6,5) * t243 + (-Icges(6,1) * t242 + Icges(6,4) * t240) * t241;
t182 = Icges(6,6) * t243 + (-Icges(6,4) * t242 + Icges(6,2) * t240) * t241;
t181 = Icges(6,3) * t243 + (-Icges(6,5) * t242 + Icges(6,6) * t240) * t241;
t180 = rSges(7,3) * t243 + (-rSges(7,1) * t234 + rSges(7,2) * t233) * t241;
t179 = Icges(7,5) * t243 + (-Icges(7,1) * t234 + Icges(7,4) * t233) * t241;
t178 = Icges(7,6) * t243 + (-Icges(7,4) * t234 + Icges(7,2) * t233) * t241;
t177 = Icges(7,3) * t243 + (-Icges(7,5) * t234 + Icges(7,6) * t233) * t241;
t176 = pkin(9) * t243 - t241 * t290;
t175 = V_base(5) * rSges(2,3) - t219 * t230 + t272;
t174 = t222 * t230 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t173 = t219 * V_base(4) - t222 * V_base(5) + V_base(3);
t169 = -rSges(4,1) * t194 - rSges(4,2) * t195;
t168 = -rSges(4,1) * t195 + rSges(4,2) * t194;
t167 = -Icges(4,1) * t194 - t191;
t166 = -Icges(4,1) * t195 + t286;
t165 = -Icges(4,2) * t195 - t286;
t164 = Icges(4,2) * t194 - t191;
t157 = -t194 * t276 + t281;
t156 = t194 * t277 + t195 * t242;
t155 = -t195 * t276 - t283;
t154 = -t194 * t242 + t195 * t277;
t151 = -t194 * t278 + t195 * t233;
t150 = t194 * t279 + t195 * t234;
t149 = -t194 * t233 - t195 * t278;
t148 = -t194 * t234 + t195 * t279;
t146 = V_base(5) * rSges(3,2) + (-t217 - t218) * t230 + t265;
t145 = t230 * t221 + (-rSges(3,2) - pkin(6)) * V_base(4) + t255;
t144 = rSges(5,3) * t195 + t194 * t260;
t143 = -rSges(5,3) * t194 + t195 * t260;
t136 = t218 * V_base(4) + (-t220 - t221) * V_base(5) + t273;
t135 = -t194 * t274 + t153;
t134 = -t195 * t274 + t152;
t132 = (-qJ(3) - rSges(4,3)) * V_base(5) + (-t168 + t264) * t230 + t265;
t131 = t230 * t169 + (rSges(4,3) - pkin(6)) * V_base(4) + t251;
t130 = pkin(5) * t281 + t194 * t253;
t129 = -pkin(5) * t283 + t195 * t253;
t128 = rSges(6,1) * t157 + rSges(6,2) * t156 - rSges(6,3) * t282;
t127 = rSges(6,1) * t155 + rSges(6,2) * t154 - rSges(6,3) * t280;
t126 = Icges(6,1) * t157 + Icges(6,4) * t156 - Icges(6,5) * t282;
t125 = Icges(6,1) * t155 + Icges(6,4) * t154 - Icges(6,5) * t280;
t124 = Icges(6,4) * t157 + Icges(6,2) * t156 - Icges(6,6) * t282;
t123 = Icges(6,4) * t155 + Icges(6,2) * t154 - Icges(6,6) * t280;
t122 = Icges(6,5) * t157 + Icges(6,6) * t156 - Icges(6,3) * t282;
t121 = Icges(6,5) * t155 + Icges(6,6) * t154 - Icges(6,3) * t280;
t120 = rSges(7,1) * t151 + rSges(7,2) * t150 - rSges(7,3) * t282;
t119 = rSges(7,1) * t149 + rSges(7,2) * t148 - rSges(7,3) * t280;
t118 = Icges(7,1) * t151 + Icges(7,4) * t150 - Icges(7,5) * t282;
t117 = Icges(7,1) * t149 + Icges(7,4) * t148 - Icges(7,5) * t280;
t116 = Icges(7,4) * t151 + Icges(7,2) * t150 - Icges(7,6) * t282;
t115 = Icges(7,4) * t149 + Icges(7,2) * t148 - Icges(7,6) * t280;
t114 = Icges(7,5) * t151 + Icges(7,6) * t150 - Icges(7,3) * t282;
t113 = Icges(7,5) * t149 + Icges(7,6) * t148 - Icges(7,3) * t280;
t112 = V_base(4) * t168 + (-t169 + t263) * V_base(5) + t261;
t111 = t185 * t216 + (-t143 + t256) * t230 + t254;
t110 = t230 * t144 - t186 * t216 + t250;
t109 = t186 * t143 - t185 * t144 + t249;
t108 = -t215 * t127 + t152 * t184 + t247;
t107 = t215 * t128 - t153 * t184 + t248;
t106 = t153 * t127 - t152 * t128 + t246;
t105 = -t193 * t119 - t215 * t129 + t134 * t180 + t152 * t176 + t247;
t104 = t193 * t120 + t215 * t130 - t135 * t180 - t153 * t176 + t248;
t103 = t135 * t119 - t134 * t120 + t153 * t129 - t152 * t130 + t246;
t1 = m(1) * (t197 ^ 2 + t198 ^ 2 + t199 ^ 2) / 0.2e1 + m(2) * (t173 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(3) * (t136 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + t186 * (t194 * t245 + t195 * t252) / 0.2e1 + t185 * (-t194 * t252 + t245 * t195) / 0.2e1 + t153 * ((-t122 * t282 + t156 * t124 + t157 * t126) * t153 + (-t121 * t282 + t123 * t156 + t125 * t157) * t152 + (t156 * t182 + t157 * t183 - t181 * t282) * t215) / 0.2e1 + t135 * ((-t114 * t282 + t150 * t116 + t151 * t118) * t135 + (-t113 * t282 + t115 * t150 + t117 * t151) * t134 + (t150 * t178 + t151 * t179 - t177 * t282) * t193) / 0.2e1 + t152 * ((-t122 * t280 + t124 * t154 + t126 * t155) * t153 + (-t121 * t280 + t154 * t123 + t155 * t125) * t152 + (t154 * t182 + t155 * t183 - t181 * t280) * t215) / 0.2e1 + t134 * ((-t114 * t280 + t116 * t148 + t118 * t149) * t135 + (-t113 * t280 + t148 * t115 + t149 * t117) * t134 + (t148 * t178 + t149 * t179 - t177 * t280) * t193) / 0.2e1 + t215 * ((t121 * t152 + t122 * t153 + t181 * t215) * t243 + ((t124 * t240 - t126 * t242) * t153 + (t123 * t240 - t125 * t242) * t152 + (t182 * t240 - t183 * t242) * t215) * t241) / 0.2e1 + t193 * ((t113 * t134 + t114 * t135 + t177 * t193) * t243 + ((t116 * t233 - t118 * t234) * t135 + (t115 * t233 - t117 * t234) * t134 + (t178 * t233 - t179 * t234) * t193) * t241) / 0.2e1 + ((-t140 * t243 - t142 * t241) * t186 + (-t139 * t243 - t141 * t241) * t185 + (-t243 * t205 - t241 * t210 + Icges(3,2) + Icges(2,3) + Icges(4,3)) * t230) * t230 / 0.2e1 + ((-t164 * t195 - t166 * t194 + t291 * t299 + t297 * t292 + Icges(1,4)) * V_base(5) + (-t195 * t165 - t194 * t167 + t298 * t291 + t296 * t292 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t194 * t164 - t195 * t166 + t297 * t291 - t299 * t292 + Icges(1,2)) * V_base(5) + (t165 * t194 - t167 * t195 + t291 * t296 - t292 * t298 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t230 * (Icges(4,5) * t195 - Icges(4,6) * t194 + t304 * t291 + t302 * t292) + V_base(4) * t230 * (Icges(4,5) * t194 + Icges(4,6) * t195 - t302 * t291 + t304 * t292) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
