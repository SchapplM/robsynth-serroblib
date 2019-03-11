% Calculate kinetic energy for
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:31
% EndTime: 2019-03-09 05:22:34
% DurationCPUTime: 3.99s
% Computational Cost: add. (1588->353), mult. (2104->499), div. (0->0), fcn. (1998->10), ass. (0->167)
t320 = Icges(2,4) + Icges(3,6);
t319 = Icges(2,1) + Icges(3,2);
t318 = -Icges(3,4) + Icges(2,5);
t317 = Icges(3,5) - Icges(2,6);
t316 = Icges(2,2) + Icges(3,3);
t315 = Icges(6,3) + Icges(5,3);
t245 = cos(qJ(1));
t314 = t320 * t245;
t242 = sin(qJ(1));
t313 = t320 * t242;
t238 = qJ(4) + pkin(10);
t225 = cos(t238);
t241 = sin(qJ(3));
t224 = sin(t238);
t282 = t242 * t224;
t157 = t225 * t245 - t241 * t282;
t281 = t242 * t225;
t158 = t224 * t245 + t241 * t281;
t243 = cos(qJ(4));
t277 = t243 * t245;
t240 = sin(qJ(4));
t280 = t242 * t240;
t177 = -t241 * t280 + t277;
t279 = t242 * t243;
t286 = t240 * t245;
t178 = t241 * t279 + t286;
t244 = cos(qJ(3));
t278 = t242 * t244;
t312 = Icges(5,5) * t178 + Icges(6,5) * t158 + Icges(5,6) * t177 + Icges(6,6) * t157 - t278 * t315;
t285 = t241 * t245;
t159 = t224 * t285 + t281;
t160 = -t225 * t285 + t282;
t179 = t240 * t285 + t279;
t180 = -t241 * t277 + t280;
t276 = t244 * t245;
t311 = Icges(5,5) * t180 + Icges(6,5) * t160 + Icges(5,6) * t179 + Icges(6,6) * t159 + t276 * t315;
t310 = (Icges(5,5) * t243 + Icges(6,5) * t225 - Icges(5,6) * t240 - Icges(6,6) * t224) * t244 + t315 * t241;
t309 = -t245 * t316 - t313;
t308 = t242 * t316 - t314;
t307 = t242 * t319 + t314;
t306 = t245 * t319 - t313;
t275 = pkin(5) * t225;
t303 = -pkin(9) * t244 + t241 * t275;
t293 = t243 * pkin(4);
t302 = -qJ(5) * t244 + t241 * t293;
t290 = Icges(4,4) * t241;
t260 = Icges(4,2) * t244 + t290;
t166 = Icges(4,6) * t245 + t242 * t260;
t167 = Icges(4,6) * t242 - t245 * t260;
t289 = Icges(4,4) * t244;
t261 = Icges(4,1) * t241 + t289;
t169 = Icges(4,5) * t245 + t242 * t261;
t170 = Icges(4,5) * t242 - t245 * t261;
t198 = -Icges(4,2) * t241 + t289;
t203 = Icges(4,1) * t244 - t290;
t216 = qJD(3) * t242 + V_base(5);
t217 = qJD(3) * t245 + V_base(4);
t227 = V_base(6) + qJD(1);
t301 = (t166 * t244 + t169 * t241) * t217 + (t167 * t244 + t170 * t241) * t216 + (t198 * t244 + t203 * t241) * t227;
t295 = pkin(7) * t242;
t294 = pkin(7) * t245;
t226 = qJ(6) + t238;
t221 = sin(t226);
t284 = t242 * t221;
t222 = cos(t226);
t283 = t242 * t222;
t273 = qJD(4) * t244;
t272 = qJD(5) * t244;
t207 = t242 * pkin(1) - qJ(2) * t245;
t271 = V_base(4) * t207 + V_base(3);
t270 = V_base(5) * pkin(6) + V_base(1);
t267 = pkin(5) * t224;
t266 = -t207 - t295;
t265 = qJD(2) * t242 + t270;
t175 = t245 * t273 + t216;
t206 = qJD(4) * t241 + t227;
t264 = V_base(5) * pkin(2) + t265;
t263 = pkin(3) * t241 - pkin(8) * t244;
t262 = rSges(4,1) * t241 + rSges(4,2) * t244;
t259 = Icges(4,5) * t241 + Icges(4,6) * t244;
t211 = pkin(1) * t245 + t242 * qJ(2);
t255 = -qJD(2) * t245 + t227 * t211 + V_base(2);
t254 = (Icges(4,3) * t245 + t242 * t259) * t217 + (Icges(4,3) * t242 - t245 * t259) * t216 + (Icges(4,5) * t244 - Icges(4,6) * t241) * t227;
t253 = V_base(4) * t295 + (-t211 - t294) * V_base(5) + t271;
t252 = t227 * t294 + (-pkin(2) - pkin(6)) * V_base(4) + t255;
t183 = t263 * t245;
t214 = pkin(3) * t244 + pkin(8) * t241;
t251 = t216 * t214 + (t183 + t266) * t227 + t264;
t182 = t263 * t242;
t250 = -t216 * t182 - t217 * t183 + t253;
t249 = t227 * t182 - t217 * t214 + t252;
t133 = pkin(4) * t280 - t302 * t245;
t176 = -t242 * t273 + t217;
t248 = qJD(5) * t241 + t176 * t133 + t250;
t132 = pkin(4) * t286 + t242 * t302;
t247 = t206 * t132 + t245 * t272 + t249;
t145 = qJ(5) * t241 + t244 * t293;
t246 = t175 * t145 - t242 * t272 + t251;
t213 = rSges(2,1) * t245 - t242 * rSges(2,2);
t212 = -rSges(3,2) * t245 + t242 * rSges(3,3);
t210 = rSges(4,1) * t244 - rSges(4,2) * t241;
t209 = t242 * rSges(2,1) + rSges(2,2) * t245;
t208 = -t242 * rSges(3,2) - rSges(3,3) * t245;
t189 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t188 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t187 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t184 = qJD(6) * t241 + t206;
t173 = t242 * rSges(4,3) - t245 * t262;
t172 = rSges(5,3) * t241 + (rSges(5,1) * t243 - rSges(5,2) * t240) * t244;
t171 = rSges(4,3) * t245 + t242 * t262;
t168 = Icges(5,5) * t241 + (Icges(5,1) * t243 - Icges(5,4) * t240) * t244;
t165 = Icges(5,6) * t241 + (Icges(5,4) * t243 - Icges(5,2) * t240) * t244;
t156 = -t222 * t285 + t284;
t155 = t221 * t285 + t283;
t154 = t221 * t245 + t241 * t283;
t153 = t222 * t245 - t241 * t284;
t152 = (-qJD(4) - qJD(6)) * t278 + t217;
t151 = qJD(6) * t276 + t175;
t149 = rSges(6,3) * t241 + (rSges(6,1) * t225 - rSges(6,2) * t224) * t244;
t148 = Icges(6,5) * t241 + (Icges(6,1) * t225 - Icges(6,4) * t224) * t244;
t147 = Icges(6,6) * t241 + (Icges(6,4) * t225 - Icges(6,2) * t224) * t244;
t144 = rSges(7,3) * t241 + (rSges(7,1) * t222 - rSges(7,2) * t221) * t244;
t143 = V_base(5) * rSges(2,3) - t209 * t227 + t270;
t142 = t213 * t227 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t141 = Icges(7,5) * t241 + (Icges(7,1) * t222 - Icges(7,4) * t221) * t244;
t140 = Icges(7,6) * t241 + (Icges(7,4) * t222 - Icges(7,2) * t221) * t244;
t139 = Icges(7,3) * t241 + (Icges(7,5) * t222 - Icges(7,6) * t221) * t244;
t138 = t209 * V_base(4) - t213 * V_base(5) + V_base(3);
t137 = pkin(9) * t241 + t244 * t275;
t135 = t180 * rSges(5,1) + t179 * rSges(5,2) + rSges(5,3) * t276;
t134 = rSges(5,1) * t178 + rSges(5,2) * t177 - rSges(5,3) * t278;
t131 = Icges(5,1) * t180 + Icges(5,4) * t179 + Icges(5,5) * t276;
t130 = Icges(5,1) * t178 + Icges(5,4) * t177 - Icges(5,5) * t278;
t129 = Icges(5,4) * t180 + Icges(5,2) * t179 + Icges(5,6) * t276;
t128 = Icges(5,4) * t178 + Icges(5,2) * t177 - Icges(5,6) * t278;
t125 = V_base(5) * rSges(3,1) + (-t207 - t208) * t227 + t265;
t124 = t227 * t212 + (-rSges(3,1) - pkin(6)) * V_base(4) + t255;
t123 = t160 * rSges(6,1) + t159 * rSges(6,2) + rSges(6,3) * t276;
t122 = rSges(6,1) * t158 + rSges(6,2) * t157 - rSges(6,3) * t278;
t121 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t276;
t120 = Icges(6,1) * t158 + Icges(6,4) * t157 - Icges(6,5) * t278;
t119 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t276;
t118 = Icges(6,4) * t158 + Icges(6,2) * t157 - Icges(6,6) * t278;
t114 = t208 * V_base(4) + (-t211 - t212) * V_base(5) + t271;
t113 = t156 * rSges(7,1) + t155 * rSges(7,2) + rSges(7,3) * t276;
t112 = rSges(7,1) * t154 + rSges(7,2) * t153 - rSges(7,3) * t278;
t111 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t276;
t110 = Icges(7,1) * t154 + Icges(7,4) * t153 - Icges(7,5) * t278;
t109 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t276;
t108 = Icges(7,4) * t154 + Icges(7,2) * t153 - Icges(7,6) * t278;
t107 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t276;
t106 = Icges(7,5) * t154 + Icges(7,6) * t153 - Icges(7,3) * t278;
t104 = t267 * t242 - t303 * t245;
t103 = t242 * t303 + t267 * t245;
t102 = t210 * t216 + (-t173 + t266) * t227 + t264;
t101 = t227 * t171 - t217 * t210 + t252;
t100 = -t216 * t171 + t217 * t173 + t253;
t99 = -t135 * t206 + t172 * t175 + t251;
t98 = t206 * t134 - t176 * t172 + t249;
t97 = -t175 * t134 + t176 * t135 + t250;
t96 = t149 * t175 + (-t123 - t133) * t206 + t246;
t95 = t206 * t122 + (-t145 - t149) * t176 + t247;
t94 = t176 * t123 + (-t122 - t132) * t175 + t248;
t93 = -t113 * t184 + t137 * t175 + t144 * t151 + (-t104 - t133) * t206 + t246;
t92 = t206 * t103 + t184 * t112 - t152 * t144 + (-t137 - t145) * t176 + t247;
t91 = t176 * t104 - t151 * t112 + t152 * t113 + (-t103 - t132) * t175 + t248;
t1 = t217 * (t301 * t242 + t254 * t245) / 0.2e1 + t216 * (t254 * t242 - t301 * t245) / 0.2e1 + t152 * ((-t106 * t278 + t153 * t108 + t154 * t110) * t152 + (-t107 * t278 + t109 * t153 + t111 * t154) * t151 + (-t139 * t278 + t140 * t153 + t141 * t154) * t184) / 0.2e1 + t151 * ((t106 * t276 + t155 * t108 + t156 * t110) * t152 + (t107 * t276 + t155 * t109 + t156 * t111) * t151 + (t139 * t276 + t155 * t140 + t156 * t141) * t184) / 0.2e1 + m(2) * (t138 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(3) * (t114 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(7) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(6) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(5) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(1) * (t187 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + t184 * ((t106 * t152 + t107 * t151 + t139 * t184) * t241 + ((-t108 * t221 + t110 * t222) * t152 + (-t109 * t221 + t111 * t222) * t151 + (-t140 * t221 + t141 * t222) * t184) * t244) / 0.2e1 + ((t159 * t147 + t160 * t148 + t179 * t165 + t180 * t168 + t276 * t310) * t206 + (t159 * t118 + t160 * t120 + t179 * t128 + t180 * t130 + t276 * t312) * t176 + (t159 * t119 + t160 * t121 + t179 * t129 + t180 * t131 + t311 * t276) * t175) * t175 / 0.2e1 + ((t147 * t157 + t148 * t158 + t165 * t177 + t168 * t178 - t278 * t310) * t206 + (t157 * t118 + t158 * t120 + t177 * t128 + t178 * t130 - t312 * t278) * t176 + (t119 * t157 + t121 * t158 + t129 * t177 + t131 * t178 - t278 * t311) * t175) * t176 / 0.2e1 + (((-t147 * t224 + t148 * t225 - t165 * t240 + t168 * t243) * t206 + (-t118 * t224 + t120 * t225 - t128 * t240 + t130 * t243) * t176 + (-t119 * t224 + t121 * t225 - t129 * t240 + t131 * t243) * t175) * t244 + (t311 * t175 + t176 * t312 + t310 * t206) * t241) * t206 / 0.2e1 + ((-t166 * t241 + t169 * t244) * t217 + (-t167 * t241 + t170 * t244) * t216 + (-t241 * t198 + t244 * t203 + Icges(3,1) + Icges(2,3)) * t227) * t227 / 0.2e1 + ((t242 * t309 + t307 * t245 + Icges(1,4)) * V_base(5) + (t242 * t308 + t306 * t245 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t242 * t307 - t309 * t245 + Icges(1,2)) * V_base(5) + (t242 * t306 - t308 * t245 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t227 * (t242 * t318 - t245 * t317) + V_base(4) * t227 * (t242 * t317 + t318 * t245) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
