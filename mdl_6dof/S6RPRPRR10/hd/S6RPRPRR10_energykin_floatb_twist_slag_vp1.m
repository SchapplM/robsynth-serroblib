% Calculate kinetic energy for
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:33
% EndTime: 2019-03-09 04:07:37
% DurationCPUTime: 4.05s
% Computational Cost: add. (1561->361), mult. (2059->521), div. (0->0), fcn. (1953->10), ass. (0->169)
t311 = Icges(2,4) + Icges(3,6);
t310 = Icges(2,1) + Icges(3,2);
t309 = -Icges(3,4) + Icges(2,5);
t308 = Icges(3,5) - Icges(2,6);
t307 = Icges(2,2) + Icges(3,3);
t241 = cos(qJ(1));
t306 = t311 * t241;
t239 = sin(qJ(1));
t305 = t311 * t239;
t304 = -t241 * t307 - t305;
t303 = t239 * t307 - t306;
t302 = t239 * t310 + t306;
t301 = t241 * t310 - t305;
t238 = sin(qJ(3));
t240 = cos(qJ(3));
t234 = pkin(10) + qJ(5);
t221 = cos(t234);
t272 = pkin(5) * t221;
t298 = -pkin(9) * t240 + t238 * t272;
t236 = cos(pkin(10));
t289 = t236 * pkin(4);
t297 = -pkin(8) * t240 + t238 * t289;
t287 = Icges(4,4) * t238;
t256 = Icges(4,2) * t240 + t287;
t164 = Icges(4,6) * t241 + t239 * t256;
t165 = Icges(4,6) * t239 - t241 * t256;
t286 = Icges(4,4) * t240;
t257 = Icges(4,1) * t238 + t286;
t166 = Icges(4,5) * t241 + t239 * t257;
t167 = Icges(4,5) * t239 - t241 * t257;
t194 = -Icges(4,2) * t238 + t286;
t199 = Icges(4,1) * t240 - t287;
t212 = qJD(3) * t239 + V_base(5);
t213 = qJD(3) * t241 + V_base(4);
t223 = V_base(6) + qJD(1);
t296 = (t164 * t240 + t166 * t238) * t213 + (t165 * t240 + t167 * t238) * t212 + (t194 * t240 + t199 * t238) * t223;
t291 = pkin(7) * t239;
t290 = pkin(7) * t241;
t235 = sin(pkin(10));
t283 = t235 * t241;
t282 = t238 * t241;
t222 = qJ(6) + t234;
t217 = sin(t222);
t281 = t239 * t217;
t218 = cos(t222);
t280 = t239 * t218;
t220 = sin(t234);
t279 = t239 * t220;
t278 = t239 * t221;
t277 = t239 * t235;
t276 = t239 * t236;
t275 = t239 * t240;
t274 = t240 * t241;
t270 = qJD(4) * t240;
t269 = qJD(5) * t240;
t203 = t239 * pkin(1) - qJ(2) * t241;
t268 = V_base(4) * t203 + V_base(3);
t267 = V_base(5) * pkin(6) + V_base(1);
t264 = pkin(5) * t220;
t263 = -t203 - t291;
t262 = qJD(2) * t239 + t267;
t175 = t241 * t269 + t212;
t202 = qJD(5) * t238 + t223;
t258 = pkin(3) * t238 - qJ(4) * t240;
t179 = t258 * t241;
t261 = t179 + t263;
t260 = V_base(5) * pkin(2) + t262;
t259 = rSges(4,1) * t238 + rSges(4,2) * t240;
t255 = Icges(4,5) * t238 + Icges(4,6) * t240;
t208 = pkin(1) * t241 + t239 * qJ(2);
t251 = -qJD(2) * t241 + t223 * t208 + V_base(2);
t250 = (Icges(4,3) * t241 + t239 * t255) * t213 + (Icges(4,3) * t239 - t241 * t255) * t212 + (Icges(4,5) * t240 - Icges(4,6) * t238) * t223;
t206 = pkin(3) * t240 + qJ(4) * t238;
t249 = t212 * t206 - t239 * t270 + t260;
t248 = V_base(4) * t291 + (-t208 - t290) * V_base(5) + t268;
t247 = t223 * t290 + (-pkin(2) - pkin(6)) * V_base(4) + t251;
t246 = qJD(4) * t238 - t213 * t179 + t248;
t178 = t258 * t239;
t245 = t223 * t178 + t241 * t270 + t247;
t130 = pkin(4) * t283 + t239 * t297;
t131 = pkin(4) * t277 - t241 * t297;
t244 = t213 * t131 + (-t130 - t178) * t212 + t246;
t141 = pkin(8) * t238 + t240 * t289;
t243 = t212 * t141 + (-t131 + t261) * t223 + t249;
t242 = t223 * t130 + (-t141 - t206) * t213 + t245;
t210 = rSges(2,1) * t241 - t239 * rSges(2,2);
t209 = -rSges(3,2) * t241 + t239 * rSges(3,3);
t207 = rSges(4,1) * t240 - rSges(4,2) * t238;
t205 = t239 * rSges(2,1) + rSges(2,2) * t241;
t204 = -t239 * rSges(3,2) - rSges(3,3) * t241;
t186 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t185 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t184 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t180 = qJD(6) * t238 + t202;
t176 = -t239 * t269 + t213;
t174 = -t236 * t282 + t277;
t173 = t235 * t282 + t276;
t172 = t238 * t276 + t283;
t171 = t236 * t241 - t238 * t277;
t169 = t239 * rSges(4,3) - t241 * t259;
t168 = rSges(4,3) * t241 + t239 * t259;
t161 = -t221 * t282 + t279;
t160 = t220 * t282 + t278;
t159 = t220 * t241 + t238 * t278;
t158 = t221 * t241 - t238 * t279;
t157 = rSges(5,3) * t238 + (rSges(5,1) * t236 - rSges(5,2) * t235) * t240;
t155 = Icges(5,5) * t238 + (Icges(5,1) * t236 - Icges(5,4) * t235) * t240;
t154 = Icges(5,6) * t238 + (Icges(5,4) * t236 - Icges(5,2) * t235) * t240;
t153 = Icges(5,3) * t238 + (Icges(5,5) * t236 - Icges(5,6) * t235) * t240;
t152 = -t218 * t282 + t281;
t151 = t217 * t282 + t280;
t150 = t217 * t241 + t238 * t280;
t149 = t218 * t241 - t238 * t281;
t148 = (-qJD(5) - qJD(6)) * t275 + t213;
t147 = qJD(6) * t274 + t175;
t145 = rSges(6,3) * t238 + (rSges(6,1) * t221 - rSges(6,2) * t220) * t240;
t144 = Icges(6,5) * t238 + (Icges(6,1) * t221 - Icges(6,4) * t220) * t240;
t143 = Icges(6,6) * t238 + (Icges(6,4) * t221 - Icges(6,2) * t220) * t240;
t142 = Icges(6,3) * t238 + (Icges(6,5) * t221 - Icges(6,6) * t220) * t240;
t140 = rSges(7,3) * t238 + (rSges(7,1) * t218 - rSges(7,2) * t217) * t240;
t139 = V_base(5) * rSges(2,3) - t205 * t223 + t267;
t138 = t210 * t223 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t137 = Icges(7,5) * t238 + (Icges(7,1) * t218 - Icges(7,4) * t217) * t240;
t136 = Icges(7,6) * t238 + (Icges(7,4) * t218 - Icges(7,2) * t217) * t240;
t135 = Icges(7,3) * t238 + (Icges(7,5) * t218 - Icges(7,6) * t217) * t240;
t134 = t205 * V_base(4) - t210 * V_base(5) + V_base(3);
t132 = pkin(9) * t238 + t240 * t272;
t129 = t174 * rSges(5,1) + t173 * rSges(5,2) + rSges(5,3) * t274;
t128 = rSges(5,1) * t172 + rSges(5,2) * t171 - rSges(5,3) * t275;
t127 = Icges(5,1) * t174 + Icges(5,4) * t173 + Icges(5,5) * t274;
t126 = Icges(5,1) * t172 + Icges(5,4) * t171 - Icges(5,5) * t275;
t125 = Icges(5,4) * t174 + Icges(5,2) * t173 + Icges(5,6) * t274;
t124 = Icges(5,4) * t172 + Icges(5,2) * t171 - Icges(5,6) * t275;
t123 = Icges(5,5) * t174 + Icges(5,6) * t173 + Icges(5,3) * t274;
t122 = Icges(5,5) * t172 + Icges(5,6) * t171 - Icges(5,3) * t275;
t121 = V_base(5) * rSges(3,1) + (-t203 - t204) * t223 + t262;
t120 = t223 * t209 + (-rSges(3,1) - pkin(6)) * V_base(4) + t251;
t118 = t161 * rSges(6,1) + t160 * rSges(6,2) + rSges(6,3) * t274;
t117 = rSges(6,1) * t159 + rSges(6,2) * t158 - rSges(6,3) * t275;
t116 = Icges(6,1) * t161 + Icges(6,4) * t160 + Icges(6,5) * t274;
t115 = Icges(6,1) * t159 + Icges(6,4) * t158 - Icges(6,5) * t275;
t114 = Icges(6,4) * t161 + Icges(6,2) * t160 + Icges(6,6) * t274;
t113 = Icges(6,4) * t159 + Icges(6,2) * t158 - Icges(6,6) * t275;
t112 = Icges(6,5) * t161 + Icges(6,6) * t160 + Icges(6,3) * t274;
t111 = Icges(6,5) * t159 + Icges(6,6) * t158 - Icges(6,3) * t275;
t109 = t204 * V_base(4) + (-t208 - t209) * V_base(5) + t268;
t108 = t152 * rSges(7,1) + t151 * rSges(7,2) + rSges(7,3) * t274;
t107 = rSges(7,1) * t150 + rSges(7,2) * t149 - rSges(7,3) * t275;
t106 = Icges(7,1) * t152 + Icges(7,4) * t151 + Icges(7,5) * t274;
t105 = Icges(7,1) * t150 + Icges(7,4) * t149 - Icges(7,5) * t275;
t104 = Icges(7,4) * t152 + Icges(7,2) * t151 + Icges(7,6) * t274;
t103 = Icges(7,4) * t150 + Icges(7,2) * t149 - Icges(7,6) * t275;
t102 = Icges(7,5) * t152 + Icges(7,6) * t151 + Icges(7,3) * t274;
t101 = Icges(7,5) * t150 + Icges(7,6) * t149 - Icges(7,3) * t275;
t100 = t264 * t239 - t241 * t298;
t99 = t239 * t298 + t264 * t241;
t98 = t207 * t212 + (-t169 + t263) * t223 + t260;
t97 = t223 * t168 - t213 * t207 + t247;
t96 = -t212 * t168 + t213 * t169 + t248;
t95 = t157 * t212 + (-t129 + t261) * t223 + t249;
t94 = t223 * t128 + (-t157 - t206) * t213 + t245;
t93 = t213 * t129 + (-t128 - t178) * t212 + t246;
t92 = -t118 * t202 + t145 * t175 + t243;
t91 = t202 * t117 - t176 * t145 + t242;
t90 = -t175 * t117 + t176 * t118 + t244;
t89 = -t100 * t202 - t108 * t180 + t132 * t175 + t140 * t147 + t243;
t88 = t180 * t107 - t176 * t132 - t148 * t140 + t202 * t99 + t242;
t87 = t176 * t100 - t147 * t107 + t148 * t108 - t175 * t99 + t244;
t1 = t176 * ((-t111 * t275 + t158 * t113 + t159 * t115) * t176 + (-t112 * t275 + t114 * t158 + t116 * t159) * t175 + (-t142 * t275 + t143 * t158 + t144 * t159) * t202) / 0.2e1 + t148 * ((-t101 * t275 + t149 * t103 + t150 * t105) * t148 + (-t102 * t275 + t104 * t149 + t106 * t150) * t147 + (-t135 * t275 + t136 * t149 + t137 * t150) * t180) / 0.2e1 + t175 * ((t111 * t274 + t160 * t113 + t161 * t115) * t176 + (t112 * t274 + t160 * t114 + t161 * t116) * t175 + (t142 * t274 + t160 * t143 + t161 * t144) * t202) / 0.2e1 + t147 * ((t101 * t274 + t151 * t103 + t152 * t105) * t148 + (t102 * t274 + t151 * t104 + t152 * t106) * t147 + (t135 * t274 + t151 * t136 + t152 * t137) * t180) / 0.2e1 + m(3) * (t109 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(5) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(4) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(7) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(6) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + t202 * ((t111 * t176 + t112 * t175 + t142 * t202) * t238 + ((-t113 * t220 + t115 * t221) * t176 + (-t114 * t220 + t116 * t221) * t175 + (-t143 * t220 + t144 * t221) * t202) * t240) / 0.2e1 + t180 * ((t101 * t148 + t102 * t147 + t135 * t180) * t238 + ((-t103 * t217 + t105 * t218) * t148 + (-t104 * t217 + t106 * t218) * t147 + (-t136 * t217 + t137 * t218) * t180) * t240) / 0.2e1 + m(2) * (t134 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(1) * (t184 ^ 2 + t185 ^ 2 + t186 ^ 2) / 0.2e1 + (t250 * t239 - t296 * t241 + (t122 * t274 + t173 * t124 + t174 * t126) * t213 + (t123 * t274 + t173 * t125 + t174 * t127) * t212 + (t153 * t274 + t173 * t154 + t174 * t155) * t223) * t212 / 0.2e1 + ((-t122 * t275 + t124 * t171 + t126 * t172) * t213 + (-t123 * t275 + t125 * t171 + t127 * t172) * t212 + (-t153 * t275 + t154 * t171 + t155 * t172) * t223 + t250 * t241 + t296 * t239) * t213 / 0.2e1 + ((t239 * t304 + t302 * t241 + Icges(1,4)) * V_base(5) + (t239 * t303 + t301 * t241 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t302 * t239 - t304 * t241 + Icges(1,2)) * V_base(5) + (t239 * t301 - t241 * t303 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t122 * t213 + t123 * t212) * t238 + ((-t124 * t235 + t126 * t236) * t213 + (-t125 * t235 + t127 * t236) * t212) * t240 + (-t164 * t238 + t166 * t240) * t213 + (-t165 * t238 + t167 * t240) * t212 + (Icges(2,3) + Icges(3,1) + (-t154 * t235 + t155 * t236 + t199) * t240 + (t153 - t194) * t238) * t223) * t223 / 0.2e1 + t223 * V_base(5) * (t239 * t309 - t241 * t308) + t223 * V_base(4) * (t239 * t308 + t309 * t241) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
