% Calculate kinetic energy for
% S6RPPRRR8
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
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:20
% EndTime: 2019-03-09 02:35:23
% DurationCPUTime: 3.06s
% Computational Cost: add. (1539->338), mult. (1736->479), div. (0->0), fcn. (1574->10), ass. (0->171)
t308 = Icges(2,4) + Icges(3,6);
t307 = Icges(2,1) + Icges(3,2);
t306 = -Icges(3,4) + Icges(2,5);
t305 = Icges(3,5) - Icges(2,6);
t304 = Icges(2,2) + Icges(3,3);
t230 = cos(qJ(1));
t303 = t308 * t230;
t228 = sin(qJ(1));
t302 = t308 * t228;
t301 = t307 * t228 + t303;
t300 = t307 * t230 - t302;
t299 = t306 * t228 - t305 * t230;
t298 = t305 * t228 + t306 * t230;
t225 = cos(pkin(10));
t224 = sin(pkin(10));
t286 = Icges(4,4) * t224;
t252 = Icges(4,2) * t225 + t286;
t152 = Icges(4,6) * t228 - t230 * t252;
t285 = Icges(4,4) * t225;
t254 = Icges(4,1) * t224 + t285;
t154 = Icges(4,5) * t228 - t230 * t254;
t297 = t152 * t225 + t154 * t224 - t304 * t230 - t302;
t151 = Icges(4,6) * t230 + t228 * t252;
t153 = Icges(4,5) * t230 + t228 * t254;
t296 = t151 * t225 + t153 * t224 + t304 * t228 - t303;
t222 = pkin(10) + qJ(4);
t209 = sin(t222);
t210 = cos(t222);
t229 = cos(qJ(5));
t289 = pkin(5) * t229;
t295 = -pkin(9) * t210 + t209 * t289;
t284 = Icges(5,4) * t209;
t251 = Icges(5,2) * t210 + t284;
t136 = Icges(5,6) * t230 + t228 * t251;
t137 = Icges(5,6) * t228 - t230 * t251;
t283 = Icges(5,4) * t210;
t253 = Icges(5,1) * t209 + t283;
t138 = Icges(5,5) * t230 + t228 * t253;
t139 = Icges(5,5) * t228 - t230 * t253;
t170 = -Icges(5,2) * t209 + t283;
t171 = Icges(5,1) * t210 - t284;
t203 = qJD(4) * t228 + V_base(5);
t204 = qJD(4) * t230 + V_base(4);
t211 = V_base(6) + qJD(1);
t294 = (t136 * t210 + t138 * t209) * t204 + (t137 * t210 + t139 * t209) * t203 + (t170 * t210 + t171 * t209) * t211;
t293 = -pkin(2) - pkin(6);
t291 = pkin(3) * t224;
t290 = pkin(3) * t225;
t280 = qJ(3) * t228;
t279 = qJ(3) * t230;
t278 = t210 * t228;
t277 = t210 * t230;
t223 = qJ(5) + qJ(6);
t217 = sin(t223);
t276 = t217 * t228;
t275 = t217 * t230;
t218 = cos(t223);
t274 = t218 * t228;
t273 = t218 * t230;
t227 = sin(qJ(5));
t272 = t227 * t228;
t271 = t227 * t230;
t270 = t228 * t229;
t269 = t229 * t230;
t267 = qJD(5) * t210;
t196 = pkin(1) * t228 - qJ(2) * t230;
t266 = V_base(4) * t196 + V_base(3);
t265 = V_base(5) * pkin(6) + V_base(1);
t262 = V_base(4) * t280 + t266;
t261 = qJD(2) * t228 + t265;
t260 = -t196 - t280;
t199 = pkin(1) * t230 + qJ(2) * t228;
t259 = -t199 - t279;
t162 = t230 * t267 + t203;
t175 = qJD(5) * t209 + t211;
t258 = pkin(4) * t209 - pkin(8) * t210;
t160 = pkin(7) * t228 - t230 * t291;
t257 = -t160 + t260;
t256 = rSges(4,1) * t224 + rSges(4,2) * t225;
t255 = rSges(5,1) * t209 + rSges(5,2) * t210;
t250 = Icges(4,5) * t224 + Icges(4,6) * t225;
t249 = Icges(5,5) * t209 + Icges(5,6) * t210;
t180 = -Icges(4,2) * t224 + t285;
t181 = Icges(4,1) * t225 - t286;
t243 = t180 * t225 + t181 * t224;
t242 = -qJD(2) * t230 + t211 * t199 + V_base(2);
t241 = V_base(5) * pkin(2) + qJD(3) * t230 + t261;
t240 = V_base(5) * t290 + t241;
t239 = qJD(3) * t228 + t211 * t279 + t242;
t238 = (Icges(5,3) * t230 + t228 * t249) * t204 + (Icges(5,3) * t228 - t230 * t249) * t203 + (Icges(5,5) * t210 - Icges(5,6) * t209) * t211;
t237 = (Icges(4,3) * t230 + t228 * t250) * V_base(4) + (Icges(4,3) * t228 - t230 * t250) * V_base(5) + (Icges(4,5) * t225 - Icges(4,6) * t224) * t211;
t161 = pkin(7) * t230 + t228 * t291;
t236 = V_base(4) * t160 + (-t161 + t259) * V_base(5) + t262;
t235 = t211 * t161 + (-t290 + t293) * V_base(4) + t239;
t158 = t258 * t230;
t173 = t210 * pkin(4) + t209 * pkin(8);
t234 = t203 * t173 + (t158 + t257) * t211 + t240;
t157 = t258 * t228;
t233 = -t157 * t203 - t204 * t158 + t236;
t232 = t211 * t157 - t173 * t204 + t235;
t201 = rSges(2,1) * t230 - rSges(2,2) * t228;
t200 = -rSges(3,2) * t230 + rSges(3,3) * t228;
t198 = rSges(2,1) * t228 + rSges(2,2) * t230;
t197 = -rSges(3,2) * t228 - rSges(3,3) * t230;
t182 = rSges(4,1) * t225 - rSges(4,2) * t224;
t178 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t177 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t176 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t172 = rSges(5,1) * t210 - rSges(5,2) * t209;
t167 = -t209 * t269 + t272;
t166 = t209 * t271 + t270;
t165 = t209 * t270 + t271;
t164 = -t209 * t272 + t269;
t163 = -t228 * t267 + t204;
t159 = qJD(6) * t209 + t175;
t156 = rSges(4,3) * t228 - t230 * t256;
t155 = rSges(4,3) * t230 + t228 * t256;
t148 = -t209 * t273 + t276;
t147 = t209 * t275 + t274;
t146 = t209 * t274 + t275;
t145 = -t209 * t276 + t273;
t142 = rSges(5,3) * t228 - t230 * t255;
t141 = rSges(5,3) * t230 + t228 * t255;
t132 = rSges(6,3) * t209 + (rSges(6,1) * t229 - rSges(6,2) * t227) * t210;
t131 = V_base(5) * rSges(2,3) - t198 * t211 + t265;
t130 = t201 * t211 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t129 = Icges(6,5) * t209 + (Icges(6,1) * t229 - Icges(6,4) * t227) * t210;
t128 = Icges(6,6) * t209 + (Icges(6,4) * t229 - Icges(6,2) * t227) * t210;
t127 = Icges(6,3) * t209 + (Icges(6,5) * t229 - Icges(6,6) * t227) * t210;
t126 = (-qJD(5) - qJD(6)) * t278 + t204;
t125 = qJD(6) * t277 + t162;
t123 = t198 * V_base(4) - t201 * V_base(5) + V_base(3);
t122 = rSges(7,3) * t209 + (rSges(7,1) * t218 - rSges(7,2) * t217) * t210;
t121 = Icges(7,5) * t209 + (Icges(7,1) * t218 - Icges(7,4) * t217) * t210;
t120 = Icges(7,6) * t209 + (Icges(7,4) * t218 - Icges(7,2) * t217) * t210;
t119 = Icges(7,3) * t209 + (Icges(7,5) * t218 - Icges(7,6) * t217) * t210;
t118 = pkin(9) * t209 + t210 * t289;
t117 = V_base(5) * rSges(3,1) + (-t196 - t197) * t211 + t261;
t116 = t200 * t211 + (-rSges(3,1) - pkin(6)) * V_base(4) + t242;
t115 = rSges(6,1) * t167 + rSges(6,2) * t166 + rSges(6,3) * t277;
t114 = rSges(6,1) * t165 + rSges(6,2) * t164 - rSges(6,3) * t278;
t113 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t277;
t112 = Icges(6,1) * t165 + Icges(6,4) * t164 - Icges(6,5) * t278;
t111 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t277;
t110 = Icges(6,4) * t165 + Icges(6,2) * t164 - Icges(6,6) * t278;
t109 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t277;
t108 = Icges(6,5) * t165 + Icges(6,6) * t164 - Icges(6,3) * t278;
t107 = pkin(5) * t272 - t230 * t295;
t106 = pkin(5) * t271 + t228 * t295;
t105 = t197 * V_base(4) + (-t199 - t200) * V_base(5) + t266;
t104 = rSges(7,1) * t148 + rSges(7,2) * t147 + rSges(7,3) * t277;
t103 = rSges(7,1) * t146 + rSges(7,2) * t145 - rSges(7,3) * t278;
t102 = Icges(7,1) * t148 + Icges(7,4) * t147 + Icges(7,5) * t277;
t101 = Icges(7,1) * t146 + Icges(7,4) * t145 - Icges(7,5) * t278;
t100 = Icges(7,4) * t148 + Icges(7,2) * t147 + Icges(7,6) * t277;
t99 = Icges(7,4) * t146 + Icges(7,2) * t145 - Icges(7,6) * t278;
t98 = Icges(7,5) * t148 + Icges(7,6) * t147 + Icges(7,3) * t277;
t97 = Icges(7,5) * t146 + Icges(7,6) * t145 - Icges(7,3) * t278;
t96 = t182 * V_base(5) + (-t156 + t260) * t211 + t241;
t95 = t155 * t211 + (-t182 + t293) * V_base(4) + t239;
t94 = t156 * V_base(4) + (-t155 + t259) * V_base(5) + t262;
t93 = t172 * t203 + (-t142 + t257) * t211 + t240;
t92 = t141 * t211 - t172 * t204 + t235;
t91 = -t141 * t203 + t142 * t204 + t236;
t90 = -t115 * t175 + t132 * t162 + t234;
t89 = t114 * t175 - t132 * t163 + t232;
t88 = -t114 * t162 + t115 * t163 + t233;
t87 = -t104 * t159 - t107 * t175 + t118 * t162 + t122 * t125 + t234;
t86 = t103 * t159 + t106 * t175 - t118 * t163 - t122 * t126 + t232;
t85 = -t103 * t125 + t104 * t126 - t106 * t162 + t107 * t163 + t233;
t1 = t175 * ((t108 * t163 + t109 * t162 + t127 * t175) * t209 + ((-t110 * t227 + t112 * t229) * t163 + (-t111 * t227 + t113 * t229) * t162 + (-t128 * t227 + t129 * t229) * t175) * t210) / 0.2e1 + t159 * ((t119 * t159 + t98 * t125 + t97 * t126) * t209 + ((t101 * t218 - t217 * t99) * t126 + (-t100 * t217 + t102 * t218) * t125 + (-t120 * t217 + t121 * t218) * t159) * t210) / 0.2e1 + m(1) * (t176 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(2) * (t123 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(3) * (t105 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + t204 * (t228 * t294 + t238 * t230) / 0.2e1 + t203 * (t238 * t228 - t294 * t230) / 0.2e1 + t162 * ((t108 * t277 + t110 * t166 + t112 * t167) * t163 + (t109 * t277 + t111 * t166 + t113 * t167) * t162 + (t127 * t277 + t128 * t166 + t129 * t167) * t175) / 0.2e1 + t125 * ((t101 * t148 + t147 * t99 + t277 * t97) * t126 + (t100 * t147 + t102 * t148 + t98 * t277) * t125 + (t119 * t277 + t120 * t147 + t121 * t148) * t159) / 0.2e1 + t163 * ((-t108 * t278 + t110 * t164 + t112 * t165) * t163 + (-t109 * t278 + t111 * t164 + t113 * t165) * t162 + (-t127 * t278 + t128 * t164 + t129 * t165) * t175) / 0.2e1 + t126 * ((t101 * t146 + t145 * t99 - t97 * t278) * t126 + (t100 * t145 + t102 * t146 - t278 * t98) * t125 + (-t119 * t278 + t120 * t145 + t121 * t146) * t159) / 0.2e1 + ((-t136 * t209 + t138 * t210) * t204 + (-t137 * t209 + t139 * t210) * t203 + (-t152 * t224 + t154 * t225 + t299) * V_base(5) + (-t151 * t224 + t153 * t225 + t298) * V_base(4) + (-t170 * t209 + t171 * t210 - t180 * t224 + t181 * t225 + Icges(3,1) + Icges(2,3)) * t211) * t211 / 0.2e1 + (t237 * t230 + (t243 * t228 + t298) * t211 + (t297 * t228 + t230 * t301 + Icges(1,4)) * V_base(5) + (t228 * t296 + t230 * t300 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t237 * t228 + (-t243 * t230 + t299) * t211 + (t228 * t301 - t297 * t230 + Icges(1,2)) * V_base(5) + (t228 * t300 - t230 * t296 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
