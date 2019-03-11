% Calculate kinetic energy for
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:43:48
% EndTime: 2019-03-09 01:43:51
% DurationCPUTime: 2.50s
% Computational Cost: add. (1502->292), mult. (1268->385), div. (0->0), fcn. (1048->10), ass. (0->150)
t291 = Icges(3,4) + Icges(4,6);
t290 = Icges(3,1) + Icges(4,2);
t289 = -Icges(4,4) + Icges(3,5);
t288 = Icges(4,5) - Icges(3,6);
t287 = Icges(3,2) + Icges(4,3);
t286 = Icges(5,3) + Icges(6,3);
t201 = qJ(4) + pkin(10);
t192 = sin(t201);
t194 = cos(t201);
t205 = sin(qJ(4));
t208 = cos(qJ(4));
t285 = Icges(5,5) * t205 + Icges(6,5) * t192 + Icges(5,6) * t208 + Icges(6,6) * t194;
t202 = qJ(1) + pkin(9);
t195 = cos(t202);
t284 = t291 * t195;
t193 = sin(t202);
t283 = t291 * t193;
t256 = Icges(6,4) * t192;
t229 = Icges(6,2) * t194 + t256;
t106 = Icges(6,6) * t195 + t193 * t229;
t107 = Icges(6,6) * t193 - t195 * t229;
t255 = Icges(6,4) * t194;
t231 = Icges(6,1) * t192 + t255;
t108 = Icges(6,5) * t195 + t193 * t231;
t109 = Icges(6,5) * t193 - t195 * t231;
t258 = Icges(5,4) * t205;
t230 = Icges(5,2) * t208 + t258;
t118 = Icges(5,6) * t195 + t193 * t230;
t119 = Icges(5,6) * t193 - t195 * t230;
t257 = Icges(5,4) * t208;
t232 = Icges(5,1) * t205 + t257;
t121 = Icges(5,5) * t195 + t193 * t232;
t122 = Icges(5,5) * t193 - t195 * t232;
t148 = -Icges(6,2) * t192 + t255;
t153 = Icges(6,1) * t194 - t256;
t169 = qJD(4) * t193 + V_base(5);
t170 = qJD(4) * t195 + V_base(4);
t174 = -Icges(5,2) * t205 + t257;
t177 = Icges(5,1) * t208 - t258;
t196 = V_base(6) + qJD(1);
t282 = t169 * (t107 * t194 + t109 * t192 + t119 * t208 + t122 * t205) + t170 * (t106 * t194 + t108 * t192 + t118 * t208 + t121 * t205) + t196 * (t148 * t194 + t153 * t192 + t174 * t208 + t177 * t205);
t281 = -t195 * t287 - t283;
t280 = t193 * t287 - t284;
t279 = t193 * t290 + t284;
t278 = t195 * t290 - t283;
t277 = (Icges(5,5) * t208 + Icges(6,5) * t194 - Icges(5,6) * t205 - Icges(6,6) * t192) * t196 + (t193 * t285 + t195 * t286) * t170 + (t193 * t286 - t195 * t285) * t169;
t206 = sin(qJ(1));
t268 = pkin(1) * t206;
t209 = cos(qJ(1));
t267 = pkin(1) * t209;
t266 = pkin(4) * t205;
t265 = pkin(4) * t208;
t264 = pkin(7) * t193;
t263 = pkin(7) * t195;
t262 = -pkin(6) - qJ(2);
t260 = Icges(2,4) * t206;
t252 = t193 * t194;
t204 = sin(qJ(6));
t251 = t193 * t204;
t207 = cos(qJ(6));
t250 = t193 * t207;
t249 = t194 * t195;
t248 = t195 * t204;
t247 = t195 * t207;
t246 = qJD(6) * t194;
t245 = t196 * t267 + V_base(2);
t244 = V_base(5) * pkin(6) + V_base(1);
t156 = pkin(2) * t193 - qJ(3) * t195;
t241 = -t156 - t268;
t160 = pkin(2) * t195 + qJ(3) * t193;
t240 = -t160 - t267;
t239 = V_base(5) * qJ(2) + t244;
t238 = t268 * V_base(4) + qJD(2) + V_base(3);
t237 = qJD(3) * t193 + t239;
t236 = pkin(5) * t192 - pkin(8) * t194;
t235 = t156 * V_base(4) + t238;
t234 = rSges(5,1) * t205 + rSges(5,2) * t208;
t233 = rSges(6,1) * t192 + rSges(6,2) * t194;
t220 = V_base(5) * pkin(3) + t237;
t219 = t241 - t264;
t218 = -qJD(3) * t195 + t160 * t196 + t245;
t128 = qJ(5) * t193 - t195 * t266;
t217 = -t128 + t219;
t216 = qJD(5) * t195 + t169 * t265 + t220;
t213 = t196 * t263 + (-pkin(3) + t262) * V_base(4) + t218;
t212 = V_base(4) * t264 + (t240 - t263) * V_base(5) + t235;
t129 = qJ(5) * t195 + t193 * t266;
t211 = qJD(5) * t193 + t129 * t196 + t213;
t210 = t128 * t170 + t212;
t198 = Icges(2,4) * t209;
t182 = rSges(2,1) * t209 - rSges(2,2) * t206;
t181 = rSges(5,1) * t208 - rSges(5,2) * t205;
t180 = rSges(2,1) * t206 + rSges(2,2) * t209;
t179 = Icges(2,1) * t209 - t260;
t178 = Icges(2,1) * t206 + t198;
t176 = -Icges(2,2) * t206 + t198;
t175 = Icges(2,2) * t209 + t260;
t167 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t166 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t165 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t164 = qJD(6) * t192 + t196;
t163 = pkin(5) * t194 + pkin(8) * t192;
t162 = rSges(3,1) * t195 - rSges(3,2) * t193;
t161 = -rSges(4,2) * t195 + rSges(4,3) * t193;
t159 = rSges(6,1) * t194 - rSges(6,2) * t192;
t158 = rSges(3,1) * t193 + rSges(3,2) * t195;
t157 = -rSges(4,2) * t193 - rSges(4,3) * t195;
t137 = -t192 * t247 + t251;
t136 = t192 * t248 + t250;
t135 = t192 * t250 + t248;
t134 = -t192 * t251 + t247;
t133 = -t193 * t246 + t170;
t132 = t195 * t246 + t169;
t131 = t236 * t195;
t130 = t236 * t193;
t127 = rSges(5,3) * t193 - t195 * t234;
t126 = rSges(7,3) * t192 + (rSges(7,1) * t207 - rSges(7,2) * t204) * t194;
t125 = rSges(5,3) * t195 + t193 * t234;
t124 = V_base(5) * rSges(2,3) - t180 * t196 + t244;
t123 = t182 * t196 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t120 = Icges(7,5) * t192 + (Icges(7,1) * t207 - Icges(7,4) * t204) * t194;
t117 = Icges(7,6) * t192 + (Icges(7,4) * t207 - Icges(7,2) * t204) * t194;
t114 = Icges(7,3) * t192 + (Icges(7,5) * t207 - Icges(7,6) * t204) * t194;
t113 = t180 * V_base(4) - t182 * V_base(5) + V_base(3);
t111 = rSges(6,3) * t193 - t195 * t233;
t110 = rSges(6,3) * t195 + t193 * t233;
t102 = V_base(5) * rSges(3,3) + (-t158 - t268) * t196 + t239;
t101 = t162 * t196 + (-rSges(3,3) + t262) * V_base(4) + t245;
t100 = V_base(4) * t158 + (-t162 - t267) * V_base(5) + t238;
t99 = rSges(7,1) * t137 + rSges(7,2) * t136 + rSges(7,3) * t249;
t98 = rSges(7,1) * t135 + rSges(7,2) * t134 - rSges(7,3) * t252;
t97 = Icges(7,1) * t137 + Icges(7,4) * t136 + Icges(7,5) * t249;
t96 = Icges(7,1) * t135 + Icges(7,4) * t134 - Icges(7,5) * t252;
t95 = Icges(7,4) * t137 + Icges(7,2) * t136 + Icges(7,6) * t249;
t94 = Icges(7,4) * t135 + Icges(7,2) * t134 - Icges(7,6) * t252;
t93 = Icges(7,5) * t137 + Icges(7,6) * t136 + Icges(7,3) * t249;
t92 = Icges(7,5) * t135 + Icges(7,6) * t134 - Icges(7,3) * t252;
t91 = V_base(5) * rSges(4,1) + (-t157 + t241) * t196 + t237;
t90 = t161 * t196 + (-rSges(4,1) + t262) * V_base(4) + t218;
t89 = V_base(4) * t157 + (-t161 + t240) * V_base(5) + t235;
t88 = t169 * t181 + (-t127 + t219) * t196 + t220;
t87 = t125 * t196 - t170 * t181 + t213;
t86 = -t125 * t169 + t127 * t170 + t212;
t85 = t159 * t169 + (-t111 + t217) * t196 + t216;
t84 = t110 * t196 + (-t159 - t265) * t170 + t211;
t83 = t170 * t111 + (-t110 - t129) * t169 + t210;
t82 = t126 * t132 + t163 * t169 - t164 * t99 + (t131 + t217) * t196 + t216;
t81 = -t126 * t133 + t130 * t196 + t164 * t98 + (-t163 - t265) * t170 + t211;
t80 = -t170 * t131 - t132 * t98 + t133 * t99 + (-t129 - t130) * t169 + t210;
t1 = t133 * ((t134 * t94 + t135 * t96 - t92 * t252) * t133 + (t134 * t95 + t135 * t97 - t252 * t93) * t132 + (-t114 * t252 + t117 * t134 + t120 * t135) * t164) / 0.2e1 + t132 * ((t136 * t94 + t137 * t96 + t249 * t92) * t133 + (t136 * t95 + t137 * t97 + t93 * t249) * t132 + (t114 * t249 + t117 * t136 + t120 * t137) * t164) / 0.2e1 + m(1) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + t164 * ((t114 * t164 + t132 * t93 + t133 * t92) * t192 + ((-t204 * t94 + t207 * t96) * t133 + (-t204 * t95 + t207 * t97) * t132 + (-t117 * t204 + t120 * t207) * t164) * t194) / 0.2e1 + m(2) * (t113 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(7) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + (t277 * t193 - t282 * t195) * t169 / 0.2e1 + (t282 * t193 + t277 * t195) * t170 / 0.2e1 + ((-t206 * t175 + t178 * t209 + t193 * t281 + t195 * t279 + Icges(1,4)) * V_base(5) + (-t206 * t176 + t209 * t179 + t193 * t280 + t195 * t278 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t209 * t175 + t206 * t178 + t279 * t193 - t195 * t281 + Icges(1,2)) * V_base(5) + (t176 * t209 + t179 * t206 + t193 * t278 - t195 * t280 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t106 * t192 + t108 * t194 - t118 * t205 + t121 * t208) * t170 + (-t107 * t192 + t109 * t194 - t119 * t205 + t122 * t208) * t169 + (-t192 * t148 + t194 * t153 - t205 * t174 + t208 * t177 + Icges(4,1) + Icges(2,3) + Icges(3,3)) * t196) * t196 / 0.2e1 + t196 * V_base(5) * (Icges(2,5) * t206 + Icges(2,6) * t209 + t193 * t289 - t195 * t288) + t196 * V_base(4) * (Icges(2,5) * t209 - Icges(2,6) * t206 + t193 * t288 + t195 * t289) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
