% Calculate kinetic energy for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:33
% EndTime: 2019-12-05 15:12:35
% DurationCPUTime: 2.09s
% Computational Cost: add. (1343->286), mult. (1368->416), div. (0->0), fcn. (1206->10), ass. (0->149)
t196 = sin(pkin(8));
t198 = cos(pkin(8));
t258 = Icges(2,5) * t198 - Icges(2,6) * t196 + Icges(1,5);
t257 = Icges(2,5) * t196 + Icges(2,6) * t198 + Icges(1,6);
t195 = sin(pkin(9));
t256 = pkin(2) * t195;
t194 = pkin(9) + qJ(3);
t187 = sin(t194);
t255 = pkin(3) * t187;
t197 = cos(pkin(9));
t254 = t197 * pkin(2);
t253 = Icges(2,4) * t196;
t252 = Icges(3,4) * t195;
t251 = Icges(3,4) * t197;
t250 = Icges(4,4) * t187;
t188 = cos(t194);
t249 = Icges(4,4) * t188;
t190 = qJ(4) + t194;
t181 = sin(t190);
t248 = Icges(5,4) * t181;
t182 = cos(t190);
t247 = Icges(5,4) * t182;
t246 = t181 * t196;
t245 = t181 * t198;
t200 = sin(qJ(5));
t244 = t196 * t200;
t201 = cos(qJ(5));
t243 = t196 * t201;
t242 = t198 * t200;
t241 = t198 * t201;
t117 = -pkin(5) * t198 + t254 * t196;
t174 = pkin(1) * t196 - qJ(2) * t198;
t239 = -t117 - t174;
t238 = pkin(3) * t188;
t236 = qJD(5) * t181;
t235 = V_base(5) * qJ(1) + V_base(1);
t231 = qJD(1) + V_base(3);
t100 = -pkin(6) * t198 + t238 * t196;
t230 = -t100 + t239;
t179 = qJD(3) * t196 + V_base(4);
t229 = qJD(2) * t196 + t235;
t228 = V_base(4) * t174 + t231;
t156 = qJD(4) * t196 + t179;
t227 = V_base(5) * t256 + t229;
t226 = pkin(4) * t182 + pkin(7) * t181;
t225 = rSges(3,1) * t197 - rSges(3,2) * t195;
t224 = rSges(4,1) * t188 - rSges(4,2) * t187;
t223 = rSges(5,1) * t182 - rSges(5,2) * t181;
t222 = Icges(3,1) * t197 - t252;
t221 = Icges(4,1) * t188 - t250;
t220 = Icges(5,1) * t182 - t248;
t219 = -Icges(3,2) * t195 + t251;
t218 = -Icges(4,2) * t187 + t249;
t217 = -Icges(5,2) * t181 + t247;
t216 = Icges(3,5) * t197 - Icges(3,6) * t195;
t215 = Icges(4,5) * t188 - Icges(4,6) * t187;
t214 = Icges(5,5) * t182 - Icges(5,6) * t181;
t176 = pkin(1) * t198 + qJ(2) * t196;
t213 = -qJD(2) * t198 + V_base(6) * t176 + V_base(2);
t178 = -qJD(3) * t198 + V_base(5);
t212 = t178 * t255 + t227;
t155 = V_base(5) + (-qJD(3) - qJD(4)) * t198;
t211 = (-Icges(5,3) * t198 + t214 * t196) * t155 + (Icges(5,3) * t196 + t214 * t198) * t156 + (Icges(5,5) * t181 + Icges(5,6) * t182) * V_base(6);
t210 = (-Icges(4,3) * t198 + t215 * t196) * t178 + (Icges(4,3) * t196 + t215 * t198) * t179 + (Icges(4,5) * t187 + Icges(4,6) * t188) * V_base(6);
t118 = pkin(5) * t196 + t254 * t198;
t209 = V_base(4) * t117 + (-t118 - t176) * V_base(5) + t228;
t208 = (-Icges(3,3) * t198 + t216 * t196) * V_base(5) + (Icges(3,3) * t196 + t216 * t198) * V_base(4) + (Icges(3,5) * t195 + Icges(3,6) * t197) * V_base(6);
t207 = V_base(6) * t118 + (-qJ(1) - t256) * V_base(4) + t213;
t101 = pkin(6) * t196 + t238 * t198;
t206 = t179 * t100 - t101 * t178 + t209;
t205 = V_base(6) * t101 - t179 * t255 + t207;
t111 = -Icges(5,6) * t198 + t217 * t196;
t112 = Icges(5,6) * t196 + t217 * t198;
t113 = -Icges(5,5) * t198 + t220 * t196;
t114 = Icges(5,5) * t196 + t220 * t198;
t146 = Icges(5,2) * t182 + t248;
t147 = Icges(5,1) * t181 + t247;
t204 = (-t112 * t181 + t114 * t182) * t156 + (-t111 * t181 + t113 * t182) * t155 + (-t146 * t181 + t147 * t182) * V_base(6);
t121 = -Icges(4,6) * t198 + t218 * t196;
t122 = Icges(4,6) * t196 + t218 * t198;
t123 = -Icges(4,5) * t198 + t221 * t196;
t124 = Icges(4,5) * t196 + t221 * t198;
t152 = Icges(4,2) * t188 + t250;
t153 = Icges(4,1) * t187 + t249;
t203 = (-t122 * t187 + t124 * t188) * t179 + (-t121 * t187 + t123 * t188) * t178 + (-t152 * t187 + t153 * t188) * V_base(6);
t135 = -Icges(3,6) * t198 + t219 * t196;
t136 = Icges(3,6) * t196 + t219 * t198;
t137 = -Icges(3,5) * t198 + t222 * t196;
t138 = Icges(3,5) * t196 + t222 * t198;
t167 = Icges(3,2) * t197 + t252;
t170 = Icges(3,1) * t195 + t251;
t202 = (-t136 * t195 + t138 * t197) * V_base(4) + (-t135 * t195 + t137 * t197) * V_base(5) + (-t167 * t195 + t170 * t197) * V_base(6);
t189 = Icges(2,4) * t198;
t177 = rSges(2,1) * t198 - rSges(2,2) * t196;
t175 = rSges(2,1) * t196 + rSges(2,2) * t198;
t173 = rSges(3,1) * t195 + rSges(3,2) * t197;
t172 = Icges(2,1) * t198 - t253;
t171 = Icges(2,1) * t196 + t189;
t169 = -Icges(2,2) * t196 + t189;
t168 = Icges(2,2) * t198 + t253;
t163 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t162 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t161 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t160 = -qJD(5) * t182 + V_base(6);
t154 = rSges(4,1) * t187 + rSges(4,2) * t188;
t149 = pkin(4) * t181 - pkin(7) * t182;
t148 = rSges(5,1) * t181 + rSges(5,2) * t182;
t144 = t182 * t241 + t244;
t143 = -t182 * t242 + t243;
t142 = t182 * t243 - t242;
t141 = -t182 * t244 - t241;
t140 = rSges(3,3) * t196 + t225 * t198;
t139 = -rSges(3,3) * t198 + t225 * t196;
t132 = t226 * t198;
t131 = t226 * t196;
t130 = rSges(4,3) * t196 + t224 * t198;
t129 = -rSges(4,3) * t198 + t224 * t196;
t128 = t198 * t236 + t156;
t127 = t196 * t236 + t155;
t126 = V_base(5) * rSges(2,3) - t175 * V_base(6) + t235;
t125 = t177 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t116 = rSges(5,3) * t196 + t223 * t198;
t115 = -rSges(5,3) * t198 + t223 * t196;
t108 = -t182 * rSges(6,3) + (rSges(6,1) * t201 - rSges(6,2) * t200) * t181;
t105 = -Icges(6,5) * t182 + (Icges(6,1) * t201 - Icges(6,4) * t200) * t181;
t104 = -Icges(6,6) * t182 + (Icges(6,4) * t201 - Icges(6,2) * t200) * t181;
t103 = -Icges(6,3) * t182 + (Icges(6,5) * t201 - Icges(6,6) * t200) * t181;
t102 = t175 * V_base(4) - t177 * V_base(5) + t231;
t97 = rSges(6,1) * t144 + rSges(6,2) * t143 + rSges(6,3) * t245;
t96 = rSges(6,1) * t142 + rSges(6,2) * t141 + rSges(6,3) * t246;
t95 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t245;
t94 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t246;
t93 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t245;
t92 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t246;
t91 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t245;
t90 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t246;
t89 = t173 * V_base(5) + (-t139 - t174) * V_base(6) + t229;
t88 = t140 * V_base(6) + (-qJ(1) - t173) * V_base(4) + t213;
t87 = t139 * V_base(4) + (-t140 - t176) * V_base(5) + t228;
t86 = t154 * t178 + (-t129 + t239) * V_base(6) + t227;
t85 = t130 * V_base(6) - t154 * t179 + t207;
t84 = t129 * t179 - t130 * t178 + t209;
t83 = t148 * t155 + (-t115 + t230) * V_base(6) + t212;
t82 = t116 * V_base(6) - t148 * t156 + t205;
t81 = t115 * t156 - t116 * t155 + t206;
t80 = t108 * t127 + t149 * t155 - t160 * t96 + (-t131 + t230) * V_base(6) + t212;
t79 = -t108 * t128 + t132 * V_base(6) - t149 * t156 + t160 * t97 + t205;
t78 = -t127 * t97 + t128 * t96 + t131 * t156 - t132 * t155 + t206;
t1 = m(1) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(2) * (t102 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + t179 * (t210 * t196 + t203 * t198) / 0.2e1 + t178 * (t203 * t196 - t210 * t198) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t156 * (t211 * t196 + t204 * t198) / 0.2e1 + t155 * (t204 * t196 - t211 * t198) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t128 * ((t143 * t93 + t144 * t95 + t91 * t245) * t128 + (t143 * t92 + t144 * t94 + t90 * t245) * t127 + (t103 * t245 + t104 * t143 + t105 * t144) * t160) / 0.2e1 + t127 * ((t141 * t93 + t142 * t95 + t91 * t246) * t128 + (t141 * t92 + t142 * t94 + t90 * t246) * t127 + (t103 * t246 + t104 * t141 + t105 * t142) * t160) / 0.2e1 + t160 * ((-t103 * t160 - t90 * t127 - t91 * t128) * t182 + ((-t200 * t93 + t201 * t95) * t128 + (-t200 * t92 + t201 * t94) * t127 + (-t104 * t200 + t105 * t201) * t160) * t181) / 0.2e1 + (t208 * t196 + t202 * t198 + t258 * V_base(6) + (-t168 * t196 + t171 * t198 + Icges(1,4)) * V_base(5) + (-t196 * t169 + t198 * t172 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t202 * t196 - t208 * t198 + t257 * V_base(6) + (t198 * t168 + t196 * t171 + Icges(1,2)) * V_base(5) + (t169 * t198 + t172 * t196 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t122 * t188 + t124 * t187) * t179 + (t121 * t188 + t123 * t187) * t178 + (t112 * t182 + t114 * t181) * t156 + (t111 * t182 + t113 * t181) * t155 + (t135 * t197 + t137 * t195 + t257) * V_base(5) + (t136 * t197 + t138 * t195 + t258) * V_base(4) + (t182 * t146 + t181 * t147 + t188 * t152 + t187 * t153 + t197 * t167 + t195 * t170 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1;
T = t1;
