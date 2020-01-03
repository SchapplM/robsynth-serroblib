% Calculate kinetic energy for
% S5RPRRR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:15
% EndTime: 2019-12-31 19:07:17
% DurationCPUTime: 1.96s
% Computational Cost: add. (1388->286), mult. (1368->421), div. (0->0), fcn. (1206->10), ass. (0->149)
t196 = sin(pkin(9));
t255 = pkin(2) * t196;
t195 = pkin(9) + qJ(3);
t184 = sin(t195);
t254 = pkin(3) * t184;
t197 = cos(pkin(9));
t253 = t197 * pkin(2);
t200 = sin(qJ(1));
t252 = Icges(2,4) * t200;
t251 = Icges(3,4) * t196;
t250 = Icges(3,4) * t197;
t249 = Icges(4,4) * t184;
t185 = cos(t195);
t248 = Icges(4,4) * t185;
t186 = qJ(4) + t195;
t181 = sin(t186);
t247 = Icges(5,4) * t181;
t182 = cos(t186);
t246 = Icges(5,4) * t182;
t245 = t181 * t200;
t202 = cos(qJ(1));
t244 = t181 * t202;
t199 = sin(qJ(5));
t243 = t199 * t202;
t242 = t200 * t199;
t201 = cos(qJ(5));
t241 = t200 * t201;
t240 = t201 * t202;
t119 = -pkin(6) * t202 + t200 * t253;
t174 = t200 * pkin(1) - qJ(2) * t202;
t238 = -t119 - t174;
t237 = pkin(3) * t185;
t235 = qJD(5) * t181;
t234 = V_base(4) * t174 + V_base(3);
t233 = V_base(5) * pkin(5) + V_base(1);
t100 = -pkin(7) * t202 + t200 * t237;
t230 = -t100 + t238;
t179 = qJD(3) * t200 + V_base(4);
t187 = V_base(6) + qJD(1);
t229 = qJD(2) * t200 + t233;
t158 = qJD(4) * t200 + t179;
t228 = V_base(5) * t255 + t229;
t227 = pkin(4) * t182 + pkin(8) * t181;
t226 = rSges(3,1) * t197 - rSges(3,2) * t196;
t225 = rSges(4,1) * t185 - rSges(4,2) * t184;
t224 = rSges(5,1) * t182 - rSges(5,2) * t181;
t223 = Icges(3,1) * t197 - t251;
t222 = Icges(4,1) * t185 - t249;
t221 = Icges(5,1) * t182 - t247;
t220 = -Icges(3,2) * t196 + t250;
t219 = -Icges(4,2) * t184 + t248;
t218 = -Icges(5,2) * t181 + t246;
t217 = Icges(3,5) * t197 - Icges(3,6) * t196;
t216 = Icges(4,5) * t185 - Icges(4,6) * t184;
t215 = Icges(5,5) * t182 - Icges(5,6) * t181;
t176 = pkin(1) * t202 + t200 * qJ(2);
t214 = -qJD(2) * t202 + t187 * t176 + V_base(2);
t178 = -qJD(3) * t202 + V_base(5);
t213 = t178 * t254 + t228;
t157 = V_base(5) + (-qJD(3) - qJD(4)) * t202;
t212 = (-Icges(5,3) * t202 + t200 * t215) * t157 + (Icges(5,3) * t200 + t202 * t215) * t158 + (Icges(5,5) * t181 + Icges(5,6) * t182) * t187;
t211 = (-Icges(4,3) * t202 + t200 * t216) * t178 + (Icges(4,3) * t200 + t202 * t216) * t179 + (Icges(4,5) * t184 + Icges(4,6) * t185) * t187;
t120 = pkin(6) * t200 + t202 * t253;
t210 = V_base(4) * t119 + (-t120 - t176) * V_base(5) + t234;
t209 = (-Icges(3,3) * t202 + t200 * t217) * V_base(5) + (Icges(3,3) * t200 + t202 * t217) * V_base(4) + (Icges(3,5) * t196 + Icges(3,6) * t197) * t187;
t101 = pkin(7) * t200 + t202 * t237;
t208 = t179 * t100 - t101 * t178 + t210;
t207 = t187 * t120 + (-pkin(5) - t255) * V_base(4) + t214;
t206 = t187 * t101 - t179 * t254 + t207;
t111 = -Icges(5,6) * t202 + t200 * t218;
t112 = Icges(5,6) * t200 + t202 * t218;
t113 = -Icges(5,5) * t202 + t200 * t221;
t114 = Icges(5,5) * t200 + t202 * t221;
t146 = Icges(5,2) * t182 + t247;
t147 = Icges(5,1) * t181 + t246;
t205 = (-t112 * t181 + t114 * t182) * t158 + (-t111 * t181 + t113 * t182) * t157 + (-t146 * t181 + t147 * t182) * t187;
t123 = -Icges(4,6) * t202 + t200 * t219;
t124 = Icges(4,6) * t200 + t202 * t219;
t125 = -Icges(4,5) * t202 + t200 * t222;
t126 = Icges(4,5) * t200 + t202 * t222;
t153 = Icges(4,2) * t185 + t249;
t154 = Icges(4,1) * t184 + t248;
t204 = (-t124 * t184 + t126 * t185) * t179 + (-t123 * t184 + t125 * t185) * t178 + (-t153 * t184 + t154 * t185) * t187;
t135 = -Icges(3,6) * t202 + t200 * t220;
t136 = Icges(3,6) * t200 + t202 * t220;
t137 = -Icges(3,5) * t202 + t200 * t223;
t138 = Icges(3,5) * t200 + t202 * t223;
t165 = Icges(3,2) * t197 + t251;
t166 = Icges(3,1) * t196 + t250;
t203 = (-t136 * t196 + t138 * t197) * V_base(4) + (-t135 * t196 + t137 * t197) * V_base(5) + (-t165 * t196 + t166 * t197) * t187;
t192 = Icges(2,4) * t202;
t177 = rSges(2,1) * t202 - t200 * rSges(2,2);
t175 = t200 * rSges(2,1) + rSges(2,2) * t202;
t173 = Icges(2,1) * t202 - t252;
t172 = Icges(2,1) * t200 + t192;
t171 = -Icges(2,2) * t200 + t192;
t170 = Icges(2,2) * t202 + t252;
t169 = Icges(2,5) * t202 - Icges(2,6) * t200;
t168 = Icges(2,5) * t200 + Icges(2,6) * t202;
t167 = rSges(3,1) * t196 + rSges(3,2) * t197;
t163 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t162 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t161 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t156 = -qJD(5) * t182 + t187;
t155 = rSges(4,1) * t184 + rSges(4,2) * t185;
t149 = pkin(4) * t181 - pkin(8) * t182;
t148 = rSges(5,1) * t181 + rSges(5,2) * t182;
t144 = t182 * t240 + t242;
t143 = -t182 * t243 + t241;
t142 = t182 * t241 - t243;
t141 = -t182 * t242 - t240;
t140 = t200 * rSges(3,3) + t202 * t226;
t139 = -rSges(3,3) * t202 + t200 * t226;
t132 = t227 * t202;
t131 = t227 * t200;
t130 = t200 * rSges(4,3) + t202 * t225;
t129 = -rSges(4,3) * t202 + t200 * t225;
t128 = t202 * t235 + t158;
t127 = t200 * t235 + t157;
t118 = t200 * rSges(5,3) + t202 * t224;
t117 = -rSges(5,3) * t202 + t200 * t224;
t116 = V_base(5) * rSges(2,3) - t175 * t187 + t233;
t115 = t177 * t187 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t108 = t175 * V_base(4) - t177 * V_base(5) + V_base(3);
t106 = -rSges(6,3) * t182 + (rSges(6,1) * t201 - rSges(6,2) * t199) * t181;
t105 = -Icges(6,5) * t182 + (Icges(6,1) * t201 - Icges(6,4) * t199) * t181;
t104 = -Icges(6,6) * t182 + (Icges(6,4) * t201 - Icges(6,2) * t199) * t181;
t103 = -Icges(6,3) * t182 + (Icges(6,5) * t201 - Icges(6,6) * t199) * t181;
t97 = t144 * rSges(6,1) + t143 * rSges(6,2) + rSges(6,3) * t244;
t96 = rSges(6,1) * t142 + rSges(6,2) * t141 + rSges(6,3) * t245;
t95 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t244;
t94 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t245;
t93 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t244;
t92 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t245;
t91 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t244;
t90 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t245;
t89 = t167 * V_base(5) + (-t139 - t174) * t187 + t229;
t88 = t187 * t140 + (-pkin(5) - t167) * V_base(4) + t214;
t87 = t139 * V_base(4) + (-t140 - t176) * V_base(5) + t234;
t86 = t155 * t178 + (-t129 + t238) * t187 + t228;
t85 = t187 * t130 - t179 * t155 + t207;
t84 = t129 * t179 - t130 * t178 + t210;
t83 = t148 * t157 + (-t117 + t230) * t187 + t213;
t82 = t187 * t118 - t158 * t148 + t206;
t81 = t117 * t158 - t118 * t157 + t208;
t80 = t106 * t127 + t149 * t157 - t156 * t96 + (-t131 + t230) * t187 + t213;
t79 = -t128 * t106 + t187 * t132 - t158 * t149 + t156 * t97 + t206;
t78 = -t127 * t97 + t128 * t96 + t131 * t158 - t132 * t157 + t208;
t1 = m(1) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(2) * (t108 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(3) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + t179 * (t211 * t200 + t204 * t202) / 0.2e1 + t178 * (t204 * t200 - t211 * t202) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t158 * (t212 * t200 + t205 * t202) / 0.2e1 + t157 * (t205 * t200 - t212 * t202) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t128 * ((t143 * t93 + t144 * t95 + t91 * t244) * t128 + (t143 * t92 + t144 * t94 + t244 * t90) * t127 + (t103 * t244 + t143 * t104 + t144 * t105) * t156) / 0.2e1 + t127 * ((t141 * t93 + t142 * t95 + t245 * t91) * t128 + (t141 * t92 + t142 * t94 + t90 * t245) * t127 + (t103 * t245 + t104 * t141 + t105 * t142) * t156) / 0.2e1 + t156 * ((-t103 * t156 - t90 * t127 - t91 * t128) * t182 + ((-t199 * t93 + t201 * t95) * t128 + (-t199 * t92 + t201 * t94) * t127 + (-t104 * t199 + t105 * t201) * t156) * t181) / 0.2e1 + (t169 * t187 + t209 * t200 + t203 * t202 + (-t200 * t170 + t172 * t202 + Icges(1,4)) * V_base(5) + (-t200 * t171 + t202 * t173 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t168 * t187 + t203 * t200 - t209 * t202 + (t202 * t170 + t200 * t172 + Icges(1,2)) * V_base(5) + (t171 * t202 + t200 * t173 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t124 * t185 + t126 * t184) * t179 + (t123 * t185 + t125 * t184) * t178 + (t112 * t182 + t114 * t181) * t158 + (t111 * t182 + t113 * t181) * t157 + (t135 * t197 + t137 * t196 + t168) * V_base(5) + (t136 * t197 + t138 * t196 + t169) * V_base(4) + (t182 * t146 + t181 * t147 + t185 * t153 + t184 * t154 + t197 * t165 + t196 * t166 + Icges(2,3)) * t187) * t187 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
