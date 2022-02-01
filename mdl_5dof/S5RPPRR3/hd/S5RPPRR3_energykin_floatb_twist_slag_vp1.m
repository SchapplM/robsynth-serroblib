% Calculate kinetic energy for
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:05
% EndTime: 2022-01-23 09:14:06
% DurationCPUTime: 1.70s
% Computational Cost: add. (1226->249), mult. (928->337), div. (0->0), fcn. (708->10), ass. (0->132)
t170 = qJ(1) + pkin(8);
t159 = sin(t170);
t161 = cos(t170);
t174 = sin(qJ(1));
t175 = cos(qJ(1));
t228 = Icges(2,5) * t174 + Icges(3,5) * t159 + Icges(2,6) * t175 + Icges(3,6) * t161;
t227 = Icges(2,5) * t175 + Icges(3,5) * t161 - Icges(2,6) * t174 - Icges(3,6) * t159;
t225 = pkin(1) * t174;
t224 = pkin(1) * t175;
t171 = sin(pkin(9));
t223 = pkin(3) * t171;
t169 = pkin(9) + qJ(4);
t158 = sin(t169);
t222 = pkin(4) * t158;
t172 = cos(pkin(9));
t221 = t172 * pkin(3);
t220 = -pkin(5) - qJ(2);
t219 = Icges(2,4) * t174;
t218 = Icges(3,4) * t159;
t217 = Icges(4,4) * t171;
t216 = Icges(4,4) * t172;
t215 = Icges(5,4) * t158;
t160 = cos(t169);
t214 = Icges(5,4) * t160;
t162 = qJ(5) + t169;
t155 = sin(t162);
t213 = Icges(6,4) * t155;
t156 = cos(t162);
t212 = Icges(6,4) * t156;
t210 = pkin(4) * t160;
t163 = V_base(6) + qJD(1);
t208 = t163 * t224 + V_base(2);
t207 = V_base(5) * pkin(5) + V_base(1);
t139 = qJD(4) * t159 + V_base(4);
t127 = pkin(2) * t159 - qJ(3) * t161;
t204 = -t127 - t225;
t129 = pkin(2) * t161 + qJ(3) * t159;
t203 = -t129 - t224;
t202 = V_base(5) * qJ(2) + t207;
t201 = V_base(4) * t225 + qJD(2) + V_base(3);
t80 = -pkin(6) * t161 + t221 * t159;
t200 = t204 - t80;
t199 = qJD(3) * t159 + t202;
t198 = V_base(4) * t127 + t201;
t197 = rSges(4,1) * t172 - rSges(4,2) * t171;
t196 = rSges(5,1) * t160 - rSges(5,2) * t158;
t195 = rSges(6,1) * t156 - rSges(6,2) * t155;
t194 = Icges(4,1) * t172 - t217;
t193 = Icges(5,1) * t160 - t215;
t192 = Icges(6,1) * t156 - t213;
t191 = -Icges(4,2) * t171 + t216;
t190 = -Icges(5,2) * t158 + t214;
t189 = -Icges(6,2) * t155 + t212;
t188 = Icges(4,5) * t172 - Icges(4,6) * t171;
t187 = Icges(5,5) * t160 - Icges(5,6) * t158;
t186 = Icges(6,5) * t156 - Icges(6,6) * t155;
t185 = V_base(5) * t223 + t199;
t184 = -qJD(3) * t161 + t163 * t129 + t208;
t110 = V_base(5) + (-qJD(4) - qJD(5)) * t161;
t111 = qJD(5) * t159 + t139;
t183 = t110 * (-Icges(6,3) * t161 + t186 * t159) + t111 * (Icges(6,3) * t159 + t186 * t161) + (Icges(6,5) * t155 + Icges(6,6) * t156) * t163;
t138 = -qJD(4) * t161 + V_base(5);
t182 = (Icges(5,5) * t158 + Icges(5,6) * t160) * t163 + t138 * (-Icges(5,3) * t161 + t187 * t159) + t139 * (Icges(5,3) * t159 + t187 * t161);
t181 = (Icges(4,3) * t159 + t188 * t161) * V_base(4) + (Icges(4,5) * t171 + Icges(4,6) * t172) * t163 + (-Icges(4,3) * t161 + t188 * t159) * V_base(5);
t81 = pkin(6) * t159 + t221 * t161;
t180 = V_base(4) * t80 + (t203 - t81) * V_base(5) + t198;
t179 = t163 * t81 + (t220 - t223) * V_base(4) + t184;
t114 = Icges(6,2) * t156 + t213;
t115 = Icges(6,1) * t155 + t212;
t84 = -Icges(6,6) * t161 + t189 * t159;
t85 = Icges(6,6) * t159 + t189 * t161;
t86 = -Icges(6,5) * t161 + t192 * t159;
t87 = Icges(6,5) * t159 + t192 * t161;
t178 = (-t155 * t85 + t156 * t87) * t111 + (-t155 * t84 + t156 * t86) * t110 + (-t114 * t155 + t115 * t156) * t163;
t120 = Icges(5,2) * t160 + t215;
t123 = Icges(5,1) * t158 + t214;
t92 = -Icges(5,6) * t161 + t190 * t159;
t93 = Icges(5,6) * t159 + t190 * t161;
t94 = -Icges(5,5) * t161 + t193 * t159;
t95 = Icges(5,5) * t159 + t193 * t161;
t177 = (-t158 * t93 + t160 * t95) * t139 + (-t158 * t92 + t160 * t94) * t138 + (-t120 * t158 + t123 * t160) * t163;
t101 = -Icges(4,6) * t161 + t191 * t159;
t102 = Icges(4,6) * t159 + t191 * t161;
t103 = -Icges(4,5) * t161 + t194 * t159;
t104 = Icges(4,5) * t159 + t194 * t161;
t136 = Icges(4,2) * t172 + t217;
t137 = Icges(4,1) * t171 + t216;
t176 = (-t102 * t171 + t104 * t172) * V_base(4) + (-t101 * t171 + t103 * t172) * V_base(5) + (-t136 * t171 + t137 * t172) * t163;
t166 = Icges(2,4) * t175;
t154 = Icges(3,4) * t161;
t148 = rSges(2,1) * t175 - t174 * rSges(2,2);
t147 = t174 * rSges(2,1) + rSges(2,2) * t175;
t146 = Icges(2,1) * t175 - t219;
t145 = Icges(2,1) * t174 + t166;
t144 = -Icges(2,2) * t174 + t166;
t143 = Icges(2,2) * t175 + t219;
t140 = rSges(4,1) * t171 + rSges(4,2) * t172;
t134 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t133 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t132 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t130 = rSges(3,1) * t161 - rSges(3,2) * t159;
t128 = rSges(3,1) * t159 + rSges(3,2) * t161;
t126 = rSges(5,1) * t158 + rSges(5,2) * t160;
t125 = Icges(3,1) * t161 - t218;
t124 = Icges(3,1) * t159 + t154;
t122 = -Icges(3,2) * t159 + t154;
t121 = Icges(3,2) * t161 + t218;
t116 = rSges(6,1) * t155 + rSges(6,2) * t156;
t108 = V_base(5) * rSges(2,3) - t147 * t163 + t207;
t107 = t148 * t163 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t106 = rSges(4,3) * t159 + t197 * t161;
t105 = -rSges(4,3) * t161 + t197 * t159;
t98 = t147 * V_base(4) - t148 * V_base(5) + V_base(3);
t97 = rSges(5,3) * t159 + t196 * t161;
t96 = -rSges(5,3) * t161 + t196 * t159;
t89 = rSges(6,3) * t159 + t195 * t161;
t88 = -rSges(6,3) * t161 + t195 * t159;
t77 = V_base(5) * rSges(3,3) + (-t128 - t225) * t163 + t202;
t76 = t130 * t163 + (-rSges(3,3) + t220) * V_base(4) + t208;
t75 = pkin(7) * t159 + t210 * t161;
t74 = -pkin(7) * t161 + t210 * t159;
t73 = V_base(4) * t128 + (-t130 - t224) * V_base(5) + t201;
t72 = t140 * V_base(5) + (-t105 + t204) * t163 + t199;
t71 = t106 * t163 + (-t140 + t220) * V_base(4) + t184;
t70 = V_base(4) * t105 + (-t106 + t203) * V_base(5) + t198;
t69 = t126 * t138 + (t200 - t96) * t163 + t185;
t68 = -t126 * t139 + t163 * t97 + t179;
t67 = -t138 * t97 + t139 * t96 + t180;
t66 = t138 * t222 + t110 * t116 + (t200 - t74 - t88) * t163 + t185;
t65 = -t139 * t222 - t111 * t116 + (t75 + t89) * t163 + t179;
t64 = -t110 * t89 + t111 * t88 - t138 * t75 + t139 * t74 + t180;
t1 = m(1) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(2) * (t107 ^ 2 + t108 ^ 2 + t98 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + t139 * (t182 * t159 + t177 * t161) / 0.2e1 + t138 * (t177 * t159 - t182 * t161) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + t111 * (t183 * t159 + t178 * t161) / 0.2e1 + t110 * (t178 * t159 - t183 * t161) / 0.2e1 + (t181 * t159 + t176 * t161 + t227 * t163 + (-t121 * t159 + t124 * t161 - t174 * t143 + t145 * t175 + Icges(1,4)) * V_base(5) + (-t122 * t159 + t125 * t161 - t174 * t144 + t146 * t175 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t176 * t159 - t181 * t161 + t228 * t163 + (t121 * t161 + t124 * t159 + t143 * t175 + t174 * t145 + Icges(1,2)) * V_base(5) + (t122 * t161 + t125 * t159 + t144 * t175 + t174 * t146 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t158 * t95 + t160 * t93) * t139 + (t158 * t94 + t160 * t92) * t138 + (t155 * t87 + t156 * t85) * t111 + (t155 * t86 + t156 * t84) * t110 + (t101 * t172 + t103 * t171 + t228) * V_base(5) + (t102 * t172 + t104 * t171 + t227) * V_base(4) + (t114 * t156 + t115 * t155 + t120 * t160 + t123 * t158 + t136 * t172 + t137 * t171 + Icges(2,3) + Icges(3,3)) * t163) * t163 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
