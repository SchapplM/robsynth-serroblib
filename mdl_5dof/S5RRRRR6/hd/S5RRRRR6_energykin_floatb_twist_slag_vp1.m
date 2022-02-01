% Calculate kinetic energy for
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:36
% EndTime: 2022-01-20 12:07:37
% DurationCPUTime: 1.42s
% Computational Cost: add. (1342->245), mult. (970->351), div. (0->0), fcn. (750->10), ass. (0->131)
t227 = -pkin(5) - pkin(6);
t175 = sin(qJ(1));
t225 = pkin(1) * t175;
t177 = cos(qJ(1));
t224 = pkin(1) * t177;
t174 = sin(qJ(3));
t223 = pkin(3) * t174;
t172 = qJ(3) + qJ(4);
t162 = sin(t172);
t222 = pkin(4) * t162;
t176 = cos(qJ(3));
t221 = t176 * pkin(3);
t173 = qJ(1) + qJ(2);
t163 = sin(t173);
t165 = cos(t173);
t132 = t163 * pkin(2) - t165 * pkin(7);
t80 = -pkin(8) * t165 + t221 * t163;
t219 = -t132 - t80;
t218 = Icges(2,4) * t175;
t217 = Icges(3,4) * t163;
t216 = Icges(4,4) * t174;
t215 = Icges(4,4) * t176;
t214 = Icges(5,4) * t162;
t164 = cos(t172);
t213 = Icges(5,4) * t164;
t167 = qJ(5) + t172;
t158 = sin(t167);
t212 = Icges(6,4) * t158;
t159 = cos(t167);
t211 = Icges(6,4) * t159;
t210 = pkin(4) * t164;
t208 = -qJD(3) - qJD(4);
t161 = V_base(6) + qJD(1);
t207 = t161 * t224 + V_base(2);
t206 = V_base(4) * t225 + V_base(3);
t205 = V_base(5) * pkin(5) + V_base(1);
t139 = qJD(3) * t163 + V_base(4);
t113 = qJD(4) * t163 + t139;
t202 = rSges(4,1) * t176 - rSges(4,2) * t174;
t201 = rSges(5,1) * t164 - rSges(5,2) * t162;
t200 = rSges(6,1) * t159 - rSges(6,2) * t158;
t199 = Icges(4,1) * t176 - t216;
t198 = Icges(5,1) * t164 - t214;
t197 = Icges(6,1) * t159 - t212;
t196 = -Icges(4,2) * t174 + t215;
t195 = -Icges(5,2) * t162 + t213;
t194 = -Icges(6,2) * t158 + t211;
t193 = Icges(4,5) * t176 - Icges(4,6) * t174;
t192 = Icges(5,5) * t164 - Icges(5,6) * t162;
t191 = Icges(6,5) * t159 - Icges(6,6) * t158;
t190 = V_base(5) * pkin(6) - t161 * t225 + t205;
t109 = V_base(5) + (-qJD(5) + t208) * t165;
t110 = qJD(5) * t163 + t113;
t157 = qJD(2) + t161;
t189 = t109 * (-Icges(6,3) * t165 + t191 * t163) + t110 * (Icges(6,3) * t163 + t191 * t165) + (Icges(6,5) * t158 + Icges(6,6) * t159) * t157;
t112 = t208 * t165 + V_base(5);
t188 = t112 * (-Icges(5,3) * t165 + t192 * t163) + t113 * (Icges(5,3) * t163 + t192 * t165) + (Icges(5,5) * t162 + Icges(5,6) * t164) * t157;
t138 = -qJD(3) * t165 + V_base(5);
t187 = (-Icges(4,3) * t165 + t193 * t163) * t138 + (Icges(4,3) * t163 + t193 * t165) * t139 + (Icges(4,5) * t174 + Icges(4,6) * t176) * t157;
t186 = t138 * t223 + t190;
t133 = t165 * pkin(2) + t163 * pkin(7);
t185 = t157 * t133 + t227 * V_base(4) + t207;
t184 = V_base(4) * t132 + (-t133 - t224) * V_base(5) + t206;
t81 = pkin(8) * t163 + t221 * t165;
t183 = -t139 * t223 + t157 * t81 + t185;
t182 = -t138 * t81 + t139 * t80 + t184;
t116 = Icges(6,2) * t159 + t212;
t117 = Icges(6,1) * t158 + t211;
t84 = -Icges(6,6) * t165 + t194 * t163;
t85 = Icges(6,6) * t163 + t194 * t165;
t86 = -Icges(6,5) * t165 + t197 * t163;
t87 = Icges(6,5) * t163 + t197 * t165;
t181 = (-t158 * t85 + t159 * t87) * t110 + (-t158 * t84 + t159 * t86) * t109 + (-t116 * t158 + t117 * t159) * t157;
t123 = Icges(5,2) * t164 + t214;
t126 = Icges(5,1) * t162 + t213;
t92 = -Icges(5,6) * t165 + t195 * t163;
t93 = Icges(5,6) * t163 + t195 * t165;
t94 = -Icges(5,5) * t165 + t198 * t163;
t95 = Icges(5,5) * t163 + t198 * t165;
t180 = (-t162 * t93 + t164 * t95) * t113 + (-t162 * t92 + t164 * t94) * t112 + (-t123 * t162 + t126 * t164) * t157;
t103 = -Icges(4,6) * t165 + t196 * t163;
t104 = Icges(4,6) * t163 + t196 * t165;
t105 = -Icges(4,5) * t165 + t199 * t163;
t106 = Icges(4,5) * t163 + t199 * t165;
t143 = Icges(4,2) * t176 + t216;
t146 = Icges(4,1) * t174 + t215;
t179 = (-t104 * t174 + t106 * t176) * t139 + (-t103 * t174 + t105 * t176) * t138 + (-t143 * t174 + t146 * t176) * t157;
t166 = Icges(2,4) * t177;
t156 = Icges(3,4) * t165;
t151 = rSges(2,1) * t177 - rSges(2,2) * t175;
t150 = rSges(2,1) * t175 + rSges(2,2) * t177;
t149 = rSges(4,1) * t174 + rSges(4,2) * t176;
t148 = Icges(2,1) * t177 - t218;
t147 = Icges(2,1) * t175 + t166;
t145 = -Icges(2,2) * t175 + t166;
t144 = Icges(2,2) * t177 + t218;
t137 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t136 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t135 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t131 = rSges(3,1) * t165 - rSges(3,2) * t163;
t130 = rSges(3,1) * t163 + rSges(3,2) * t165;
t129 = rSges(5,1) * t162 + rSges(5,2) * t164;
t128 = Icges(3,1) * t165 - t217;
t127 = Icges(3,1) * t163 + t156;
t125 = -Icges(3,2) * t163 + t156;
t124 = Icges(3,2) * t165 + t217;
t118 = rSges(6,1) * t158 + rSges(6,2) * t159;
t108 = rSges(4,3) * t163 + t202 * t165;
t107 = -rSges(4,3) * t165 + t202 * t163;
t100 = V_base(5) * rSges(2,3) - t150 * t161 + t205;
t99 = t151 * t161 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t98 = t150 * V_base(4) - t151 * V_base(5) + V_base(3);
t97 = rSges(5,3) * t163 + t201 * t165;
t96 = -rSges(5,3) * t165 + t201 * t163;
t89 = rSges(6,3) * t163 + t200 * t165;
t88 = -rSges(6,3) * t165 + t200 * t163;
t77 = V_base(5) * rSges(3,3) - t130 * t157 + t190;
t76 = t131 * t157 + (-rSges(3,3) + t227) * V_base(4) + t207;
t75 = pkin(9) * t163 + t210 * t165;
t74 = -pkin(9) * t165 + t210 * t163;
t73 = t130 * V_base(4) + (-t131 - t224) * V_base(5) + t206;
t72 = t138 * t149 + (-t107 - t132) * t157 + t190;
t71 = t108 * t157 - t139 * t149 + t185;
t70 = t107 * t139 - t108 * t138 + t184;
t69 = t112 * t129 + (-t96 + t219) * t157 + t186;
t68 = -t113 * t129 + t157 * t97 + t183;
t67 = -t112 * t97 + t113 * t96 + t182;
t66 = t112 * t222 + t109 * t118 + (-t74 - t88 + t219) * t157 + t186;
t65 = -t113 * t222 - t110 * t118 + (t75 + t89) * t157 + t183;
t64 = -t109 * t89 + t110 * t88 - t112 * t75 + t113 * t74 + t182;
t1 = m(1) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + t139 * (t187 * t163 + t179 * t165) / 0.2e1 + t138 * (t179 * t163 - t187 * t165) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + t113 * (t188 * t163 + t180 * t165) / 0.2e1 + t112 * (t180 * t163 - t188 * t165) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + t110 * (t189 * t163 + t181 * t165) / 0.2e1 + t109 * (t181 * t163 - t189 * t165) / 0.2e1 + ((-t124 * t163 + t127 * t165 - t144 * t175 + t147 * t177 + Icges(1,4)) * V_base(5) + (-t125 * t163 + t128 * t165 - t145 * t175 + t148 * t177 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t124 * t165 + t127 * t163 + t144 * t177 + t147 * t175 + Icges(1,2)) * V_base(5) + (t125 * t165 + t128 * t163 + t145 * t177 + t148 * t175 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t104 * t176 + t106 * t174) * t139 + (t103 * t176 + t105 * t174) * t138 + (t162 * t95 + t164 * t93) * t113 + (t162 * t94 + t164 * t92) * t112 + (t158 * t87 + t159 * t85) * t110 + (t158 * t86 + t159 * t84) * t109 + (t116 * t159 + t117 * t158 + t123 * t164 + t126 * t162 + t143 * t176 + t146 * t174 + Icges(3,3)) * t157) * t157 / 0.2e1 + V_base(4) * t157 * (Icges(3,5) * t165 - Icges(3,6) * t163) + V_base(5) * t157 * (Icges(3,5) * t163 + Icges(3,6) * t165) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t175 + Icges(2,6) * t177) * V_base(5) + (Icges(2,5) * t177 - Icges(2,6) * t175) * V_base(4) + Icges(2,3) * t161 / 0.2e1) * t161;
T = t1;
