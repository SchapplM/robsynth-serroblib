% Calculate kinetic energy for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:00
% EndTime: 2022-01-20 10:48:01
% DurationCPUTime: 1.26s
% Computational Cost: add. (1186->229), mult. (760->307), div. (0->0), fcn. (538->10), ass. (0->119)
t205 = -pkin(5) - pkin(6);
t160 = sin(qJ(1));
t203 = pkin(1) * t160;
t162 = cos(qJ(1));
t202 = pkin(1) * t162;
t158 = qJ(1) + qJ(2);
t150 = sin(t158);
t201 = pkin(2) * t150;
t152 = cos(t158);
t200 = pkin(2) * t152;
t159 = sin(qJ(4));
t199 = pkin(4) * t159;
t161 = cos(qJ(4));
t198 = pkin(4) * t161;
t196 = Icges(2,4) * t160;
t195 = Icges(3,4) * t150;
t147 = pkin(9) + t158;
t143 = sin(t147);
t194 = Icges(4,4) * t143;
t193 = Icges(5,4) * t159;
t192 = Icges(5,4) * t161;
t157 = qJ(4) + qJ(5);
t149 = sin(t157);
t191 = Icges(6,4) * t149;
t151 = cos(t157);
t190 = Icges(6,4) * t151;
t189 = -qJ(3) + t205;
t148 = V_base(6) + qJD(1);
t188 = t148 * t202 + V_base(2);
t187 = V_base(4) * t203 + V_base(3);
t186 = V_base(5) * pkin(5) + V_base(1);
t121 = qJD(4) * t143 + V_base(4);
t144 = cos(t147);
t105 = t143 * pkin(3) - t144 * pkin(7);
t183 = -t105 - t201;
t145 = qJD(2) + t148;
t182 = t145 * t200 + t188;
t181 = -t200 - t202;
t180 = V_base(4) * t201 + qJD(3) + t187;
t179 = rSges(5,1) * t161 - rSges(5,2) * t159;
t178 = rSges(6,1) * t151 - rSges(6,2) * t149;
t177 = Icges(5,1) * t161 - t193;
t176 = Icges(6,1) * t151 - t191;
t175 = -Icges(5,2) * t159 + t192;
t174 = -Icges(6,2) * t149 + t190;
t173 = Icges(5,5) * t161 - Icges(5,6) * t159;
t172 = Icges(6,5) * t151 - Icges(6,6) * t149;
t171 = V_base(5) * pkin(6) - t148 * t203 + t186;
t94 = V_base(5) + (-qJD(4) - qJD(5)) * t144;
t95 = qJD(5) * t143 + t121;
t170 = (Icges(6,5) * t149 + Icges(6,6) * t151) * t145 + (-Icges(6,3) * t144 + t143 * t172) * t94 + (Icges(6,3) * t143 + t144 * t172) * t95;
t120 = -qJD(4) * t144 + V_base(5);
t169 = t120 * (-Icges(5,3) * t144 + t143 * t173) + t121 * (Icges(5,3) * t143 + t144 * t173) + (Icges(5,5) * t159 + Icges(5,6) * t161) * t145;
t168 = V_base(5) * qJ(3) + t171;
t106 = t144 * pkin(3) + t143 * pkin(7);
t167 = t145 * t106 + t189 * V_base(4) + t182;
t166 = V_base(4) * t105 + (-t106 + t181) * V_base(5) + t180;
t110 = Icges(6,2) * t151 + t191;
t113 = Icges(6,1) * t149 + t190;
t76 = -Icges(6,6) * t144 + t143 * t174;
t77 = Icges(6,6) * t143 + t144 * t174;
t78 = -Icges(6,5) * t144 + t143 * t176;
t79 = Icges(6,5) * t143 + t144 * t176;
t165 = (-t149 * t77 + t151 * t79) * t95 + (-t149 * t76 + t151 * t78) * t94 + (-t110 * t149 + t113 * t151) * t145;
t128 = Icges(5,2) * t161 + t193;
t131 = Icges(5,1) * t159 + t192;
t84 = -Icges(5,6) * t144 + t143 * t175;
t85 = Icges(5,6) * t143 + t144 * t175;
t86 = -Icges(5,5) * t144 + t143 * t177;
t87 = Icges(5,5) * t143 + t144 * t177;
t164 = (-t159 * t85 + t161 * t87) * t121 + (-t159 * t84 + t161 * t86) * t120 + (-t128 * t159 + t131 * t161) * t145;
t154 = Icges(2,4) * t162;
t142 = Icges(3,4) * t152;
t140 = Icges(4,4) * t144;
t136 = rSges(2,1) * t162 - rSges(2,2) * t160;
t135 = rSges(2,1) * t160 + rSges(2,2) * t162;
t134 = rSges(5,1) * t159 + rSges(5,2) * t161;
t133 = Icges(2,1) * t162 - t196;
t132 = Icges(2,1) * t160 + t154;
t130 = -Icges(2,2) * t160 + t154;
t129 = Icges(2,2) * t162 + t196;
t124 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t123 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t122 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t118 = rSges(3,1) * t152 - rSges(3,2) * t150;
t117 = rSges(3,1) * t150 + rSges(3,2) * t152;
t116 = rSges(6,1) * t149 + rSges(6,2) * t151;
t115 = Icges(3,1) * t152 - t195;
t114 = Icges(3,1) * t150 + t142;
t112 = -Icges(3,2) * t150 + t142;
t111 = Icges(3,2) * t152 + t195;
t104 = rSges(4,1) * t144 - rSges(4,2) * t143;
t103 = rSges(4,1) * t143 + rSges(4,2) * t144;
t102 = Icges(4,1) * t144 - t194;
t101 = Icges(4,1) * t143 + t140;
t100 = -Icges(4,2) * t143 + t140;
t99 = Icges(4,2) * t144 + t194;
t92 = V_base(5) * rSges(2,3) - t135 * t148 + t186;
t91 = t136 * t148 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t90 = t135 * V_base(4) - t136 * V_base(5) + V_base(3);
t89 = rSges(5,3) * t143 + t144 * t179;
t88 = -rSges(5,3) * t144 + t143 * t179;
t81 = rSges(6,3) * t143 + t144 * t178;
t80 = -rSges(6,3) * t144 + t143 * t178;
t73 = pkin(8) * t143 + t144 * t198;
t72 = -pkin(8) * t144 + t143 * t198;
t71 = V_base(5) * rSges(3,3) - t117 * t145 + t171;
t70 = t118 * t145 + (-rSges(3,3) + t205) * V_base(4) + t188;
t69 = t117 * V_base(4) + (-t118 - t202) * V_base(5) + t187;
t68 = V_base(5) * rSges(4,3) + (-t103 - t201) * t145 + t168;
t67 = t104 * t145 + (-rSges(4,3) + t189) * V_base(4) + t182;
t66 = t103 * V_base(4) + (-t104 + t181) * V_base(5) + t180;
t65 = t120 * t134 + (t183 - t88) * t145 + t168;
t64 = -t121 * t134 + t145 * t89 + t167;
t63 = -t120 * t89 + t121 * t88 + t166;
t62 = t120 * t199 + t116 * t94 + (t183 - t72 - t80) * t145 + t168;
t61 = -t121 * t199 - t116 * t95 + (t73 + t81) * t145 + t167;
t60 = -t120 * t73 + t121 * t72 + t80 * t95 - t81 * t94 + t166;
t1 = m(1) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(2) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t121 * (t169 * t143 + t164 * t144) / 0.2e1 + t120 * (t164 * t143 - t169 * t144) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t95 * (t170 * t143 + t165 * t144) / 0.2e1 + t94 * (t165 * t143 - t170 * t144) / 0.2e1 + ((t159 * t87 + t161 * t85) * t121 + (t159 * t86 + t161 * t84) * t120 + (t149 * t79 + t151 * t77) * t95 + (t149 * t78 + t151 * t76) * t94 + (t110 * t151 + t113 * t149 + t128 * t161 + t131 * t159 + Icges(3,3) + Icges(4,3)) * t145) * t145 / 0.2e1 + ((t101 * t144 - t111 * t150 + t114 * t152 - t129 * t160 + t132 * t162 - t143 * t99 + Icges(1,4)) * V_base(5) + (-t100 * t143 + t102 * t144 - t112 * t150 + t115 * t152 - t130 * t160 + t133 * t162 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t101 * t143 + t111 * t152 + t114 * t150 + t129 * t162 + t132 * t160 + t144 * t99 + Icges(1,2)) * V_base(5) + (t100 * t144 + t102 * t143 + t112 * t152 + t115 * t150 + t130 * t162 + t133 * t160 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t145 * (Icges(3,5) * t150 + Icges(4,5) * t143 + Icges(3,6) * t152 + Icges(4,6) * t144) + V_base(4) * t145 * (Icges(3,5) * t152 + Icges(4,5) * t144 - Icges(3,6) * t150 - Icges(4,6) * t143) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t160 + Icges(2,6) * t162) * V_base(5) + (Icges(2,5) * t162 - Icges(2,6) * t160) * V_base(4) + Icges(2,3) * t148 / 0.2e1) * t148;
T = t1;
