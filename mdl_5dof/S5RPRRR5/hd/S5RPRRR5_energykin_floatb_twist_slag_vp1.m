% Calculate kinetic energy for
% S5RPRRR5
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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:31
% EndTime: 2022-01-20 09:48:32
% DurationCPUTime: 1.59s
% Computational Cost: add. (1175->229), mult. (760->308), div. (0->0), fcn. (538->10), ass. (0->118)
t160 = sin(qJ(1));
t203 = pkin(1) * t160;
t162 = cos(qJ(1));
t202 = pkin(1) * t162;
t157 = qJ(1) + pkin(9);
t147 = sin(t157);
t201 = pkin(2) * t147;
t148 = cos(t157);
t200 = pkin(2) * t148;
t159 = sin(qJ(4));
t199 = pkin(4) * t159;
t161 = cos(qJ(4));
t198 = pkin(4) * t161;
t197 = -pkin(5) - qJ(2);
t195 = Icges(2,4) * t160;
t194 = Icges(3,4) * t147;
t149 = qJ(3) + t157;
t143 = sin(t149);
t193 = Icges(4,4) * t143;
t192 = Icges(5,4) * t159;
t191 = Icges(5,4) * t161;
t158 = qJ(4) + qJ(5);
t151 = sin(t158);
t190 = Icges(6,4) * t151;
t152 = cos(t158);
t189 = Icges(6,4) * t152;
t188 = -pkin(6) + t197;
t150 = V_base(6) + qJD(1);
t187 = t150 * t202 + V_base(2);
t186 = V_base(5) * pkin(5) + V_base(1);
t120 = qJD(4) * t143 + V_base(4);
t183 = t150 * t200 + t187;
t182 = V_base(5) * qJ(2) + t186;
t181 = V_base(4) * t203 + qJD(2) + V_base(3);
t180 = -t200 - t202;
t179 = V_base(4) * t201 + t181;
t178 = rSges(5,1) * t161 - rSges(5,2) * t159;
t177 = rSges(6,1) * t152 - rSges(6,2) * t151;
t176 = Icges(5,1) * t161 - t192;
t175 = Icges(6,1) * t152 - t190;
t174 = -Icges(5,2) * t159 + t191;
t173 = -Icges(6,2) * t151 + t189;
t172 = Icges(5,5) * t161 - Icges(5,6) * t159;
t171 = Icges(6,5) * t152 - Icges(6,6) * t151;
t144 = cos(t149);
t145 = qJD(3) + t150;
t94 = V_base(5) + (-qJD(4) - qJD(5)) * t144;
t95 = qJD(5) * t143 + t120;
t170 = (Icges(6,5) * t151 + Icges(6,6) * t152) * t145 + (-Icges(6,3) * t144 + t143 * t171) * t94 + (Icges(6,3) * t143 + t144 * t171) * t95;
t119 = -qJD(4) * t144 + V_base(5);
t169 = (-Icges(5,3) * t144 + t143 * t172) * t119 + (Icges(5,3) * t143 + t144 * t172) * t120 + (Icges(5,5) * t159 + Icges(5,6) * t161) * t145;
t106 = t144 * pkin(3) + t143 * pkin(7);
t168 = t145 * t106 + t188 * V_base(4) + t183;
t167 = V_base(5) * pkin(6) + (-t201 - t203) * t150 + t182;
t105 = t143 * pkin(3) - t144 * pkin(7);
t166 = V_base(4) * t105 + (-t106 + t180) * V_base(5) + t179;
t116 = Icges(6,2) * t152 + t190;
t117 = Icges(6,1) * t151 + t189;
t76 = -Icges(6,6) * t144 + t143 * t173;
t77 = Icges(6,6) * t143 + t144 * t173;
t78 = -Icges(6,5) * t144 + t143 * t175;
t79 = Icges(6,5) * t143 + t144 * t175;
t165 = (-t151 * t77 + t152 * t79) * t95 + (-t151 * t76 + t152 * t78) * t94 + (-t116 * t151 + t117 * t152) * t145;
t128 = Icges(5,2) * t161 + t192;
t131 = Icges(5,1) * t159 + t191;
t84 = -Icges(5,6) * t144 + t143 * t174;
t85 = Icges(5,6) * t143 + t144 * t174;
t86 = -Icges(5,5) * t144 + t143 * t176;
t87 = Icges(5,5) * t143 + t144 * t176;
t164 = (-t159 * t85 + t161 * t87) * t120 + (-t159 * t84 + t161 * t86) * t119 + (-t128 * t159 + t131 * t161) * t145;
t154 = Icges(2,4) * t162;
t142 = Icges(3,4) * t148;
t140 = Icges(4,4) * t144;
t136 = rSges(2,1) * t162 - rSges(2,2) * t160;
t135 = rSges(2,1) * t160 + rSges(2,2) * t162;
t134 = rSges(5,1) * t159 + rSges(5,2) * t161;
t133 = Icges(2,1) * t162 - t195;
t132 = Icges(2,1) * t160 + t154;
t130 = -Icges(2,2) * t160 + t154;
t129 = Icges(2,2) * t162 + t195;
t123 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t122 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t121 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t118 = rSges(6,1) * t151 + rSges(6,2) * t152;
t114 = rSges(3,1) * t148 - rSges(3,2) * t147;
t113 = rSges(3,1) * t147 + rSges(3,2) * t148;
t112 = Icges(3,1) * t148 - t194;
t111 = Icges(3,1) * t147 + t142;
t110 = -Icges(3,2) * t147 + t142;
t109 = Icges(3,2) * t148 + t194;
t104 = rSges(4,1) * t144 - rSges(4,2) * t143;
t103 = rSges(4,1) * t143 + rSges(4,2) * t144;
t102 = Icges(4,1) * t144 - t193;
t101 = Icges(4,1) * t143 + t140;
t100 = -Icges(4,2) * t143 + t140;
t99 = Icges(4,2) * t144 + t193;
t92 = V_base(5) * rSges(2,3) - t135 * t150 + t186;
t91 = t136 * t150 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t90 = t135 * V_base(4) - t136 * V_base(5) + V_base(3);
t89 = rSges(5,3) * t143 + t144 * t178;
t88 = -rSges(5,3) * t144 + t143 * t178;
t81 = rSges(6,3) * t143 + t144 * t177;
t80 = -rSges(6,3) * t144 + t143 * t177;
t73 = pkin(8) * t143 + t144 * t198;
t72 = -pkin(8) * t144 + t143 * t198;
t71 = V_base(5) * rSges(3,3) + (-t113 - t203) * t150 + t182;
t70 = t114 * t150 + (-rSges(3,3) + t197) * V_base(4) + t187;
t69 = t113 * V_base(4) + (-t114 - t202) * V_base(5) + t181;
t68 = V_base(5) * rSges(4,3) - t103 * t145 + t167;
t67 = t104 * t145 + (-rSges(4,3) + t188) * V_base(4) + t183;
t66 = t103 * V_base(4) + (-t104 + t180) * V_base(5) + t179;
t65 = t119 * t134 + (-t105 - t88) * t145 + t167;
t64 = -t120 * t134 + t145 * t89 + t168;
t63 = -t119 * t89 + t120 * t88 + t166;
t62 = t119 * t199 + t118 * t94 + (-t105 - t72 - t80) * t145 + t167;
t61 = -t120 * t199 - t118 * t95 + (t73 + t81) * t145 + t168;
t60 = -t119 * t73 + t120 * t72 + t80 * t95 - t81 * t94 + t166;
t1 = m(1) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(2) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t120 * (t169 * t143 + t164 * t144) / 0.2e1 + t119 * (t164 * t143 - t169 * t144) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t95 * (t143 * t170 + t144 * t165) / 0.2e1 + t94 * (t143 * t165 - t170 * t144) / 0.2e1 + ((t159 * t87 + t161 * t85) * t120 + (t159 * t86 + t161 * t84) * t119 + (t151 * t79 + t152 * t77) * t95 + (t151 * t78 + t152 * t76) * t94 + (t152 * t116 + t151 * t117 + t161 * t128 + t159 * t131 + Icges(4,3)) * t145) * t145 / 0.2e1 + ((t101 * t144 - t109 * t147 + t111 * t148 - t129 * t160 + t132 * t162 - t143 * t99 + Icges(1,4)) * V_base(5) + (-t143 * t100 + t144 * t102 - t147 * t110 + t148 * t112 - t160 * t130 + t162 * t133 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t143 * t101 + t148 * t109 + t147 * t111 + t162 * t129 + t160 * t132 + t144 * t99 + Icges(1,2)) * V_base(5) + (t100 * t144 + t102 * t143 + t110 * t148 + t112 * t147 + t130 * t162 + t133 * t160 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t145 * (Icges(4,5) * t144 - Icges(4,6) * t143) + V_base(5) * t145 * (Icges(4,5) * t143 + Icges(4,6) * t144) + ((Icges(2,5) * t162 + Icges(3,5) * t148 - Icges(2,6) * t160 - Icges(3,6) * t147) * V_base(4) + (Icges(2,5) * t160 + Icges(3,5) * t147 + Icges(2,6) * t162 + Icges(3,6) * t148) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * t150) * t150 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
