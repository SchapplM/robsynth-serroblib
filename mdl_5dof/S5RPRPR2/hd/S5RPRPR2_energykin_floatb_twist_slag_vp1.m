% Calculate kinetic energy for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:33:34
% EndTime: 2020-01-03 11:33:36
% DurationCPUTime: 1.66s
% Computational Cost: add. (1103->234), mult. (736->300), div. (0->0), fcn. (514->10), ass. (0->119)
t152 = qJ(1) + pkin(8);
t145 = qJ(3) + t152;
t137 = sin(t145);
t135 = Icges(4,4) * t137;
t138 = cos(t145);
t153 = sin(pkin(9));
t154 = cos(pkin(9));
t192 = Icges(5,4) * t154;
t169 = -Icges(5,2) * t153 + t192;
t82 = -Icges(5,6) * t137 - t138 * t169;
t193 = Icges(5,4) * t153;
t171 = Icges(5,1) * t154 - t193;
t84 = -Icges(5,5) * t137 - t138 * t171;
t209 = Icges(4,1) * t138 + t153 * t82 - t154 * t84 - t135;
t194 = Icges(4,4) * t138;
t81 = -Icges(5,6) * t138 + t137 * t169;
t83 = -Icges(5,5) * t138 + t137 * t171;
t208 = -Icges(4,1) * t137 + t153 * t81 - t154 * t83 - t194;
t151 = pkin(9) + qJ(5);
t143 = cos(t151);
t141 = sin(t151);
t191 = Icges(6,4) * t141;
t105 = Icges(6,2) * t143 + t191;
t190 = Icges(6,4) * t143;
t108 = Icges(6,1) * t141 + t190;
t114 = -qJD(5) * t137 + V_base(6);
t115 = -qJD(5) * t138 + V_base(5);
t146 = V_base(4) + qJD(1);
t139 = qJD(3) + t146;
t168 = -Icges(6,2) * t141 + t190;
t73 = -Icges(6,6) * t138 + t137 * t168;
t74 = -Icges(6,6) * t137 - t138 * t168;
t170 = Icges(6,1) * t143 - t191;
t75 = -Icges(6,5) * t138 + t137 * t170;
t76 = -Icges(6,5) * t137 - t138 * t170;
t205 = (t105 * t141 - t108 * t143) * t139 + (t141 * t73 - t143 * t75) * t115 + (t141 * t74 - t143 * t76) * t114;
t156 = sin(qJ(1));
t202 = pkin(1) * t156;
t157 = cos(qJ(1));
t201 = pkin(1) * t157;
t200 = pkin(2) * t146;
t199 = pkin(4) * t153;
t198 = pkin(4) * t154;
t197 = -pkin(5) - qJ(2);
t196 = Icges(2,4) * t157;
t144 = cos(t152);
t195 = Icges(3,4) * t144;
t188 = -pkin(6) + t197;
t187 = t146 * t202 + V_base(3);
t186 = V_base(6) * pkin(5) + V_base(2);
t183 = qJD(2) + V_base(1);
t100 = -pkin(3) * t138 - qJ(4) * t137;
t182 = V_base(5) * t100 + t183;
t142 = sin(t152);
t181 = t142 * t200 + t187;
t180 = V_base(6) * qJ(2) + t146 * t201 + t186;
t179 = -pkin(2) * t142 - t202;
t178 = -pkin(2) * t144 - t201;
t177 = rSges(5,1) * t154 - rSges(5,2) * t153;
t176 = rSges(6,1) * t143 - rSges(6,2) * t141;
t167 = Icges(5,5) * t154 - Icges(5,6) * t153;
t166 = Icges(6,5) * t143 - Icges(6,6) * t141;
t122 = Icges(5,2) * t154 + t193;
t123 = Icges(5,1) * t153 + t192;
t164 = t122 * t153 - t123 * t154;
t98 = pkin(3) * t137 - qJ(4) * t138;
t163 = t179 - t98;
t162 = V_base(6) * pkin(6) + t144 * t200 + t180;
t161 = -qJD(4) * t137 + t139 * t98 + t181;
t160 = -(Icges(6,5) * t141 + Icges(6,6) * t143) * t139 - (-Icges(6,3) * t137 - t138 * t166) * t114 - (-Icges(6,3) * t138 + t137 * t166) * t115;
t159 = -qJD(4) * t138 + t162;
t158 = -(Icges(5,5) * t153 + Icges(5,6) * t154) * t139 - (-Icges(5,3) * t138 + t137 * t167) * V_base(5) - (-Icges(5,3) * t137 - t138 * t167) * V_base(6);
t148 = Icges(2,4) * t156;
t136 = Icges(3,4) * t142;
t132 = -rSges(2,1) * t157 + t156 * rSges(2,2);
t131 = t156 * rSges(2,1) + rSges(2,2) * t157;
t130 = -Icges(2,1) * t157 + t148;
t129 = Icges(2,1) * t156 + t196;
t128 = Icges(2,2) * t156 - t196;
t127 = Icges(2,2) * t157 + t148;
t124 = rSges(5,1) * t153 + rSges(5,2) * t154;
t118 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t117 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t116 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t113 = -rSges(3,1) * t144 + rSges(3,2) * t142;
t112 = rSges(3,1) * t142 + rSges(3,2) * t144;
t111 = rSges(6,1) * t141 + rSges(6,2) * t143;
t110 = -Icges(3,1) * t144 + t136;
t109 = Icges(3,1) * t142 + t195;
t107 = Icges(3,2) * t142 - t195;
t106 = Icges(3,2) * t144 + t136;
t101 = -rSges(4,1) * t138 + rSges(4,2) * t137;
t99 = rSges(4,1) * t137 + rSges(4,2) * t138;
t95 = Icges(4,2) * t137 - t194;
t94 = Icges(4,2) * t138 + t135;
t93 = -Icges(4,5) * t138 + Icges(4,6) * t137;
t92 = Icges(4,5) * t137 + Icges(4,6) * t138;
t89 = V_base(6) * rSges(2,3) - t132 * t146 + t186;
t88 = t131 * t146 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t87 = -t131 * V_base(6) + t132 * V_base(5) + V_base(1);
t86 = -rSges(5,3) * t137 - t138 * t177;
t85 = -rSges(5,3) * t138 + t137 * t177;
t78 = -rSges(6,3) * t137 - t138 * t176;
t77 = -rSges(6,3) * t138 + t137 * t176;
t70 = -pkin(7) * t137 - t138 * t198;
t69 = -pkin(7) * t138 + t137 * t198;
t68 = V_base(6) * rSges(3,3) - t113 * t146 + t180;
t67 = t112 * t146 + (-rSges(3,3) + t197) * V_base(5) + t187;
t66 = -V_base(6) * t112 + V_base(5) * t113 + (-t156 * V_base(6) - t157 * V_base(5)) * pkin(1) + t183;
t65 = V_base(6) * rSges(4,3) - t101 * t139 + t162;
t64 = t139 * t99 + (-rSges(4,3) + t188) * V_base(5) + t181;
t63 = (t179 - t99) * V_base(6) + (t101 + t178) * V_base(5) + t183;
t62 = t124 * V_base(6) + (-t100 - t86) * t139 + t159;
t61 = t139 * t85 + (-t124 + t188) * V_base(5) + t161;
t60 = (t178 + t86) * V_base(5) + (t163 - t85) * V_base(6) + t182;
t59 = V_base(6) * t199 + t111 * t114 + (-t100 - t70 - t78) * t139 + t159;
t58 = -t111 * t115 + (t69 + t77) * t139 + (t188 - t199) * V_base(5) + t161;
t57 = -t114 * t77 + t115 * t78 + (t178 + t70) * V_base(5) + (t163 - t69) * V_base(6) + t182;
t1 = m(1) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(2) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(3) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t115 * (-t205 * t137 + t160 * t138) / 0.2e1 + t114 * (t160 * t137 + t205 * t138) / 0.2e1 + ((t141 * t75 + t143 * t73) * t115 + (t141 * t76 + t143 * t74) * t114 + (t153 * t84 + t154 * t82 + t93) * V_base(6) + (t153 * t83 + t154 * t81 + t92) * V_base(5) + (t105 * t143 + t108 * t141 + t122 * t154 + t123 * t153 + Icges(4,3)) * t139) * t139 / 0.2e1 + (t158 * t138 + (-t164 * t137 + t92) * t139 + (t107 * t144 + t110 * t142 + t128 * t157 + t156 * t130 - t137 * t209 + t138 * t95 + Icges(1,6)) * V_base(6) + (t106 * t144 + t109 * t142 + t127 * t157 + t156 * t129 - t208 * t137 + t138 * t94 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t158 * t137 + (t164 * t138 + t93) * t139 + (t107 * t142 - t110 * t144 + t156 * t128 - t130 * t157 + t137 * t95 + t138 * t209 + Icges(1,3)) * V_base(6) + (t106 * t142 - t109 * t144 + t156 * t127 - t129 * t157 + t137 * t94 + t138 * t208 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((-Icges(2,5) * t157 - Icges(3,5) * t144 + Icges(2,6) * t156 + Icges(3,6) * t142) * V_base(6) + (Icges(2,5) * t156 + Icges(3,5) * t142 + Icges(2,6) * t157 + Icges(3,6) * t144) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * t146) * t146 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
