% Calculate kinetic energy for
% S5RRPRR3
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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:16
% EndTime: 2019-12-05 18:30:17
% DurationCPUTime: 1.52s
% Computational Cost: add. (1001->210), mult. (596->272), div. (0->0), fcn. (372->10), ass. (0->106)
t142 = qJ(1) + qJ(2);
t134 = pkin(9) + t142;
t133 = qJ(4) + t134;
t126 = sin(t133);
t104 = qJD(5) * t126 + V_base(6);
t127 = cos(t133);
t105 = qJD(5) * t127 + V_base(5);
t145 = cos(qJ(5));
t143 = sin(qJ(5));
t169 = Icges(6,4) * t143;
t114 = Icges(6,2) * t145 + t169;
t168 = Icges(6,4) * t145;
t117 = Icges(6,1) * t143 + t168;
t135 = V_base(4) + qJD(1);
t132 = qJD(2) + t135;
t125 = qJD(4) + t132;
t156 = -Icges(6,2) * t143 + t168;
t69 = Icges(6,6) * t127 - t126 * t156;
t70 = Icges(6,6) * t126 + t127 * t156;
t157 = Icges(6,1) * t145 - t169;
t71 = Icges(6,5) * t127 - t126 * t157;
t72 = Icges(6,5) * t126 + t127 * t157;
t187 = (t114 * t143 - t117 * t145) * t125 + (t143 * t69 - t145 * t71) * t105 + (t143 * t70 - t145 * t72) * t104;
t186 = -pkin(5) - pkin(6);
t144 = sin(qJ(1));
t183 = pkin(1) * t144;
t146 = cos(qJ(1));
t182 = pkin(1) * t146;
t136 = sin(t142);
t181 = pkin(2) * t136;
t137 = cos(t142);
t180 = pkin(2) * t137;
t130 = sin(t134);
t179 = pkin(3) * t130;
t131 = cos(t134);
t178 = pkin(3) * t131;
t177 = Icges(2,4) * t144;
t176 = Icges(2,4) * t146;
t175 = Icges(3,4) * t136;
t174 = Icges(3,4) * t137;
t173 = Icges(4,4) * t130;
t172 = Icges(4,4) * t131;
t171 = Icges(5,4) * t126;
t170 = Icges(5,4) * t127;
t167 = -qJ(3) + t186;
t166 = V_base(6) * pkin(5) + V_base(2);
t163 = -pkin(7) + t167;
t162 = V_base(5) * t182 + V_base(6) * t183 + V_base(1);
t161 = rSges(6,1) * t145 - rSges(6,2) * t143;
t160 = -t135 * t183 + V_base(3);
t155 = Icges(6,5) * t145 - Icges(6,6) * t143;
t153 = V_base(5) * t180 + V_base(6) * t181 + qJD(3) + t162;
t152 = V_base(6) * pkin(6) - t135 * t182 + t166;
t151 = t104 * (Icges(6,3) * t126 + t127 * t155) + t105 * (Icges(6,3) * t127 - t126 * t155) + (Icges(6,5) * t143 + Icges(6,6) * t145) * t125;
t150 = V_base(6) * qJ(3) + t152;
t149 = V_base(5) * t178 + V_base(6) * t179 + t153;
t148 = (-t179 - t181) * t132 + t160;
t147 = V_base(6) * pkin(7) + (-t178 - t180) * t132 + t150;
t122 = rSges(2,1) * t146 - t144 * rSges(2,2);
t121 = -t144 * rSges(2,1) - rSges(2,2) * t146;
t120 = rSges(6,1) * t143 + rSges(6,2) * t145;
t119 = Icges(2,1) * t146 - t177;
t118 = -Icges(2,1) * t144 - t176;
t116 = -Icges(2,2) * t144 + t176;
t115 = -Icges(2,2) * t146 - t177;
t108 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t107 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t106 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t103 = rSges(3,1) * t137 - rSges(3,2) * t136;
t102 = -rSges(3,1) * t136 - rSges(3,2) * t137;
t101 = Icges(3,1) * t137 - t175;
t100 = -Icges(3,1) * t136 - t174;
t99 = -Icges(3,2) * t136 + t174;
t98 = -Icges(3,2) * t137 - t175;
t95 = rSges(4,1) * t131 - rSges(4,2) * t130;
t94 = -rSges(4,1) * t130 - rSges(4,2) * t131;
t93 = Icges(4,1) * t131 - t173;
t92 = -Icges(4,1) * t130 - t172;
t91 = -Icges(4,2) * t130 + t172;
t90 = -Icges(4,2) * t131 - t173;
t87 = pkin(4) * t127 + pkin(8) * t126;
t86 = -pkin(4) * t126 + pkin(8) * t127;
t85 = rSges(5,1) * t127 - rSges(5,2) * t126;
t84 = -rSges(5,1) * t126 - rSges(5,2) * t127;
t83 = Icges(5,1) * t127 - t171;
t82 = -Icges(5,1) * t126 - t170;
t81 = -Icges(5,2) * t126 + t170;
t80 = -Icges(5,2) * t127 - t171;
t77 = V_base(6) * rSges(2,3) - t122 * t135 + t166;
t76 = t121 * t135 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t75 = -t121 * V_base(6) + t122 * V_base(5) + V_base(1);
t74 = rSges(6,3) * t126 + t127 * t161;
t73 = rSges(6,3) * t127 - t126 * t161;
t66 = V_base(6) * rSges(3,3) - t132 * t103 + t152;
t65 = t102 * t132 + (-rSges(3,3) + t186) * V_base(5) + t160;
t64 = -t102 * V_base(6) + t103 * V_base(5) + t162;
t63 = V_base(6) * rSges(4,3) + (-t95 - t180) * t132 + t150;
t62 = (t94 - t181) * t132 + (-rSges(4,3) + t167) * V_base(5) + t160;
t61 = -t94 * V_base(6) + t95 * V_base(5) + t153;
t60 = V_base(6) * rSges(5,3) - t125 * t85 + t147;
t59 = t125 * t84 + (-rSges(5,3) + t163) * V_base(5) + t148;
t58 = -t84 * V_base(6) + t85 * V_base(5) + t149;
t57 = t104 * t120 + (-t74 - t87) * t125 + t147;
t56 = -t105 * t120 + (t73 + t86) * t125 + t163 * V_base(5) + t148;
t55 = -t104 * t73 + t105 * t74 - t86 * V_base(6) + t87 * V_base(5) + t149;
t1 = m(1) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(2) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(3) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(4) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(5) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + t105 * (t187 * t126 + t151 * t127) / 0.2e1 + t104 * (t151 * t126 - t187 * t127) / 0.2e1 + ((t143 * t71 + t145 * t69) * t105 + (t143 * t72 + t145 * t70) * t104 + (t145 * t114 + t143 * t117 + Icges(5,3)) * t125) * t125 / 0.2e1 + ((-t101 * t136 - t116 * t146 - t144 * t119 - t126 * t83 - t127 * t81 - t130 * t93 - t131 * t91 - t137 * t99 + Icges(1,6)) * V_base(6) + (-t136 * t100 - t146 * t115 - t144 * t118 - t126 * t82 - t127 * t80 - t130 * t92 - t131 * t90 - t137 * t98 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((t137 * t101 - t144 * t116 + t146 * t119 - t126 * t81 + t127 * t83 - t130 * t91 + t131 * t93 - t136 * t99 + Icges(1,3)) * V_base(6) + (t100 * t137 - t144 * t115 + t118 * t146 - t126 * t80 + t127 * t82 - t130 * t90 + t131 * t92 - t136 * t98 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + t125 * V_base(6) * (Icges(5,5) * t127 - Icges(5,6) * t126) + t125 * V_base(5) * (-Icges(5,5) * t126 - Icges(5,6) * t127) + ((Icges(3,5) * t137 + Icges(4,5) * t131 - Icges(3,6) * t136 - Icges(4,6) * t130) * V_base(6) + (-Icges(3,5) * t136 - Icges(4,5) * t130 - Icges(3,6) * t137 - Icges(4,6) * t131) * V_base(5) + (Icges(3,3) / 0.2e1 + Icges(4,3) / 0.2e1) * t132) * t132 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((-Icges(2,5) * t144 - Icges(2,6) * t146) * V_base(5) + (Icges(2,5) * t146 - Icges(2,6) * t144) * V_base(6) + Icges(2,3) * t135 / 0.2e1) * t135;
T = t1;
