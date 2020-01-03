% Calculate kinetic energy for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:14
% EndTime: 2019-12-31 17:24:15
% DurationCPUTime: 1.02s
% Computational Cost: add. (794->207), mult. (862->307), div. (0->0), fcn. (694->8), ass. (0->110)
t140 = sin(qJ(2));
t185 = pkin(2) * t140;
t139 = qJ(2) + qJ(3);
t132 = sin(t139);
t184 = pkin(3) * t132;
t142 = cos(qJ(2));
t183 = t142 * pkin(2);
t141 = sin(qJ(1));
t143 = cos(qJ(1));
t122 = t141 * pkin(1) - t143 * pkin(5);
t73 = -pkin(6) * t143 + t183 * t141;
t181 = -t122 - t73;
t180 = Icges(2,4) * t141;
t179 = Icges(3,4) * t140;
t178 = Icges(3,4) * t142;
t177 = Icges(4,4) * t132;
t133 = cos(t139);
t176 = Icges(4,4) * t133;
t135 = qJ(4) + t139;
t126 = sin(t135);
t175 = Icges(5,4) * t126;
t127 = cos(t135);
t174 = Icges(5,4) * t127;
t173 = pkin(3) * t133;
t171 = -qJD(2) - qJD(3);
t170 = V_base(5) * pkin(4) + V_base(1);
t125 = qJD(2) * t141 + V_base(4);
t124 = -qJD(2) * t143 + V_base(5);
t167 = t124 * t185 + t170;
t103 = qJD(3) * t141 + t125;
t166 = rSges(3,1) * t142 - rSges(3,2) * t140;
t165 = rSges(4,1) * t133 - rSges(4,2) * t132;
t164 = rSges(5,1) * t127 - rSges(5,2) * t126;
t163 = Icges(3,1) * t142 - t179;
t162 = Icges(4,1) * t133 - t177;
t161 = Icges(5,1) * t127 - t175;
t160 = -Icges(3,2) * t140 + t178;
t159 = -Icges(4,2) * t132 + t176;
t158 = -Icges(5,2) * t126 + t174;
t157 = Icges(3,5) * t142 - Icges(3,6) * t140;
t156 = Icges(4,5) * t133 - Icges(4,6) * t132;
t155 = Icges(5,5) * t127 - Icges(5,6) * t126;
t123 = t143 * pkin(1) + t141 * pkin(5);
t129 = V_base(6) + qJD(1);
t154 = -V_base(4) * pkin(4) + t129 * t123 + V_base(2);
t153 = V_base(4) * t122 - t123 * V_base(5) + V_base(3);
t91 = V_base(5) + (-qJD(4) + t171) * t143;
t92 = qJD(4) * t141 + t103;
t152 = (Icges(5,5) * t126 + Icges(5,6) * t127) * t129 + (-Icges(5,3) * t143 + t155 * t141) * t91 + (Icges(5,3) * t141 + t155 * t143) * t92;
t102 = t171 * t143 + V_base(5);
t151 = (-Icges(4,3) * t143 + t156 * t141) * t102 + (Icges(4,3) * t141 + t156 * t143) * t103 + (Icges(4,5) * t132 + Icges(4,6) * t133) * t129;
t150 = (Icges(3,5) * t140 + Icges(3,6) * t142) * t129 + (-Icges(3,3) * t143 + t157 * t141) * t124 + (Icges(3,3) * t141 + t157 * t143) * t125;
t74 = pkin(6) * t141 + t183 * t143;
t149 = -t124 * t74 + t125 * t73 + t153;
t148 = -t125 * t185 + t129 * t74 + t154;
t67 = -Icges(5,6) * t143 + t158 * t141;
t68 = Icges(5,6) * t141 + t158 * t143;
t69 = -Icges(5,5) * t143 + t161 * t141;
t70 = Icges(5,5) * t141 + t161 * t143;
t94 = Icges(5,2) * t127 + t175;
t95 = Icges(5,1) * t126 + t174;
t147 = (-t126 * t68 + t127 * t70) * t92 + (-t126 * t67 + t127 * t69) * t91 + (-t126 * t94 + t127 * t95) * t129;
t100 = Icges(4,1) * t132 + t176;
t77 = -Icges(4,6) * t143 + t159 * t141;
t78 = Icges(4,6) * t141 + t159 * t143;
t79 = -Icges(4,5) * t143 + t162 * t141;
t80 = Icges(4,5) * t141 + t162 * t143;
t99 = Icges(4,2) * t133 + t177;
t146 = (-t132 * t78 + t133 * t80) * t103 + (-t132 * t77 + t133 * t79) * t102 + (t100 * t133 - t132 * t99) * t129;
t113 = Icges(3,2) * t142 + t179;
t116 = Icges(3,1) * t140 + t178;
t85 = -Icges(3,6) * t143 + t160 * t141;
t86 = Icges(3,6) * t141 + t160 * t143;
t87 = -Icges(3,5) * t143 + t163 * t141;
t88 = Icges(3,5) * t141 + t163 * t143;
t145 = (-t140 * t86 + t142 * t88) * t125 + (-t140 * t85 + t142 * t87) * t124 + (-t113 * t140 + t116 * t142) * t129;
t134 = Icges(2,4) * t143;
t121 = rSges(2,1) * t143 - rSges(2,2) * t141;
t120 = rSges(2,1) * t141 + rSges(2,2) * t143;
t119 = rSges(3,1) * t140 + rSges(3,2) * t142;
t118 = Icges(2,1) * t143 - t180;
t117 = Icges(2,1) * t141 + t134;
t115 = -Icges(2,2) * t141 + t134;
t114 = Icges(2,2) * t143 + t180;
t109 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t108 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t107 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t101 = rSges(4,1) * t132 + rSges(4,2) * t133;
t96 = rSges(5,1) * t126 + rSges(5,2) * t127;
t90 = rSges(3,3) * t141 + t166 * t143;
t89 = -rSges(3,3) * t143 + t166 * t141;
t82 = rSges(4,3) * t141 + t165 * t143;
t81 = -rSges(4,3) * t143 + t165 * t141;
t72 = rSges(5,3) * t141 + t164 * t143;
t71 = -rSges(5,3) * t143 + t164 * t141;
t64 = V_base(5) * rSges(2,3) - t120 * t129 + t170;
t63 = t121 * t129 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t62 = t120 * V_base(4) - t121 * V_base(5) + V_base(3);
t59 = pkin(7) * t141 + t173 * t143;
t58 = -pkin(7) * t143 + t173 * t141;
t57 = t119 * t124 + (-t122 - t89) * t129 + t170;
t56 = -t119 * t125 + t129 * t90 + t154;
t55 = -t124 * t90 + t125 * t89 + t153;
t54 = t101 * t102 + (-t81 + t181) * t129 + t167;
t53 = -t101 * t103 + t129 * t82 + t148;
t52 = -t102 * t82 + t103 * t81 + t149;
t51 = t102 * t184 + t91 * t96 + (-t58 - t71 + t181) * t129 + t167;
t50 = -t103 * t184 - t92 * t96 + (t59 + t72) * t129 + t148;
t49 = -t102 * t59 + t103 * t58 + t71 * t92 - t72 * t91 + t149;
t1 = m(1) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(2) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(3) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + t125 * (t150 * t141 + t145 * t143) / 0.2e1 + t124 * (t145 * t141 - t150 * t143) / 0.2e1 + m(4) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + t103 * (t151 * t141 + t146 * t143) / 0.2e1 + t102 * (t146 * t141 - t151 * t143) / 0.2e1 + m(5) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + t92 * (t152 * t141 + t147 * t143) / 0.2e1 + t91 * (t147 * t141 - t152 * t143) / 0.2e1 + ((-t114 * t141 + t117 * t143 + Icges(1,4)) * V_base(5) + (-t115 * t141 + t118 * t143 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t114 * t143 + t117 * t141 + Icges(1,2)) * V_base(5) + (t143 * t115 + t118 * t141 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t140 * t88 + t142 * t86) * t125 + (t140 * t87 + t142 * t85) * t124 + (t132 * t80 + t133 * t78) * t103 + (t132 * t79 + t133 * t77) * t102 + (t126 * t70 + t127 * t68) * t92 + (t126 * t69 + t127 * t67) * t91 + (t100 * t132 + t113 * t142 + t116 * t140 + t126 * t95 + t127 * t94 + t133 * t99 + Icges(2,3)) * t129) * t129 / 0.2e1 + V_base(4) * t129 * (Icges(2,5) * t143 - Icges(2,6) * t141) + V_base(5) * t129 * (Icges(2,5) * t141 + Icges(2,6) * t143) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
