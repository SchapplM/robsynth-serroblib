% Calculate kinetic energy for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:02:49
% EndTime: 2019-12-05 17:02:50
% DurationCPUTime: 1.67s
% Computational Cost: add. (768->229), mult. (1006->341), div. (0->0), fcn. (876->8), ass. (0->113)
t178 = Icges(1,5) + Icges(2,5);
t137 = cos(qJ(2));
t121 = -qJD(3) * t137 + V_base(6);
t175 = pkin(2) * t121;
t134 = sin(qJ(2));
t174 = Icges(3,4) * t134;
t133 = sin(qJ(3));
t173 = Icges(4,4) * t133;
t136 = cos(qJ(3));
t172 = Icges(4,4) * t136;
t131 = qJ(3) + qJ(4);
t126 = sin(t131);
t171 = Icges(5,4) * t126;
t127 = cos(t131);
t170 = Icges(5,4) * t127;
t169 = t126 * t134;
t168 = t126 * t137;
t132 = sin(qJ(5));
t167 = t132 * t137;
t166 = t134 * t132;
t135 = cos(qJ(5));
t165 = t134 * t135;
t164 = t135 * t137;
t163 = t136 * t137;
t162 = qJD(5) * t126;
t161 = V_base(5) * qJ(1) + V_base(1);
t160 = V_base(6) * pkin(1) + V_base(2);
t159 = pkin(2) * t134 * t136;
t155 = qJD(1) + V_base(3);
t122 = qJD(3) * t134 + V_base(4);
t123 = V_base(5) - qJD(2);
t102 = qJD(4) * t134 + t122;
t154 = t123 * t159 + t133 * t175 + t161;
t153 = rSges(4,1) * t136 - rSges(4,2) * t133;
t152 = rSges(5,1) * t127 - rSges(5,2) * t126;
t151 = Icges(4,1) * t136 - t173;
t150 = Icges(5,1) * t127 - t171;
t149 = -Icges(4,2) * t133 + t172;
t148 = -Icges(5,2) * t126 + t170;
t147 = Icges(4,5) * t136 - Icges(4,6) * t133;
t146 = Icges(5,5) * t127 - Icges(5,6) * t126;
t145 = -V_base(5) * pkin(1) + t155;
t144 = -V_base(4) * qJ(1) + t160;
t101 = V_base(6) + (-qJD(3) - qJD(4)) * t137;
t143 = (-Icges(5,3) * t137 + t134 * t146) * t101 + (Icges(5,3) * t134 + t137 * t146) * t102 + (-Icges(5,5) * t126 - Icges(5,6) * t127) * t123;
t142 = (-Icges(4,5) * t133 - Icges(4,6) * t136) * t123 + (-Icges(4,3) * t137 + t134 * t147) * t121 + (Icges(4,3) * t134 + t137 * t147) * t122;
t141 = -t122 * t159 + t163 * t175 + t144;
t140 = (-t133 * t122 - t123 * t163) * pkin(2) + t145;
t73 = -Icges(5,6) * t137 + t134 * t148;
t74 = Icges(5,6) * t134 + t137 * t148;
t75 = -Icges(5,5) * t137 + t134 * t150;
t76 = Icges(5,5) * t134 + t137 * t150;
t97 = -Icges(5,2) * t127 - t171;
t98 = -Icges(5,1) * t126 - t170;
t139 = (-t126 * t74 + t127 * t76) * t102 + (-t126 * t97 + t127 * t98) * t123 + (-t126 * t73 + t127 * t75) * t101;
t112 = -Icges(4,2) * t136 - t173;
t115 = -Icges(4,1) * t133 - t172;
t83 = -Icges(4,6) * t137 + t134 * t149;
t84 = Icges(4,6) * t134 + t137 * t149;
t85 = -Icges(4,5) * t137 + t134 * t151;
t86 = Icges(4,5) * t134 + t137 * t151;
t138 = (-t133 * t84 + t136 * t86) * t122 + (-t112 * t133 + t115 * t136) * t123 + (-t133 * t83 + t136 * t85) * t121;
t129 = Icges(3,4) * t137;
t120 = rSges(3,1) * t137 - t134 * rSges(3,2);
t119 = t134 * rSges(3,1) + rSges(3,2) * t137;
t118 = -rSges(4,1) * t133 - rSges(4,2) * t136;
t117 = Icges(3,1) * t137 - t174;
t116 = Icges(3,1) * t134 + t129;
t114 = -Icges(3,2) * t134 + t129;
t113 = Icges(3,2) * t137 + t174;
t108 = qJD(5) * t127 + t123;
t107 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t106 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t105 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t100 = -V_base(5) * rSges(2,1) + V_base(4) * rSges(2,2) + t155;
t99 = -rSges(5,1) * t126 - rSges(5,2) * t127;
t94 = V_base(6) * rSges(2,1) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t93 = -V_base(6) * rSges(2,2) + V_base(5) * rSges(2,3) + t161;
t92 = t127 * t164 + t166;
t91 = -t127 * t167 + t165;
t90 = t127 * t165 - t167;
t89 = -t127 * t166 - t164;
t88 = t134 * rSges(4,3) + t137 * t153;
t87 = -rSges(4,3) * t137 + t134 * t153;
t80 = t137 * t162 + t102;
t79 = t134 * t162 + t101;
t78 = t134 * rSges(5,3) + t137 * t152;
t77 = -rSges(5,3) * t137 + t134 * t152;
t70 = rSges(6,3) * t127 + (-rSges(6,1) * t135 + rSges(6,2) * t132) * t126;
t69 = Icges(6,5) * t127 + (-Icges(6,1) * t135 + Icges(6,4) * t132) * t126;
t68 = Icges(6,6) * t127 + (-Icges(6,4) * t135 + Icges(6,2) * t132) * t126;
t67 = Icges(6,3) * t127 + (-Icges(6,5) * t135 + Icges(6,6) * t132) * t126;
t66 = V_base(6) * rSges(3,3) + t119 * t123 + t161;
t65 = -V_base(4) * rSges(3,3) - t120 * t123 + t145;
t64 = t120 * V_base(6) + (-qJ(1) - t119) * V_base(4) + t160;
t63 = t92 * rSges(6,1) + t91 * rSges(6,2) + rSges(6,3) * t168;
t62 = rSges(6,1) * t90 + rSges(6,2) * t89 + rSges(6,3) * t169;
t61 = Icges(6,1) * t92 + Icges(6,4) * t91 + Icges(6,5) * t168;
t60 = Icges(6,1) * t90 + Icges(6,4) * t89 + Icges(6,5) * t169;
t59 = Icges(6,4) * t92 + Icges(6,2) * t91 + Icges(6,6) * t168;
t58 = Icges(6,4) * t90 + Icges(6,2) * t89 + Icges(6,6) * t169;
t57 = Icges(6,5) * t92 + Icges(6,6) * t91 + Icges(6,3) * t168;
t56 = Icges(6,5) * t90 + Icges(6,6) * t89 + Icges(6,3) * t169;
t55 = -t118 * t121 + t123 * t87 + t161;
t54 = t118 * t122 - t123 * t88 + t145;
t53 = t121 * t88 - t122 * t87 + t144;
t52 = -t101 * t99 + t123 * t77 + t154;
t51 = t102 * t99 - t123 * t78 + t140;
t50 = t101 * t78 - t102 * t77 + t141;
t49 = t108 * t62 - t70 * t79 + t154;
t48 = -t108 * t63 + t80 * t70 + t140;
t47 = -t62 * t80 + t63 * t79 + t141;
t1 = m(1) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(3) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(4) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + t122 * (t134 * t142 + t137 * t138) / 0.2e1 + t121 * (t134 * t138 - t137 * t142) / 0.2e1 + m(5) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + t102 * (t134 * t143 + t137 * t139) / 0.2e1 + t101 * (t134 * t139 - t137 * t143) / 0.2e1 + m(6) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + t80 * ((t57 * t168 + t91 * t59 + t92 * t61) * t80 + (t168 * t67 + t91 * t68 + t92 * t69) * t108 + (t168 * t56 + t91 * t58 + t92 * t60) * t79) / 0.2e1 + t108 * ((t67 * t108 + t56 * t79 + t57 * t80) * t127 + ((t132 * t59 - t135 * t61) * t80 + (t132 * t68 - t135 * t69) * t108 + (t132 * t58 - t135 * t60) * t79) * t126) / 0.2e1 + t79 * ((t169 * t57 + t59 * t89 + t61 * t90) * t80 + (t67 * t169 + t68 * t89 + t69 * t90) * t108 + (t56 * t169 + t89 * t58 + t90 * t60) * t79) / 0.2e1 + ((-t133 * t86 - t136 * t84) * t122 + (-t133 * t85 - t136 * t83) * t121 + (-t126 * t76 - t127 * t74) * t102 + (-t126 * t75 - t127 * t73) * t101 + (-t136 * t112 - t133 * t115 - t126 * t98 - t127 * t97 + Icges(3,3)) * t123) * t123 / 0.2e1 + ((-t134 * t113 + t116 * t137 + t178) * V_base(6) + (-t134 * t114 + t117 * t137 + Icges(1,1) + Icges(2,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t113 * t137 + t134 * t116 + Icges(1,3) + Icges(2,3)) * V_base(6) + (t114 * t137 + t134 * t117 + t178) * V_base(4)) * V_base(6) / 0.2e1 + V_base(4) * t123 * (-Icges(3,5) * t137 + Icges(3,6) * t134) + t123 * V_base(6) * (-Icges(3,5) * t134 - Icges(3,6) * t137) + ((Icges(1,6) + Icges(2,6)) * V_base(6) + (Icges(1,4) + Icges(2,4)) * V_base(4) + (Icges(1,2) / 0.2e1 + Icges(2,2) / 0.2e1) * V_base(5)) * V_base(5);
T = t1;
