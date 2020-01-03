% Calculate kinetic energy for
% S4RPRR6
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR6_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR6_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:27
% EndTime: 2019-12-31 16:52:28
% DurationCPUTime: 1.09s
% Computational Cost: add. (752->210), mult. (820->299), div. (0->0), fcn. (652->8), ass. (0->110)
t138 = sin(pkin(7));
t183 = pkin(2) * t138;
t137 = pkin(7) + qJ(3);
t127 = sin(t137);
t182 = pkin(3) * t127;
t139 = cos(pkin(7));
t181 = t139 * pkin(2);
t141 = sin(qJ(1));
t142 = cos(qJ(1));
t117 = t141 * pkin(1) - qJ(2) * t142;
t73 = -pkin(5) * t142 + t181 * t141;
t180 = -t117 - t73;
t179 = Icges(2,4) * t141;
t178 = Icges(3,4) * t138;
t177 = Icges(3,4) * t139;
t176 = Icges(4,4) * t127;
t128 = cos(t137);
t175 = Icges(4,4) * t128;
t129 = qJ(4) + t137;
t124 = sin(t129);
t174 = Icges(5,4) * t124;
t125 = cos(t129);
t173 = Icges(5,4) * t125;
t171 = pkin(3) * t128;
t169 = V_base(4) * t117 + V_base(3);
t168 = V_base(5) * pkin(4) + V_base(1);
t122 = qJD(3) * t141 + V_base(4);
t165 = qJD(2) * t141 + t168;
t164 = V_base(5) * t183 + t165;
t163 = rSges(3,1) * t139 - rSges(3,2) * t138;
t162 = rSges(4,1) * t128 - rSges(4,2) * t127;
t161 = rSges(5,1) * t125 - rSges(5,2) * t124;
t160 = Icges(3,1) * t139 - t178;
t159 = Icges(4,1) * t128 - t176;
t158 = Icges(5,1) * t125 - t174;
t157 = -Icges(3,2) * t138 + t177;
t156 = -Icges(4,2) * t127 + t175;
t155 = -Icges(5,2) * t124 + t173;
t154 = Icges(3,5) * t139 - Icges(3,6) * t138;
t153 = Icges(4,5) * t128 - Icges(4,6) * t127;
t152 = Icges(5,5) * t125 - Icges(5,6) * t124;
t119 = pkin(1) * t142 + t141 * qJ(2);
t130 = V_base(6) + qJD(1);
t151 = -qJD(2) * t142 + t130 * t119 + V_base(2);
t100 = V_base(5) + (-qJD(3) - qJD(4)) * t142;
t101 = qJD(4) * t141 + t122;
t150 = (-Icges(5,3) * t142 + t152 * t141) * t100 + (Icges(5,3) * t141 + t152 * t142) * t101 + (Icges(5,5) * t124 + Icges(5,6) * t125) * t130;
t121 = -qJD(3) * t142 + V_base(5);
t149 = (-Icges(4,3) * t142 + t153 * t141) * t121 + (Icges(4,3) * t141 + t153 * t142) * t122 + (Icges(4,5) * t127 + Icges(4,6) * t128) * t130;
t74 = pkin(5) * t141 + t181 * t142;
t148 = V_base(4) * t73 + (-t119 - t74) * V_base(5) + t169;
t147 = (Icges(3,5) * t138 + Icges(3,6) * t139) * t130 + (-Icges(3,3) * t142 + t154 * t141) * V_base(5) + (Icges(3,3) * t141 + t154 * t142) * V_base(4);
t146 = t130 * t74 + (-pkin(4) - t183) * V_base(4) + t151;
t65 = -Icges(5,6) * t142 + t155 * t141;
t66 = Icges(5,6) * t141 + t155 * t142;
t67 = -Icges(5,5) * t142 + t158 * t141;
t68 = Icges(5,5) * t141 + t158 * t142;
t92 = Icges(5,2) * t125 + t174;
t93 = Icges(5,1) * t124 + t173;
t145 = (-t124 * t66 + t125 * t68) * t101 + (-t124 * t65 + t125 * t67) * t100 + (-t124 * t92 + t125 * t93) * t130;
t77 = -Icges(4,6) * t142 + t156 * t141;
t78 = Icges(4,6) * t141 + t156 * t142;
t79 = -Icges(4,5) * t142 + t159 * t141;
t80 = Icges(4,5) * t141 + t159 * t142;
t97 = Icges(4,2) * t128 + t176;
t98 = Icges(4,1) * t127 + t175;
t144 = (-t127 * t78 + t128 * t80) * t122 + (-t127 * t77 + t128 * t79) * t121 + (-t127 * t97 + t128 * t98) * t130;
t108 = Icges(3,2) * t139 + t178;
t109 = Icges(3,1) * t138 + t177;
t85 = -Icges(3,6) * t142 + t157 * t141;
t86 = Icges(3,6) * t141 + t157 * t142;
t87 = -Icges(3,5) * t142 + t160 * t141;
t88 = Icges(3,5) * t141 + t160 * t142;
t143 = (-t138 * t86 + t139 * t88) * V_base(4) + (-t138 * t85 + t139 * t87) * V_base(5) + (-t108 * t138 + t109 * t139) * t130;
t134 = Icges(2,4) * t142;
t120 = rSges(2,1) * t142 - t141 * rSges(2,2);
t118 = t141 * rSges(2,1) + rSges(2,2) * t142;
t116 = Icges(2,1) * t142 - t179;
t115 = Icges(2,1) * t141 + t134;
t114 = -Icges(2,2) * t141 + t134;
t113 = Icges(2,2) * t142 + t179;
t112 = Icges(2,5) * t142 - Icges(2,6) * t141;
t111 = Icges(2,5) * t141 + Icges(2,6) * t142;
t110 = rSges(3,1) * t138 + rSges(3,2) * t139;
t106 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t105 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t104 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t99 = rSges(4,1) * t127 + rSges(4,2) * t128;
t94 = rSges(5,1) * t124 + rSges(5,2) * t125;
t90 = t141 * rSges(3,3) + t163 * t142;
t89 = -rSges(3,3) * t142 + t163 * t141;
t82 = t141 * rSges(4,3) + t162 * t142;
t81 = -rSges(4,3) * t142 + t162 * t141;
t72 = t141 * rSges(5,3) + t161 * t142;
t71 = -rSges(5,3) * t142 + t161 * t141;
t70 = V_base(5) * rSges(2,3) - t118 * t130 + t168;
t69 = t120 * t130 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t62 = t118 * V_base(4) - t120 * V_base(5) + V_base(3);
t59 = pkin(6) * t141 + t171 * t142;
t58 = -pkin(6) * t142 + t171 * t141;
t57 = t110 * V_base(5) + (-t117 - t89) * t130 + t165;
t56 = t130 * t90 + (-pkin(4) - t110) * V_base(4) + t151;
t55 = t89 * V_base(4) + (-t119 - t90) * V_base(5) + t169;
t54 = t121 * t99 + (-t81 + t180) * t130 + t164;
t53 = -t122 * t99 + t130 * t82 + t146;
t52 = -t121 * t82 + t122 * t81 + t148;
t51 = t121 * t182 + t100 * t94 + (-t58 - t71 + t180) * t130 + t164;
t50 = -t122 * t182 - t101 * t94 + (t59 + t72) * t130 + t146;
t49 = -t100 * t72 + t101 * t71 - t121 * t59 + t122 * t58 + t148;
t1 = m(1) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(2) * (t62 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(3) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(4) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + t122 * (t149 * t141 + t144 * t142) / 0.2e1 + t121 * (t144 * t141 - t149 * t142) / 0.2e1 + m(5) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + t101 * (t150 * t141 + t145 * t142) / 0.2e1 + t100 * (t145 * t141 - t150 * t142) / 0.2e1 + (t112 * t130 + t141 * t147 + t142 * t143 + (-t141 * t113 + t115 * t142 + Icges(1,4)) * V_base(5) + (-t141 * t114 + t116 * t142 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t111 * t130 + t141 * t143 - t142 * t147 + (t113 * t142 + t141 * t115 + Icges(1,2)) * V_base(5) + (t114 * t142 + t141 * t116 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t127 * t80 + t128 * t78) * t122 + (t127 * t79 + t128 * t77) * t121 + (t124 * t68 + t125 * t66) * t101 + (t124 * t67 + t125 * t65) * t100 + (t138 * t87 + t139 * t85 + t111) * V_base(5) + (t138 * t88 + t139 * t86 + t112) * V_base(4) + (t108 * t139 + t109 * t138 + t124 * t93 + t125 * t92 + t127 * t98 + t128 * t97 + Icges(2,3)) * t130) * t130 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
