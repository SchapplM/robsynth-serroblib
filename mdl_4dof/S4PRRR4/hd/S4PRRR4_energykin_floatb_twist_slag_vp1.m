% Calculate kinetic energy for
% S4PRRR4
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:28
% EndTime: 2019-12-31 16:32:29
% DurationCPUTime: 0.90s
% Computational Cost: add. (741->190), mult. (652->266), div. (0->0), fcn. (482->8), ass. (0->97)
t125 = cos(pkin(7));
t162 = pkin(1) * t125;
t126 = sin(qJ(3));
t161 = pkin(3) * t126;
t127 = cos(qJ(3));
t160 = pkin(3) * t127;
t159 = -pkin(4) - qJ(1);
t124 = sin(pkin(7));
t157 = Icges(2,4) * t124;
t122 = pkin(7) + qJ(2);
t114 = sin(t122);
t156 = Icges(3,4) * t114;
t155 = Icges(4,4) * t126;
t154 = Icges(4,4) * t127;
t123 = qJ(3) + qJ(4);
t118 = sin(t123);
t153 = Icges(5,4) * t118;
t119 = cos(t123);
t152 = Icges(5,4) * t119;
t145 = pkin(1) * V_base(6);
t151 = t125 * t145 + V_base(2);
t150 = V_base(5) * qJ(1) + V_base(1);
t146 = qJD(1) + V_base(3);
t102 = qJD(3) * t114 + V_base(4);
t144 = V_base(4) * t124 * pkin(1) + t146;
t143 = rSges(4,1) * t127 - rSges(4,2) * t126;
t142 = rSges(5,1) * t119 - rSges(5,2) * t118;
t141 = Icges(4,1) * t127 - t155;
t140 = Icges(5,1) * t119 - t153;
t139 = -Icges(4,2) * t126 + t154;
t138 = -Icges(5,2) * t118 + t152;
t137 = Icges(4,5) * t127 - Icges(4,6) * t126;
t136 = Icges(5,5) * t119 - Icges(5,6) * t118;
t115 = cos(t122);
t117 = V_base(6) + qJD(2);
t75 = V_base(5) + (-qJD(3) - qJD(4)) * t115;
t76 = qJD(4) * t114 + t102;
t135 = (Icges(5,5) * t118 + Icges(5,6) * t119) * t117 + (-Icges(5,3) * t115 + t114 * t136) * t75 + (Icges(5,3) * t114 + t115 * t136) * t76;
t101 = -qJD(3) * t115 + V_base(5);
t134 = (-Icges(4,3) * t115 + t114 * t137) * t101 + (Icges(4,3) * t114 + t115 * t137) * t102 + (Icges(4,5) * t126 + Icges(4,6) * t127) * t117;
t133 = V_base(5) * pkin(4) - t124 * t145 + t150;
t87 = t115 * pkin(2) + pkin(5) * t114;
t132 = t117 * t87 + t159 * V_base(4) + t151;
t86 = t114 * pkin(2) - t115 * pkin(5);
t131 = V_base(4) * t86 + (-t87 - t162) * V_base(5) + t144;
t58 = -Icges(5,6) * t115 + t114 * t138;
t59 = Icges(5,6) * t114 + t115 * t138;
t60 = -Icges(5,5) * t115 + t114 * t140;
t61 = Icges(5,5) * t114 + t115 * t140;
t89 = Icges(5,2) * t119 + t153;
t90 = Icges(5,1) * t118 + t152;
t130 = (-t118 * t59 + t119 * t61) * t76 + (-t118 * t58 + t119 * t60) * t75 + (-t118 * t89 + t119 * t90) * t117;
t106 = Icges(4,2) * t127 + t155;
t107 = Icges(4,1) * t126 + t154;
t66 = -Icges(4,6) * t115 + t114 * t139;
t67 = Icges(4,6) * t114 + t115 * t139;
t68 = -Icges(4,5) * t115 + t114 * t141;
t69 = Icges(4,5) * t114 + t115 * t141;
t129 = (-t126 * t67 + t127 * t69) * t102 + (-t126 * t66 + t127 * t68) * t101 + (-t106 * t126 + t107 * t127) * t117;
t116 = Icges(2,4) * t125;
t112 = Icges(3,4) * t115;
t108 = rSges(4,1) * t126 + rSges(4,2) * t127;
t104 = rSges(2,1) * t125 - rSges(2,2) * t124;
t103 = rSges(2,1) * t124 + rSges(2,2) * t125;
t100 = Icges(2,1) * t125 - t157;
t99 = Icges(2,1) * t124 + t116;
t98 = -Icges(2,2) * t124 + t116;
t97 = Icges(2,2) * t125 + t157;
t94 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t93 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t92 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t91 = rSges(5,1) * t118 + rSges(5,2) * t119;
t85 = rSges(3,1) * t115 - rSges(3,2) * t114;
t84 = rSges(3,1) * t114 + rSges(3,2) * t115;
t83 = Icges(3,1) * t115 - t156;
t82 = Icges(3,1) * t114 + t112;
t81 = -Icges(3,2) * t114 + t112;
t80 = Icges(3,2) * t115 + t156;
t73 = V_base(5) * rSges(2,3) - t103 * V_base(6) + t150;
t72 = t104 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t71 = rSges(4,3) * t114 + t115 * t143;
t70 = -rSges(4,3) * t115 + t114 * t143;
t63 = rSges(5,3) * t114 + t115 * t142;
t62 = -rSges(5,3) * t115 + t114 * t142;
t55 = t103 * V_base(4) - t104 * V_base(5) + t146;
t54 = pkin(6) * t114 + t115 * t160;
t53 = -pkin(6) * t115 + t114 * t160;
t52 = V_base(5) * rSges(3,3) - t117 * t84 + t133;
t51 = t117 * t85 + (-rSges(3,3) + t159) * V_base(4) + t151;
t50 = t84 * V_base(4) + (-t85 - t162) * V_base(5) + t144;
t49 = t101 * t108 + (-t70 - t86) * t117 + t133;
t48 = -t102 * t108 + t117 * t71 + t132;
t47 = -t101 * t71 + t102 * t70 + t131;
t46 = t101 * t161 + t75 * t91 + (-t53 - t62 - t86) * t117 + t133;
t45 = -t102 * t161 - t76 * t91 + (t54 + t63) * t117 + t132;
t44 = -t101 * t54 + t102 * t53 + t62 * t76 - t63 * t75 + t131;
t1 = m(1) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(2) * (t55 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(3) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + t102 * (t134 * t114 + t129 * t115) / 0.2e1 + t101 * (t129 * t114 - t134 * t115) / 0.2e1 + m(5) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + t76 * (t135 * t114 + t130 * t115) / 0.2e1 + t75 * (t130 * t114 - t135 * t115) / 0.2e1 + ((t126 * t69 + t127 * t67) * t102 + (t126 * t68 + t127 * t66) * t101 + (t118 * t61 + t119 * t59) * t76 + (t118 * t60 + t119 * t58) * t75 + (t127 * t106 + t126 * t107 + t118 * t90 + t119 * t89 + Icges(3,3)) * t117) * t117 / 0.2e1 + ((-t114 * t80 + t115 * t82 - t124 * t97 + t125 * t99 + Icges(1,4)) * V_base(5) + (t125 * t100 - t114 * t81 + t115 * t83 - t124 * t98 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t114 * t82 + t115 * t80 + t124 * t99 + t125 * t97 + Icges(1,2)) * V_base(5) + (t100 * t124 + t114 * t83 + t115 * t81 + t125 * t98 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t117 * (Icges(3,5) * t115 - Icges(3,6) * t114) + V_base(5) * t117 * (Icges(3,5) * t114 + Icges(3,6) * t115) + ((Icges(2,5) * t124 + Icges(2,6) * t125 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t125 - Icges(2,6) * t124 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
