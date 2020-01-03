% Calculate kinetic energy for
% S4RPPR3
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
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:48
% EndTime: 2019-12-31 16:37:49
% DurationCPUTime: 0.83s
% Computational Cost: add. (704->194), mult. (628->255), div. (0->0), fcn. (458->8), ass. (0->99)
t121 = qJ(1) + pkin(6);
t113 = sin(t121);
t115 = cos(t121);
t125 = sin(qJ(1));
t126 = cos(qJ(1));
t164 = Icges(2,5) * t125 + Icges(3,5) * t113 + Icges(2,6) * t126 + Icges(3,6) * t115;
t163 = Icges(2,5) * t126 + Icges(3,5) * t115 - Icges(2,6) * t125 - Icges(3,6) * t113;
t161 = pkin(1) * t125;
t160 = pkin(1) * t126;
t122 = sin(pkin(7));
t159 = pkin(3) * t122;
t123 = cos(pkin(7));
t158 = pkin(3) * t123;
t157 = -pkin(4) - qJ(2);
t156 = Icges(2,4) * t125;
t155 = Icges(3,4) * t113;
t154 = Icges(4,4) * t122;
t153 = Icges(4,4) * t123;
t120 = pkin(7) + qJ(4);
t112 = sin(t120);
t152 = Icges(5,4) * t112;
t114 = cos(t120);
t151 = Icges(5,4) * t114;
t116 = V_base(6) + qJD(1);
t149 = t116 * t160 + V_base(2);
t148 = V_base(5) * pkin(4) + V_base(1);
t86 = pkin(2) * t113 - qJ(3) * t115;
t145 = -t86 - t161;
t88 = pkin(2) * t115 + qJ(3) * t113;
t144 = -t88 - t160;
t143 = V_base(5) * qJ(2) + t148;
t142 = V_base(4) * t161 + qJD(2) + V_base(3);
t141 = V_base(4) * t86 + t142;
t140 = qJD(3) * t113 + t143;
t139 = rSges(4,1) * t123 - rSges(4,2) * t122;
t138 = rSges(5,1) * t114 - rSges(5,2) * t112;
t137 = Icges(4,1) * t123 - t154;
t136 = Icges(5,1) * t114 - t152;
t135 = -Icges(4,2) * t122 + t153;
t134 = -Icges(5,2) * t112 + t151;
t133 = Icges(4,5) * t123 - Icges(4,6) * t122;
t132 = Icges(5,5) * t114 - Icges(5,6) * t112;
t131 = -qJD(3) * t115 + t116 * t88 + t149;
t96 = -qJD(4) * t115 + V_base(5);
t97 = qJD(4) * t113 + V_base(4);
t130 = t116 * (Icges(5,5) * t112 + Icges(5,6) * t114) + (-Icges(5,3) * t115 + t132 * t113) * t96 + (Icges(5,3) * t113 + t132 * t115) * t97;
t129 = t116 * (Icges(4,5) * t122 + Icges(4,6) * t123) + (-Icges(4,3) * t115 + t133 * t113) * V_base(5) + (Icges(4,3) * t113 + t133 * t115) * V_base(4);
t57 = -Icges(5,6) * t115 + t134 * t113;
t58 = Icges(5,6) * t113 + t134 * t115;
t59 = -Icges(5,5) * t115 + t136 * t113;
t60 = Icges(5,5) * t113 + t136 * t115;
t79 = Icges(5,2) * t114 + t152;
t82 = Icges(5,1) * t112 + t151;
t128 = (-t112 * t58 + t114 * t60) * t97 + (-t112 * t57 + t114 * t59) * t96 + (-t112 * t79 + t114 * t82) * t116;
t66 = -Icges(4,6) * t115 + t135 * t113;
t67 = Icges(4,6) * t113 + t135 * t115;
t68 = -Icges(4,5) * t115 + t137 * t113;
t69 = Icges(4,5) * t113 + t137 * t115;
t94 = Icges(4,2) * t123 + t154;
t95 = Icges(4,1) * t122 + t153;
t127 = (-t122 * t67 + t123 * t69) * V_base(4) + (-t122 * t66 + t123 * t68) * V_base(5) + (-t122 * t94 + t123 * t95) * t116;
t118 = Icges(2,4) * t126;
t110 = Icges(3,4) * t115;
t106 = rSges(2,1) * t126 - t125 * rSges(2,2);
t105 = t125 * rSges(2,1) + rSges(2,2) * t126;
t104 = Icges(2,1) * t126 - t156;
t103 = Icges(2,1) * t125 + t118;
t102 = -Icges(2,2) * t125 + t118;
t101 = Icges(2,2) * t126 + t156;
t98 = rSges(4,1) * t122 + rSges(4,2) * t123;
t92 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t91 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t90 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t89 = rSges(3,1) * t115 - rSges(3,2) * t113;
t87 = rSges(3,1) * t113 + rSges(3,2) * t115;
t85 = rSges(5,1) * t112 + rSges(5,2) * t114;
t84 = Icges(3,1) * t115 - t155;
t83 = Icges(3,1) * t113 + t110;
t81 = -Icges(3,2) * t113 + t110;
t80 = Icges(3,2) * t115 + t155;
t73 = V_base(5) * rSges(2,3) - t105 * t116 + t148;
t72 = t106 * t116 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t71 = rSges(4,3) * t113 + t139 * t115;
t70 = -rSges(4,3) * t115 + t139 * t113;
t63 = t105 * V_base(4) - t106 * V_base(5) + V_base(3);
t62 = rSges(5,3) * t113 + t138 * t115;
t61 = -rSges(5,3) * t115 + t138 * t113;
t54 = pkin(5) * t113 + t158 * t115;
t53 = -pkin(5) * t115 + t158 * t113;
t52 = V_base(5) * rSges(3,3) + (-t87 - t161) * t116 + t143;
t51 = t116 * t89 + (-rSges(3,3) + t157) * V_base(4) + t149;
t50 = V_base(4) * t87 + (-t89 - t160) * V_base(5) + t142;
t49 = t98 * V_base(5) + (t145 - t70) * t116 + t140;
t48 = t116 * t71 + (-t98 + t157) * V_base(4) + t131;
t47 = V_base(4) * t70 + (t144 - t71) * V_base(5) + t141;
t46 = V_base(5) * t159 + t85 * t96 + (t145 - t53 - t61) * t116 + t140;
t45 = -t85 * t97 + (t54 + t62) * t116 + (t157 - t159) * V_base(4) + t131;
t44 = V_base(4) * t53 + t97 * t61 - t96 * t62 + (t144 - t54) * V_base(5) + t141;
t1 = m(1) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(2) * (t63 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(3) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + m(5) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + t97 * (t130 * t113 + t128 * t115) / 0.2e1 + t96 * (t128 * t113 - t130 * t115) / 0.2e1 + ((t112 * t60 + t114 * t58) * t97 + (t112 * t59 + t114 * t57) * t96 + (t122 * t68 + t123 * t66 + t164) * V_base(5) + (t122 * t69 + t123 * t67 + t163) * V_base(4) + (t112 * t82 + t114 * t79 + t122 * t95 + t123 * t94 + Icges(2,3) + Icges(3,3)) * t116) * t116 / 0.2e1 + (t129 * t113 + t127 * t115 + t163 * t116 + (-t125 * t101 + t103 * t126 - t113 * t80 + t115 * t83 + Icges(1,4)) * V_base(5) + (-t125 * t102 + t126 * t104 - t113 * t81 + t115 * t84 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t127 * t113 - t129 * t115 + t164 * t116 + (t126 * t101 + t125 * t103 + t113 * t83 + t115 * t80 + Icges(1,2)) * V_base(5) + (t102 * t126 + t125 * t104 + t113 * t84 + t115 * t81 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
