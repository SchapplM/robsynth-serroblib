% Calculate kinetic energy for
% S4PRPR3
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:46
% EndTime: 2019-12-31 16:20:47
% DurationCPUTime: 0.94s
% Computational Cost: add. (693->193), mult. (628->258), div. (0->0), fcn. (458->8), ass. (0->99)
t125 = cos(pkin(6));
t162 = pkin(1) * t125;
t122 = sin(pkin(7));
t161 = pkin(3) * t122;
t124 = cos(pkin(7));
t160 = pkin(3) * t124;
t159 = -pkin(4) - qJ(1);
t123 = sin(pkin(6));
t158 = Icges(2,4) * t123;
t121 = pkin(6) + qJ(2);
t113 = sin(t121);
t157 = Icges(3,4) * t113;
t156 = Icges(4,4) * t122;
t155 = Icges(4,4) * t124;
t120 = pkin(7) + qJ(4);
t112 = sin(t120);
t154 = Icges(5,4) * t112;
t114 = cos(t120);
t153 = Icges(5,4) * t114;
t145 = pkin(1) * V_base(6);
t151 = t125 * t145 + V_base(2);
t150 = V_base(5) * qJ(1) + V_base(1);
t146 = qJD(1) + V_base(3);
t115 = cos(t121);
t88 = t115 * pkin(2) + t113 * qJ(3);
t144 = -t88 - t162;
t143 = V_base(4) * t123 * pkin(1) + t146;
t86 = t113 * pkin(2) - t115 * qJ(3);
t142 = V_base(4) * t86 + t143;
t141 = rSges(4,1) * t124 - rSges(4,2) * t122;
t140 = rSges(5,1) * t114 - rSges(5,2) * t112;
t139 = Icges(4,1) * t124 - t156;
t138 = Icges(5,1) * t114 - t154;
t137 = -Icges(4,2) * t122 + t155;
t136 = -Icges(5,2) * t112 + t153;
t135 = Icges(4,5) * t124 - Icges(4,6) * t122;
t134 = Icges(5,5) * t114 - Icges(5,6) * t112;
t117 = V_base(6) + qJD(2);
t133 = -qJD(3) * t115 + t117 * t88 + t151;
t102 = -qJD(4) * t115 + V_base(5);
t103 = qJD(4) * t113 + V_base(4);
t132 = t102 * (-Icges(5,3) * t115 + t113 * t134) + t103 * (Icges(5,3) * t113 + t115 * t134) + t117 * (Icges(5,5) * t112 + Icges(5,6) * t114);
t131 = V_base(5) * pkin(4) - t123 * t145 + t150;
t130 = qJD(3) * t113 + t131;
t129 = t117 * (Icges(4,5) * t122 + Icges(4,6) * t124) + (-Icges(4,3) * t115 + t113 * t135) * V_base(5) + (Icges(4,3) * t113 + t115 * t135) * V_base(4);
t58 = -Icges(5,6) * t115 + t113 * t136;
t59 = Icges(5,6) * t113 + t115 * t136;
t60 = -Icges(5,5) * t115 + t113 * t138;
t61 = Icges(5,5) * t113 + t115 * t138;
t79 = Icges(5,2) * t114 + t154;
t82 = Icges(5,1) * t112 + t153;
t128 = (-t112 * t59 + t114 * t61) * t103 + (-t112 * t58 + t114 * t60) * t102 + (-t112 * t79 + t114 * t82) * t117;
t66 = -Icges(4,6) * t115 + t113 * t137;
t67 = Icges(4,6) * t113 + t115 * t137;
t68 = -Icges(4,5) * t115 + t113 * t139;
t69 = Icges(4,5) * t113 + t115 * t139;
t96 = Icges(4,2) * t124 + t156;
t99 = Icges(4,1) * t122 + t155;
t127 = (-t122 * t67 + t124 * t69) * V_base(4) + (-t122 * t66 + t124 * t68) * V_base(5) + (-t122 * t96 + t124 * t99) * t117;
t116 = Icges(2,4) * t125;
t110 = Icges(3,4) * t115;
t106 = rSges(2,1) * t125 - rSges(2,2) * t123;
t105 = rSges(2,1) * t123 + rSges(2,2) * t125;
t104 = rSges(4,1) * t122 + rSges(4,2) * t124;
t101 = Icges(2,1) * t125 - t158;
t100 = Icges(2,1) * t123 + t116;
t98 = -Icges(2,2) * t123 + t116;
t97 = Icges(2,2) * t125 + t158;
t92 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t91 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t90 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t89 = rSges(3,1) * t115 - rSges(3,2) * t113;
t87 = rSges(3,1) * t113 + rSges(3,2) * t115;
t85 = rSges(5,1) * t112 + rSges(5,2) * t114;
t84 = Icges(3,1) * t115 - t157;
t83 = Icges(3,1) * t113 + t110;
t81 = -Icges(3,2) * t113 + t110;
t80 = Icges(3,2) * t115 + t157;
t78 = Icges(3,5) * t115 - Icges(3,6) * t113;
t77 = Icges(3,5) * t113 + Icges(3,6) * t115;
t73 = V_base(5) * rSges(2,3) - t105 * V_base(6) + t150;
t72 = t106 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t71 = rSges(4,3) * t113 + t115 * t141;
t70 = -rSges(4,3) * t115 + t113 * t141;
t63 = rSges(5,3) * t113 + t115 * t140;
t62 = -rSges(5,3) * t115 + t113 * t140;
t55 = t105 * V_base(4) - t106 * V_base(5) + t146;
t54 = pkin(5) * t113 + t115 * t160;
t53 = -pkin(5) * t115 + t113 * t160;
t52 = V_base(5) * rSges(3,3) - t117 * t87 + t131;
t51 = t117 * t89 + (-rSges(3,3) + t159) * V_base(4) + t151;
t50 = t87 * V_base(4) + (-t89 - t162) * V_base(5) + t143;
t49 = t104 * V_base(5) + (-t70 - t86) * t117 + t130;
t48 = t117 * t71 + (-t104 + t159) * V_base(4) + t133;
t47 = t70 * V_base(4) + (t144 - t71) * V_base(5) + t142;
t46 = V_base(5) * t161 + t102 * t85 + (-t53 - t62 - t86) * t117 + t130;
t45 = -t103 * t85 + (t54 + t63) * t117 + (t159 - t161) * V_base(4) + t133;
t44 = -t102 * t63 + t103 * t62 + t53 * V_base(4) + (t144 - t54) * V_base(5) + t142;
t1 = m(1) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(2) * (t55 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(3) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + m(5) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + t103 * (t132 * t113 + t128 * t115) / 0.2e1 + t102 * (t128 * t113 - t132 * t115) / 0.2e1 + ((t112 * t61 + t114 * t59) * t103 + (t112 * t60 + t114 * t58) * t102 + (t122 * t68 + t124 * t66 + t77) * V_base(5) + (t122 * t69 + t124 * t67 + t78) * V_base(4) + (t112 * t82 + t114 * t79 + t122 * t99 + t124 * t96 + Icges(3,3)) * t117) * t117 / 0.2e1 + (t129 * t113 + t127 * t115 + t78 * t117 + (t100 * t125 - t113 * t80 + t115 * t83 - t123 * t97 + Icges(1,4)) * V_base(5) + (t125 * t101 - t113 * t81 + t115 * t84 - t123 * t98 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t127 * t113 - t129 * t115 + t77 * t117 + (t123 * t100 + t113 * t83 + t115 * t80 + t125 * t97 + Icges(1,2)) * V_base(5) + (t101 * t123 + t113 * t84 + t115 * t81 + t125 * t98 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t123 + Icges(2,6) * t125 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t125 - Icges(2,6) * t123 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
