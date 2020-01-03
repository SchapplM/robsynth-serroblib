% Calculate kinetic energy for
% S4RPRR3
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:03
% EndTime: 2019-12-31 16:49:04
% DurationCPUTime: 0.85s
% Computational Cost: add. (752->191), mult. (652->263), div. (0->0), fcn. (482->8), ass. (0->97)
t125 = sin(qJ(1));
t161 = pkin(1) * t125;
t127 = cos(qJ(1));
t160 = pkin(1) * t127;
t124 = sin(qJ(3));
t159 = pkin(3) * t124;
t126 = cos(qJ(3));
t158 = pkin(3) * t126;
t157 = -pkin(4) - qJ(2);
t155 = Icges(2,4) * t125;
t122 = qJ(1) + pkin(7);
t114 = sin(t122);
t154 = Icges(3,4) * t114;
t153 = Icges(4,4) * t124;
t152 = Icges(4,4) * t126;
t123 = qJ(3) + qJ(4);
t117 = sin(t123);
t151 = Icges(5,4) * t117;
t118 = cos(t123);
t150 = Icges(5,4) * t118;
t116 = V_base(6) + qJD(1);
t149 = t116 * t160 + V_base(2);
t148 = V_base(5) * pkin(4) + V_base(1);
t96 = qJD(3) * t114 + V_base(4);
t115 = cos(t122);
t86 = t114 * pkin(2) - t115 * pkin(5);
t145 = -t86 - t161;
t144 = V_base(5) * qJ(2) + t148;
t143 = V_base(4) * t161 + qJD(2) + V_base(3);
t142 = rSges(4,1) * t126 - rSges(4,2) * t124;
t141 = rSges(5,1) * t118 - rSges(5,2) * t117;
t140 = Icges(4,1) * t126 - t153;
t139 = Icges(5,1) * t118 - t151;
t138 = -Icges(4,2) * t124 + t152;
t137 = -Icges(5,2) * t117 + t150;
t136 = Icges(4,5) * t126 - Icges(4,6) * t124;
t135 = Icges(5,5) * t118 - Icges(5,6) * t117;
t75 = V_base(5) + (-qJD(3) - qJD(4)) * t115;
t76 = qJD(4) * t114 + t96;
t134 = (Icges(5,5) * t117 + Icges(5,6) * t118) * t116 + (-Icges(5,3) * t115 + t114 * t135) * t75 + (Icges(5,3) * t114 + t115 * t135) * t76;
t95 = -qJD(3) * t115 + V_base(5);
t133 = (Icges(4,5) * t124 + Icges(4,6) * t126) * t116 + (-Icges(4,3) * t115 + t114 * t136) * t95 + (Icges(4,3) * t114 + t115 * t136) * t96;
t87 = t115 * pkin(2) + t114 * pkin(5);
t132 = t116 * t87 + t157 * V_base(4) + t149;
t131 = V_base(4) * t86 + (-t87 - t160) * V_base(5) + t143;
t57 = -Icges(5,6) * t115 + t114 * t137;
t58 = Icges(5,6) * t114 + t115 * t137;
t59 = -Icges(5,5) * t115 + t114 * t139;
t60 = Icges(5,5) * t114 + t115 * t139;
t89 = Icges(5,2) * t118 + t151;
t90 = Icges(5,1) * t117 + t150;
t130 = (-t117 * t58 + t118 * t60) * t76 + (-t117 * t57 + t118 * t59) * t75 + (-t117 * t89 + t118 * t90) * t116;
t100 = Icges(4,2) * t126 + t153;
t103 = Icges(4,1) * t124 + t152;
t66 = -Icges(4,6) * t115 + t114 * t138;
t67 = Icges(4,6) * t114 + t115 * t138;
t68 = -Icges(4,5) * t115 + t114 * t140;
t69 = Icges(4,5) * t114 + t115 * t140;
t129 = (-t124 * t67 + t126 * t69) * t96 + (-t124 * t66 + t126 * t68) * t95 + (-t100 * t124 + t103 * t126) * t116;
t120 = Icges(2,4) * t127;
t112 = Icges(3,4) * t115;
t108 = rSges(2,1) * t127 - rSges(2,2) * t125;
t107 = rSges(2,1) * t125 + rSges(2,2) * t127;
t106 = rSges(4,1) * t124 + rSges(4,2) * t126;
t105 = Icges(2,1) * t127 - t155;
t104 = Icges(2,1) * t125 + t120;
t102 = -Icges(2,2) * t125 + t120;
t101 = Icges(2,2) * t127 + t155;
t94 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t93 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t92 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t91 = rSges(5,1) * t117 + rSges(5,2) * t118;
t85 = rSges(3,1) * t115 - rSges(3,2) * t114;
t84 = rSges(3,1) * t114 + rSges(3,2) * t115;
t83 = Icges(3,1) * t115 - t154;
t82 = Icges(3,1) * t114 + t112;
t81 = -Icges(3,2) * t114 + t112;
t80 = Icges(3,2) * t115 + t154;
t73 = rSges(4,3) * t114 + t115 * t142;
t72 = -rSges(4,3) * t115 + t114 * t142;
t71 = V_base(5) * rSges(2,3) - t107 * t116 + t148;
t70 = t108 * t116 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t63 = t107 * V_base(4) - t108 * V_base(5) + V_base(3);
t62 = rSges(5,3) * t114 + t115 * t141;
t61 = -rSges(5,3) * t115 + t114 * t141;
t54 = pkin(6) * t114 + t115 * t158;
t53 = -pkin(6) * t115 + t114 * t158;
t52 = V_base(5) * rSges(3,3) + (-t84 - t161) * t116 + t144;
t51 = t116 * t85 + (-rSges(3,3) + t157) * V_base(4) + t149;
t50 = t84 * V_base(4) + (-t85 - t160) * V_base(5) + t143;
t49 = t106 * t95 + (t145 - t72) * t116 + t144;
t48 = -t106 * t96 + t116 * t73 + t132;
t47 = t72 * t96 - t73 * t95 + t131;
t46 = t95 * t159 + t75 * t91 + (t145 - t53 - t61) * t116 + t144;
t45 = -t96 * t159 - t76 * t91 + (t54 + t62) * t116 + t132;
t44 = t53 * t96 - t54 * t95 + t61 * t76 - t62 * t75 + t131;
t1 = m(1) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(2) * (t63 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(3) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + t96 * (t133 * t114 + t129 * t115) / 0.2e1 + t95 * (t129 * t114 - t133 * t115) / 0.2e1 + m(5) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + t76 * (t134 * t114 + t130 * t115) / 0.2e1 + t75 * (t130 * t114 - t134 * t115) / 0.2e1 + ((-t101 * t125 + t104 * t127 - t114 * t80 + t115 * t82 + Icges(1,4)) * V_base(5) + (-t125 * t102 + t127 * t105 - t114 * t81 + t115 * t83 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t127 * t101 + t125 * t104 + t114 * t82 + t115 * t80 + Icges(1,2)) * V_base(5) + (t102 * t127 + t105 * t125 + t114 * t83 + t115 * t81 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t124 * t69 + t126 * t67) * t96 + (t124 * t68 + t126 * t66) * t95 + (t117 * t60 + t118 * t58) * t76 + (t117 * t59 + t118 * t57) * t75 + (t126 * t100 + t124 * t103 + t117 * t90 + t118 * t89 + Icges(2,3) + Icges(3,3)) * t116) * t116 / 0.2e1 + t116 * V_base(5) * (Icges(2,5) * t125 + Icges(3,5) * t114 + Icges(2,6) * t127 + Icges(3,6) * t115) + t116 * V_base(4) * (Icges(2,5) * t127 + Icges(3,5) * t115 - Icges(2,6) * t125 - Icges(3,6) * t114) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
