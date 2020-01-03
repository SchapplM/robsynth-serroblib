% Calculate kinetic energy for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:34
% EndTime: 2019-12-31 16:41:35
% DurationCPUTime: 1.05s
% Computational Cost: add. (457->183), mult. (632->236), div. (0->0), fcn. (464->6), ass. (0->95)
t173 = Icges(2,4) + Icges(3,6);
t172 = Icges(2,1) + Icges(3,2);
t171 = -Icges(3,4) + Icges(2,5);
t170 = Icges(3,5) - Icges(2,6);
t169 = Icges(2,2) + Icges(3,3);
t118 = cos(qJ(1));
t168 = t173 * t118;
t117 = sin(qJ(1));
t167 = t173 * t117;
t166 = t172 * t117 + t168;
t165 = t172 * t118 - t167;
t164 = t171 * t117 - t170 * t118;
t163 = t170 * t117 + t171 * t118;
t114 = sin(pkin(6));
t115 = cos(pkin(6));
t154 = Icges(4,4) * t114;
t127 = Icges(4,2) * t115 + t154;
t62 = Icges(4,6) * t117 - t127 * t118;
t153 = Icges(4,4) * t115;
t129 = Icges(4,1) * t114 + t153;
t64 = Icges(4,5) * t117 - t129 * t118;
t162 = t114 * t64 + t115 * t62 - t169 * t118 - t167;
t61 = Icges(4,6) * t118 + t127 * t117;
t63 = Icges(4,5) * t118 + t129 * t117;
t161 = t114 * t63 + t115 * t61 + t169 * t117 - t168;
t101 = qJD(4) * t117 + V_base(5);
t102 = qJD(4) * t118 + V_base(4);
t113 = pkin(6) + qJ(4);
t104 = sin(t113);
t105 = cos(t113);
t106 = V_base(6) + qJD(1);
t152 = Icges(5,4) * t104;
t126 = Icges(5,2) * t105 + t152;
t53 = Icges(5,6) * t118 + t126 * t117;
t54 = Icges(5,6) * t117 - t126 * t118;
t151 = Icges(5,4) * t105;
t128 = Icges(5,1) * t104 + t151;
t55 = Icges(5,5) * t118 + t128 * t117;
t56 = Icges(5,5) * t117 - t128 * t118;
t71 = -Icges(5,2) * t104 + t151;
t72 = Icges(5,1) * t105 - t152;
t160 = (t104 * t55 + t105 * t53) * t102 + (t104 * t56 + t105 * t54) * t101 + (t104 * t72 + t105 * t71) * t106;
t159 = -pkin(2) - pkin(4);
t157 = pkin(3) * t114;
t156 = pkin(3) * t115;
t148 = qJ(3) * t118;
t147 = t117 * qJ(3);
t94 = t117 * pkin(1) - qJ(2) * t118;
t145 = V_base(4) * t94 + V_base(3);
t144 = V_base(5) * pkin(4) + V_base(1);
t141 = V_base(4) * t147 + t145;
t140 = -t94 - t147;
t97 = pkin(1) * t118 + t117 * qJ(2);
t139 = -t97 - t148;
t138 = qJD(2) * t117 + t144;
t137 = rSges(4,1) * t114 + rSges(4,2) * t115;
t136 = rSges(5,1) * t104 + rSges(5,2) * t105;
t79 = -Icges(4,2) * t114 + t153;
t80 = Icges(4,1) * t115 - t154;
t130 = t114 * t80 + t115 * t79;
t125 = Icges(4,5) * t114 + Icges(4,6) * t115;
t124 = Icges(5,5) * t104 + Icges(5,6) * t105;
t123 = -qJD(2) * t118 + t106 * t97 + V_base(2);
t122 = V_base(5) * pkin(2) + qJD(3) * t118 + t138;
t121 = qJD(3) * t117 + t106 * t148 + t123;
t120 = (Icges(5,3) * t117 - t124 * t118) * t101 + (Icges(5,3) * t118 + t124 * t117) * t102 + t106 * (Icges(5,5) * t105 - Icges(5,6) * t104);
t119 = t106 * (Icges(4,5) * t115 - Icges(4,6) * t114) + (Icges(4,3) * t118 + t125 * t117) * V_base(4) + (Icges(4,3) * t117 - t125 * t118) * V_base(5);
t99 = rSges(2,1) * t118 - t117 * rSges(2,2);
t98 = -rSges(3,2) * t118 + t117 * rSges(3,3);
t96 = t117 * rSges(2,1) + rSges(2,2) * t118;
t95 = -t117 * rSges(3,2) - rSges(3,3) * t118;
t81 = rSges(4,1) * t115 - rSges(4,2) * t114;
t77 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t76 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t75 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t73 = rSges(5,1) * t105 - rSges(5,2) * t104;
t68 = pkin(5) * t118 + t117 * t157;
t67 = pkin(5) * t117 - t118 * t157;
t66 = t117 * rSges(4,3) - t137 * t118;
t65 = rSges(4,3) * t118 + t137 * t117;
t58 = t117 * rSges(5,3) - t136 * t118;
t57 = rSges(5,3) * t118 + t136 * t117;
t50 = V_base(5) * rSges(2,3) - t106 * t96 + t144;
t49 = t106 * t99 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t48 = t96 * V_base(4) - t99 * V_base(5) + V_base(3);
t47 = V_base(5) * rSges(3,1) + (-t94 - t95) * t106 + t138;
t46 = t106 * t98 + (-rSges(3,1) - pkin(4)) * V_base(4) + t123;
t45 = t95 * V_base(4) + (-t97 - t98) * V_base(5) + t145;
t44 = t81 * V_base(5) + (t140 - t66) * t106 + t122;
t43 = t106 * t65 + (-t81 + t159) * V_base(4) + t121;
t42 = V_base(4) * t66 + (t139 - t65) * V_base(5) + t141;
t41 = V_base(5) * t156 + t101 * t73 + (t140 - t58 - t67) * t106 + t122;
t40 = -t102 * t73 + (t57 + t68) * t106 + (-t156 + t159) * V_base(4) + t121;
t39 = -t101 * t57 + t102 * t58 + V_base(4) * t67 + (t139 - t68) * V_base(5) + t141;
t1 = m(1) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(2) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + m(3) * (t45 ^ 2 + t46 ^ 2 + t47 ^ 2) / 0.2e1 + m(4) * (t42 ^ 2 + t43 ^ 2 + t44 ^ 2) / 0.2e1 + m(5) * (t39 ^ 2 + t40 ^ 2 + t41 ^ 2) / 0.2e1 + t102 * (t160 * t117 + t120 * t118) / 0.2e1 + t101 * (t120 * t117 - t160 * t118) / 0.2e1 + ((-t104 * t53 + t105 * t55) * t102 + (-t104 * t54 + t105 * t56) * t101 + (-t114 * t62 + t115 * t64 + t164) * V_base(5) + (-t114 * t61 + t115 * t63 + t163) * V_base(4) + (-t104 * t71 + t105 * t72 - t114 * t79 + t115 * t80 + Icges(3,1) + Icges(2,3)) * t106) * t106 / 0.2e1 + (t119 * t118 + (t130 * t117 + t163) * t106 + (t162 * t117 + t166 * t118 + Icges(1,4)) * V_base(5) + (t161 * t117 + t165 * t118 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t119 * t117 + (-t130 * t118 + t164) * t106 + (t166 * t117 - t162 * t118 + Icges(1,2)) * V_base(5) + (t165 * t117 - t161 * t118 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
