% Calculate kinetic energy for
% S4RPPR5
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR5_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:41
% EndTime: 2019-12-31 16:39:42
% DurationCPUTime: 1.00s
% Computational Cost: add. (460->166), mult. (776->209), div. (0->0), fcn. (748->6), ass. (0->81)
t162 = Icges(2,4) - Icges(3,5);
t161 = Icges(2,1) + Icges(3,1);
t160 = Icges(3,4) + Icges(2,5);
t159 = Icges(2,2) + Icges(3,3);
t158 = Icges(2,6) - Icges(3,6);
t148 = sin(qJ(1));
t157 = t162 * t148;
t122 = cos(qJ(1));
t156 = t162 * t122;
t155 = -t122 * t159 - t157;
t154 = t148 * t159 - t156;
t153 = t148 * t161 + t156;
t152 = t122 * t161 - t157;
t147 = pkin(2) * t122;
t144 = sin(pkin(6));
t145 = cos(pkin(6));
t81 = -t122 * t145 - t144 * t148;
t146 = Icges(4,4) * t81;
t120 = sin(qJ(4));
t143 = Icges(5,4) * t120;
t121 = cos(qJ(4));
t142 = Icges(5,4) * t121;
t103 = pkin(1) * t148 - qJ(2) * t122;
t140 = t103 * V_base(4) + V_base(3);
t139 = V_base(5) * pkin(4) + V_base(1);
t136 = t148 * pkin(2);
t106 = pkin(1) * t122 + qJ(2) * t148;
t134 = -t106 - t147;
t133 = qJD(2) * t148 + t139;
t132 = -t103 - t136;
t131 = t136 * V_base(4) - qJD(3) + t140;
t130 = -rSges(5,1) * t121 + rSges(5,2) * t120;
t129 = -Icges(5,1) * t121 + t143;
t128 = Icges(5,2) * t120 - t142;
t127 = -Icges(5,5) * t121 + Icges(5,6) * t120;
t114 = V_base(6) + qJD(1);
t126 = -qJD(2) * t122 + t106 * t114 + V_base(2);
t77 = -qJD(4) * t81 + V_base(5);
t82 = t122 * t144 - t145 * t148;
t78 = qJD(4) * t82 + V_base(4);
t125 = t114 * (-Icges(5,5) * t120 - Icges(5,6) * t121) + (-Icges(5,3) * t81 + t127 * t82) * t77 + (Icges(5,3) * t82 + t127 * t81) * t78;
t124 = V_base(4) * qJ(3) + t114 * t147 + t126;
t56 = -Icges(5,6) * t81 + t128 * t82;
t57 = Icges(5,6) * t82 + t128 * t81;
t58 = -Icges(5,5) * t81 + t129 * t82;
t59 = Icges(5,5) * t82 + t129 * t81;
t92 = -Icges(5,2) * t121 - t143;
t97 = -Icges(5,1) * t120 - t142;
t123 = (t120 * t57 - t121 * t59) * t78 + (t120 * t56 - t121 * t58) * t77 + (t120 * t92 - t121 * t97) * t114;
t108 = rSges(2,1) * t122 - rSges(2,2) * t148;
t107 = rSges(3,1) * t122 + rSges(3,3) * t148;
t105 = rSges(2,1) * t148 + rSges(2,2) * t122;
t104 = rSges(3,1) * t148 - rSges(3,3) * t122;
t102 = -rSges(5,1) * t120 - rSges(5,2) * t121;
t86 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t85 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t84 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t79 = Icges(4,4) * t82;
t76 = V_base(5) * rSges(2,3) - t105 * t114 + t139;
t75 = t108 * t114 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t74 = t105 * V_base(4) - t108 * V_base(5) + V_base(3);
t73 = -pkin(3) * t81 + pkin(5) * t82;
t72 = -pkin(3) * t82 - pkin(5) * t81;
t71 = -rSges(4,1) * t81 - rSges(4,2) * t82;
t70 = -rSges(4,1) * t82 + rSges(4,2) * t81;
t69 = -Icges(4,1) * t81 - t79;
t68 = -Icges(4,1) * t82 + t146;
t67 = -Icges(4,2) * t82 - t146;
t66 = Icges(4,2) * t81 - t79;
t63 = V_base(5) * rSges(3,2) + (-t103 - t104) * t114 + t133;
t62 = t114 * t107 + (-rSges(3,2) - pkin(4)) * V_base(4) + t126;
t61 = rSges(5,3) * t82 + t130 * t81;
t60 = -rSges(5,3) * t81 + t130 * t82;
t53 = t104 * V_base(4) + (-t106 - t107) * V_base(5) + t140;
t52 = (-qJ(3) - rSges(4,3)) * V_base(5) + (-t70 + t132) * t114 + t133;
t51 = t114 * t71 + (rSges(4,3) - pkin(4)) * V_base(4) + t124;
t50 = V_base(4) * t70 + (t134 - t71) * V_base(5) + t131;
t49 = -V_base(5) * qJ(3) + t77 * t102 + (-t60 - t72 + t132) * t114 + t133;
t48 = -V_base(4) * pkin(4) - t78 * t102 + (t61 + t73) * t114 + t124;
t47 = t78 * t60 - t77 * t61 + V_base(4) * t72 + (t134 - t73) * V_base(5) + t131;
t1 = m(1) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(2) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(3) * (t53 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(4) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(5) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + t78 * (t123 * t81 + t125 * t82) / 0.2e1 + t77 * (t123 * t82 - t125 * t81) / 0.2e1 + ((-t120 * t59 - t121 * t57) * t78 + (-t120 * t58 - t121 * t56) * t77 + (-t120 * t97 - t121 * t92 + Icges(3,2) + Icges(2,3) + Icges(4,3)) * t114) * t114 / 0.2e1 + ((t122 * t153 + t148 * t155 - t66 * t82 - t68 * t81 + Icges(1,4)) * V_base(5) + (t152 * t122 + t154 * t148 - t82 * t67 - t81 * t69 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t155 * t122 + t153 * t148 + t81 * t66 - t82 * t68 + Icges(1,2)) * V_base(5) + (-t122 * t154 + t148 * t152 + t67 * t81 - t69 * t82 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t114 * (Icges(4,5) * t82 - Icges(4,6) * t81 + t122 * t158 + t148 * t160) + V_base(4) * t114 * (Icges(4,5) * t81 + Icges(4,6) * t82 + t122 * t160 - t148 * t158) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
