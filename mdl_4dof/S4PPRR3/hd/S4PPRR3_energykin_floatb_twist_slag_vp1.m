% Calculate kinetic energy for
% S4PPRR3
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:20
% EndTime: 2019-12-31 16:17:21
% DurationCPUTime: 1.09s
% Computational Cost: add. (440->165), mult. (776->214), div. (0->0), fcn. (748->6), ass. (0->82)
t164 = Icges(2,4) - Icges(3,5);
t163 = Icges(2,1) + Icges(3,1);
t162 = Icges(3,4) + Icges(2,5);
t161 = Icges(2,2) + Icges(3,3);
t160 = Icges(2,6) - Icges(3,6);
t146 = sin(pkin(6));
t159 = t164 * t146;
t120 = cos(pkin(6));
t158 = t164 * t120;
t157 = -t161 * t120 - t159;
t156 = t161 * t146 - t158;
t155 = t163 * t146 + t158;
t154 = t163 * t120 - t159;
t150 = cos(qJ(3));
t149 = sin(qJ(3));
t148 = pkin(2) * t120;
t82 = -t120 * t150 - t146 * t149;
t147 = Icges(4,4) * t82;
t121 = sin(qJ(4));
t145 = Icges(5,4) * t121;
t122 = cos(qJ(4));
t144 = Icges(5,4) * t122;
t142 = V_base(5) * qJ(1) + V_base(1);
t138 = qJD(1) + V_base(3);
t137 = t146 * pkin(2);
t102 = t120 * pkin(1) + qJ(2) * t146;
t136 = -t102 - t148;
t99 = pkin(1) * t146 - t120 * qJ(2);
t134 = V_base(4) * t99 + t138;
t133 = qJD(2) * t146 + t142;
t132 = V_base(4) * t137 + t134;
t131 = -rSges(5,1) * t122 + rSges(5,2) * t121;
t130 = -Icges(5,1) * t122 + t145;
t129 = Icges(5,2) * t121 - t144;
t128 = -Icges(5,5) * t122 + Icges(5,6) * t121;
t127 = -qJD(2) * t120 + V_base(6) * t102 + V_base(2);
t117 = V_base(6) - qJD(3);
t77 = -qJD(4) * t82 + V_base(5);
t83 = t120 * t149 - t146 * t150;
t78 = qJD(4) * t83 + V_base(4);
t126 = (-Icges(5,5) * t121 - Icges(5,6) * t122) * t117 + (-Icges(5,3) * t82 + t128 * t83) * t77 + (Icges(5,3) * t83 + t128 * t82) * t78;
t125 = V_base(4) * pkin(4) + V_base(6) * t148 + t127;
t124 = (-t137 - t99) * V_base(6) + t133;
t106 = -Icges(5,2) * t122 - t145;
t107 = -Icges(5,1) * t121 - t144;
t56 = -Icges(5,6) * t82 + t129 * t83;
t57 = Icges(5,6) * t83 + t129 * t82;
t58 = -Icges(5,5) * t82 + t130 * t83;
t59 = Icges(5,5) * t83 + t130 * t82;
t123 = (t121 * t57 - t122 * t59) * t78 + (t121 * t56 - t122 * t58) * t77 + (t106 * t121 - t107 * t122) * t117;
t108 = -t121 * rSges(5,1) - rSges(5,2) * t122;
t104 = t120 * rSges(2,1) - rSges(2,2) * t146;
t103 = t120 * rSges(3,1) + rSges(3,3) * t146;
t101 = rSges(2,1) * t146 + t120 * rSges(2,2);
t100 = rSges(3,1) * t146 - t120 * rSges(3,3);
t86 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t85 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t84 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t79 = Icges(4,4) * t83;
t76 = V_base(5) * rSges(2,3) - t101 * V_base(6) + t142;
t75 = t104 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t74 = -pkin(3) * t82 + pkin(5) * t83;
t73 = -pkin(3) * t83 - pkin(5) * t82;
t72 = -rSges(4,1) * t82 - rSges(4,2) * t83;
t71 = -rSges(4,1) * t83 + rSges(4,2) * t82;
t70 = t101 * V_base(4) - t104 * V_base(5) + t138;
t69 = -Icges(4,1) * t82 - t79;
t68 = -Icges(4,1) * t83 + t147;
t67 = -Icges(4,2) * t83 - t147;
t66 = Icges(4,2) * t82 - t79;
t63 = V_base(5) * rSges(3,2) + (-t100 - t99) * V_base(6) + t133;
t62 = t103 * V_base(6) + (-rSges(3,2) - qJ(1)) * V_base(4) + t127;
t61 = t83 * rSges(5,3) + t131 * t82;
t60 = -t82 * rSges(5,3) + t131 * t83;
t53 = t100 * V_base(4) + (-t102 - t103) * V_base(5) + t134;
t52 = -t117 * t71 + (-pkin(4) - rSges(4,3)) * V_base(5) + t124;
t51 = t117 * t72 + (rSges(4,3) - qJ(1)) * V_base(4) + t125;
t50 = t71 * V_base(4) + (t136 - t72) * V_base(5) + t132;
t49 = -V_base(5) * pkin(4) + t77 * t108 + (-t60 - t73) * t117 + t124;
t48 = -V_base(4) * qJ(1) - t108 * t78 + (t61 + t74) * t117 + t125;
t47 = t60 * t78 - t61 * t77 + t73 * V_base(4) + (t136 - t74) * V_base(5) + t132;
t1 = m(1) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(2) * (t70 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(3) * (t53 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(4) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(5) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + t78 * (t123 * t82 + t126 * t83) / 0.2e1 + t77 * (t123 * t83 - t126 * t82) / 0.2e1 + ((-t121 * t59 - t122 * t57) * t78 + (-t121 * t58 - t122 * t56) * t77 + (-t122 * t106 - t121 * t107 + Icges(4,3)) * t117) * t117 / 0.2e1 + ((t155 * t120 + t146 * t157 - t66 * t83 - t68 * t82 + Icges(1,4)) * V_base(5) + (t154 * t120 + t156 * t146 - t83 * t67 - t82 * t69 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t157 * t120 + t155 * t146 + t82 * t66 - t83 * t68 + Icges(1,2)) * V_base(5) + (-t120 * t156 + t146 * t154 + t67 * t82 - t69 * t83 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t117 * (Icges(4,5) * t82 + Icges(4,6) * t83) + V_base(5) * t117 * (Icges(4,5) * t83 - Icges(4,6) * t82) + (Icges(1,6) * V_base(5) + Icges(1,5) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * V_base(6) + (-t160 * V_base(4) + t162 * V_base(5)) * t146 + (t160 * V_base(5) + t162 * V_base(4)) * t120) * V_base(6);
T = t1;
