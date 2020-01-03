% Calculate kinetic energy for
% S4RPRR5
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR5_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR5_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:30
% EndTime: 2019-12-31 16:51:31
% DurationCPUTime: 1.18s
% Computational Cost: add. (472->165), mult. (776->216), div. (0->0), fcn. (748->6), ass. (0->82)
t164 = Icges(2,4) - Icges(3,5);
t163 = Icges(2,1) + Icges(3,1);
t162 = Icges(3,4) + Icges(2,5);
t161 = Icges(2,2) + Icges(3,3);
t160 = Icges(2,6) - Icges(3,6);
t148 = sin(qJ(1));
t159 = t164 * t148;
t123 = cos(qJ(1));
t158 = t164 * t123;
t157 = -t161 * t123 - t159;
t156 = t161 * t148 - t158;
t153 = t163 * t148 + t158;
t152 = t163 * t123 - t159;
t149 = cos(qJ(3));
t147 = sin(qJ(3));
t146 = pkin(2) * t123;
t82 = -t123 * t149 - t148 * t147;
t145 = Icges(4,4) * t82;
t121 = sin(qJ(4));
t144 = Icges(5,4) * t121;
t122 = cos(qJ(4));
t143 = Icges(5,4) * t122;
t103 = t148 * pkin(1) - t123 * qJ(2);
t141 = V_base(4) * t103 + V_base(3);
t140 = V_base(5) * pkin(4) + V_base(1);
t137 = t148 * pkin(2);
t115 = V_base(6) + qJD(1);
t106 = t123 * pkin(1) + t148 * qJ(2);
t135 = -t106 - t146;
t134 = V_base(4) * t137 + t141;
t133 = qJD(2) * t148 + t140;
t132 = -rSges(5,1) * t122 + rSges(5,2) * t121;
t131 = -Icges(5,1) * t122 + t144;
t130 = Icges(5,2) * t121 - t143;
t129 = -Icges(5,5) * t122 + Icges(5,6) * t121;
t128 = -qJD(2) * t123 + t115 * t106 + V_base(2);
t111 = -qJD(3) + t115;
t77 = -qJD(4) * t82 + V_base(5);
t83 = t123 * t147 - t148 * t149;
t78 = qJD(4) * t83 + V_base(4);
t127 = (-Icges(5,5) * t121 - Icges(5,6) * t122) * t111 + (-Icges(5,3) * t82 + t129 * t83) * t77 + (Icges(5,3) * t83 + t129 * t82) * t78;
t126 = V_base(4) * pkin(5) + t115 * t146 + t128;
t125 = (-t137 - t103) * t115 + t133;
t56 = -Icges(5,6) * t82 + t130 * t83;
t57 = Icges(5,6) * t83 + t130 * t82;
t58 = -Icges(5,5) * t82 + t131 * t83;
t59 = Icges(5,5) * t83 + t131 * t82;
t92 = -Icges(5,2) * t122 - t144;
t97 = -Icges(5,1) * t121 - t143;
t124 = (t121 * t57 - t122 * t59) * t78 + (t121 * t56 - t122 * t58) * t77 + (t121 * t92 - t122 * t97) * t111;
t108 = t123 * rSges(2,1) - t148 * rSges(2,2);
t107 = t123 * rSges(3,1) + t148 * rSges(3,3);
t105 = t148 * rSges(2,1) + t123 * rSges(2,2);
t104 = t148 * rSges(3,1) - t123 * rSges(3,3);
t102 = -rSges(5,1) * t121 - rSges(5,2) * t122;
t86 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t85 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t84 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t80 = Icges(4,4) * t83;
t76 = V_base(5) * rSges(2,3) - t105 * t115 + t140;
t75 = t108 * t115 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t74 = t105 * V_base(4) - t108 * V_base(5) + V_base(3);
t73 = -pkin(3) * t82 + pkin(6) * t83;
t72 = -pkin(3) * t83 - pkin(6) * t82;
t71 = -rSges(4,1) * t82 - rSges(4,2) * t83;
t70 = -rSges(4,1) * t83 + rSges(4,2) * t82;
t69 = -Icges(4,1) * t82 - t80;
t68 = -Icges(4,1) * t83 + t145;
t67 = -Icges(4,2) * t83 - t145;
t66 = Icges(4,2) * t82 - t80;
t63 = V_base(5) * rSges(3,2) + (-t103 - t104) * t115 + t133;
t62 = t115 * t107 + (-rSges(3,2) - pkin(4)) * V_base(4) + t128;
t61 = rSges(5,3) * t83 + t132 * t82;
t60 = -rSges(5,3) * t82 + t132 * t83;
t53 = t104 * V_base(4) + (-t106 - t107) * V_base(5) + t141;
t52 = -t111 * t70 + (-pkin(5) - rSges(4,3)) * V_base(5) + t125;
t51 = t111 * t71 + (rSges(4,3) - pkin(4)) * V_base(4) + t126;
t50 = V_base(4) * t70 + (t135 - t71) * V_base(5) + t134;
t49 = -V_base(5) * pkin(5) + t77 * t102 + (-t60 - t72) * t111 + t125;
t48 = -V_base(4) * pkin(4) - t78 * t102 + (t61 + t73) * t111 + t126;
t47 = t78 * t60 - t77 * t61 + V_base(4) * t72 + (t135 - t73) * V_base(5) + t134;
t1 = m(1) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(2) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(3) * (t53 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(4) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(5) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + t78 * (t124 * t82 + t127 * t83) / 0.2e1 + t77 * (t124 * t83 - t127 * t82) / 0.2e1 + ((-t121 * t59 - t122 * t57) * t78 + (-t121 * t58 - t122 * t56) * t77 + (-t121 * t97 - t122 * t92 + Icges(4,3)) * t111) * t111 / 0.2e1 + ((t153 * t123 + t148 * t157 - t66 * t83 - t68 * t82 + Icges(1,4)) * V_base(5) + (t123 * t152 + t148 * t156 - t83 * t67 - t82 * t69 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t123 * t157 + t153 * t148 + t82 * t66 - t83 * t68 + Icges(1,2)) * V_base(5) + (-t123 * t156 + t148 * t152 + t67 * t82 - t69 * t83 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t111 * (Icges(4,5) * t82 + Icges(4,6) * t83) + V_base(5) * t111 * (Icges(4,5) * t83 - Icges(4,6) * t82) + ((Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * t115 + (-t160 * V_base(4) + t162 * V_base(5)) * t148 + (t160 * V_base(5) + t162 * V_base(4)) * t123) * t115 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
