% Calculate kinetic energy for
% S4RPRR9
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
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR9_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR9_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:06
% EndTime: 2019-12-31 16:56:07
% DurationCPUTime: 1.53s
% Computational Cost: add. (478->201), mult. (870->287), div. (0->0), fcn. (758->6), ass. (0->99)
t181 = Icges(2,4) + Icges(3,6);
t180 = Icges(2,1) + Icges(3,2);
t179 = -Icges(3,4) + Icges(2,5);
t178 = Icges(3,5) - Icges(2,6);
t177 = Icges(2,2) + Icges(3,3);
t133 = cos(qJ(1));
t176 = t181 * t133;
t130 = sin(qJ(1));
t175 = t181 * t130;
t174 = t177 * t133 + t175;
t173 = -t177 * t130 + t176;
t172 = t180 * t130 + t176;
t171 = t180 * t133 - t175;
t129 = sin(qJ(3));
t132 = cos(qJ(3));
t162 = Icges(4,4) * t132;
t100 = -Icges(4,2) * t129 + t162;
t163 = Icges(4,4) * t129;
t105 = Icges(4,1) * t132 - t163;
t118 = qJD(3) * t130 + V_base(5);
t119 = qJD(3) * t133 + V_base(4);
t121 = V_base(6) + qJD(1);
t140 = Icges(4,2) * t132 + t163;
t72 = Icges(4,6) * t133 + t130 * t140;
t73 = Icges(4,6) * t130 - t133 * t140;
t141 = Icges(4,1) * t129 + t162;
t75 = Icges(4,5) * t133 + t130 * t141;
t76 = Icges(4,5) * t130 - t133 * t141;
t168 = (t129 * t75 + t132 * t72) * t119 + (t129 * t76 + t132 * t73) * t118 + (t100 * t132 + t105 * t129) * t121;
t166 = pkin(5) * t130;
t165 = pkin(5) * t133;
t128 = sin(qJ(4));
t159 = t128 * t133;
t158 = t130 * t128;
t131 = cos(qJ(4));
t157 = t130 * t131;
t156 = t130 * t132;
t155 = t131 * t133;
t154 = t132 * t133;
t153 = qJD(4) * t132;
t109 = t130 * pkin(1) - qJ(2) * t133;
t152 = V_base(4) * t109 + V_base(3);
t151 = V_base(5) * pkin(4) + V_base(1);
t148 = -t109 - t166;
t147 = qJD(2) * t130 + t151;
t146 = V_base(5) * pkin(2) + t147;
t145 = pkin(3) * t129 - pkin(6) * t132;
t144 = rSges(4,1) * t129 + rSges(4,2) * t132;
t139 = Icges(4,5) * t129 + Icges(4,6) * t132;
t113 = pkin(1) * t133 + t130 * qJ(2);
t137 = -qJD(2) * t133 + t121 * t113 + V_base(2);
t136 = (Icges(4,3) * t130 - t133 * t139) * t118 + (Icges(4,3) * t133 + t130 * t139) * t119 + t121 * (Icges(4,5) * t132 - Icges(4,6) * t129);
t135 = V_base(4) * t166 + (-t113 - t165) * V_base(5) + t152;
t134 = t121 * t165 + (-pkin(2) - pkin(4)) * V_base(4) + t137;
t116 = pkin(3) * t132 + pkin(6) * t129;
t115 = rSges(2,1) * t133 - t130 * rSges(2,2);
t114 = -rSges(3,2) * t133 + t130 * rSges(3,3);
t112 = rSges(4,1) * t132 - rSges(4,2) * t129;
t111 = t130 * rSges(2,1) + rSges(2,2) * t133;
t110 = -t130 * rSges(3,2) - rSges(3,3) * t133;
t108 = qJD(4) * t129 + t121;
t92 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t91 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t90 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t88 = t145 * t133;
t87 = t145 * t130;
t85 = -t129 * t155 + t158;
t84 = t129 * t159 + t157;
t83 = t129 * t157 + t159;
t82 = -t129 * t158 + t155;
t81 = -t130 * t153 + t119;
t80 = t133 * t153 + t118;
t79 = t130 * rSges(4,3) - t133 * t144;
t78 = rSges(5,3) * t129 + (rSges(5,1) * t131 - rSges(5,2) * t128) * t132;
t77 = rSges(4,3) * t133 + t130 * t144;
t74 = Icges(5,5) * t129 + (Icges(5,1) * t131 - Icges(5,4) * t128) * t132;
t71 = Icges(5,6) * t129 + (Icges(5,4) * t131 - Icges(5,2) * t128) * t132;
t68 = Icges(5,3) * t129 + (Icges(5,5) * t131 - Icges(5,6) * t128) * t132;
t67 = V_base(5) * rSges(2,3) - t111 * t121 + t151;
t66 = t115 * t121 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t65 = t111 * V_base(4) - t115 * V_base(5) + V_base(3);
t64 = t85 * rSges(5,1) + t84 * rSges(5,2) + rSges(5,3) * t154;
t63 = rSges(5,1) * t83 + rSges(5,2) * t82 - rSges(5,3) * t156;
t62 = Icges(5,1) * t85 + Icges(5,4) * t84 + Icges(5,5) * t154;
t61 = Icges(5,1) * t83 + Icges(5,4) * t82 - Icges(5,5) * t156;
t60 = Icges(5,4) * t85 + Icges(5,2) * t84 + Icges(5,6) * t154;
t59 = Icges(5,4) * t83 + Icges(5,2) * t82 - Icges(5,6) * t156;
t58 = Icges(5,5) * t85 + Icges(5,6) * t84 + Icges(5,3) * t154;
t57 = Icges(5,5) * t83 + Icges(5,6) * t82 - Icges(5,3) * t156;
t56 = V_base(5) * rSges(3,1) + (-t109 - t110) * t121 + t147;
t55 = t121 * t114 + (-rSges(3,1) - pkin(4)) * V_base(4) + t137;
t54 = t110 * V_base(4) + (-t113 - t114) * V_base(5) + t152;
t53 = t112 * t118 + (t148 - t79) * t121 + t146;
t52 = -t119 * t112 + t121 * t77 + t134;
t51 = -t118 * t77 + t119 * t79 + t135;
t50 = -t108 * t64 + t116 * t118 + t78 * t80 + (t148 + t88) * t121 + t146;
t49 = t108 * t63 - t119 * t116 + t121 * t87 - t81 * t78 + t134;
t48 = -t118 * t87 - t119 * t88 - t80 * t63 + t81 * t64 + t135;
t1 = m(1) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(2) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(3) * (t54 ^ 2 + t55 ^ 2 + t56 ^ 2) / 0.2e1 + m(4) * (t51 ^ 2 + t52 ^ 2 + t53 ^ 2) / 0.2e1 + t119 * (t168 * t130 + t136 * t133) / 0.2e1 + t118 * (t136 * t130 - t168 * t133) / 0.2e1 + m(5) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + t81 * ((-t57 * t156 + t82 * t59 + t83 * t61) * t81 + (-t156 * t58 + t60 * t82 + t62 * t83) * t80 + (-t156 * t68 + t71 * t82 + t74 * t83) * t108) / 0.2e1 + t80 * ((t154 * t57 + t84 * t59 + t85 * t61) * t81 + (t58 * t154 + t84 * t60 + t85 * t62) * t80 + (t154 * t68 + t84 * t71 + t85 * t74) * t108) / 0.2e1 + t108 * ((t68 * t108 + t57 * t81 + t58 * t80) * t129 + ((-t128 * t59 + t131 * t61) * t81 + (-t128 * t60 + t131 * t62) * t80 + (-t128 * t71 + t131 * t74) * t108) * t132) / 0.2e1 + ((-t129 * t72 + t132 * t75) * t119 + (-t129 * t73 + t132 * t76) * t118 + (-t129 * t100 + t132 * t105 + Icges(3,1) + Icges(2,3)) * t121) * t121 / 0.2e1 + ((-t130 * t174 + t172 * t133 + Icges(1,4)) * V_base(5) + (-t173 * t130 + t171 * t133 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t172 * t130 + t174 * t133 + Icges(1,2)) * V_base(5) + (t130 * t171 + t133 * t173 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t121 * (t179 * t130 - t178 * t133) + V_base(4) * t121 * (t178 * t130 + t179 * t133) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
