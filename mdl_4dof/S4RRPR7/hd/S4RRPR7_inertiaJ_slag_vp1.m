% Calculate joint inertia matrix for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR7_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR7_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:58
% EndTime: 2019-12-31 17:06:01
% DurationCPUTime: 1.16s
% Computational Cost: add. (1727->195), mult. (2341->302), div. (0->0), fcn. (2442->8), ass. (0->108)
t88 = qJ(2) + pkin(7);
t84 = cos(t88);
t97 = cos(qJ(1));
t136 = t84 * t97;
t157 = Icges(4,6) * t84;
t96 = cos(qJ(2));
t153 = Icges(3,6) * t96;
t83 = sin(t88);
t154 = Icges(4,5) * t83;
t93 = sin(qJ(2));
t155 = Icges(3,5) * t93;
t156 = t153 / 0.2e1 + t154 / 0.2e1 + t155 / 0.2e1;
t152 = Icges(3,5) * t96 + Icges(4,5) * t84 - Icges(3,6) * t93 - Icges(4,6) * t83;
t151 = Icges(3,3) + Icges(4,3);
t94 = sin(qJ(1));
t147 = t151 * t97 - t152 * t94;
t146 = t151 * t94 + t152 * t97;
t137 = t83 * t97;
t92 = sin(qJ(4));
t130 = t97 * t92;
t95 = cos(qJ(4));
t133 = t94 * t95;
t61 = -t84 * t130 + t133;
t129 = t97 * t95;
t134 = t94 * t92;
t62 = t84 * t129 + t134;
t30 = t62 * rSges(5,1) + t61 * rSges(5,2) + rSges(5,3) * t137;
t145 = pkin(3) * t136 + pkin(6) * t137 + t30;
t89 = t94 ^ 2;
t90 = t97 ^ 2;
t144 = -t84 / 0.2e1;
t143 = t94 / 0.2e1;
t142 = pkin(2) * t93;
t141 = pkin(3) * t84;
t140 = rSges(3,1) * t96;
t139 = rSges(3,2) * t93;
t138 = t83 * t94;
t40 = -Icges(5,6) * t84 + (Icges(5,4) * t95 - Icges(5,2) * t92) * t83;
t135 = t92 * t40;
t132 = t97 * rSges(3,3);
t91 = -qJ(3) - pkin(5);
t131 = t97 * t91;
t82 = t96 * pkin(2) + pkin(1);
t79 = t97 * t82;
t87 = t97 * pkin(5);
t128 = t94 * (t131 + t87 + (-pkin(1) + t82) * t94) + t97 * (-t97 * pkin(1) + t79 + (-pkin(5) - t91) * t94);
t126 = t94 * rSges(3,3) + t97 * t140;
t125 = t89 + t90;
t121 = Icges(4,4) * t84;
t120 = Icges(5,5) * t83;
t119 = Icges(5,6) * t83;
t118 = Icges(5,3) * t83;
t59 = -t84 * t134 - t129;
t60 = t84 * t133 - t130;
t23 = Icges(5,5) * t60 + Icges(5,6) * t59 + t94 * t118;
t25 = Icges(5,4) * t60 + Icges(5,2) * t59 + t94 * t119;
t27 = Icges(5,1) * t60 + Icges(5,4) * t59 + t94 * t120;
t10 = -t84 * t23 + (-t25 * t92 + t27 * t95) * t83;
t39 = -Icges(5,3) * t84 + (Icges(5,5) * t95 - Icges(5,6) * t92) * t83;
t41 = -Icges(5,5) * t84 + (Icges(5,1) * t95 - Icges(5,4) * t92) * t83;
t12 = t39 * t138 + t59 * t40 + t60 * t41;
t117 = t10 / 0.2e1 + t12 / 0.2e1;
t24 = Icges(5,5) * t62 + Icges(5,6) * t61 + t97 * t118;
t26 = Icges(5,4) * t62 + Icges(5,2) * t61 + t97 * t119;
t28 = Icges(5,1) * t62 + Icges(5,4) * t61 + t97 * t120;
t11 = -t84 * t24 + (-t26 * t92 + t28 * t95) * t83;
t13 = t39 * t137 + t61 * t40 + t62 * t41;
t116 = t11 / 0.2e1 + t13 / 0.2e1;
t115 = -t83 * rSges(4,1) - t84 * rSges(4,2) - t142;
t114 = -t94 * t91 + t79;
t42 = -t84 * rSges(5,3) + (rSges(5,1) * t95 - rSges(5,2) * t92) * t83;
t113 = -t83 * pkin(3) + t84 * pkin(6) - t142 - t42;
t111 = -t139 + t140;
t110 = rSges(4,1) * t84 - rSges(4,2) * t83;
t109 = -t60 * rSges(5,1) - t59 * rSges(5,2);
t101 = -Icges(4,2) * t83 + t121;
t98 = rSges(4,1) * t136 - rSges(4,2) * t137 + t94 * rSges(4,3);
t74 = t97 * rSges(2,1) - t94 * rSges(2,2);
t73 = -t94 * rSges(2,1) - t97 * rSges(2,2);
t72 = t93 * rSges(3,1) + t96 * rSges(3,2);
t44 = t115 * t97;
t43 = t115 * t94;
t38 = t94 * pkin(5) + (pkin(1) - t139) * t97 + t126;
t37 = t132 + t87 + (-pkin(1) - t111) * t94;
t34 = t83 * t95 * t41;
t33 = t114 + t98;
t32 = (rSges(4,3) - t91) * t97 + (-t110 - t82) * t94;
t31 = t97 * (-t97 * t139 + t126) + (t111 * t94 - t132) * t94;
t29 = rSges(5,3) * t138 - t109;
t22 = t113 * t97;
t21 = t113 * t94;
t20 = t114 + t145;
t19 = -t131 + (-t141 - t82 + (-rSges(5,3) - pkin(6)) * t83) * t94 + t109;
t18 = -t42 * t137 - t84 * t30;
t17 = t42 * t138 + t84 * t29;
t16 = -t83 * t135 - t84 * t39 + t34;
t15 = t97 * t98 + (-t97 * rSges(4,3) + t110 * t94) * t94 + t128;
t14 = (t29 * t97 - t30 * t94) * t83;
t9 = t24 * t137 + t61 * t26 + t62 * t28;
t8 = t23 * t137 + t61 * t25 + t62 * t27;
t7 = t24 * t138 + t59 * t26 + t60 * t28;
t6 = t23 * t138 + t59 * t25 + t60 * t27;
t5 = t145 * t97 + (t29 + (pkin(6) * t83 + t141) * t94) * t94 + t128;
t4 = -t8 * t97 + t9 * t94;
t3 = -t6 * t97 + t7 * t94;
t2 = -t13 * t84 + (t8 * t94 + t9 * t97) * t83;
t1 = -t12 * t84 + (t6 * t94 + t7 * t97) * t83;
t35 = [Icges(3,2) * t96 ^ 2 + Icges(2,3) + t34 + (Icges(4,4) * t83 + Icges(4,2) * t84 - t39) * t84 + (Icges(4,1) * t83 + t121 - t135) * t83 + m(2) * (t73 ^ 2 + t74 ^ 2) + m(3) * (t37 ^ 2 + t38 ^ 2) + m(4) * (t32 ^ 2 + t33 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2) + (Icges(3,1) * t93 + 0.2e1 * Icges(3,4) * t96) * t93; (t101 * t94 * t144 - t117 + (-Icges(4,6) * t144 + t156) * t97) * t97 + (t101 * t136 / 0.2e1 + t116 + (t157 / 0.2e1 + t156) * t94) * t94 + m(4) * (t44 * t32 + t43 * t33) + m(5) * (t22 * t19 + t21 * t20) + m(3) * (-t37 * t97 - t38 * t94) * t72 + (t153 + t154 + t155 + t157) * (t89 / 0.2e1 + t90 / 0.2e1); m(5) * (t21 ^ 2 + t22 ^ 2 + t5 ^ 2) + m(4) * (t15 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(3) * (t125 * t72 ^ 2 + t31 ^ 2) + (t146 * t89 + t4) * t94 + (-t3 + t147 * t90 + (t146 * t97 + t147 * t94) * t94) * t97; m(4) * (t94 * t32 - t97 * t33) + m(5) * (t94 * t19 - t97 * t20); m(5) * (-t97 * t21 + t94 * t22) + m(4) * (-t97 * t43 + t94 * t44); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t125; -t16 * t84 + m(5) * (t17 * t19 + t18 * t20) + (t116 * t97 + t117 * t94) * t83; m(5) * (t14 * t5 + t17 * t22 + t18 * t21) - t97 * t1 / 0.2e1 + (-t10 * t97 + t11 * t94) * t144 + t2 * t143 + (t97 * t4 / 0.2e1 + t3 * t143) * t83; m(5) * (t17 * t94 - t18 * t97); m(5) * (t14 ^ 2 + t17 ^ 2 + t18 ^ 2) + t84 ^ 2 * t16 + (t97 * t2 + t94 * t1 - t84 * (t10 * t94 + t11 * t97)) * t83;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t35(1), t35(2), t35(4), t35(7); t35(2), t35(3), t35(5), t35(8); t35(4), t35(5), t35(6), t35(9); t35(7), t35(8), t35(9), t35(10);];
Mq = res;
