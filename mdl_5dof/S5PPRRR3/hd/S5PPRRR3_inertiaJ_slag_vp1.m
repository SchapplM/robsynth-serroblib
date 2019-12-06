% Calculate joint inertia matrix for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:23
% EndTime: 2019-12-05 15:16:27
% DurationCPUTime: 1.29s
% Computational Cost: add. (5784->289), mult. (10866->447), div. (0->0), fcn. (13220->10), ass. (0->147)
t130 = cos(pkin(9));
t131 = cos(pkin(8));
t135 = cos(qJ(3));
t142 = t131 * t135;
t129 = sin(pkin(8));
t133 = sin(qJ(3));
t146 = t129 * t133;
t115 = t130 * t146 + t142;
t162 = t115 / 0.2e1;
t143 = t131 * t133;
t145 = t129 * t135;
t117 = t130 * t143 - t145;
t161 = t117 / 0.2e1;
t160 = -t130 / 0.2e1;
t134 = cos(qJ(4));
t159 = t134 * pkin(4);
t116 = t130 * t145 - t143;
t127 = qJ(4) + qJ(5);
t125 = sin(t127);
t126 = cos(t127);
t128 = sin(pkin(9));
t152 = t128 * t129;
t100 = -t116 * t125 + t126 * t152;
t101 = t116 * t126 + t125 * t152;
t64 = t101 * rSges(6,1) + t100 * rSges(6,2) + t115 * rSges(6,3);
t132 = sin(qJ(4));
t150 = t128 * t132;
t139 = t129 * t150;
t66 = pkin(4) * t139 + pkin(7) * t115 + t159 * t116;
t157 = t64 + t66;
t118 = t130 * t142 + t146;
t151 = t128 * t131;
t102 = -t118 * t125 + t126 * t151;
t103 = t118 * t126 + t125 * t151;
t65 = t103 * rSges(6,1) + t102 * rSges(6,2) + t117 * rSges(6,3);
t138 = t131 * t150;
t67 = pkin(4) * t138 + pkin(7) * t117 + t159 * t118;
t156 = -t65 - t67;
t148 = t128 * t134;
t106 = -t118 * t132 + t131 * t148;
t107 = t118 * t134 + t138;
t75 = t107 * rSges(5,1) + t106 * rSges(5,2) + t117 * rSges(5,3);
t99 = t118 * pkin(3) + t117 * pkin(6);
t155 = -t75 - t99;
t147 = t128 * t135;
t113 = -t125 * t147 - t130 * t126;
t114 = -t130 * t125 + t126 * t147;
t149 = t128 * t133;
t82 = t114 * rSges(6,1) + t113 * rSges(6,2) + rSges(6,3) * t149;
t144 = t130 * t132;
t95 = -pkin(4) * t144 + (pkin(7) * t133 + t159 * t135) * t128;
t154 = t82 + t95;
t121 = (pkin(3) * t135 + pkin(6) * t133) * t128;
t98 = t116 * pkin(3) + t115 * pkin(6);
t153 = t121 * t152 + t130 * t98;
t58 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t115;
t60 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t115;
t62 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t115;
t36 = t113 * t60 + t114 * t62 + t58 * t149;
t59 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t117;
t61 = Icges(6,4) * t103 + Icges(6,2) * t102 + Icges(6,6) * t117;
t63 = Icges(6,1) * t103 + Icges(6,4) * t102 + Icges(6,5) * t117;
t37 = t113 * t61 + t114 * t63 + t59 * t149;
t79 = Icges(6,5) * t114 + Icges(6,6) * t113 + Icges(6,3) * t149;
t80 = Icges(6,4) * t114 + Icges(6,2) * t113 + Icges(6,6) * t149;
t81 = Icges(6,1) * t114 + Icges(6,4) * t113 + Icges(6,5) * t149;
t46 = t113 * t80 + t114 * t81 + t79 * t149;
t15 = t36 * t115 + t37 * t117 + t46 * t149;
t26 = t100 * t60 + t101 * t62 + t115 * t58;
t27 = t100 * t61 + t101 * t63 + t115 * t59;
t41 = t100 * t80 + t101 * t81 + t115 * t79;
t5 = t26 * t115 + t27 * t117 + t41 * t149;
t28 = t102 * t60 + t103 * t62 + t117 * t58;
t29 = t102 * t61 + t103 * t63 + t117 * t59;
t42 = t102 * t80 + t103 * t81 + t117 * t79;
t6 = t28 * t115 + t29 * t117 + t42 * t149;
t141 = t115 * t5 + t117 * t6 + t15 * t149;
t140 = -t99 + t156;
t10 = -t42 * t130 + (t129 * t28 + t131 * t29) * t128;
t19 = -t46 * t130 + (t129 * t36 + t131 * t37) * t128;
t9 = -t41 * t130 + (t129 * t26 + t131 * t27) * t128;
t137 = t15 * t160 + t19 * t149 / 0.2e1 + t5 * t152 / 0.2e1 + t6 * t151 / 0.2e1 + t9 * t162 + t10 * t161;
t120 = t134 * t147 - t144;
t119 = -t130 * t134 - t132 * t147;
t112 = -t130 * rSges(4,3) + (rSges(4,1) * t135 - rSges(4,2) * t133) * t128;
t111 = -Icges(4,5) * t130 + (Icges(4,1) * t135 - Icges(4,4) * t133) * t128;
t110 = -Icges(4,6) * t130 + (Icges(4,4) * t135 - Icges(4,2) * t133) * t128;
t109 = -Icges(4,3) * t130 + (Icges(4,5) * t135 - Icges(4,6) * t133) * t128;
t105 = t116 * t134 + t139;
t104 = -t116 * t132 + t129 * t148;
t96 = t98 * t151;
t94 = t120 * rSges(5,1) + t119 * rSges(5,2) + rSges(5,3) * t149;
t93 = Icges(5,1) * t120 + Icges(5,4) * t119 + Icges(5,5) * t149;
t92 = Icges(5,4) * t120 + Icges(5,2) * t119 + Icges(5,6) * t149;
t91 = Icges(5,5) * t120 + Icges(5,6) * t119 + Icges(5,3) * t149;
t90 = t118 * rSges(4,1) - t117 * rSges(4,2) + rSges(4,3) * t151;
t89 = t116 * rSges(4,1) - t115 * rSges(4,2) + rSges(4,3) * t152;
t88 = Icges(4,1) * t118 - Icges(4,4) * t117 + Icges(4,5) * t151;
t87 = Icges(4,1) * t116 - Icges(4,4) * t115 + Icges(4,5) * t152;
t86 = Icges(4,4) * t118 - Icges(4,2) * t117 + Icges(4,6) * t151;
t85 = Icges(4,4) * t116 - Icges(4,2) * t115 + Icges(4,6) * t152;
t84 = Icges(4,5) * t118 - Icges(4,6) * t117 + Icges(4,3) * t151;
t83 = Icges(4,5) * t116 - Icges(4,6) * t115 + Icges(4,3) * t152;
t78 = t115 * t82;
t77 = -t112 * t151 - t130 * t90;
t76 = t112 * t152 + t130 * t89;
t74 = t105 * rSges(5,1) + t104 * rSges(5,2) + t115 * rSges(5,3);
t73 = Icges(5,1) * t107 + Icges(5,4) * t106 + Icges(5,5) * t117;
t72 = Icges(5,1) * t105 + Icges(5,4) * t104 + Icges(5,5) * t115;
t71 = Icges(5,4) * t107 + Icges(5,2) * t106 + Icges(5,6) * t117;
t70 = Icges(5,4) * t105 + Icges(5,2) * t104 + Icges(5,6) * t115;
t69 = Icges(5,5) * t107 + Icges(5,6) * t106 + Icges(5,3) * t117;
t68 = Icges(5,5) * t105 + Icges(5,6) * t104 + Icges(5,3) * t115;
t57 = t65 * t149;
t56 = (-t129 * t90 + t131 * t89) * t128;
t55 = t117 * t64;
t54 = -t117 * t94 + t75 * t149;
t53 = t115 * t94 - t74 * t149;
t52 = -t117 * t82 + t57;
t51 = -t64 * t149 + t78;
t50 = t119 * t92 + t120 * t93 + t91 * t149;
t49 = -t115 * t75 + t117 * t74;
t48 = t155 * t130 + (-t121 - t94) * t151;
t47 = t130 * t74 + t94 * t152 + t153;
t45 = t106 * t92 + t107 * t93 + t117 * t91;
t44 = t104 * t92 + t105 * t93 + t115 * t91;
t43 = -t115 * t65 + t55;
t40 = t96 + (t155 * t129 + t131 * t74) * t128;
t39 = t119 * t71 + t120 * t73 + t69 * t149;
t38 = t119 * t70 + t120 * t72 + t68 * t149;
t35 = -t154 * t117 + t67 * t149 + t57;
t34 = t115 * t95 - t157 * t149 + t78;
t33 = t106 * t71 + t107 * t73 + t117 * t69;
t32 = t106 * t70 + t107 * t72 + t117 * t68;
t31 = t104 * t71 + t105 * t73 + t115 * t69;
t30 = t104 * t70 + t105 * t72 + t115 * t68;
t25 = t140 * t130 + (-t121 - t154) * t151;
t24 = t157 * t130 + t154 * t152 + t153;
t23 = t156 * t115 + t117 * t66 + t55;
t22 = t96 + (t140 * t129 + t157 * t131) * t128;
t21 = -t50 * t130 + (t129 * t38 + t131 * t39) * t128;
t20 = t38 * t115 + t39 * t117 + t50 * t149;
t17 = -t45 * t130 + (t129 * t32 + t131 * t33) * t128;
t16 = -t44 * t130 + (t129 * t30 + t131 * t31) * t128;
t12 = t32 * t115 + t33 * t117 + t45 * t149;
t11 = t30 * t115 + t31 * t117 + t44 * t149;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t129 ^ 2 + t131 ^ 2); m(4) * t56 + m(5) * t40 + m(6) * t22; m(4) * (t76 * t129 - t77 * t131) + m(5) * (t47 * t129 - t48 * t131) + m(6) * (t24 * t129 - t25 * t131); m(6) * (t22 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t40 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(4) * (t56 ^ 2 + t76 ^ 2 + t77 ^ 2) + (t9 + t16 + (-t115 * t85 + t116 * t87 + t83 * t152) * t152) * t152 + (t10 + t17 + (-t117 * t86 + t118 * t88 + t84 * t151) * t151 + (-t115 * t86 + t116 * t88 - t117 * t85 + t118 * t87 + t83 * t151 + t84 * t152) * t152) * t151 + ((-t109 * t152 + t115 * t110 - t116 * t111) * t152 + (-t109 * t151 + t117 * t110 - t118 * t111) * t151 - t19 - t21 - ((-t133 * t86 + t135 * t88) * t131 + (-t133 * t85 + t135 * t87) * t129) * t128 ^ 2 + (-(t110 * t133 - t111 * t135 - t129 * t83 - t131 * t84) * t128 - t130 * t109) * t130) * t130; m(5) * t49 + m(6) * t23; m(5) * (t53 * t129 - t54 * t131) + m(6) * (t34 * t129 - t35 * t131); t20 * t160 + t16 * t162 + t17 * t161 + m(6) * (t23 * t22 + t34 * t24 + t35 * t25) + m(5) * (t49 * t40 + t53 * t47 + t54 * t48) + (t131 * t12 / 0.2e1 + t129 * t11 / 0.2e1 + t133 * t21 / 0.2e1) * t128 + t137; t20 * t149 + t115 * t11 + t117 * t12 + m(6) * (t23 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(5) * (t49 ^ 2 + t53 ^ 2 + t54 ^ 2) + t141; m(6) * t43; m(6) * (t51 * t129 - t52 * t131); m(6) * (t43 * t22 + t51 * t24 + t52 * t25) + t137; m(6) * (t43 * t23 + t51 * t34 + t52 * t35) + t141; m(6) * (t43 ^ 2 + t51 ^ 2 + t52 ^ 2) + t141;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
