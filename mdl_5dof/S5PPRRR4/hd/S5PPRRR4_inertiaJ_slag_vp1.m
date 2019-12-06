% Calculate joint inertia matrix for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:39
% EndTime: 2019-12-05 15:18:48
% DurationCPUTime: 2.56s
% Computational Cost: add. (17824->323), mult. (49512->492), div. (0->0), fcn. (65422->14), ass. (0->152)
t131 = sin(pkin(5));
t133 = cos(pkin(5));
t136 = sin(qJ(3));
t152 = cos(pkin(11));
t153 = cos(pkin(6));
t143 = t153 * t152;
t150 = sin(pkin(11));
t151 = sin(pkin(6));
t158 = cos(qJ(3));
t116 = t133 * t151 * t136 + (t136 * t143 + t150 * t158) * t131;
t147 = t131 * t151;
t123 = t133 * t153 - t152 * t147;
t135 = sin(qJ(4));
t157 = cos(qJ(4));
t111 = t116 * t135 - t123 * t157;
t130 = sin(pkin(10));
t132 = cos(pkin(10));
t145 = t133 * t150;
t124 = t130 * t152 + t132 * t145;
t146 = t133 * t152;
t141 = t130 * t150 - t132 * t146;
t139 = t141 * t153;
t144 = t158 * t151;
t142 = t131 * t144;
t107 = t124 * t136 + t132 * t142 + t158 * t139;
t134 = sin(qJ(5));
t137 = cos(qJ(5));
t108 = t124 * t158 + (-t132 * t147 - t139) * t136;
t148 = t131 * t153;
t117 = -t132 * t148 + t141 * t151;
t97 = t108 * t157 + t117 * t135;
t74 = t107 * t137 - t97 * t134;
t75 = t107 * t134 + t97 * t137;
t96 = t108 * t135 - t117 * t157;
t49 = Icges(6,5) * t75 + Icges(6,6) * t74 + Icges(6,3) * t96;
t51 = Icges(6,4) * t75 + Icges(6,2) * t74 + Icges(6,6) * t96;
t53 = Icges(6,1) * t75 + Icges(6,4) * t74 + Icges(6,5) * t96;
t16 = t96 * t49 + t74 * t51 + t75 * t53;
t125 = -t130 * t145 + t132 * t152;
t140 = t130 * t146 + t132 * t150;
t138 = t140 * t153;
t109 = t125 * t136 - t130 * t142 + t158 * t138;
t110 = t125 * t158 + (t130 * t147 - t138) * t136;
t118 = t130 * t148 + t140 * t151;
t99 = t110 * t157 + t118 * t135;
t76 = t109 * t137 - t99 * t134;
t77 = t109 * t134 + t99 * t137;
t98 = t110 * t135 - t118 * t157;
t50 = Icges(6,5) * t77 + Icges(6,6) * t76 + Icges(6,3) * t98;
t52 = Icges(6,4) * t77 + Icges(6,2) * t76 + Icges(6,6) * t98;
t54 = Icges(6,1) * t77 + Icges(6,4) * t76 + Icges(6,5) * t98;
t17 = t96 * t50 + t74 * t52 + t75 * t54;
t112 = t116 * t157 + t123 * t135;
t115 = -t133 * t144 + (t136 * t150 - t143 * t158) * t131;
t100 = -t112 * t134 + t115 * t137;
t101 = t112 * t137 + t115 * t134;
t68 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t111;
t69 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t111;
t70 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t111;
t26 = t96 * t68 + t74 * t69 + t75 * t70;
t1 = t26 * t111 + t16 * t96 + t17 * t98;
t164 = t1 / 0.2e1;
t18 = t98 * t49 + t76 * t51 + t77 * t53;
t19 = t98 * t50 + t76 * t52 + t77 * t54;
t27 = t98 * t68 + t76 * t69 + t77 * t70;
t2 = t27 * t111 + t18 * t96 + t19 * t98;
t163 = t2 / 0.2e1;
t21 = t100 * t51 + t101 * t53 + t111 * t49;
t22 = t100 * t52 + t101 * t54 + t111 * t50;
t35 = t100 * t69 + t101 * t70 + t111 * t68;
t7 = t35 * t111 + t21 * t96 + t22 * t98;
t162 = t7 / 0.2e1;
t161 = t96 / 0.2e1;
t160 = t98 / 0.2e1;
t159 = t111 / 0.2e1;
t55 = t75 * rSges(6,1) + t74 * rSges(6,2) + t96 * rSges(6,3);
t156 = t97 * pkin(4) + t96 * pkin(9) + t55;
t56 = t77 * rSges(6,1) + t76 * rSges(6,2) + t98 * rSges(6,3);
t155 = t99 * pkin(4) + t98 * pkin(9) + t56;
t71 = t101 * rSges(6,1) + t100 * rSges(6,2) + t111 * rSges(6,3);
t154 = t112 * pkin(4) + t111 * pkin(9) + t71;
t149 = m(3) + m(4) + m(5) + m(6);
t106 = t116 * pkin(3) + t115 * pkin(8);
t105 = t116 * rSges(4,1) - t115 * rSges(4,2) + t123 * rSges(4,3);
t104 = Icges(4,1) * t116 - Icges(4,4) * t115 + Icges(4,5) * t123;
t103 = Icges(4,4) * t116 - Icges(4,2) * t115 + Icges(4,6) * t123;
t102 = Icges(4,5) * t116 - Icges(4,6) * t115 + Icges(4,3) * t123;
t95 = t117 * t106;
t93 = t110 * pkin(3) + t109 * pkin(8);
t92 = t108 * pkin(3) + t107 * pkin(8);
t91 = t123 * t93;
t90 = t118 * t92;
t89 = t112 * rSges(5,1) - t111 * rSges(5,2) + t115 * rSges(5,3);
t88 = Icges(5,1) * t112 - Icges(5,4) * t111 + Icges(5,5) * t115;
t87 = Icges(5,4) * t112 - Icges(5,2) * t111 + Icges(5,6) * t115;
t86 = Icges(5,5) * t112 - Icges(5,6) * t111 + Icges(5,3) * t115;
t85 = t110 * rSges(4,1) - t109 * rSges(4,2) + t118 * rSges(4,3);
t84 = t108 * rSges(4,1) - t107 * rSges(4,2) + t117 * rSges(4,3);
t83 = Icges(4,1) * t110 - Icges(4,4) * t109 + Icges(4,5) * t118;
t82 = Icges(4,1) * t108 - Icges(4,4) * t107 + Icges(4,5) * t117;
t81 = Icges(4,4) * t110 - Icges(4,2) * t109 + Icges(4,6) * t118;
t80 = Icges(4,4) * t108 - Icges(4,2) * t107 + Icges(4,6) * t117;
t79 = Icges(4,5) * t110 - Icges(4,6) * t109 + Icges(4,3) * t118;
t78 = Icges(4,5) * t108 - Icges(4,6) * t107 + Icges(4,3) * t117;
t67 = t99 * rSges(5,1) - t98 * rSges(5,2) + t109 * rSges(5,3);
t66 = t97 * rSges(5,1) - t96 * rSges(5,2) + t107 * rSges(5,3);
t65 = Icges(5,1) * t99 - Icges(5,4) * t98 + Icges(5,5) * t109;
t64 = Icges(5,1) * t97 - Icges(5,4) * t96 + Icges(5,5) * t107;
t63 = Icges(5,4) * t99 - Icges(5,2) * t98 + Icges(5,6) * t109;
t62 = Icges(5,4) * t97 - Icges(5,2) * t96 + Icges(5,6) * t107;
t61 = Icges(5,5) * t99 - Icges(5,6) * t98 + Icges(5,3) * t109;
t60 = Icges(5,5) * t97 - Icges(5,6) * t96 + Icges(5,3) * t107;
t59 = -t118 * t105 + t123 * t85;
t58 = t117 * t105 - t123 * t84;
t57 = -t117 * t85 + t118 * t84;
t48 = -t109 * t89 + t115 * t67;
t47 = t107 * t89 - t115 * t66;
t46 = -t111 * t87 + t112 * t88 + t115 * t86;
t45 = -t107 * t67 + t109 * t66;
t44 = t123 * t67 + t91 + (-t106 - t89) * t118;
t43 = t117 * t89 + t95 + (-t66 - t92) * t123;
t42 = t109 * t86 - t98 * t87 + t99 * t88;
t41 = t107 * t86 - t96 * t87 + t97 * t88;
t40 = t111 * t56 - t98 * t71;
t39 = -t111 * t55 + t96 * t71;
t38 = t118 * t66 + t90 + (-t67 - t93) * t117;
t37 = -t111 * t63 + t112 * t65 + t115 * t61;
t36 = -t111 * t62 + t112 * t64 + t115 * t60;
t34 = t109 * t61 - t98 * t63 + t99 * t65;
t33 = t109 * t60 - t98 * t62 + t99 * t64;
t32 = t107 * t61 - t96 * t63 + t97 * t65;
t31 = t107 * t60 - t96 * t62 + t97 * t64;
t30 = t98 * t55 - t96 * t56;
t29 = -t154 * t109 + t155 * t115;
t28 = t154 * t107 - t156 * t115;
t25 = t91 + t155 * t123 + (-t106 - t154) * t118;
t24 = t95 + t154 * t117 + (-t92 - t156) * t123;
t23 = -t155 * t107 + t156 * t109;
t20 = t90 + t156 * t118 + (-t93 - t155) * t117;
t15 = t36 * t117 + t37 * t118 + t46 * t123;
t14 = t36 * t107 + t37 * t109 + t46 * t115;
t13 = t33 * t117 + t34 * t118 + t42 * t123;
t12 = t31 * t117 + t32 * t118 + t41 * t123;
t11 = t33 * t107 + t34 * t109 + t42 * t115;
t10 = t31 * t107 + t32 * t109 + t41 * t115;
t9 = t21 * t117 + t22 * t118 + t35 * t123;
t8 = t21 * t107 + t22 * t109 + t35 * t115;
t6 = t18 * t117 + t19 * t118 + t27 * t123;
t5 = t16 * t117 + t17 * t118 + t26 * t123;
t4 = t18 * t107 + t19 * t109 + t27 * t115;
t3 = t16 * t107 + t17 * t109 + t26 * t115;
t72 = [m(2) + t149; t149 * t133; 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (t133 ^ 2 + (t130 ^ 2 + t132 ^ 2) * t131 ^ 2); m(4) * t57 + m(5) * t38 + m(6) * t20; m(4) * (t57 * t133 + (t130 * t58 - t132 * t59) * t131) + m(5) * (t38 * t133 + (t130 * t43 - t132 * t44) * t131) + m(6) * (t20 * t133 + (t130 * t24 - t132 * t25) * t131); (t9 + t15 + (t123 * t102 - t115 * t103 + t116 * t104) * t123) * t123 + (t6 + t13 + (-t109 * t81 + t110 * t83 + t118 * t79) * t118 + (t118 * t102 - t109 * t103 + t110 * t104 - t115 * t81 + t116 * t83 + t123 * t79) * t123) * t118 + (t5 + t12 + (-t107 * t80 + t108 * t82 + t117 * t78) * t117 + (t117 * t102 - t107 * t103 + t108 * t104 - t115 * t80 + t116 * t82 + t123 * t78) * t123 + (-t107 * t81 + t108 * t83 - t109 * t80 + t110 * t82 + t117 * t79 + t118 * t78) * t118) * t117 + m(6) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t38 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(4) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2); m(5) * t45 + m(6) * t23; m(5) * (t45 * t133 + (t130 * t47 - t132 * t48) * t131) + m(6) * (t23 * t133 + (t130 * t28 - t132 * t29) * t131); (t8 / 0.2e1 + t14 / 0.2e1) * t123 + (t4 / 0.2e1 + t11 / 0.2e1) * t118 + (t3 / 0.2e1 + t10 / 0.2e1) * t117 + (t9 / 0.2e1 + t15 / 0.2e1) * t115 + (t6 / 0.2e1 + t13 / 0.2e1) * t109 + (t5 / 0.2e1 + t12 / 0.2e1) * t107 + m(6) * (t23 * t20 + t28 * t24 + t29 * t25) + m(5) * (t45 * t38 + t47 * t43 + t48 * t44); (t8 + t14) * t115 + (t4 + t11) * t109 + (t3 + t10) * t107 + m(6) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(5) * (t45 ^ 2 + t47 ^ 2 + t48 ^ 2); m(6) * t30; m(6) * (t30 * t133 + (t130 * t39 - t132 * t40) * t131); m(6) * (t30 * t20 + t39 * t24 + t40 * t25) + t6 * t160 + t117 * t164 + t123 * t162 + t9 * t159 + t5 * t161 + t118 * t163; t109 * t163 + t107 * t164 + m(6) * (t30 * t23 + t39 * t28 + t40 * t29) + t8 * t159 + t115 * t162 + t4 * t160 + t3 * t161; m(6) * (t30 ^ 2 + t39 ^ 2 + t40 ^ 2) + t98 * t2 + t96 * t1 + t111 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t72(1), t72(2), t72(4), t72(7), t72(11); t72(2), t72(3), t72(5), t72(8), t72(12); t72(4), t72(5), t72(6), t72(9), t72(13); t72(7), t72(8), t72(9), t72(10), t72(14); t72(11), t72(12), t72(13), t72(14), t72(15);];
Mq = res;
