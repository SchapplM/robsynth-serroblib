% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:28:45
% EndTime: 2019-05-05 11:28:48
% DurationCPUTime: 1.07s
% Computational Cost: add. (1207->152), mult. (3071->296), div. (0->0), fcn. (3784->14), ass. (0->112)
t111 = cos(qJ(5));
t71 = sin(pkin(7));
t78 = sin(qJ(3));
t107 = t71 * t78;
t73 = cos(pkin(7));
t77 = sin(qJ(4));
t81 = cos(qJ(4));
t47 = t77 * t107 - t81 * t73;
t48 = t81 * t107 + t77 * t73;
t76 = sin(qJ(5));
t30 = t111 * t47 + t76 * t48;
t120 = -0.2e1 * t30;
t67 = -t81 * pkin(4) - pkin(3);
t119 = 0.2e1 * t67;
t118 = 0.2e1 * t81;
t117 = pkin(10) + pkin(11);
t116 = pkin(2) * t78;
t82 = cos(qJ(3));
t115 = pkin(2) * t82;
t106 = t71 * t82;
t91 = pkin(9) * t106;
t42 = t91 + (pkin(10) + t116) * t73;
t43 = (-pkin(3) * t82 - pkin(10) * t78 - pkin(2)) * t71;
t26 = -t77 * t42 + t81 * t43;
t92 = pkin(4) * t106;
t16 = -t48 * pkin(11) + t26 - t92;
t27 = t81 * t42 + t77 * t43;
t19 = -t47 * pkin(11) + t27;
t8 = t111 * t16 - t76 * t19;
t6 = pkin(5) * t106 - t8;
t80 = cos(qJ(6));
t114 = t6 * t80;
t113 = t76 * pkin(4);
t88 = t111 * pkin(4);
t66 = -t88 - pkin(5);
t112 = pkin(5) - t66;
t83 = cos(qJ(2));
t102 = t73 * t83;
t72 = sin(pkin(6));
t74 = cos(pkin(6));
t79 = sin(qJ(2));
t36 = t74 * t107 + (t78 * t102 + t79 * t82) * t72;
t104 = t72 * t83;
t46 = -t71 * t104 + t74 * t73;
t22 = -t36 * t77 + t46 * t81;
t23 = t36 * t81 + t46 * t77;
t12 = -t111 * t22 + t76 * t23;
t110 = t12 * t80;
t31 = t111 * t48 - t76 * t47;
t75 = sin(qJ(6));
t25 = -t75 * t106 + t80 * t31;
t21 = t25 * t75;
t59 = t117 * t81;
t86 = t111 * t77;
t38 = t117 * t86 + t76 * t59;
t109 = t38 * t80;
t68 = t71 ^ 2;
t108 = t68 * t82;
t105 = t72 * t79;
t103 = t73 * t78;
t28 = t75 * t30;
t55 = t76 * t81 + t86;
t101 = t75 * t55;
t65 = pkin(12) + t113;
t100 = t75 * t65;
t99 = t75 * t80;
t98 = t76 * t77;
t29 = t80 * t30;
t97 = t80 * t55;
t96 = t80 * t65;
t54 = -t111 * t81 + t98;
t95 = -0.2e1 * t55 * t54;
t94 = -0.2e1 * t106;
t93 = 0.2e1 * t106;
t90 = t77 * t106;
t89 = t81 * t106;
t87 = t111 * t19;
t85 = -pkin(5) * t55 - pkin(12) * t54;
t84 = -t54 * t65 + t55 * t66;
t61 = pkin(9) * t107;
t41 = t61 + (-pkin(3) - t115) * t73;
t9 = t76 * t16 + t87;
t32 = t47 * pkin(4) + t41;
t70 = t80 ^ 2;
t69 = t75 ^ 2;
t63 = t68 * t82 ^ 2;
t62 = 0.2e1 * t99;
t53 = t55 ^ 2;
t52 = pkin(2) * t103 + t91;
t51 = t73 * t115 - t61;
t50 = t80 * t54;
t49 = t75 * t54;
t45 = t75 * t97;
t39 = t111 * t59 - t117 * t98;
t37 = t38 * t75;
t35 = -t72 * t82 * t102 + t78 * t105 - t74 * t106;
t34 = (-t69 + t70) * t55;
t33 = t54 * pkin(5) - t55 * pkin(12) + t67;
t24 = t80 * t106 + t75 * t31;
t18 = t75 * t33 + t80 * t39;
t17 = t80 * t33 - t75 * t39;
t14 = -t75 * t24 + t25 * t80;
t13 = t111 * t23 + t76 * t22;
t11 = t30 * pkin(5) - t31 * pkin(12) + t32;
t10 = t12 * t75;
t7 = -pkin(12) * t106 + t9;
t5 = t6 * t75;
t4 = t80 * t13 + t35 * t75;
t3 = -t75 * t13 + t35 * t80;
t2 = t75 * t11 + t80 * t7;
t1 = t80 * t11 - t75 * t7;
t15 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t104, -t105, 0, 0, 0, 0, 0, -t46 * t106 - t35 * t73, t46 * t107 - t36 * t73, 0, 0, 0, 0, 0, -t22 * t106 + t35 * t47, t23 * t106 + t35 * t48, 0, 0, 0, 0, 0, t12 * t106 + t35 * t30, t13 * t106 + t35 * t31, 0, 0, 0, 0, 0, t12 * t24 + t3 * t30, t12 * t25 - t4 * t30; 0, 1, 0, 0, t68 * t78 ^ 2, 0.2e1 * t78 * t108, 0.2e1 * t71 * t103, t73 * t93, t73 ^ 2, 0.2e1 * pkin(2) * t108 + 0.2e1 * t51 * t73, -0.2e1 * t68 * t116 - 0.2e1 * t52 * t73, t48 ^ 2, -0.2e1 * t48 * t47, t48 * t94, t47 * t93, t63, -0.2e1 * t26 * t106 + 0.2e1 * t41 * t47, 0.2e1 * t27 * t106 + 0.2e1 * t41 * t48, t31 ^ 2, t31 * t120, t31 * t94, t30 * t93, t63, -0.2e1 * t8 * t106 + 0.2e1 * t32 * t30, 0.2e1 * t9 * t106 + 0.2e1 * t32 * t31, t25 ^ 2, -0.2e1 * t25 * t24, 0.2e1 * t25 * t30, t24 * t120, t30 ^ 2, 0.2e1 * t1 * t30 + 0.2e1 * t6 * t24, -0.2e1 * t2 * t30 + 0.2e1 * t6 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, 0, 0, 0, 0, -t35 * t81, t35 * t77, 0, 0, 0, 0, 0, t35 * t54, t35 * t55, 0, 0, 0, 0, 0, t12 * t101 + t3 * t54, t12 * t97 - t4 * t54; 0, 0, 0, 0, 0, 0, t107, t106, t73, t51, -t52, t48 * t77, -t77 * t47 + t48 * t81, -t90, -t89, 0, -pkin(3) * t47 + pkin(10) * t90 - t41 * t81, -pkin(3) * t48 + pkin(10) * t89 + t41 * t77, t31 * t55, -t55 * t30 - t31 * t54, -t55 * t106, t54 * t106, 0, t38 * t106 + t67 * t30 + t32 * t54, t39 * t106 + t67 * t31 + t32 * t55, t25 * t97 (-t24 * t80 - t21) * t55, t25 * t54 + t30 * t97, -t30 * t101 - t24 * t54, t30 * t54, t1 * t54 + t6 * t101 + t17 * t30 + t38 * t24, -t18 * t30 - t2 * t54 + t38 * t25 + t6 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t77 ^ 2, t77 * t118, 0, 0, 0, pkin(3) * t118, -0.2e1 * pkin(3) * t77, t53, t95, 0, 0, 0, t54 * t119, t55 * t119, t70 * t53, -0.2e1 * t53 * t99, 0.2e1 * t54 * t97, t75 * t95, t54 ^ 2, 0.2e1 * t38 * t101 + 0.2e1 * t17 * t54, -0.2e1 * t18 * t54 + 0.2e1 * t38 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t23, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, -t110, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t47, -t106, t26, -t27, 0, 0, t31, -t30, -t106, -t88 * t106 + t8, -t87 + (-t16 + t92) * t76, t21, t14, t28, t29, 0, -t30 * t100 + t66 * t24 - t114, t66 * t25 - t30 * t96 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t81, 0, -t77 * pkin(10), -t81 * pkin(10), 0, 0, t55, -t54, 0, -t38, -t39, t45, t34, t49, t50, 0, t84 * t75 - t109, t80 * t84 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t88, -0.2e1 * t113, t69, t62, 0, 0, 0, -0.2e1 * t66 * t80, 0.2e1 * t66 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, -t110, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, -t106, t8, -t9, t21, t14, t28, t29, 0, -pkin(5) * t24 - pkin(12) * t28 - t114, -pkin(5) * t25 - pkin(12) * t29 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t54, 0, -t38, -t39, t45, t34, t49, t50, 0, t85 * t75 - t109, t80 * t85 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t88, -t113, t69, t62, 0, 0, 0, t112 * t80, -t112 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t69, t62, 0, 0, 0, 0.2e1 * pkin(5) * t80, -0.2e1 * pkin(5) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, t30, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, -t101, t54, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t80, 0, -t100, -t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t80, 0, -t75 * pkin(12), -t80 * pkin(12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t15;
