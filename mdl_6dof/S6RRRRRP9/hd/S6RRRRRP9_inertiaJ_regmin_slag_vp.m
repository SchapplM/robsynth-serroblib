% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t76 = sin(qJ(5));
t78 = sin(qJ(3));
t108 = cos(qJ(5));
t80 = cos(qJ(4));
t83 = t108 * t80;
t77 = sin(qJ(4));
t98 = t77 * t78;
t44 = -t76 * t98 + t78 * t83;
t122 = -0.2e1 * t44;
t74 = sin(pkin(6));
t79 = sin(qJ(2));
t102 = t74 * t79;
t75 = cos(pkin(6));
t81 = cos(qJ(3));
t46 = t78 * t102 - t75 * t81;
t121 = -0.2e1 * t46;
t120 = 0.2e1 * t46;
t47 = t81 * t102 + t75 * t78;
t119 = -0.2e1 * t47;
t66 = -t80 * pkin(4) - pkin(3);
t118 = 0.2e1 * t66;
t117 = -0.2e1 * t78;
t116 = 0.2e1 * t81;
t115 = pkin(10) + pkin(11);
t114 = pkin(1) * t79;
t82 = cos(qJ(2));
t113 = pkin(1) * t82;
t112 = pkin(3) * t80;
t111 = pkin(9) * t77;
t110 = t46 * pkin(4);
t109 = t76 * pkin(4);
t101 = t74 * t82;
t88 = pkin(8) * t101;
t40 = t88 + (pkin(9) + t114) * t75;
t41 = (-pkin(2) * t82 - pkin(9) * t79 - pkin(1)) * t74;
t22 = -t78 * t40 + t81 * t41;
t20 = pkin(3) * t101 - t22;
t107 = t20 * t77;
t106 = t20 * t80;
t31 = -t77 * t101 + t47 * t80;
t105 = t31 * t77;
t104 = t46 * t81;
t69 = t74 ^ 2;
t103 = t69 * t82;
t100 = t75 * t79;
t99 = t77 * t46;
t97 = t77 * t80;
t96 = t77 * t81;
t95 = t78 * t46;
t94 = t80 * t46;
t93 = t80 * t78;
t92 = t80 * t81;
t67 = t78 * pkin(9);
t54 = pkin(4) * t98 + t67;
t91 = 0.2e1 * t101;
t90 = t78 * t116;
t89 = pkin(9) * t92;
t87 = t78 * t101;
t86 = t81 * t101;
t68 = t108 * pkin(4);
t60 = pkin(8) * t102;
t39 = t60 + (-pkin(2) - t113) * t75;
t19 = t46 * pkin(3) - t47 * pkin(10) + t39;
t23 = t81 * t40 + t78 * t41;
t21 = -pkin(10) * t101 + t23;
t10 = t77 * t19 + t80 * t21;
t30 = t80 * t101 + t47 * t77;
t8 = -t30 * pkin(11) + t10;
t85 = t108 * t8;
t9 = t80 * t19 - t77 * t21;
t6 = -t31 * pkin(11) + t110 + t9;
t3 = t108 * t6 - t76 * t8;
t55 = -t81 * pkin(3) - t78 * pkin(10) - pkin(2);
t32 = t89 + (-pkin(11) * t78 + t55) * t77;
t84 = t108 * t32;
t50 = t80 * t55;
t28 = -pkin(11) * t93 + t50 + (-pkin(4) - t111) * t81;
t14 = t108 * t28 - t76 * t32;
t56 = t115 * t77;
t57 = t115 * t80;
t33 = -t108 * t56 - t76 * t57;
t4 = t76 * t6 + t85;
t15 = t76 * t28 + t84;
t34 = t108 * t57 - t76 * t56;
t52 = t108 * t77 + t76 * t80;
t13 = t30 * pkin(4) + t20;
t73 = t81 ^ 2;
t72 = t80 ^ 2;
t71 = t78 ^ 2;
t70 = t77 ^ 2;
t65 = t68 + pkin(5);
t51 = t76 * t77 - t83;
t49 = pkin(1) * t100 + t88;
t48 = t75 * t113 - t60;
t45 = t46 ^ 2;
t43 = t52 * t78;
t38 = t51 * pkin(5) + t66;
t37 = t77 * t55 + t89;
t36 = -pkin(9) * t96 + t50;
t29 = t43 * pkin(5) + t54;
t25 = -t51 * qJ(6) + t34;
t24 = -t52 * qJ(6) + t33;
t17 = t108 * t31 - t76 * t30;
t16 = t108 * t30 + t76 * t31;
t12 = -t43 * qJ(6) + t15;
t11 = -t81 * pkin(5) - t44 * qJ(6) + t14;
t7 = t16 * pkin(5) + t13;
t2 = -t16 * qJ(6) + t4;
t1 = t46 * pkin(5) - t17 * qJ(6) + t3;
t5 = [1, 0, 0, t69 * t79 ^ 2, 0.2e1 * t79 * t103, 0.2e1 * t74 * t100, t75 * t91, t75 ^ 2, 0.2e1 * pkin(1) * t103 + 0.2e1 * t48 * t75, -0.2e1 * t69 * t114 - 0.2e1 * t49 * t75, t47 ^ 2, t46 * t119, t101 * t119, t46 * t91, t69 * t82 ^ 2, -0.2e1 * t22 * t101 + 0.2e1 * t39 * t46, 0.2e1 * t23 * t101 + 0.2e1 * t39 * t47, t31 ^ 2, -0.2e1 * t31 * t30, t31 * t120, t30 * t121, t45, 0.2e1 * t20 * t30 + 0.2e1 * t9 * t46, -0.2e1 * t10 * t46 + 0.2e1 * t20 * t31, t17 ^ 2, -0.2e1 * t17 * t16, t17 * t120, t16 * t121, t45, 0.2e1 * t13 * t16 + 0.2e1 * t3 * t46, 0.2e1 * t13 * t17 - 0.2e1 * t4 * t46, -0.2e1 * t1 * t17 - 0.2e1 * t2 * t16, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, t102, t101, t75, t48, -t49, t47 * t78, t47 * t81 - t95, -t87, -t86, 0, -pkin(2) * t46 + pkin(9) * t87 - t39 * t81, -pkin(2) * t47 + pkin(9) * t86 + t39 * t78, t31 * t93 (-t30 * t80 - t105) * t78, -t31 * t81 + t46 * t93, t30 * t81 - t77 * t95, -t104, t36 * t46 - t9 * t81 + (pkin(9) * t30 + t107) * t78, t10 * t81 - t37 * t46 + (pkin(9) * t31 + t106) * t78, t17 * t44, -t44 * t16 - t17 * t43, -t17 * t81 + t44 * t46, t16 * t81 - t43 * t46, -t104, t13 * t43 + t14 * t46 + t54 * t16 - t3 * t81, t13 * t44 - t15 * t46 + t54 * t17 + t4 * t81, -t1 * t44 - t11 * t17 - t12 * t16 - t2 * t43, t1 * t11 + t2 * t12 + t7 * t29; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t71, t90, 0, 0, 0, pkin(2) * t116, pkin(2) * t117, t72 * t71, -0.2e1 * t71 * t97, t92 * t117, t77 * t90, t73, 0.2e1 * t111 * t71 - 0.2e1 * t36 * t81, 0.2e1 * t71 * pkin(9) * t80 + 0.2e1 * t37 * t81, t44 ^ 2, t43 * t122, t81 * t122, t43 * t116, t73, -0.2e1 * t14 * t81 + 0.2e1 * t54 * t43, 0.2e1 * t15 * t81 + 0.2e1 * t54 * t44, -0.2e1 * t11 * t44 - 0.2e1 * t12 * t43, t11 ^ 2 + t12 ^ 2 + t29 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t46, -t101, t22, -t23, t105, -t77 * t30 + t31 * t80, t99, t94, 0, -pkin(3) * t30 - pkin(10) * t99 - t106, -pkin(3) * t31 - pkin(10) * t94 + t107, t17 * t52, -t52 * t16 - t17 * t51, t52 * t46, -t51 * t46, 0, t13 * t51 + t66 * t16 + t33 * t46, t13 * t52 + t66 * t17 - t34 * t46, -t1 * t52 - t25 * t16 - t24 * t17 - t2 * t51, t1 * t24 + t2 * t25 + t7 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t81, 0, -t67, -t81 * pkin(9), t77 * t93 (-t70 + t72) * t78, -t96, -t92, 0, -pkin(9) * t93 + (-pkin(3) * t78 + pkin(10) * t81) * t77, pkin(10) * t92 + (t111 - t112) * t78, t44 * t52, -t43 * t52 - t44 * t51, -t52 * t81, t51 * t81, 0, -t33 * t81 + t43 * t66 + t51 * t54, t34 * t81 + t44 * t66 + t52 * t54, -t11 * t52 - t12 * t51 - t24 * t44 - t25 * t43, t11 * t24 + t12 * t25 + t29 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t70, 0.2e1 * t97, 0, 0, 0, 0.2e1 * t112, -0.2e1 * pkin(3) * t77, t52 ^ 2, -0.2e1 * t52 * t51, 0, 0, 0, t51 * t118, t52 * t118, -0.2e1 * t24 * t52 - 0.2e1 * t25 * t51, t24 ^ 2 + t25 ^ 2 + t38 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, t46, t9, -t10, 0, 0, t17, -t16, t46, t46 * t68 + t3, -t85 + (-t6 - t110) * t76, -t109 * t16 - t17 * t65, t1 * t65 + t109 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, -t98, -t81, t36, -t37, 0, 0, t44, -t43, -t81, -t68 * t81 + t14, -t84 + (pkin(4) * t81 - t28) * t76, -t109 * t43 - t44 * t65, t109 * t12 + t11 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t80, 0, -t77 * pkin(10), -t80 * pkin(10), 0, 0, t52, -t51, 0, t33, -t34, -t109 * t51 - t52 * t65, t109 * t25 + t24 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t68, -0.2e1 * t109, 0, pkin(4) ^ 2 * t76 ^ 2 + t65 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t46, t3, -t4, -pkin(5) * t17, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43, -t81, t14, -t15, -pkin(5) * t44, t11 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t51, 0, t33, -t34, -pkin(5) * t52, t24 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t68, -t109, 0, t65 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t5;
