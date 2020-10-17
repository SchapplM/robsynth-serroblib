% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:26:51
% EndTime: 2019-05-05 23:26:55
% DurationCPUTime: 1.23s
% Computational Cost: add. (3048->152), mult. (8176->327), div. (0->0), fcn. (9662->14), ass. (0->100)
t77 = sin(pkin(6));
t76 = sin(pkin(7));
t96 = cos(pkin(6));
t91 = t96 * t76;
t79 = cos(pkin(12));
t95 = cos(pkin(7));
t92 = t79 * t95;
t119 = t77 * t92 + t91;
t75 = sin(pkin(12));
t94 = pkin(1) * t96;
t97 = qJ(2) * t77;
t54 = t75 * t94 + t79 * t97;
t40 = pkin(9) * t119 + t54;
t82 = sin(qJ(3));
t85 = cos(qJ(3));
t111 = t75 * t77;
t65 = t79 * t94;
t43 = t96 * pkin(2) + t65 + (-pkin(9) * t95 - qJ(2)) * t111;
t48 = (-pkin(9) * t75 * t76 - pkin(2) * t79 - pkin(1)) * t77;
t88 = t43 * t95 + t48 * t76;
t24 = -t40 * t82 + t85 * t88;
t42 = t82 * t91 + (t75 * t85 + t82 * t92) * t77;
t108 = t77 * t79;
t50 = t108 * t76 - t95 * t96;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t30 = t42 * t81 + t50 * t84;
t31 = t42 * t84 - t50 * t81;
t74 = sin(pkin(13));
t78 = cos(pkin(13));
t23 = -t30 * t74 + t31 * t78;
t41 = t82 * t111 - t119 * t85;
t80 = sin(qJ(6));
t83 = cos(qJ(6));
t15 = t23 * t80 - t41 * t83;
t118 = -0.2e1 * t15;
t117 = -0.2e1 * t30;
t116 = 0.2e1 * t41;
t115 = -0.2e1 * t42;
t114 = 0.2e1 * t84;
t71 = t77 ^ 2;
t113 = pkin(1) * t71;
t28 = -t43 * t76 + t48 * t95;
t18 = pkin(3) * t41 - pkin(10) * t42 + t28;
t25 = t85 * t40 + t82 * t88;
t20 = -t50 * pkin(10) + t25;
t13 = t18 * t81 + t20 * t84;
t11 = -qJ(5) * t30 + t13;
t12 = t18 * t84 - t20 * t81;
t9 = pkin(4) * t41 - qJ(5) * t31 + t12;
t6 = t11 * t78 + t74 * t9;
t16 = t23 * t83 + t41 * t80;
t112 = t16 * t80;
t110 = t76 * t82;
t109 = t76 * t85;
t22 = t30 * t78 + t31 * t74;
t107 = t80 * t22;
t57 = t74 * t81 - t78 * t84;
t106 = t80 * t57;
t58 = t74 * t84 + t78 * t81;
t105 = t80 * t58;
t68 = pkin(4) * t74 + pkin(11);
t104 = t80 * t68;
t103 = t80 * t83;
t102 = t81 * t41;
t101 = t83 * t58;
t100 = t83 * t68;
t99 = t84 * t41;
t98 = -qJ(5) - pkin(10);
t70 = -pkin(4) * t84 - pkin(3);
t93 = t98 * t81;
t5 = -t11 * t74 + t78 * t9;
t69 = -pkin(4) * t78 - pkin(5);
t89 = -t57 * t68 + t58 * t69;
t87 = -t110 * t81 + t84 * t95;
t19 = t50 * pkin(3) - t24;
t14 = t30 * pkin(4) + t19;
t73 = t83 ^ 2;
t72 = t80 ^ 2;
t61 = t98 * t84;
t56 = t58 ^ 2;
t55 = t110 * t84 + t81 * t95;
t53 = -t75 * t97 + t65;
t52 = t83 * t57;
t46 = -t78 * t61 + t74 * t93;
t44 = -t61 * t74 - t78 * t93;
t38 = pkin(5) * t57 - pkin(11) * t58 + t70;
t36 = t55 * t78 + t74 * t87;
t34 = t55 * t74 - t78 * t87;
t33 = -t109 * t80 + t36 * t83;
t32 = -t109 * t83 - t36 * t80;
t27 = t38 * t80 + t46 * t83;
t26 = t38 * t83 - t46 * t80;
t21 = t83 * t22;
t7 = t22 * pkin(5) - t23 * pkin(11) + t14;
t4 = pkin(11) * t41 + t6;
t3 = -pkin(5) * t41 - t5;
t2 = t4 * t83 + t7 * t80;
t1 = -t4 * t80 + t7 * t83;
t8 = [1, 0, 0, 0.2e1 * t113 * t79 + 0.2e1 * t53 * t96, -0.2e1 * t113 * t75 - 0.2e1 * t54 * t96, 0.2e1 * (-t53 * t75 + t54 * t79) * t77, pkin(1) ^ 2 * t71 + t53 ^ 2 + t54 ^ 2, t42 ^ 2, t41 * t115, t50 * t115, t50 * t116, t50 ^ 2, -0.2e1 * t24 * t50 + 0.2e1 * t28 * t41, 0.2e1 * t25 * t50 + 0.2e1 * t28 * t42, t31 ^ 2, t31 * t117, t31 * t116, t41 * t117, t41 ^ 2, 0.2e1 * t12 * t41 + 0.2e1 * t19 * t30, -0.2e1 * t13 * t41 + 0.2e1 * t19 * t31, -0.2e1 * t22 * t6 - 0.2e1 * t23 * t5, t14 ^ 2 + t5 ^ 2 + t6 ^ 2, t16 ^ 2, t16 * t118, 0.2e1 * t16 * t22, t22 * t118, t22 ^ 2, 0.2e1 * t1 * t22 + 0.2e1 * t15 * t3, 0.2e1 * t16 * t3 - 0.2e1 * t2 * t22; 0, 0, 0, -t108, t111, 0, -t77 * pkin(1), 0, 0, 0, 0, 0, -t109 * t50 + t41 * t95, t110 * t50 + t42 * t95, 0, 0, 0, 0, 0, -t109 * t30 + t41 * t87, -t109 * t31 - t41 * t55, -t22 * t36 + t23 * t34, -t109 * t14 - t34 * t5 + t36 * t6, 0, 0, 0, 0, 0, t15 * t34 + t22 * t32, t16 * t34 - t22 * t33; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76 ^ 2 * t85 ^ 2 + t34 ^ 2 + t36 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t41, -t50, t24, -t25, t31 * t81, -t30 * t81 + t31 * t84, t102, t99, 0, -pkin(3) * t30 - pkin(10) * t102 - t19 * t84, -pkin(3) * t31 - pkin(10) * t99 + t19 * t81, -t22 * t46 + t23 * t44 - t5 * t58 - t57 * t6, t14 * t70 - t44 * t5 + t46 * t6, t16 * t101 (-t15 * t83 - t112) * t58, t101 * t22 + t16 * t57, -t105 * t22 - t15 * t57, t22 * t57, t1 * t57 + t105 * t3 + t15 * t44 + t22 * t26, t101 * t3 + t16 * t44 - t2 * t57 - t22 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, -t110, 0, 0, 0, 0, 0, t84 * t109, -t81 * t109, t34 * t58 - t36 * t57, -t109 * t70 + t34 * t44 + t36 * t46, 0, 0, 0, 0, 0, t105 * t34 + t32 * t57, t101 * t34 - t33 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t81 ^ 2, t81 * t114, 0, 0, 0, pkin(3) * t114, -0.2e1 * pkin(3) * t81, 0.2e1 * t44 * t58 - 0.2e1 * t46 * t57, t44 ^ 2 + t46 ^ 2 + t70 ^ 2, t73 * t56, -0.2e1 * t56 * t103, 0.2e1 * t57 * t101, -0.2e1 * t57 * t105, t57 ^ 2, 0.2e1 * t105 * t44 + 0.2e1 * t26 * t57, 0.2e1 * t101 * t44 - 0.2e1 * t27 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, t41, t12, -t13 (-t22 * t74 - t23 * t78) * pkin(4) (t5 * t78 + t6 * t74) * pkin(4), t112, -t15 * t80 + t16 * t83, t107, t21, 0, -t104 * t22 + t15 * t69 - t3 * t83, -t100 * t22 + t16 * t69 + t3 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t55, 0 (-t34 * t78 + t36 * t74) * pkin(4), 0, 0, 0, 0, 0, -t34 * t83, t34 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t84, 0, -t81 * pkin(10), -t84 * pkin(10) (-t57 * t74 - t58 * t78) * pkin(4) (-t44 * t78 + t46 * t74) * pkin(4), t80 * t101 (-t72 + t73) * t58, t106, t52, 0, -t44 * t83 + t80 * t89, t44 * t80 + t83 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t74 ^ 2 + t78 ^ 2) * pkin(4) ^ 2, t72, 0.2e1 * t103, 0, 0, 0, -0.2e1 * t69 * t83, 0.2e1 * t69 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, 0, 0, 0, t21, -t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, 0, 0, 0, t52, -t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t105, t57, t26, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t83, 0, -t104, -t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
