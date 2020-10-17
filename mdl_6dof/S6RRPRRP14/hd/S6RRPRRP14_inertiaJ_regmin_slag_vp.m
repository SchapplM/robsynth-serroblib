% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP14_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:24:33
% EndTime: 2019-05-06 19:24:38
% DurationCPUTime: 1.38s
% Computational Cost: add. (1266->162), mult. (2806->314), div. (0->0), fcn. (3018->8), ass. (0->109)
t120 = -2 * pkin(2);
t119 = 2 * pkin(5);
t56 = sin(pkin(6));
t63 = cos(qJ(2));
t102 = t56 * t63;
t57 = cos(pkin(6));
t59 = sin(qJ(4));
t62 = cos(qJ(4));
t28 = t62 * t102 + t57 * t59;
t118 = -0.2e1 * t28;
t117 = 0.2e1 * t56;
t58 = sin(qJ(5));
t116 = -0.2e1 * t58;
t61 = cos(qJ(5));
t115 = 0.2e1 * t61;
t114 = 2 * qJ(3);
t64 = -pkin(2) - pkin(9);
t60 = sin(qJ(2));
t113 = pkin(1) * t60;
t112 = pkin(1) * t63;
t111 = pkin(10) * t59;
t110 = t28 * pkin(5);
t109 = t58 * pkin(10);
t108 = t61 * pkin(10);
t46 = t56 * t60;
t41 = pkin(8) * t46;
t78 = -pkin(2) - t112;
t14 = pkin(3) * t46 + t41 + (-pkin(9) + t78) * t57;
t75 = -qJ(3) * t60 - pkin(1);
t20 = (t64 * t63 + t75) * t56;
t9 = t62 * t14 - t59 * t20;
t7 = -pkin(4) * t46 - t9;
t107 = t7 * t58;
t106 = t7 * t61;
t31 = pkin(8) * t102 + t57 * t113;
t50 = t57 * qJ(3);
t25 = -t50 - t31;
t19 = pkin(3) * t102 - t25;
t29 = -t59 * t102 + t57 * t62;
t13 = t28 * pkin(4) - t29 * pkin(10) + t19;
t10 = t59 * t14 + t62 * t20;
t8 = pkin(10) * t46 + t10;
t4 = t58 * t13 + t61 * t8;
t16 = t29 * t58 - t61 * t46;
t105 = t16 * t61;
t17 = t29 * t61 + t58 * t46;
t104 = t17 * t58;
t51 = t56 ^ 2;
t103 = t51 * t63;
t101 = t57 * t63;
t100 = t58 * t28;
t47 = t58 * t59;
t99 = t58 * t61;
t98 = t58 * t62;
t97 = t58 * t64;
t96 = t59 * t64;
t95 = t61 * t28;
t48 = t61 * t59;
t49 = t61 * t62;
t94 = t61 * t64;
t93 = t62 * t17;
t92 = t62 * t28;
t71 = -t61 * pkin(5) - t58 * qJ(6);
t37 = -pkin(4) + t71;
t91 = t62 * t37;
t90 = t62 * t59;
t89 = t62 * t64;
t36 = t59 * pkin(4) - t62 * pkin(10) + qJ(3);
t23 = t58 * t36 + t59 * t94;
t52 = t58 ^ 2;
t54 = t61 ^ 2;
t88 = t52 + t54;
t53 = t59 ^ 2;
t55 = t62 ^ 2;
t87 = t53 + t55;
t86 = t28 * qJ(6);
t85 = t59 * qJ(6);
t84 = 0.2e1 * t46;
t83 = -0.2e1 * t90;
t82 = pkin(10) * t100;
t81 = pkin(10) * t95;
t80 = t59 * t46;
t79 = t64 * t46;
t77 = -t61 * t13 + t58 * t8;
t76 = t88 * t59;
t74 = -pkin(4) * t62 - t111;
t1 = t86 + t4;
t2 = t77 - t110;
t73 = t1 * t61 + t2 * t58;
t72 = -t91 + t111;
t70 = -pkin(5) * t58 + t61 * qJ(6);
t69 = t104 - t105;
t18 = t85 + t23;
t33 = t61 * t36;
t21 = -t33 + (-pkin(5) + t97) * t59;
t68 = t18 * t61 + t21 * t58;
t67 = -t62 * t16 - t28 * t47;
t66 = t28 * t48 + t93;
t45 = t51 * t60 ^ 2;
t39 = t62 * t46;
t35 = t87 * t61;
t34 = t87 * t58;
t30 = pkin(1) * t101 - t41;
t27 = t78 * t57 + t41;
t26 = (-pkin(2) * t63 + t75) * t56;
t24 = (-t64 - t70) * t62;
t22 = -t58 * t96 + t33;
t5 = t16 * pkin(5) - t17 * qJ(6) + t7;
t3 = [1, 0, 0, t45, 0.2e1 * t60 * t103, t57 * t84, t101 * t117, t57 ^ 2, 0.2e1 * pkin(1) * t103 + 0.2e1 * t30 * t57, -0.2e1 * t51 * t113 - 0.2e1 * t31 * t57 (-t25 * t63 + t27 * t60) * t117, 0.2e1 * t26 * t102 + 0.2e1 * t27 * t57, -0.2e1 * t25 * t57 - 0.2e1 * t26 * t46, t25 ^ 2 + t26 ^ 2 + t27 ^ 2, t29 ^ 2, t29 * t118, t29 * t84, t46 * t118, t45, 0.2e1 * t19 * t28 + 0.2e1 * t9 * t46, -0.2e1 * t10 * t46 + 0.2e1 * t19 * t29, t17 ^ 2, -0.2e1 * t17 * t16, 0.2e1 * t17 * t28, t16 * t118, t28 ^ 2, 0.2e1 * t7 * t16 - 0.2e1 * t28 * t77, 0.2e1 * t7 * t17 - 0.2e1 * t4 * t28, 0.2e1 * t5 * t16 - 0.2e1 * t2 * t28, -0.2e1 * t1 * t16 + 0.2e1 * t2 * t17, 0.2e1 * t1 * t28 - 0.2e1 * t5 * t17, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t46, t102, t57, t30, -t31 (-pkin(2) * t60 + qJ(3) * t63) * t56, t41 + (t120 - t112) * t57, 0.2e1 * t50 + t31, -t27 * pkin(2) - t25 * qJ(3), t29 * t62, -t29 * t59 - t92, t39, -t80, 0, qJ(3) * t28 + t19 * t59 + t62 * t79, qJ(3) * t29 + t19 * t62 - t59 * t79, t61 * t93 (-t104 - t105) * t62, t17 * t59 + t61 * t92, -t16 * t59 - t58 * t92, t28 * t59, t22 * t28 - t77 * t59 + (-t16 * t64 + t107) * t62, -t23 * t28 - t4 * t59 + (-t17 * t64 + t106) * t62, t24 * t16 - t2 * t59 - t21 * t28 + t5 * t98, -t18 * t16 + t21 * t17 + (-t1 * t58 + t2 * t61) * t62, t1 * t59 - t24 * t17 + t18 * t28 - t5 * t49, t1 * t18 + t2 * t21 + t5 * t24; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t120, t114, pkin(2) ^ 2 + qJ(3) ^ 2, t55, t83, 0, 0, 0, t59 * t114, t62 * t114, t54 * t55, -0.2e1 * t55 * t99, t90 * t115, t58 * t83, t53, 0.2e1 * t22 * t59 - 0.2e1 * t55 * t97, -0.2e1 * t23 * t59 - 0.2e1 * t55 * t94, -0.2e1 * t21 * t59 + 0.2e1 * t24 * t98, 0.2e1 * (-t18 * t58 + t21 * t61) * t62, 0.2e1 * t18 * t59 - 0.2e1 * t24 * t49, t18 ^ 2 + t21 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t57, 0, t27, 0, 0, 0, 0, 0, t39, -t80, 0, 0, 0, 0, 0, t67, -t66, t67, t69 * t59, t66, -t5 * t62 + t59 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t35, -t34, 0, t35, -t24 * t62 + t59 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t88 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, t46, t9, -t10, t104, -t58 * t16 + t17 * t61, t100, t95, 0, -pkin(4) * t16 - t106 - t82, -pkin(4) * t17 + t107 - t81, t37 * t16 - t5 * t61 - t82, t69 * pkin(10) + t73, -t37 * t17 - t5 * t58 + t81, pkin(10) * t73 + t5 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t59, 0, t89, -t96, t58 * t49 (-t52 + t54) * t62, t47, t48, 0, t74 * t58 + t61 * t89, -t58 * t89 + t74 * t61, -t24 * t61 - t72 * t58, t68, -t24 * t58 + t72 * t61, pkin(10) * t68 + t24 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t59, 0, 0, 0, 0, 0, t49, -t98, t49, t76, t98, pkin(10) * t76 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t52, 0.2e1 * t99, 0, 0, 0, pkin(4) * t115, pkin(4) * t116, -0.2e1 * t37 * t61, 0.2e1 * t88 * pkin(10), t37 * t116, pkin(10) ^ 2 * t88 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t28, -t77, -t4, -t77 + 0.2e1 * t110, -t17 * pkin(5) - t16 * qJ(6), 0.2e1 * t86 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t98, t59, t22, -t23, t33 + (t119 - t97) * t59, t71 * t62, 0.2e1 * t85 + t23, -t21 * pkin(5) + t18 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t48, -t47, 0, t48, t70 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t61, 0, -t109, -t108, -t109, t70, t108, t70 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t119, 0, 0.2e1 * qJ(6) (pkin(5) ^ 2) + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t17, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t49, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
