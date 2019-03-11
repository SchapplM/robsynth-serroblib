% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t75 = sin(qJ(5));
t69 = t75 ^ 2;
t78 = cos(qJ(5));
t71 = t78 ^ 2;
t58 = t69 + t71;
t111 = cos(qJ(4));
t73 = sin(pkin(10));
t116 = t73 * pkin(2);
t74 = cos(pkin(10));
t115 = t74 * pkin(2);
t63 = pkin(3) + t115;
t76 = sin(qJ(4));
t49 = t111 * t116 + t63 * t76;
t47 = pkin(9) + t49;
t99 = t58 * t47;
t77 = sin(qJ(2));
t79 = cos(qJ(2));
t51 = -t73 * t77 + t74 * t79;
t52 = t73 * t79 + t74 * t77;
t34 = t111 * t52 + t51 * t76;
t127 = -0.2e1 * t34;
t103 = -qJ(3) - pkin(7);
t55 = t103 * t77;
t57 = t103 * t79;
t37 = t55 * t73 - t57 * t74;
t22 = pkin(8) * t51 + t37;
t36 = t55 * t74 + t57 * t73;
t83 = -pkin(8) * t52 + t36;
t13 = -t111 * t83 + t22 * t76;
t126 = t13 ^ 2;
t32 = -t111 * t51 + t52 * t76;
t30 = t32 ^ 2;
t64 = -pkin(2) * t79 - pkin(1);
t44 = -pkin(3) * t51 + t64;
t125 = 0.2e1 * t44;
t124 = 0.2e1 * t52;
t123 = -0.2e1 * t75;
t122 = -0.2e1 * t78;
t121 = 0.2e1 * t79;
t120 = pkin(5) * t32;
t119 = pkin(9) * t32;
t86 = pkin(5) * t75 - qJ(6) * t78;
t7 = t34 * t86 + t13;
t118 = t7 * t75;
t117 = t7 * t78;
t114 = t75 * pkin(9);
t113 = t78 * pkin(9);
t48 = t111 * t63 - t116 * t76;
t46 = -pkin(4) - t48;
t112 = pkin(4) - t46;
t12 = pkin(4) * t32 - pkin(9) * t34 + t44;
t15 = t111 * t22 + t76 * t83;
t6 = t12 * t75 + t15 * t78;
t110 = t13 * t78;
t109 = t32 * t47;
t108 = t34 * t69;
t25 = t75 * t32;
t107 = t75 * t34;
t106 = t75 * t47;
t105 = t75 * t78;
t28 = t78 * t32;
t29 = t78 * t34;
t104 = t78 * t47;
t87 = pkin(5) * t78 + qJ(6) * t75;
t54 = -pkin(4) - t87;
t35 = t54 - t48;
t102 = -t35 - t54;
t101 = t99 * pkin(9);
t100 = t58 * t47 ^ 2;
t98 = t58 * pkin(9) ^ 2;
t97 = t58 * pkin(9);
t70 = t77 ^ 2;
t72 = t79 ^ 2;
t96 = t70 + t72;
t95 = qJ(6) * t32;
t94 = t32 * t107;
t31 = t34 ^ 2;
t93 = t31 * t105;
t92 = -t12 * t78 + t15 * t75;
t91 = -pkin(4) * t34 - t119;
t3 = t95 + t6;
t4 = t92 - t120;
t1 = t3 * t78 + t4 * t75;
t90 = t3 * t75 - t4 * t78;
t89 = t6 * t75 - t78 * t92;
t2 = t6 * t78 + t75 * t92;
t88 = -t34 * t54 + t119;
t85 = -t34 * t35 + t109;
t84 = t34 * t46 - t109;
t61 = -0.2e1 * t105;
t60 = 0.2e1 * t105;
t53 = 0.2e1 * t97;
t27 = t71 * t34;
t26 = t71 * t31;
t24 = t69 * t31;
t23 = t75 * t29;
t20 = 0.2e1 * t99;
t19 = t97 + t99;
t18 = 0.2e1 * t32 * t29;
t17 = -t27 - t108;
t16 = -t27 + t108;
t11 = t13 * t75;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t70, t77 * t121, 0, t72, 0, 0, pkin(1) * t121, -0.2e1 * pkin(1) * t77, 0.2e1 * t96 * pkin(7), pkin(7) ^ 2 * t96 + pkin(1) ^ 2, t52 ^ 2, t51 * t124, 0, t51 ^ 2, 0, 0, -0.2e1 * t64 * t51, t64 * t124, -0.2e1 * t36 * t52 + 0.2e1 * t37 * t51, t36 ^ 2 + t37 ^ 2 + t64 ^ 2, t31, t32 * t127, 0, t30, 0, 0, t32 * t125, t34 * t125, 0.2e1 * t13 * t34 - 0.2e1 * t15 * t32, t15 ^ 2 + t44 ^ 2 + t126, t26, -0.2e1 * t93, t18, t24, -0.2e1 * t94, t30, 0.2e1 * t107 * t13 - 0.2e1 * t32 * t92, 0.2e1 * t13 * t29 - 0.2e1 * t32 * t6, t89 * t127, t6 ^ 2 + t92 ^ 2 + t126, t26, t18, 0.2e1 * t93, t30, 0.2e1 * t94, t24, 0.2e1 * t107 * t7 - 0.2e1 * t32 * t4, t90 * t127, -0.2e1 * t29 * t7 + 0.2e1 * t3 * t32, t3 ^ 2 + t4 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, t79, 0, -t77 * pkin(7), -t79 * pkin(7), 0, 0, 0, 0, t52, 0, t51, 0, t36, -t37 (t51 * t73 - t52 * t74) * pkin(2) (t36 * t74 + t37 * t73) * pkin(2), 0, 0, t34, 0, -t32, 0, -t13, -t15, -t32 * t49 - t34 * t48, -t13 * t48 + t15 * t49, t23, -t16, t25, -t23, t28, 0, t75 * t84 - t110, t78 * t84 + t11, t2, t13 * t46 + t2 * t47, t23, t25, t16, 0, -t28, -t23, -t75 * t85 - t117, t1, t78 * t85 - t118, t1 * t47 + t35 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t115, -0.2e1 * t116, 0 (t73 ^ 2 + t74 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t48, -0.2e1 * t49, 0, t48 ^ 2 + t49 ^ 2, t69, t60, 0, t71, 0, 0, t46 * t122, 0.2e1 * t46 * t75, t20, t46 ^ 2 + t100, t69, 0, t61, 0, 0, t71, t35 * t122, t20, t35 * t123, t35 ^ 2 + t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t52, 0, t64, 0, 0, 0, 0, 0, 0, t32, t34, 0, t44, 0, 0, 0, 0, 0, 0, t28, -t25, t17, t89, 0, 0, 0, 0, 0, 0, t28, t17, t25, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t32, 0, -t13, -t15, 0, 0, t23, -t16, t25, -t23, t28, 0, t75 * t91 - t110, t78 * t91 + t11, t2, -t13 * pkin(4) + pkin(9) * t2, t23, t25, t16, 0, -t28, -t23, -t75 * t88 - t117, t1, t78 * t88 - t118, pkin(9) * t1 + t7 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t48, -t49, 0, 0, t69, t60, 0, t71, 0, 0, t112 * t78, -t112 * t75, t19, -pkin(4) * t46 + t101, t69, 0, t61, 0, 0, t71, t102 * t78, t19, t102 * t75, t35 * t54 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t69, t60, 0, t71, 0, 0, 0.2e1 * pkin(4) * t78, pkin(4) * t123, t53, pkin(4) ^ 2 + t98, t69, 0, t61, 0, 0, t71, t54 * t122, t53, t54 * t123, t54 ^ 2 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t107, t32, -t92, -t6, 0, 0, 0, t29, 0, t32, t107, 0, -t92 + 0.2e1 * t120, -t87 * t34, 0.2e1 * t95 + t6, -pkin(5) * t4 + qJ(6) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, t78, 0, -t106, -t104, 0, 0, 0, t75, 0, 0, -t78, 0, -t106, -t86, t104, -t86 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t75, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, t75, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, t78, 0, -t114, -t113, 0, 0, 0, t75, 0, 0, -t78, 0, -t114, -t86, t113, -t86 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t29, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t5;
