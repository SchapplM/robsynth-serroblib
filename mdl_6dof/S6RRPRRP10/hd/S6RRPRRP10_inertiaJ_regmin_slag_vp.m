% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t112 = cos(qJ(4));
t69 = sin(pkin(11));
t71 = cos(pkin(11));
t74 = sin(qJ(4));
t54 = t112 * t69 + t74 * t71;
t124 = -0.2e1 * t54;
t70 = sin(pkin(6));
t77 = cos(qJ(2));
t106 = t70 * t77;
t75 = sin(qJ(2));
t107 = t70 * t75;
t72 = cos(pkin(6));
t45 = t71 * t107 + t69 * t72;
t58 = t69 * t107;
t89 = t71 * t72 - t58;
t30 = t112 * t45 + t74 * t89;
t73 = sin(qJ(5));
t76 = cos(qJ(5));
t23 = t76 * t106 + t30 * t73;
t123 = -0.2e1 * t23;
t122 = -0.2e1 * t30;
t63 = -pkin(3) * t71 - pkin(2);
t121 = 0.2e1 * t63;
t120 = -0.2e1 * t73;
t119 = pkin(1) * t75;
t118 = pkin(1) * t77;
t29 = -t112 * t89 + t45 * t74;
t117 = pkin(5) * t29;
t53 = -t112 * t71 + t69 * t74;
t116 = pkin(5) * t53;
t115 = pkin(10) * t53;
t114 = t73 * pkin(10);
t113 = t76 * pkin(10);
t60 = pkin(8) * t107;
t33 = t58 * pkin(3) + t60 + (t63 - t118) * t72;
t13 = t29 * pkin(4) - t30 * pkin(10) + t33;
t93 = pkin(8) * t106;
t41 = t93 + (qJ(3) + t119) * t72;
t42 = (-pkin(2) * t77 - qJ(3) * t75 - pkin(1)) * t70;
t25 = -t41 * t69 + t71 * t42;
t16 = -pkin(3) * t106 - pkin(9) * t45 + t25;
t26 = t71 * t41 + t69 * t42;
t20 = t89 * pkin(9) + t26;
t10 = t112 * t20 + t74 * t16;
t8 = -pkin(10) * t106 + t10;
t4 = t73 * t13 + t76 * t8;
t111 = t23 * t76;
t24 = -t73 * t106 + t30 * t76;
t110 = t24 * t73;
t109 = t24 * t76;
t65 = t70 ^ 2;
t108 = t65 * t77;
t105 = t72 * t75;
t27 = t73 * t29;
t46 = t73 * t53;
t104 = t73 * t54;
t103 = t73 * t76;
t28 = t76 * t29;
t47 = t76 * t53;
t48 = t76 * t54;
t102 = pkin(9) + qJ(3);
t34 = pkin(4) * t53 - pkin(10) * t54 + t63;
t56 = t102 * t69;
t57 = t102 * t71;
t37 = t112 * t57 - t74 * t56;
t18 = t73 * t34 + t76 * t37;
t101 = t69 ^ 2 + t71 ^ 2;
t67 = t73 ^ 2;
t68 = t76 ^ 2;
t100 = t67 + t68;
t99 = qJ(6) * t29;
t98 = qJ(6) * t53;
t97 = t53 * t124;
t96 = 0.2e1 * t106;
t95 = pkin(10) * t27;
t94 = pkin(10) * t28;
t92 = qJ(3) * t106;
t91 = -t76 * t13 + t73 * t8;
t90 = -t76 * t34 + t37 * t73;
t36 = t112 * t56 + t57 * t74;
t88 = -pkin(4) * t54 - t115;
t1 = t99 + t4;
t2 = t91 - t117;
t87 = t1 * t76 + t2 * t73;
t86 = t1 * t73 - t2 * t76;
t9 = t112 * t16 - t74 * t20;
t84 = pkin(5) * t76 + qJ(6) * t73;
t55 = -pkin(4) - t84;
t85 = -t54 * t55 + t115;
t83 = pkin(5) * t73 - qJ(6) * t76;
t14 = t98 + t18;
t15 = t90 - t116;
t82 = t14 * t76 + t15 * t73;
t81 = t14 * t73 - t15 * t76;
t80 = -t25 * t69 + t26 * t71;
t7 = pkin(4) * t106 - t9;
t51 = t54 ^ 2;
t50 = pkin(1) * t105 + t93;
t49 = t72 * t118 - t60;
t44 = t60 + (-pkin(2) - t118) * t72;
t22 = t73 * t23;
t21 = t83 * t54 + t36;
t5 = t23 * pkin(5) - t24 * qJ(6) + t7;
t3 = [1, 0, 0, t65 * t75 ^ 2, 0.2e1 * t75 * t108, 0.2e1 * t70 * t105, t72 * t96, t72 ^ 2, 0.2e1 * pkin(1) * t108 + 0.2e1 * t49 * t72, -0.2e1 * t65 * t119 - 0.2e1 * t50 * t72, -0.2e1 * t25 * t106 - 0.2e1 * t44 * t89, 0.2e1 * t26 * t106 + 0.2e1 * t44 * t45, -0.2e1 * t25 * t45 + 0.2e1 * t26 * t89, t25 ^ 2 + t26 ^ 2 + t44 ^ 2, t30 ^ 2, t29 * t122, t106 * t122, t29 * t96, t65 * t77 ^ 2, -0.2e1 * t9 * t106 + 0.2e1 * t29 * t33, 0.2e1 * t10 * t106 + 0.2e1 * t30 * t33, t24 ^ 2, t24 * t123, 0.2e1 * t24 * t29, t29 * t123, t29 ^ 2, 0.2e1 * t23 * t7 - 0.2e1 * t29 * t91, 0.2e1 * t24 * t7 - 0.2e1 * t29 * t4, -0.2e1 * t2 * t29 + 0.2e1 * t23 * t5, -0.2e1 * t1 * t23 + 0.2e1 * t2 * t24, 0.2e1 * t1 * t29 - 0.2e1 * t24 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t107, t106, t72, t49, -t50, pkin(2) * t89 - t44 * t71 + t69 * t92, -pkin(2) * t45 + t44 * t69 + t71 * t92 (t69 * t45 + t71 * t89) * qJ(3) + t80, -t44 * pkin(2) + t80 * qJ(3), t30 * t54, -t29 * t54 - t30 * t53, -t54 * t106, t53 * t106, 0, t36 * t106 + t29 * t63 + t33 * t53, t37 * t106 + t30 * t63 + t33 * t54, t24 * t48 (-t110 - t111) * t54, t24 * t53 + t29 * t48, -t104 * t29 - t23 * t53, t29 * t53, t104 * t7 + t23 * t36 - t29 * t90 - t53 * t91, -t18 * t29 + t24 * t36 - t4 * t53 + t48 * t7, t104 * t5 - t15 * t29 - t2 * t53 + t21 * t23, -t14 * t23 + t15 * t24 - t86 * t54, t1 * t53 + t14 * t29 - t21 * t24 - t48 * t5, t1 * t14 + t15 * t2 + t21 * t5; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t71, -0.2e1 * pkin(2) * t69, 0.2e1 * t101 * qJ(3), t101 * qJ(3) ^ 2 + pkin(2) ^ 2, t51, t97, 0, 0, 0, t53 * t121, t54 * t121, t68 * t51, -0.2e1 * t51 * t103, 0.2e1 * t53 * t48, t73 * t97, t53 ^ 2, 0.2e1 * t104 * t36 - 0.2e1 * t53 * t90, -0.2e1 * t18 * t53 + 0.2e1 * t36 * t48, 0.2e1 * t104 * t21 - 0.2e1 * t15 * t53, t81 * t124, 0.2e1 * t14 * t53 - 0.2e1 * t21 * t48, t14 ^ 2 + t15 ^ 2 + t21 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t45, 0, t44, 0, 0, 0, 0, 0, t29, t30, 0, 0, 0, 0, 0, t28, -t27, t28, -t22 - t109, t27, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t69, 0, -pkin(2), 0, 0, 0, 0, 0, t53, t54, 0, 0, 0, 0, 0, t47, -t46, t47, -t100 * t54, t46, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, -t106, t9, -t10, t110, -t22 + t109, t27, t28, 0, -pkin(4) * t23 - t7 * t76 - t95, -pkin(4) * t24 + t7 * t73 - t94, t23 * t55 - t5 * t76 - t95 (t110 - t111) * pkin(10) + t87, -t24 * t55 - t5 * t73 + t94, pkin(10) * t87 + t5 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t53, 0, -t36, -t37, t73 * t48 (-t67 + t68) * t54, t46, t47, 0, -t36 * t76 + t73 * t88, t36 * t73 + t76 * t88, -t21 * t76 - t73 * t85, t82, -t21 * t73 + t76 * t85, pkin(10) * t82 + t21 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t67, 0.2e1 * t103, 0, 0, 0, 0.2e1 * pkin(4) * t76, pkin(4) * t120, -0.2e1 * t55 * t76, 0.2e1 * t100 * pkin(10), t55 * t120, pkin(10) ^ 2 * t100 + t55 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, t29, -t91, -t4, -t91 + 0.2e1 * t117, -pkin(5) * t24 - qJ(6) * t23, 0.2e1 * t99 + t4, -pkin(5) * t2 + qJ(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t104, t53, -t90, -t18, -t90 + 0.2e1 * t116, -t84 * t54, 0.2e1 * t98 + t18, -pkin(5) * t15 + qJ(6) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t73, t76, 0, t73, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t76, 0, -t114, -t113, -t114, -t83, t113, -t83 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t24, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t48, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
