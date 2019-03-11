% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t70 = sin(pkin(9));
t72 = cos(pkin(9));
t74 = sin(qJ(2));
t75 = cos(qJ(2));
t54 = t70 * t75 + t72 * t74;
t71 = cos(pkin(10));
t73 = sin(qJ(5));
t108 = cos(qJ(5));
t69 = sin(pkin(10));
t89 = t108 * t69;
t55 = t73 * t71 + t89;
t23 = t55 * t54;
t102 = t55 * t23;
t97 = t73 * t69;
t121 = t108 * t71 - t97;
t25 = t121 * t54;
t15 = t25 * t121;
t125 = -t15 + t102;
t110 = t72 * pkin(2);
t63 = -pkin(3) - t110;
t57 = -t71 * pkin(4) + t63;
t82 = pkin(5) * t121 + t55 * qJ(6);
t21 = t57 - t82;
t124 = -0.2e1 * t21;
t123 = -0.2e1 * t54;
t122 = -t15 - t102;
t117 = t121 ^ 2;
t49 = t55 ^ 2;
t120 = t49 + t117;
t119 = t23 ^ 2;
t96 = -qJ(3) - pkin(7);
t58 = t96 * t75;
t86 = t96 * t74;
t37 = -t70 * t58 - t72 * t86;
t118 = t37 ^ 2;
t50 = t70 * t74 - t72 * t75;
t45 = t50 ^ 2;
t116 = 0.2e1 * t50;
t115 = 0.2e1 * t57;
t64 = -t75 * pkin(2) - pkin(1);
t114 = 0.2e1 * t64;
t113 = 0.2e1 * t75;
t100 = t69 * t54;
t32 = t50 * pkin(3) - t54 * qJ(4) + t64;
t39 = -t72 * t58 + t70 * t86;
t13 = t69 * t32 + t71 * t39;
t10 = -pkin(8) * t100 + t13;
t12 = t71 * t32 - t69 * t39;
t98 = t71 * t54;
t8 = t50 * pkin(4) - pkin(8) * t98 + t12;
t4 = t108 * t10 + t73 * t8;
t112 = t50 * pkin(5);
t111 = t70 * pkin(2);
t60 = qJ(4) + t111;
t109 = pkin(8) + t60;
t107 = t23 * t121;
t106 = t25 * t23;
t44 = t109 * t71;
t29 = t109 * t89 + t73 * t44;
t105 = t29 * t50;
t31 = t108 * t44 - t109 * t97;
t104 = t31 * t50;
t103 = t50 * t23;
t34 = t50 * t121;
t35 = t55 * t50;
t36 = t55 * t121;
t101 = t69 * t50;
t99 = t69 * t71;
t65 = t69 ^ 2;
t66 = t71 ^ 2;
t95 = t65 + t66;
t67 = t74 ^ 2;
t68 = t75 ^ 2;
t94 = t67 + t68;
t93 = t50 * qJ(6);
t92 = t50 * t123;
t91 = t29 ^ 2 + t31 ^ 2;
t90 = t69 * t98;
t87 = t73 * t10 - t108 * t8;
t85 = -t31 * t23 + t29 * t25;
t84 = -t121 * t29 + t31 * t55;
t20 = pkin(4) * t100 + t37;
t81 = t12 * t71 + t13 * t69;
t80 = -t12 * t69 + t13 * t71;
t79 = -t50 * t60 + t54 * t63;
t78 = 0.2e1 * t121 * t31 + 0.2e1 * t29 * t55;
t48 = t54 ^ 2;
t43 = t71 * t50;
t22 = t25 ^ 2;
t16 = t25 * t55;
t14 = t25 * t116;
t5 = t23 * pkin(5) - t25 * qJ(6) + t20;
t2 = t87 - t112;
t1 = t93 + t4;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t67, t74 * t113, 0, t68, 0, 0, pkin(1) * t113, -0.2e1 * pkin(1) * t74, 0.2e1 * t94 * pkin(7), t94 * pkin(7) ^ 2 + pkin(1) ^ 2, t48, t92, 0, t45, 0, 0, t50 * t114, t54 * t114, 0.2e1 * t37 * t54 - 0.2e1 * t39 * t50, t39 ^ 2 + t64 ^ 2 + t118, t66 * t48, -0.2e1 * t48 * t99, t98 * t116, t65 * t48, t69 * t92, t45, 0.2e1 * t37 * t100 + 0.2e1 * t12 * t50, -0.2e1 * t13 * t50 + 0.2e1 * t37 * t98, t81 * t123, t12 ^ 2 + t13 ^ 2 + t118, t22, -0.2e1 * t106, t14, t119, -0.2e1 * t103, t45, 0.2e1 * t20 * t23 - 0.2e1 * t50 * t87, 0.2e1 * t20 * t25 - 0.2e1 * t4 * t50, -0.2e1 * t4 * t23 + 0.2e1 * t25 * t87, t20 ^ 2 + t4 ^ 2 + t87 ^ 2, t22, t14, 0.2e1 * t106, t45, 0.2e1 * t103, t119, -0.2e1 * t2 * t50 + 0.2e1 * t5 * t23, -0.2e1 * t1 * t23 + 0.2e1 * t2 * t25, 0.2e1 * t1 * t50 - 0.2e1 * t5 * t25, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, t75, 0, -t74 * pkin(7), -t75 * pkin(7), 0, 0, 0, 0, t54, 0, -t50, 0, -t37, -t39 (-t50 * t70 - t54 * t72) * pkin(2) (-t37 * t72 + t39 * t70) * pkin(2), t90 (-t65 + t66) * t54, t101, -t90, t43, 0, -t37 * t71 + t79 * t69, t37 * t69 + t79 * t71, t80, t37 * t63 + t80 * t60, t16, -t125, t35, -t107, t34, 0, -t121 * t20 + t57 * t23 - t105, t20 * t55 + t57 * t25 - t104, t121 * t4 + t55 * t87 + t85, t20 * t57 + t29 * t87 + t4 * t31, t16, t35, t125, 0, -t34, -t107, -t121 * t5 + t21 * t23 - t105, t1 * t121 + t2 * t55 + t85, -t21 * t25 - t5 * t55 + t104, t1 * t31 + t2 * t29 + t5 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t110, -0.2e1 * t111, 0 (t70 ^ 2 + t72 ^ 2) * pkin(2) ^ 2, t65, 0.2e1 * t99, 0, t66, 0, 0, -0.2e1 * t63 * t71, 0.2e1 * t63 * t69, 0.2e1 * t95 * t60, t95 * t60 ^ 2 + t63 ^ 2, t49, 0.2e1 * t36, 0, t117, 0, 0, -t121 * t115, t55 * t115, t78, t57 ^ 2 + t91, t49, 0, -0.2e1 * t36, 0, 0, t117, t121 * t124, t78, t55 * t124, t21 ^ 2 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t54, 0, t64, 0, 0, 0, 0, 0, 0, t43, -t101, -t95 * t54, t81, 0, 0, 0, 0, 0, 0, t34, -t35, t122, -t121 * t87 + t4 * t55, 0, 0, 0, 0, 0, 0, t34, t122, t35, t1 * t55 - t121 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t98, 0, t37, 0, 0, 0, 0, 0, 0, t23, t25, 0, t20, 0, 0, 0, 0, 0, 0, t23, 0, -t25, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t69, 0, t63, 0, 0, 0, 0, 0, 0, -t121, t55, 0, t57, 0, 0, 0, 0, 0, 0, -t121, 0, -t55, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, t50, -t87, -t4, 0, 0, 0, t25, 0, t50, t23, 0, -t87 + 0.2e1 * t112, -pkin(5) * t25 - t23 * qJ(6), 0.2e1 * t93 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t121, 0, -t29, -t31, 0, 0, 0, t55, 0, 0, -t121, 0, -t29, -pkin(5) * t55 + qJ(6) * t121, t31, -t29 * pkin(5) + t31 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, -t55, 0, 0, 0, 0, 0, 0, 0, 0, t121, 0, t55, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t25, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
