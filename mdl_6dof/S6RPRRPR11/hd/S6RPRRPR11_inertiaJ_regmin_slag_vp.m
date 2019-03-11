% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR11
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
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t86 = cos(pkin(12));
t87 = cos(pkin(7));
t108 = t86 * t87;
t84 = sin(pkin(6));
t101 = t84 * t108;
t103 = qJ(2) * t84;
t88 = cos(pkin(6));
t119 = pkin(1) * t88;
t82 = sin(pkin(12));
t55 = t86 * t103 + t82 * t119;
t83 = sin(pkin(7));
t36 = (t83 * t88 + t101) * pkin(9) + t55;
t91 = sin(qJ(3));
t94 = cos(qJ(3));
t113 = t82 * t84;
t71 = t86 * t119;
t40 = t88 * pkin(2) + t71 + (-pkin(9) * t87 - qJ(2)) * t113;
t47 = (-pkin(9) * t82 * t83 - pkin(2) * t86 - pkin(1)) * t84;
t98 = t40 * t87 + t47 * t83;
t22 = -t91 * t36 + t98 * t94;
t112 = t83 * t91;
t39 = t88 * t112 + (t91 * t108 + t82 * t94) * t84;
t110 = t84 * t86;
t53 = -t83 * t110 + t88 * t87;
t90 = sin(qJ(4));
t93 = cos(qJ(4));
t32 = t39 * t90 - t53 * t93;
t125 = -0.2e1 * t32;
t111 = t83 * t94;
t38 = -t94 * t101 - t88 * t111 + t91 * t113;
t124 = -0.2e1 * t38;
t81 = sin(pkin(13));
t85 = cos(pkin(13));
t89 = sin(qJ(6));
t92 = cos(qJ(6));
t60 = t89 * t81 - t92 * t85;
t52 = t60 * t90;
t123 = 0.2e1 * t52;
t75 = -t85 * pkin(5) - pkin(4);
t122 = 0.2e1 * t75;
t121 = 0.2e1 * t93;
t78 = t84 ^ 2;
t120 = pkin(1) * t78;
t118 = pkin(10) * t81;
t76 = t90 * pkin(10);
t117 = t93 * pkin(10);
t20 = -t53 * pkin(3) - t22;
t33 = t39 * t93 + t53 * t90;
t15 = t32 * pkin(4) - t33 * qJ(5) + t20;
t31 = -t83 * t40 + t87 * t47;
t18 = t38 * pkin(3) - t39 * pkin(10) + t31;
t23 = t94 * t36 + t91 * t98;
t21 = t53 * pkin(10) + t23;
t12 = t90 * t18 + t93 * t21;
t9 = t38 * qJ(5) + t12;
t6 = t81 * t15 + t85 * t9;
t11 = t93 * t18 - t90 * t21;
t10 = -t38 * pkin(4) - t11;
t116 = t10 * t81;
t115 = t10 * t85;
t114 = t81 * t90;
t109 = t85 * t90;
t107 = t90 * t38;
t106 = t93 * t38;
t105 = pkin(11) + qJ(5);
t64 = -t93 * pkin(4) - t90 * qJ(5) - pkin(3);
t49 = t85 * t117 + t81 * t64;
t104 = t81 ^ 2 + t85 ^ 2;
t102 = qJ(5) * t32;
t5 = t85 * t15 - t81 * t9;
t100 = -t5 * t81 + t6 * t85;
t99 = -pkin(4) * t90 + qJ(5) * t93;
t57 = t93 * t112 + t90 * t87;
t41 = -t85 * t111 - t81 * t57;
t42 = -t81 * t111 + t85 * t57;
t97 = -t41 * t81 + t42 * t85;
t59 = t85 * t64;
t48 = -t81 * t117 + t59;
t96 = -t48 * t81 + t49 * t85;
t61 = t92 * t81 + t89 * t85;
t80 = t90 ^ 2;
t66 = t105 * t85;
t65 = t105 * t81;
t63 = pkin(5) * t114 + t76;
t56 = t90 * t112 - t93 * t87;
t54 = -t82 * t103 + t71;
t51 = t61 * t90;
t45 = -t89 * t65 + t92 * t66;
t44 = -t92 * t65 - t89 * t66;
t43 = -pkin(11) * t114 + t49;
t37 = -pkin(11) * t109 + t59 + (-pkin(5) - t118) * t93;
t30 = t89 * t41 + t92 * t42;
t29 = t92 * t41 - t89 * t42;
t28 = t89 * t37 + t92 * t43;
t27 = t92 * t37 - t89 * t43;
t26 = t33 * t85 + t38 * t81;
t25 = t33 * t81 - t38 * t85;
t17 = -t89 * t25 + t92 * t26;
t16 = t92 * t25 + t89 * t26;
t7 = t25 * pkin(5) + t10;
t4 = -t25 * pkin(11) + t6;
t3 = t32 * pkin(5) - t26 * pkin(11) + t5;
t2 = t89 * t3 + t92 * t4;
t1 = t92 * t3 - t89 * t4;
t8 = [1, 0, 0, 0.2e1 * t86 * t120 + 0.2e1 * t54 * t88, -0.2e1 * t82 * t120 - 0.2e1 * t55 * t88, 0.2e1 * (-t54 * t82 + t55 * t86) * t84, t78 * pkin(1) ^ 2 + t54 ^ 2 + t55 ^ 2, t39 ^ 2, t39 * t124, 0.2e1 * t39 * t53, t53 * t124, t53 ^ 2, 0.2e1 * t22 * t53 + 0.2e1 * t31 * t38, -0.2e1 * t23 * t53 + 0.2e1 * t31 * t39, t33 ^ 2, t33 * t125, 0.2e1 * t33 * t38, t32 * t124, t38 ^ 2, 0.2e1 * t11 * t38 + 0.2e1 * t20 * t32, -0.2e1 * t12 * t38 + 0.2e1 * t20 * t33, 0.2e1 * t10 * t25 + 0.2e1 * t5 * t32, 0.2e1 * t10 * t26 - 0.2e1 * t6 * t32, -0.2e1 * t6 * t25 - 0.2e1 * t5 * t26, t10 ^ 2 + t5 ^ 2 + t6 ^ 2, t17 ^ 2, -0.2e1 * t17 * t16, 0.2e1 * t17 * t32, t16 * t125, t32 ^ 2, 0.2e1 * t1 * t32 + 0.2e1 * t7 * t16, 0.2e1 * t7 * t17 - 0.2e1 * t2 * t32; 0, 0, 0, -t110, t113, 0, -t84 * pkin(1), 0, 0, 0, 0, 0, t53 * t111 + t87 * t38, -t53 * t112 + t87 * t39, 0, 0, 0, 0, 0, -t32 * t111 - t56 * t38, -t33 * t111 - t57 * t38, t56 * t25 + t41 * t32, t56 * t26 - t42 * t32, -t42 * t25 - t41 * t26, t10 * t56 + t5 * t41 + t6 * t42, 0, 0, 0, 0, 0, t56 * t16 + t29 * t32, t56 * t17 - t30 * t32; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 ^ 2 + t42 ^ 2 + t56 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t38, t53, t22, -t23, t33 * t90, -t90 * t32 + t33 * t93, t107, t106, 0, -pkin(3) * t32 - pkin(10) * t107 - t20 * t93, -pkin(3) * t33 - pkin(10) * t106 + t20 * t90, t48 * t32 - t5 * t93 + (pkin(10) * t25 + t116) * t90, -t49 * t32 + t6 * t93 + (pkin(10) * t26 + t115) * t90, -t49 * t25 - t48 * t26 + (-t5 * t85 - t6 * t81) * t90, t10 * t76 + t5 * t48 + t6 * t49, -t17 * t52, t52 * t16 - t17 * t51, -t17 * t93 - t52 * t32, t16 * t93 - t51 * t32, -t32 * t93, -t1 * t93 + t63 * t16 + t27 * t32 + t7 * t51, t63 * t17 + t2 * t93 - t28 * t32 - t7 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t112, 0, 0, 0, 0, 0, t93 * t111, -t90 * t111, t56 * t114 - t41 * t93, t56 * t109 + t42 * t93 (-t41 * t85 - t42 * t81) * t90, t41 * t48 + t42 * t49 + t56 * t76, 0, 0, 0, 0, 0, -t29 * t93 + t56 * t51, t30 * t93 - t56 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t80, t90 * t121, 0, 0, 0, pkin(3) * t121, -0.2e1 * pkin(3) * t90, 0.2e1 * t80 * t118 - 0.2e1 * t48 * t93, 0.2e1 * t80 * pkin(10) * t85 + 0.2e1 * t49 * t93, 0.2e1 * (-t48 * t85 - t49 * t81) * t90, t80 * pkin(10) ^ 2 + t48 ^ 2 + t49 ^ 2, t52 ^ 2, t51 * t123, t93 * t123, t51 * t121, t93 ^ 2, -0.2e1 * t27 * t93 + 0.2e1 * t63 * t51, 0.2e1 * t28 * t93 - 0.2e1 * t63 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, t38, t11, -t12, -pkin(4) * t25 - t81 * t102 - t115, -pkin(4) * t26 - t85 * t102 + t116 (-t25 * t85 + t26 * t81) * qJ(5) + t100, -t10 * pkin(4) + t100 * qJ(5), t17 * t61, -t61 * t16 - t17 * t60, t61 * t32, -t60 * t32, 0, t75 * t16 + t44 * t32 + t7 * t60, t75 * t17 - t45 * t32 + t7 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t57, -t56 * t85, t56 * t81, t97, -t56 * pkin(4) + qJ(5) * t97, 0, 0, 0, 0, 0, t56 * t60, t56 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t93, 0, -t76, -t117, -pkin(10) * t109 + t99 * t81, pkin(10) * t114 + t99 * t85, t96, -pkin(4) * t76 + qJ(5) * t96, -t52 * t61, -t61 * t51 + t52 * t60, -t61 * t93, t60 * t93, 0, -t44 * t93 + t75 * t51 + t63 * t60, t45 * t93 - t75 * t52 + t63 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t85, -0.2e1 * pkin(4) * t81, 0.2e1 * t104 * qJ(5), t104 * qJ(5) ^ 2 + pkin(4) ^ 2, t61 ^ 2, -0.2e1 * t61 * t60, 0, 0, 0, t60 * t122, t61 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t26, 0, t10, 0, 0, 0, 0, 0, t16, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t109, 0, t76, 0, 0, 0, 0, 0, t51, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t81, 0, -pkin(4), 0, 0, 0, 0, 0, t60, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t32, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t51, -t93, t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t60, 0, t44, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
