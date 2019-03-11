% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP6
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
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t63 = sin(pkin(11));
t64 = sin(pkin(6));
t65 = cos(pkin(11));
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t39 = (t63 * t72 + t65 * t69) * t64;
t66 = cos(pkin(6));
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t30 = t39 * t71 + t66 * t68;
t103 = t64 * t72;
t104 = t64 * t69;
t38 = -t65 * t103 + t63 * t104;
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t21 = t30 * t70 + t38 * t67;
t109 = t21 * t71;
t29 = t39 * t68 - t66 * t71;
t56 = t70 * t68;
t23 = t29 * t56;
t122 = t23 + t109;
t121 = -0.2e1 * t29;
t78 = -t70 * pkin(5) - t67 * qJ(6);
t47 = -pkin(4) + t78;
t120 = -0.2e1 * t47;
t119 = 0.2e1 * t68;
t49 = t66 * t72 * pkin(1);
t92 = pkin(8) + qJ(3);
t31 = t66 * pkin(2) - t92 * t104 + t49;
t118 = pkin(1) * t69;
t86 = t66 * t118;
t35 = t92 * t103 + t86;
t18 = t65 * t31 - t63 * t35;
t15 = -t66 * pkin(3) - t18;
t11 = t29 * pkin(4) - t30 * pkin(10) + t15;
t19 = t63 * t31 + t65 * t35;
t16 = t66 * pkin(9) + t19;
t45 = (-pkin(2) * t72 - pkin(1)) * t64;
t22 = t38 * pkin(3) - t39 * pkin(9) + t45;
t13 = t71 * t16 + t68 * t22;
t8 = t38 * pkin(10) + t13;
t4 = t67 * t11 + t70 * t8;
t117 = pkin(4) * t67;
t116 = pkin(4) * t70;
t115 = t29 * pkin(5);
t114 = t67 * pkin(10);
t12 = -t68 * t16 + t71 * t22;
t7 = -t38 * pkin(4) - t12;
t113 = t7 * t67;
t112 = t7 * t70;
t111 = t70 * pkin(10);
t110 = t21 * t67;
t52 = t63 * pkin(2) + pkin(9);
t108 = t52 * t67;
t107 = t52 * t70;
t58 = t64 ^ 2;
t106 = t58 * t72;
t59 = t67 ^ 2;
t105 = t59 * t68;
t102 = t67 * t29;
t101 = t67 * t68;
t100 = t67 * t70;
t99 = t67 * t71;
t98 = t68 * t29;
t97 = t68 * t38;
t96 = t68 * t52;
t95 = t70 * t29;
t57 = t70 * t71;
t20 = t30 * t67 - t38 * t70;
t94 = t71 * t20;
t93 = t71 * t52;
t53 = -t65 * pkin(2) - pkin(3);
t44 = -t71 * pkin(4) - t68 * pkin(10) + t53;
t28 = t67 * t44 + t70 * t93;
t61 = t70 ^ 2;
t91 = t59 + t61;
t90 = t29 * qJ(6);
t89 = t71 * qJ(6);
t88 = 0.2e1 * t64 * t66;
t87 = t71 * t119;
t85 = pkin(10) * t102;
t84 = pkin(10) * t95;
t83 = t21 * t101;
t82 = t67 * t98;
t81 = t91 * pkin(10);
t80 = -t70 * t11 + t67 * t8;
t1 = t90 + t4;
t2 = t80 - t115;
t79 = t1 * t70 + t2 * t67;
t77 = -pkin(5) * t67 + t70 * qJ(6);
t24 = -t89 + t28;
t41 = t70 * t44;
t25 = -t41 + (pkin(5) + t108) * t71;
t76 = t24 * t70 + t25 * t67;
t75 = -t82 - t94;
t62 = t71 ^ 2;
t60 = t68 ^ 2;
t55 = t61 * t68;
t54 = t61 * t60;
t50 = pkin(10) * t99;
t43 = pkin(8) * t103 + t86;
t42 = -pkin(8) * t104 + t49;
t36 = t71 * t38;
t34 = (t52 - t77) * t68;
t27 = -t67 * t93 + t41;
t17 = t20 * t56;
t5 = t20 * pkin(5) - t21 * qJ(6) + t7;
t3 = [1, 0, 0, t58 * t69 ^ 2, 0.2e1 * t69 * t106, t69 * t88, t72 * t88, t66 ^ 2, 0.2e1 * pkin(1) * t106 + 0.2e1 * t42 * t66, -0.2e1 * t58 * t118 - 0.2e1 * t43 * t66, -0.2e1 * t18 * t39 - 0.2e1 * t19 * t38, t18 ^ 2 + t19 ^ 2 + t45 ^ 2, t30 ^ 2, t30 * t121, 0.2e1 * t30 * t38, t38 * t121, t38 ^ 2, 0.2e1 * t12 * t38 + 0.2e1 * t15 * t29, -0.2e1 * t13 * t38 + 0.2e1 * t15 * t30, t21 ^ 2, -0.2e1 * t21 * t20, 0.2e1 * t21 * t29, t20 * t121, t29 ^ 2, 0.2e1 * t7 * t20 - 0.2e1 * t29 * t80, 0.2e1 * t7 * t21 - 0.2e1 * t4 * t29, -0.2e1 * t2 * t29 + 0.2e1 * t5 * t20, -0.2e1 * t1 * t20 + 0.2e1 * t2 * t21, 0.2e1 * t1 * t29 - 0.2e1 * t5 * t21, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t104, t103, t66, t42, -t43 (-t38 * t63 - t39 * t65) * pkin(2) (t18 * t65 + t19 * t63) * pkin(2), t30 * t68, t30 * t71 - t98, t97, t36, 0, -t15 * t71 + t53 * t29 - t38 * t96, t15 * t68 + t53 * t30 - t38 * t93, t21 * t56, -t17 - t83, t23 - t109, -t82 + t94, -t29 * t71, t27 * t29 + t80 * t71 + (t20 * t52 + t113) * t68, -t28 * t29 + t4 * t71 + (t21 * t52 + t112) * t68, t101 * t5 + t2 * t71 + t34 * t20 - t25 * t29, -t24 * t20 + t25 * t21 + (-t1 * t67 + t2 * t70) * t68, -t1 * t71 - t34 * t21 + t24 * t29 - t5 * t56, t1 * t24 + t2 * t25 + t5 * t34; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t63 ^ 2 + t65 ^ 2) * pkin(2) ^ 2, t60, t87, 0, 0, 0, -0.2e1 * t53 * t71, t53 * t119, t54, -0.2e1 * t60 * t100, -0.2e1 * t68 * t57, t67 * t87, t62, 0.2e1 * t60 * t108 - 0.2e1 * t27 * t71, 0.2e1 * t60 * t107 + 0.2e1 * t28 * t71, 0.2e1 * t101 * t34 + 0.2e1 * t25 * t71 (-t24 * t67 + t25 * t70) * t119, -0.2e1 * t24 * t71 - 0.2e1 * t34 * t56, t24 ^ 2 + t25 ^ 2 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, t36, -t97, 0, 0, 0, 0, 0, t75, -t122, t75, -t17 + t83, t122, -t5 * t71 + t68 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34 * t71 + t68 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t60 + t54 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, t38, t12, -t13, t110, -t67 * t20 + t21 * t70, t102, t95, 0, -pkin(4) * t20 - t112 - t85, -pkin(4) * t21 + t113 - t84, t47 * t20 - t5 * t70 - t85 (-t20 * t70 + t110) * pkin(10) + t79, -t47 * t21 - t5 * t67 + t84, pkin(10) * t79 + t5 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t71, 0, -t96, -t93, t67 * t56, t55 - t105, -t99, -t57, 0, t50 + (-t107 - t117) * t68, pkin(10) * t57 + (t108 - t116) * t68, t101 * t47 - t34 * t70 + t50, t76, -t34 * t67 + (-pkin(10) * t71 - t47 * t68) * t70, pkin(10) * t76 + t34 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t68, 0, 0, 0, 0, 0, t57, -t99, t57, t55 + t105, t99, -t71 * t47 + t68 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t59, 0.2e1 * t100, 0, 0, 0, 0.2e1 * t116, -0.2e1 * t117, t70 * t120, 0.2e1 * t81, t67 * t120, pkin(10) ^ 2 * t91 + t47 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, t29, -t80, -t4, -t80 + 0.2e1 * t115, -t21 * pkin(5) - t20 * qJ(6), 0.2e1 * t90 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t101, -t71, t27, -t28, t41 + (-0.2e1 * pkin(5) - t108) * t71, t78 * t68, -0.2e1 * t89 + t28, -t25 * pkin(5) + t24 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t56, -t101, 0, t56, t77 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t70, 0, -t114, -t111, -t114, t77, t111, t77 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t21, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t56, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
