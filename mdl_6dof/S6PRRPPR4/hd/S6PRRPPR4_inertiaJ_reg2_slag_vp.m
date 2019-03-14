% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPPR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t67 = sin(pkin(11));
t62 = t67 ^ 2;
t69 = cos(pkin(11));
t64 = t69 ^ 2;
t113 = t62 + t64;
t71 = sin(qJ(6));
t74 = cos(qJ(6));
t37 = t74 * t67 - t71 * t69;
t70 = cos(pkin(6));
t72 = sin(qJ(3));
t75 = cos(qJ(3));
t68 = sin(pkin(6));
t73 = sin(qJ(2));
t97 = t68 * t73;
t29 = -t70 * t75 + t72 * t97;
t112 = t29 ^ 2;
t24 = t37 * t72;
t111 = 0.2e1 * t24;
t106 = pkin(4) + pkin(5);
t86 = t67 * qJ(5) + pkin(3);
t32 = t106 * t69 + t86;
t110 = 0.2e1 * t32;
t109 = -0.2e1 * t69;
t108 = 0.2e1 * t72;
t107 = 0.2e1 * t75;
t105 = pkin(3) * t67;
t65 = t72 ^ 2;
t104 = t65 * pkin(8);
t103 = t72 * pkin(8);
t31 = t70 * t72 + t75 * t97;
t76 = cos(qJ(2));
t96 = t68 * t76;
t13 = t31 * t69 - t67 * t96;
t9 = t13 * t69;
t102 = t29 * t67;
t101 = t29 * t69;
t100 = t29 * t72;
t99 = t67 * t69;
t51 = t67 * t72;
t98 = t67 * t75;
t55 = t69 * t72;
t95 = t69 * t75;
t93 = t72 * t75;
t40 = -t75 * pkin(3) - t72 * qJ(4) - pkin(2);
t22 = pkin(8) * t95 + t67 * t40;
t91 = t113 * qJ(4) ^ 2;
t90 = qJ(4) * t75;
t58 = t67 * qJ(4);
t89 = t67 * t93;
t88 = t65 * t99;
t11 = t31 * t67 + t69 * t96;
t87 = t11 * t67 + t9;
t85 = t11 * t75 + t29 * t51;
t48 = pkin(8) * t98;
t21 = t69 * t40 - t48;
t84 = t11 ^ 2 + t13 ^ 2 + t112;
t83 = qJ(4) * t9 + t11 * t58;
t18 = -t75 * qJ(5) + t22;
t61 = t75 * pkin(4);
t19 = -t21 + t61;
t82 = t18 * t69 + t19 * t67;
t81 = -t21 * t67 + t22 * t69;
t80 = t31 * t75 + t100;
t35 = t71 * t67 + t74 * t69;
t79 = (t11 * t69 - t13 * t67) * t72;
t78 = pkin(8) ^ 2;
t66 = t75 ^ 2;
t63 = t68 ^ 2;
t60 = t65 * t78;
t54 = t64 * t65;
t53 = t63 * t76 ^ 2;
t50 = t62 * t65;
t47 = qJ(5) * t55;
t46 = t75 * t58;
t44 = t67 * t55;
t43 = t93 * t109;
t42 = (-pkin(9) + qJ(4)) * t69;
t41 = -t67 * pkin(9) + t58;
t39 = -t69 * pkin(4) - t86;
t38 = 0.2e1 * t113 * qJ(4);
t34 = (t62 - t64) * t72;
t26 = t35 * t72;
t23 = -t47 + (pkin(4) * t67 + pkin(8)) * t72;
t16 = -t47 - (-t106 * t67 - pkin(8)) * t72;
t15 = t71 * t41 + t74 * t42;
t14 = t74 * t41 - t71 * t42;
t8 = pkin(9) * t51 + t18;
t6 = t75 * pkin(5) + t48 + t61 + (-pkin(9) * t72 - t40) * t69;
t5 = t13 * t75 + t29 * t55;
t4 = t11 * t71 + t13 * t74;
t3 = t11 * t74 - t13 * t71;
t2 = t71 * t6 + t74 * t8;
t1 = t74 * t6 - t71 * t8;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t73 ^ 2 + t70 ^ 2 + t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 ^ 2 + t112 + t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, -t97, 0, 0, 0, 0, 0, 0, 0, 0, t75 * t96, -t72 * t96, t80, pkin(2) * t96 + t80 * pkin(8), 0, 0, 0, 0, 0, 0, t85, t5, t79, pkin(8) * t100 - t11 * t21 + t13 * t22, 0, 0, 0, 0, 0, 0, t85, t79, -t5, t11 * t19 + t13 * t18 + t29 * t23, 0, 0, 0, 0, 0, 0, t29 * t24 + t3 * t75, -t29 * t26 - t4 * t75, t4 * t24 - t3 * t26, t3 * t1 + t29 * t16 + t4 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t65, 0.2e1 * t93, 0, t66, 0, 0, pkin(2) * t107, -0.2e1 * pkin(2) * t72, 0.2e1 * (t65 + t66) * pkin(8), pkin(2) ^ 2 + t66 * t78 + t60, t54, -0.2e1 * t88, t43, t50, 0.2e1 * t89, t66, 0.2e1 * t67 * t104 - 0.2e1 * t21 * t75, 0.2e1 * t69 * t104 + 0.2e1 * t22 * t75 (-t21 * t69 - t22 * t67) * t108, t21 ^ 2 + t22 ^ 2 + t60, t54, t43, 0.2e1 * t88, t66, -0.2e1 * t89, t50, 0.2e1 * t19 * t75 + 0.2e1 * t23 * t51 (-t18 * t67 + t19 * t69) * t108, -0.2e1 * t18 * t75 - 0.2e1 * t23 * t55, t18 ^ 2 + t19 ^ 2 + t23 ^ 2, t26 ^ 2, t26 * t111, t26 * t107, t24 ^ 2, t75 * t111, t66, 0.2e1 * t1 * t75 + 0.2e1 * t16 * t24, -0.2e1 * t16 * t26 - 0.2e1 * t2 * t75, -0.2e1 * t1 * t26 + 0.2e1 * t2 * t24, t1 ^ 2 + t16 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t31, 0, 0, 0, 0, 0, 0, 0, 0, -t101, t102, t87, -t29 * pkin(3) + t83, 0, 0, 0, 0, 0, 0, -t101, t87, -t102, t29 * t39 + t83, 0, 0, 0, 0, 0, 0, -t29 * t35, -t29 * t37, -t3 * t37 - t4 * t35, t3 * t14 + t4 * t15 - t29 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, t75, 0, -t103, -t75 * pkin(8), 0, 0, t44, -t34, -t98, -t44, -t95, 0, t46 + (-pkin(8) * t69 - t105) * t72, pkin(8) * t51 + (-pkin(3) * t72 + t90) * t69, t81, -pkin(3) * t103 + t81 * qJ(4), t44, -t98, t34, 0, t95, -t44, -t23 * t69 + t39 * t51 + t46, t82, -t23 * t67 + (-t39 * t72 - t90) * t69, t82 * qJ(4) + t23 * t39, t26 * t37, t37 * t24 - t26 * t35, t37 * t75, -t24 * t35, -t35 * t75, 0, t14 * t75 - t16 * t35 - t32 * t24, -t15 * t75 - t16 * t37 + t32 * t26, -t1 * t37 - t14 * t26 + t15 * t24 - t2 * t35, t1 * t14 + t2 * t15 - t16 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t62, 0.2e1 * t99, 0, t64, 0, 0, 0.2e1 * pkin(3) * t69, -0.2e1 * t105, t38, pkin(3) ^ 2 + t91, t62, 0, -0.2e1 * t99, 0, 0, t64, t39 * t109, t38, -0.2e1 * t39 * t67, t39 ^ 2 + t91, t37 ^ 2, -0.2e1 * t37 * t35, 0, t35 ^ 2, 0, 0, t35 * t110, t37 * t110, -0.2e1 * t14 * t37 - 0.2e1 * t15 * t35, t14 ^ 2 + t15 ^ 2 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t55, 0, t103, 0, 0, 0, 0, 0, 0, t51, 0, -t55, t23, 0, 0, 0, 0, 0, 0, t24, -t26, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t67, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t69, 0, -t67, t39, 0, 0, 0, 0, 0, 0, -t35, -t37, 0, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t74 + t4 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t55, 0, t19, 0, 0, 0, 0, 0, 0, t74 * t75, -t71 * t75, t71 * t24 - t74 * t26, t1 * t74 + t2 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, -t71 * t35 - t74 * t37, t14 * t74 + t15 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 ^ 2 + t74 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t24, t75, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t35, 0, t14, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t71, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t7;