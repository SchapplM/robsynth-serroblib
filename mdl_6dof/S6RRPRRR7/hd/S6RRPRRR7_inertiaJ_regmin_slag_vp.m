% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t56 = sin(qJ(6));
t57 = sin(qJ(5));
t60 = cos(qJ(6));
t61 = cos(qJ(5));
t31 = t56 * t61 + t60 * t57;
t113 = t31 ^ 2;
t58 = sin(qJ(4));
t62 = cos(qJ(4));
t64 = -pkin(2) - pkin(3);
t35 = t58 * qJ(3) - t62 * t64;
t33 = pkin(4) + t35;
t99 = t61 * pkin(5);
t24 = t33 + t99;
t112 = -0.2e1 * t24;
t59 = sin(qJ(2));
t63 = cos(qJ(2));
t29 = t59 * t58 + t63 * t62;
t111 = -0.2e1 * t29;
t110 = 0.2e1 * t29;
t45 = -pkin(4) - t99;
t109 = 0.2e1 * t45;
t108 = -0.2e1 * t57;
t107 = -0.2e1 * t59;
t106 = 0.2e1 * t61;
t105 = 0.2e1 * t63;
t104 = pkin(9) + pkin(10);
t103 = t29 * pkin(5);
t102 = t56 * pkin(5);
t101 = t60 * pkin(5);
t37 = -t63 * pkin(2) - t59 * qJ(3) - pkin(1);
t25 = t63 * pkin(3) - t37;
t30 = -t63 * t58 + t59 * t62;
t11 = t29 * pkin(4) - t30 * pkin(9) + t25;
t48 = t59 * pkin(7);
t39 = -t59 * pkin(8) + t48;
t49 = t63 * pkin(7);
t41 = -t63 * pkin(8) + t49;
t19 = t58 * t39 + t62 * t41;
t84 = t61 * t19;
t5 = t84 + (-pkin(10) * t30 + t11) * t57;
t100 = t60 * t5;
t98 = pkin(4) + t33;
t36 = t62 * qJ(3) + t58 * t64;
t34 = -pkin(9) + t36;
t97 = pkin(10) - t34;
t17 = -t62 * t39 + t58 * t41;
t87 = t57 * t30;
t12 = pkin(5) * t87 + t17;
t28 = t56 * t57 - t60 * t61;
t96 = t12 * t28;
t95 = t12 * t31;
t14 = t28 * t30;
t94 = t14 * t31;
t93 = t17 * t57;
t92 = t17 * t61;
t91 = t28 * t29;
t90 = t31 * t28;
t89 = t31 * t29;
t88 = t57 * t29;
t86 = t57 * t58;
t85 = t57 * t61;
t83 = t61 * t29;
t82 = t61 * t30;
t81 = t61 * t58;
t80 = t62 * t28;
t79 = t62 * t31;
t78 = t62 * t57;
t77 = t62 * t61;
t76 = t24 - t45;
t53 = t59 ^ 2;
t75 = t63 ^ 2 + t53;
t74 = t30 * t111;
t73 = -0.2e1 * t90;
t72 = -0.2e1 * t85;
t71 = t57 * t82;
t6 = t61 * t11 - t57 * t19;
t4 = -pkin(10) * t82 + t103 + t6;
t1 = t60 * t4 - t56 * t5;
t70 = -pkin(4) * t30 - pkin(9) * t29;
t69 = -t59 * pkin(2) + t63 * qJ(3);
t13 = t31 * t30;
t68 = t31 * t13 - t14 * t28;
t67 = -t29 * t34 + t30 * t33;
t66 = -t29 * t58 - t30 * t62;
t54 = t61 ^ 2;
t52 = t57 ^ 2;
t43 = 0.2e1 * t85;
t40 = t104 * t61;
t38 = t104 * t57;
t27 = t30 ^ 2;
t26 = t29 ^ 2;
t23 = -t56 * t86 + t60 * t81;
t22 = t31 * t58;
t21 = t97 * t61;
t20 = t97 * t57;
t18 = -t56 * t38 + t60 * t40;
t16 = -t60 * t38 - t56 * t40;
t15 = (t52 - t54) * t30;
t9 = t56 * t20 - t60 * t21;
t8 = t60 * t20 + t56 * t21;
t7 = t57 * t11 + t84;
t2 = t56 * t4 + t100;
t3 = [1, 0, 0, t53, t59 * t105, 0, 0, 0, pkin(1) * t105, pkin(1) * t107, -0.2e1 * t37 * t63, 0.2e1 * t75 * pkin(7), t37 * t107, pkin(7) ^ 2 * t75 + t37 ^ 2, t27, t74, 0, 0, 0, t25 * t110, 0.2e1 * t25 * t30, t54 * t27, t27 * t72, t82 * t110, t57 * t74, t26, 0.2e1 * t17 * t87 + 0.2e1 * t6 * t29, 0.2e1 * t17 * t82 - 0.2e1 * t7 * t29, t14 ^ 2, 0.2e1 * t14 * t13, -t14 * t110, t13 * t111, t26, 0.2e1 * t1 * t29 + 0.2e1 * t12 * t13, -0.2e1 * t12 * t14 - 0.2e1 * t2 * t29; 0, 0, 0, 0, 0, t59, t63, 0, -t48, -t49, -t48, t69, t49, t69 * pkin(7), 0, 0, -t30, t29, 0, t17, t19, -t71, t15, -t88, -t83, 0, t57 * t67 + t92, t61 * t67 - t93, t94, t68, -t89, t91, 0, t24 * t13 + t8 * t29 - t96, -t24 * t14 - t9 * t29 - t95; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t35, 0.2e1 * t36, t52, t43, 0, 0, 0, t33 * t106, t33 * t108, t113, t73, 0, 0, 0, t28 * t112, t31 * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t57, t66 * t61, 0, 0, 0, 0, 0, -t62 * t13 - t22 * t29, t62 * t14 - t23 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, -t62, t58, 0, 0, 0, 0, 0, -t77, t78, 0, 0, 0, 0, 0, t80, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, -t17, -t19, t71, -t15, t88, t83, 0, t57 * t70 - t92, t61 * t70 + t93, -t94, -t68, t89, -t91, 0, t45 * t13 + t16 * t29 + t96, -t45 * t14 - t18 * t29 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t35, -t36, -t52, t72, 0, 0, 0, -t98 * t61, t98 * t57, -t113, 0.2e1 * t90, 0, 0, 0, t76 * t28, t76 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t58, 0, 0, 0, 0, 0, t77, -t78, 0, 0, 0, 0, 0, -t80, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t52, t43, 0, 0, 0, pkin(4) * t106, pkin(4) * t108, t113, t73, 0, 0, 0, t28 * t109, t31 * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t87, t29, t6, -t7, 0, 0, -t14, -t13, t29, t29 * t101 + t1, -t100 + (-t4 - t103) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t61, 0, -t57 * t34, -t61 * t34, 0, 0, -t31, t28, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t81, 0, 0, 0, 0, 0, -t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t61, 0, -t57 * pkin(9), -t61 * pkin(9), 0, 0, t31, -t28, 0, t16, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t101, -0.2e1 * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t13, t29, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t28, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t28, 0, t16, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t101, -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
