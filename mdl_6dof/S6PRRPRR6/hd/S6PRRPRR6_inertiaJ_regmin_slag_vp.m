% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t61 = sin(pkin(13));
t64 = cos(pkin(13));
t65 = cos(pkin(7));
t62 = sin(pkin(7));
t69 = sin(qJ(3));
t95 = t62 * t69;
t39 = t61 * t95 - t64 * t65;
t41 = t61 * t65 + t64 * t95;
t68 = sin(qJ(5));
t98 = cos(qJ(5));
t25 = -t39 * t68 + t41 * t98;
t67 = sin(qJ(6));
t71 = cos(qJ(6));
t72 = cos(qJ(3));
t94 = t62 * t72;
t19 = t25 * t67 + t71 * t94;
t103 = -0.2e1 * t19;
t102 = -0.2e1 * t25;
t55 = -pkin(4) * t64 - pkin(3);
t101 = 0.2e1 * t55;
t100 = pkin(2) * t69;
t99 = pkin(2) * t72;
t20 = t25 * t71 - t67 * t94;
t97 = t20 * t67;
t57 = t62 ^ 2;
t96 = t57 * t72;
t63 = sin(pkin(6));
t70 = sin(qJ(2));
t93 = t63 * t70;
t73 = cos(qJ(2));
t92 = t63 * t73;
t91 = t65 * t69;
t90 = t65 * t73;
t24 = t39 * t98 + t41 * t68;
t89 = t67 * t24;
t47 = t61 * t68 - t64 * t98;
t88 = t67 * t47;
t48 = t61 * t98 + t64 * t68;
t87 = t67 * t48;
t86 = t67 * t71;
t23 = t71 * t24;
t85 = t71 * t48;
t84 = pkin(10) + qJ(4);
t80 = pkin(9) * t94;
t35 = t80 + (qJ(4) + t100) * t65;
t36 = (-pkin(3) * t72 - qJ(4) * t69 - pkin(2)) * t62;
t22 = t35 * t64 + t36 * t61;
t83 = t61 ^ 2 + t64 ^ 2;
t82 = -0.2e1 * t48 * t47;
t81 = 0.2e1 * t94;
t79 = qJ(4) * t94;
t78 = t84 * t61;
t21 = -t35 * t61 + t36 * t64;
t77 = -pkin(5) * t48 - pkin(11) * t47;
t66 = cos(pkin(6));
t29 = t66 * t95 + (t69 * t90 + t70 * t72) * t63;
t40 = -t62 * t92 + t65 * t66;
t17 = -t29 * t61 + t40 * t64;
t18 = t29 * t64 + t40 * t61;
t76 = -t17 * t61 + t18 * t64;
t75 = -t21 * t61 + t22 * t64;
t53 = pkin(9) * t95;
t38 = t53 + (-pkin(3) - t99) * t65;
t12 = -pkin(4) * t94 - pkin(10) * t41 + t21;
t15 = -pkin(10) * t39 + t22;
t7 = t12 * t98 - t15 * t68;
t8 = t12 * t68 + t15 * t98;
t26 = t39 * pkin(4) + t38;
t60 = t71 ^ 2;
t59 = t67 ^ 2;
t50 = t84 * t64;
t45 = t48 ^ 2;
t44 = pkin(2) * t91 + t80;
t43 = t65 * t99 - t53;
t42 = t71 * t47;
t31 = t50 * t98 - t68 * t78;
t30 = t50 * t68 + t78 * t98;
t28 = -t63 * t72 * t90 - t66 * t94 + t69 * t93;
t27 = pkin(5) * t47 - pkin(11) * t48 + t55;
t14 = t27 * t67 + t31 * t71;
t13 = t27 * t71 - t31 * t67;
t11 = t17 * t68 + t18 * t98;
t10 = -t17 * t98 + t18 * t68;
t9 = t24 * pkin(5) - t25 * pkin(11) + t26;
t6 = -pkin(11) * t94 + t8;
t5 = pkin(5) * t94 - t7;
t4 = t11 * t71 + t28 * t67;
t3 = -t11 * t67 + t28 * t71;
t2 = t6 * t71 + t67 * t9;
t1 = -t6 * t67 + t71 * t9;
t16 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 ^ 2 + t18 ^ 2 + t28 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t92, -t93, 0, 0, 0, 0, 0, -t28 * t65 - t40 * t94, -t29 * t65 + t40 * t95, -t17 * t94 + t28 * t39, t18 * t94 + t28 * t41, -t17 * t41 - t18 * t39, t17 * t21 + t18 * t22 + t28 * t38, 0, 0, 0, 0, 0, t10 * t94 + t24 * t28, t11 * t94 + t25 * t28, 0, 0, 0, 0, 0, t10 * t19 + t24 * t3, t10 * t20 - t24 * t4; 0, 1, 0, 0, t57 * t69 ^ 2, 0.2e1 * t69 * t96, 0.2e1 * t62 * t91, t65 * t81, t65 ^ 2, 0.2e1 * pkin(2) * t96 + 0.2e1 * t43 * t65, -0.2e1 * t100 * t57 - 0.2e1 * t44 * t65, -0.2e1 * t21 * t94 + 0.2e1 * t38 * t39, 0.2e1 * t22 * t94 + 0.2e1 * t38 * t41, -0.2e1 * t21 * t41 - 0.2e1 * t22 * t39, t21 ^ 2 + t22 ^ 2 + t38 ^ 2, t25 ^ 2, t24 * t102, t94 * t102, t24 * t81, t57 * t72 ^ 2, 0.2e1 * t24 * t26 - 0.2e1 * t7 * t94, 0.2e1 * t25 * t26 + 0.2e1 * t8 * t94, t20 ^ 2, t20 * t103, 0.2e1 * t20 * t24, t24 * t103, t24 ^ 2, 0.2e1 * t1 * t24 + 0.2e1 * t19 * t5, -0.2e1 * t2 * t24 + 0.2e1 * t20 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29, -t28 * t64, t28 * t61, t76, -t28 * pkin(3) + qJ(4) * t76, 0, 0, 0, 0, 0, t28 * t47, t28 * t48, 0, 0, 0, 0, 0, t10 * t87 + t3 * t47, t10 * t85 - t4 * t47; 0, 0, 0, 0, 0, 0, t95, t94, t65, t43, -t44, -pkin(3) * t39 - t38 * t64 + t61 * t79, -pkin(3) * t41 + t38 * t61 + t64 * t79 (-t39 * t64 + t41 * t61) * qJ(4) + t75, -t38 * pkin(3) + qJ(4) * t75, t25 * t48, -t24 * t48 - t25 * t47, -t48 * t94, t47 * t94, 0, t24 * t55 + t26 * t47 + t30 * t94, t25 * t55 + t26 * t48 + t31 * t94, t20 * t85 (-t19 * t71 - t97) * t48, t20 * t47 + t24 * t85, -t19 * t47 - t24 * t87, t24 * t47, t1 * t47 + t13 * t24 + t19 * t30 + t5 * t87, -t14 * t24 - t2 * t47 + t20 * t30 + t5 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t64, -0.2e1 * pkin(3) * t61, 0.2e1 * t83 * qJ(4), qJ(4) ^ 2 * t83 + pkin(3) ^ 2, t45, t82, 0, 0, 0, t47 * t101, t48 * t101, t60 * t45, -0.2e1 * t45 * t86, 0.2e1 * t47 * t85, t67 * t82, t47 ^ 2, 0.2e1 * t13 * t47 + 0.2e1 * t30 * t87, -0.2e1 * t14 * t47 + 0.2e1 * t30 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t41, 0, t38, 0, 0, 0, 0, 0, t24, t25, 0, 0, 0, 0, 0, t23, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t61, 0, -pkin(3), 0, 0, 0, 0, 0, t47, t48, 0, 0, 0, 0, 0, t42, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t10 * t71, t10 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, -t94, t7, -t8, t97, -t19 * t67 + t20 * t71, t89, t23, 0, -pkin(5) * t19 - pkin(11) * t89 - t5 * t71, -pkin(5) * t20 - pkin(11) * t23 + t5 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t47, 0, -t30, -t31, t67 * t85 (-t59 + t60) * t48, t88, t42, 0, -t30 * t71 + t67 * t77, t30 * t67 + t71 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t59, 0.2e1 * t86, 0, 0, 0, 0.2e1 * pkin(5) * t71, -0.2e1 * pkin(5) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t24, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t87, t47, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t71, 0, -t67 * pkin(11), -t71 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t16;
