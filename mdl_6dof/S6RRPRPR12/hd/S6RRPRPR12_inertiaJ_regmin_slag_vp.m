% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t100 = sin(qJ(4));
t65 = sin(pkin(11));
t67 = cos(pkin(11));
t72 = cos(qJ(4));
t42 = -t65 * t100 + t67 * t72;
t43 = -t67 * t100 - t65 * t72;
t109 = (t42 * t67 - t43 * t65) * pkin(4);
t40 = t42 ^ 2;
t41 = t43 ^ 2;
t108 = t41 + t40;
t107 = -2 * pkin(2);
t68 = cos(pkin(6));
t66 = sin(pkin(6));
t73 = cos(qJ(2));
t96 = t66 * t73;
t35 = t68 * t100 + t72 * t96;
t84 = t66 * t100;
t36 = t68 * t72 - t73 * t84;
t21 = -t65 * t35 + t67 * t36;
t70 = sin(qJ(2));
t55 = t66 * t70;
t69 = sin(qJ(6));
t71 = cos(qJ(6));
t16 = t69 * t21 - t71 * t55;
t106 = -0.2e1 * t16;
t105 = -0.2e1 * t35;
t104 = 0.2e1 * t66;
t103 = 2 * qJ(3);
t74 = -pkin(2) - pkin(9);
t102 = pkin(1) * t70;
t101 = pkin(1) * t73;
t48 = pkin(8) * t55;
t85 = -pkin(2) - t101;
t24 = pkin(3) * t55 + t48 + (-pkin(9) + t85) * t68;
t82 = -qJ(3) * t70 - pkin(1);
t29 = (t74 * t73 + t82) * t66;
t15 = t100 * t24 + t72 * t29;
t11 = -t35 * qJ(5) + t15;
t14 = -t100 * t29 + t72 * t24;
t9 = pkin(4) * t55 - t36 * qJ(5) + t14;
t6 = t67 * t11 + t65 * t9;
t17 = t71 * t21 + t69 * t55;
t99 = t17 * t69;
t20 = t67 * t35 + t65 * t36;
t98 = t20 * t43;
t62 = t66 ^ 2;
t97 = t62 * t73;
t95 = t68 * t73;
t94 = t69 * t20;
t93 = t69 * t42;
t92 = t69 * t43;
t57 = t65 * pkin(4) + pkin(10);
t91 = t69 * t57;
t90 = t69 * t71;
t89 = t71 * t42;
t34 = t71 * t43;
t88 = t71 * t57;
t38 = pkin(8) * t96 + t68 * t102;
t59 = t100 * pkin(4) + qJ(3);
t87 = 0.2e1 * t55;
t47 = t72 * t55;
t60 = t68 * qJ(3);
t30 = -t60 - t38;
t83 = t100 * t74;
t81 = (-qJ(5) + t74) * t72;
t28 = pkin(3) * t96 - t30;
t80 = t70 * t84;
t5 = -t65 * t11 + t67 * t9;
t79 = t5 * t42 - t6 * t43;
t45 = -t100 * qJ(5) + t83;
t25 = t65 * t45 - t67 * t81;
t27 = t67 * t45 + t65 * t81;
t78 = t25 * t42 + t27 * t43;
t58 = -t67 * pkin(4) - pkin(5);
t77 = t42 * t58 + t43 * t57;
t18 = t35 * pkin(4) + t28;
t64 = t71 ^ 2;
t63 = t69 ^ 2;
t54 = t62 * t70 ^ 2;
t37 = pkin(1) * t95 - t48;
t32 = t85 * t68 + t48;
t31 = (-pkin(2) * t73 + t82) * t66;
t22 = -t43 * pkin(5) - t42 * pkin(10) + t59;
t19 = t71 * t20;
t13 = t69 * t22 + t71 * t27;
t12 = t71 * t22 - t69 * t27;
t7 = t20 * pkin(5) - t21 * pkin(10) + t18;
t4 = pkin(10) * t55 + t6;
t3 = -pkin(5) * t55 - t5;
t2 = t71 * t4 + t69 * t7;
t1 = -t69 * t4 + t71 * t7;
t8 = [1, 0, 0, t54, 0.2e1 * t70 * t97, t68 * t87, t95 * t104, t68 ^ 2, 0.2e1 * pkin(1) * t97 + 0.2e1 * t37 * t68, -0.2e1 * t62 * t102 - 0.2e1 * t38 * t68 (-t30 * t73 + t32 * t70) * t104, 0.2e1 * t31 * t96 + 0.2e1 * t32 * t68, -0.2e1 * t30 * t68 - 0.2e1 * t31 * t55, t30 ^ 2 + t31 ^ 2 + t32 ^ 2, t36 ^ 2, t36 * t105, t36 * t87, t55 * t105, t54, 0.2e1 * t14 * t55 + 0.2e1 * t28 * t35, -0.2e1 * t15 * t55 + 0.2e1 * t28 * t36, -0.2e1 * t6 * t20 - 0.2e1 * t5 * t21, t18 ^ 2 + t5 ^ 2 + t6 ^ 2, t17 ^ 2, t17 * t106, 0.2e1 * t17 * t20, t20 * t106, t20 ^ 2, 0.2e1 * t1 * t20 + 0.2e1 * t3 * t16, 0.2e1 * t3 * t17 - 0.2e1 * t2 * t20; 0, 0, 0, 0, 0, t55, t96, t68, t37, -t38 (-pkin(2) * t70 + qJ(3) * t73) * t66, t48 + (t107 - t101) * t68, 0.2e1 * t60 + t38, -t32 * pkin(2) - t30 * qJ(3), t36 * t72, -t36 * t100 - t72 * t35, t47, -t80, 0, qJ(3) * t35 + t28 * t100 + t74 * t47, qJ(3) * t36 + t28 * t72 - t74 * t80, -t27 * t20 + t25 * t21 - t79, t18 * t59 - t5 * t25 + t6 * t27, t17 * t89 (-t16 * t71 - t99) * t42, -t17 * t43 + t20 * t89, t16 * t43 - t20 * t93, -t98, -t1 * t43 + t12 * t20 + t25 * t16 + t3 * t93, -t13 * t20 + t25 * t17 + t2 * t43 + t3 * t89; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t107, t103, pkin(2) ^ 2 + qJ(3) ^ 2, t72 ^ 2, -0.2e1 * t72 * t100, 0, 0, 0, t100 * t103, t72 * t103, 0.2e1 * t78, t25 ^ 2 + t27 ^ 2 + t59 ^ 2, t64 * t40, -0.2e1 * t40 * t90, -0.2e1 * t42 * t34, 0.2e1 * t42 * t92, t41, -0.2e1 * t12 * t43 + 0.2e1 * t25 * t93, 0.2e1 * t13 * t43 + 0.2e1 * t25 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t68, 0, t32, 0, 0, 0, 0, 0, t47, -t80, -t42 * t21 + t98, t79, 0, 0, 0, 0, 0, -t42 * t16 + t20 * t92, -t42 * t17 + t20 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, -t108, -t78, 0, 0, 0, 0, 0, -t108 * t69, -t108 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t108, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, t55, t14, -t15 (-t20 * t65 - t21 * t67) * pkin(4) (t5 * t67 + t6 * t65) * pkin(4), t99, -t69 * t16 + t17 * t71, t94, t19, 0, t58 * t16 - t20 * t91 - t3 * t71, t58 * t17 - t20 * t88 + t3 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t100, 0, t72 * t74, -t83, -t109 (-t25 * t67 + t27 * t65) * pkin(4), t69 * t89 (-t63 + t64) * t42, -t92, -t34, 0, -t25 * t71 + t77 * t69, t25 * t69 + t77 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t100, 0, t109, 0, 0, 0, 0, 0, t89, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t65 ^ 2 + t67 ^ 2) * pkin(4) ^ 2, t63, 0.2e1 * t90, 0, 0, 0, -0.2e1 * t58 * t71, 0.2e1 * t58 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, t19, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, 0, 0, -t34, t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t20, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -t93, -t43, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t71, 0, -t91, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
