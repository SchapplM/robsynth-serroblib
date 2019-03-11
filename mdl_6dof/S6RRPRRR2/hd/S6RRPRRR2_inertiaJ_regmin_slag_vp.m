% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t69 = sin(pkin(11));
t70 = cos(pkin(11));
t74 = sin(qJ(2));
t77 = cos(qJ(2));
t49 = t69 * t77 + t70 * t74;
t73 = sin(qJ(4));
t79 = t69 * t74 - t70 * t77;
t91 = cos(qJ(4));
t33 = t73 * t49 + t91 * t79;
t105 = -0.2e1 * t33;
t104 = 0.2e1 * t33;
t62 = t70 * pkin(2) + pkin(3);
t98 = t69 * pkin(2);
t46 = t91 * t62 - t73 * t98;
t44 = -pkin(4) - t46;
t76 = cos(qJ(5));
t94 = t76 * pkin(5);
t42 = t44 - t94;
t103 = 0.2e1 * t42;
t65 = -t77 * pkin(2) - pkin(1);
t43 = t79 * pkin(3) + t65;
t102 = 0.2e1 * t43;
t64 = -pkin(4) - t94;
t101 = 0.2e1 * t64;
t100 = 0.2e1 * t77;
t99 = t33 * pkin(5);
t71 = sin(qJ(6));
t97 = t71 * pkin(5);
t75 = cos(qJ(6));
t96 = t75 * pkin(5);
t34 = t91 * t49 - t73 * t79;
t15 = t33 * pkin(4) - t34 * pkin(9) + t43;
t72 = sin(qJ(5));
t84 = -qJ(3) - pkin(7);
t56 = t84 * t74;
t57 = t84 * t77;
t35 = t70 * t56 + t69 * t57;
t26 = -t49 * pkin(8) + t35;
t36 = t69 * t56 - t70 * t57;
t27 = -t79 * pkin(8) + t36;
t17 = t73 * t26 + t91 * t27;
t87 = t76 * t17;
t5 = t87 + (-pkin(10) * t34 + t15) * t72;
t95 = t75 * t5;
t93 = t76 * pkin(9);
t92 = pkin(4) - t44;
t16 = -t91 * t26 + t73 * t27;
t90 = t16 * t76;
t55 = t71 * t76 + t75 * t72;
t24 = t55 * t33;
t29 = t72 * t33;
t89 = t72 * t34;
t88 = t72 * t76;
t86 = t76 * t34;
t47 = -t73 * t62 - t91 * t98;
t45 = pkin(9) - t47;
t85 = t76 * t45;
t83 = t42 + t64;
t82 = t34 * t105;
t6 = t76 * t15 - t72 * t17;
t4 = -pkin(10) * t86 + t6 + t99;
t1 = t75 * t4 - t71 * t5;
t81 = -pkin(4) * t34 - pkin(9) * t33;
t80 = -t33 * t45 + t34 * t44;
t54 = t71 * t72 - t75 * t76;
t68 = t76 ^ 2;
t67 = t72 ^ 2;
t66 = t76 * pkin(10);
t61 = 0.2e1 * t88;
t59 = t66 + t93;
t58 = (-pkin(9) - pkin(10)) * t72;
t50 = t55 ^ 2;
t41 = t71 * t58 + t75 * t59;
t40 = t75 * t58 - t71 * t59;
t39 = t66 + t85;
t38 = (-pkin(10) - t45) * t72;
t37 = -0.2e1 * t55 * t54;
t32 = t34 ^ 2;
t31 = t33 ^ 2;
t30 = t76 * t33;
t28 = t72 * t86;
t23 = t54 * t33;
t22 = t71 * t38 + t75 * t39;
t21 = t75 * t38 - t71 * t39;
t20 = (-t67 + t68) * t34;
t19 = t54 * t34;
t18 = t55 * t34;
t14 = t19 * t55;
t13 = t16 * t72;
t11 = pkin(5) * t89 + t16;
t10 = t11 * t55;
t9 = t11 * t54;
t8 = -t55 * t18 + t19 * t54;
t7 = t72 * t15 + t87;
t2 = t71 * t4 + t95;
t3 = [1, 0, 0, t74 ^ 2, t74 * t100, 0, 0, 0, pkin(1) * t100, -0.2e1 * pkin(1) * t74, -0.2e1 * t35 * t49 - 0.2e1 * t36 * t79, t35 ^ 2 + t36 ^ 2 + t65 ^ 2, t32, t82, 0, 0, 0, t33 * t102, t34 * t102, t68 * t32, -0.2e1 * t32 * t88, t86 * t104, t72 * t82, t31, 0.2e1 * t16 * t89 + 0.2e1 * t6 * t33, 0.2e1 * t16 * t86 - 0.2e1 * t7 * t33, t19 ^ 2, 0.2e1 * t19 * t18, -t19 * t104, t18 * t105, t31, 0.2e1 * t1 * t33 + 0.2e1 * t11 * t18, -0.2e1 * t11 * t19 - 0.2e1 * t2 * t33; 0, 0, 0, 0, 0, t74, t77, 0, -t74 * pkin(7), -t77 * pkin(7) (-t70 * t49 - t69 * t79) * pkin(2) (t35 * t70 + t36 * t69) * pkin(2), 0, 0, t34, -t33, 0, -t16, -t17, t28, t20, t29, t30, 0, t80 * t72 - t90, t80 * t76 + t13, -t14, t8, t24, -t23, 0, t42 * t18 + t21 * t33 + t9, -t42 * t19 - t22 * t33 + t10; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t69 ^ 2 + t70 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t46, 0.2e1 * t47, t67, t61, 0, 0, 0, -0.2e1 * t44 * t76, 0.2e1 * t44 * t72, t50, t37, 0, 0, 0, t54 * t103, t55 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, t33, t34, 0, 0, 0, 0, 0, t30, -t29, 0, 0, 0, 0, 0, -t23, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t33, 0, -t16, -t17, t28, t20, t29, t30, 0, t81 * t72 - t90, t81 * t76 + t13, -t14, t8, t24, -t23, 0, t64 * t18 + t40 * t33 + t9, -t64 * t19 - t41 * t33 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t46, t47, t67, t61, 0, 0, 0, t92 * t76, -t92 * t72, t50, t37, 0, 0, 0, t83 * t54, t83 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t67, t61, 0, 0, 0, 0.2e1 * pkin(4) * t76, -0.2e1 * pkin(4) * t72, t50, t37, 0, 0, 0, t54 * t101, t55 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t89, t33, t6, -t7, 0, 0, -t19, -t18, t33, t33 * t96 + t1, -t95 + (-t4 - t99) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t76, 0, -t72 * t45, -t85, 0, 0, t55, -t54, 0, t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t72, 0, 0, 0, 0, 0, -t54, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t76, 0, -t72 * pkin(9), -t93, 0, 0, t55, -t54, 0, t40, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t96, -0.2e1 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t18, t33, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t54, 0, t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t54, 0, t40, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t96, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
