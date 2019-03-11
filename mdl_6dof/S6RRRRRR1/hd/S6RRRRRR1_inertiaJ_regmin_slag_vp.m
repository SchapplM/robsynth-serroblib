% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x38]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t64 = sin(qJ(3));
t65 = sin(qJ(2));
t69 = cos(qJ(2));
t89 = cos(qJ(3));
t40 = t64 * t65 - t89 * t69;
t54 = -t69 * pkin(2) - pkin(1);
t32 = t40 * pkin(3) + t54;
t41 = t64 * t69 + t89 * t65;
t63 = sin(qJ(4));
t68 = cos(qJ(4));
t71 = t68 * t40 + t63 * t41;
t18 = t71 * pkin(4) + t32;
t99 = 0.2e1 * t18;
t98 = 0.2e1 * t32;
t97 = 0.2e1 * t54;
t96 = 0.2e1 * t69;
t95 = pkin(7) + pkin(8);
t61 = sin(qJ(6));
t94 = pkin(5) * t61;
t62 = sin(qJ(5));
t67 = cos(qJ(5));
t43 = t95 * t65;
t44 = t95 * t69;
t28 = -t89 * t43 - t64 * t44;
t24 = -t41 * pkin(9) + t28;
t29 = t64 * t43 - t89 * t44;
t25 = -t40 * pkin(9) - t29;
t10 = t68 * t24 - t63 * t25;
t27 = -t63 * t40 + t68 * t41;
t70 = -t27 * pkin(10) + t10;
t11 = -t63 * t24 - t68 * t25;
t9 = -t71 * pkin(10) - t11;
t4 = t62 * t9 - t67 * t70;
t66 = cos(qJ(6));
t93 = t4 * t66;
t92 = t62 * pkin(4);
t91 = t63 * pkin(3);
t90 = t64 * pkin(2);
t57 = t89 * pkin(2);
t53 = t57 + pkin(3);
t37 = t68 * t53 - t63 * t90;
t34 = pkin(4) + t37;
t30 = t67 * t34;
t79 = t68 * t90;
t39 = t63 * t53 + t79;
t77 = t62 * t39 - t30;
t20 = -pkin(5) + t77;
t88 = t20 * t66;
t56 = t68 * pkin(3);
t52 = t56 + pkin(4);
t45 = t67 * t52;
t76 = t62 * t91 - t45;
t33 = -pkin(5) + t76;
t87 = t33 * t66;
t55 = t67 * pkin(4);
t51 = -t55 - pkin(5);
t86 = t51 * t66;
t17 = t67 * t27 - t62 * t71;
t85 = t61 * t17;
t84 = t61 * t66;
t83 = t66 * t17;
t82 = t67 * t39;
t16 = t62 * t27 + t67 * t71;
t81 = -0.2e1 * t17 * t16;
t80 = t67 * t91;
t78 = -t39 - t91;
t75 = -pkin(5) * t17 - pkin(11) * t16;
t23 = -t62 * t34 - t82;
t21 = pkin(11) - t23;
t74 = -t16 * t21 + t17 * t20;
t38 = -t62 * t52 - t80;
t35 = pkin(11) - t38;
t73 = -t16 * t35 + t17 * t33;
t50 = pkin(11) + t92;
t72 = -t16 * t50 + t17 * t51;
t60 = t66 ^ 2;
t59 = t61 ^ 2;
t58 = pkin(5) * t66;
t48 = 0.2e1 * t84;
t47 = t51 * t61;
t31 = t33 * t61;
t19 = t20 * t61;
t15 = t17 ^ 2;
t14 = t66 * t16;
t13 = t61 * t16;
t12 = t61 * t83;
t7 = (-t59 + t60) * t17;
t6 = t16 * pkin(5) - t17 * pkin(11) + t18;
t5 = t62 * t70 + t67 * t9;
t3 = t4 * t61;
t2 = t66 * t5 + t61 * t6;
t1 = -t61 * t5 + t66 * t6;
t8 = [1, 0, 0, t65 ^ 2, t65 * t96, 0, 0, 0, pkin(1) * t96, -0.2e1 * pkin(1) * t65, t41 ^ 2, -0.2e1 * t41 * t40, 0, 0, 0, t40 * t97, t41 * t97, t27 ^ 2, -0.2e1 * t27 * t71, 0, 0, 0, t71 * t98, t27 * t98, t15, t81, 0, 0, 0, t16 * t99, t17 * t99, t60 * t15, -0.2e1 * t15 * t84, 0.2e1 * t16 * t83, t61 * t81, t16 ^ 2, 0.2e1 * t1 * t16 + 0.2e1 * t4 * t85, -0.2e1 * t2 * t16 + 0.2e1 * t4 * t83; 0, 0, 0, 0, 0, t65, t69, 0, -t65 * pkin(7), -t69 * pkin(7), 0, 0, t41, -t40, 0, t28, t29, 0, 0, t27, -t71, 0, t10, t11, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t74 * t61 - t93, t74 * t66 + t3; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t90, 0, 0, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t39, 0, 0, 0, 0, 1, -0.2e1 * t77, 0.2e1 * t23, t59, t48, 0, 0, 0, -0.2e1 * t88, 0.2e1 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, 0, t28, t29, 0, 0, t27, -t71, 0, t10, t11, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t73 * t61 - t93, t73 * t66 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t57, -t90, 0, 0, 0, 0, 1, t37 + t56, -t79 + (-pkin(3) - t53) * t63, 0, 0, 0, 0, 1, t78 * t62 + t30 + t45, t78 * t67 + (-t34 - t52) * t62, t59, t48, 0, 0, 0 (-t20 - t33) * t66, t31 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t91, 0, 0, 0, 0, 1, -0.2e1 * t76, 0.2e1 * t38, t59, t48, 0, 0, 0, -0.2e1 * t87, 0.2e1 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t71, 0, t10, t11, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t72 * t61 - t93, t72 * t66 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t37, -t39, 0, 0, 0, 0, 1, t55 - t77, -t82 + (-pkin(4) - t34) * t62, t59, t48, 0, 0, 0 (-t20 - t51) * t66, t47 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t56, -t91, 0, 0, 0, 0, 1, t55 - t76, -t80 + (-pkin(4) - t52) * t62, t59, t48, 0, 0, 0 (-t33 - t51) * t66, t47 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t55, -0.2e1 * t92, t59, t48, 0, 0, 0, -0.2e1 * t86, 0.2e1 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, -t4, -t5, t12, t7, t13, t14, 0, t75 * t61 - t93, t75 * t66 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t77, t23, t59, t48, 0, 0, 0, t58 - t88, t19 - t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t76, t38, t59, t48, 0, 0, 0, t58 - t87, t31 - t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t55, -t92, t59, t48, 0, 0, 0, t58 - t86, t47 - t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t59, t48, 0, 0, 0, 0.2e1 * t58, -0.2e1 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t85, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t66, 0, -t61 * t21, -t66 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t66, 0, -t61 * t35, -t66 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t66, 0, -t61 * t50, -t66 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t66, 0, -t61 * pkin(11), -t66 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
