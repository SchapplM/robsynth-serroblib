% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t64 = sin(qJ(4));
t65 = sin(qJ(3));
t82 = cos(qJ(4));
t83 = cos(qJ(3));
t42 = t64 * t65 - t82 * t83;
t97 = -0.2e1 * t42;
t96 = 0.2e1 * t42;
t61 = cos(pkin(11));
t51 = -t61 * pkin(1) - pkin(2);
t45 = -t83 * pkin(3) + t51;
t95 = 0.2e1 * t45;
t72 = t82 * pkin(3);
t55 = -t72 - pkin(4);
t67 = cos(qJ(5));
t86 = t67 * pkin(5);
t46 = t55 - t86;
t94 = 0.2e1 * t46;
t56 = -pkin(4) - t86;
t93 = 0.2e1 * t56;
t92 = 0.2e1 * t65;
t91 = t42 * pkin(5);
t62 = sin(qJ(6));
t90 = t62 * pkin(5);
t89 = t64 * pkin(3);
t66 = cos(qJ(6));
t88 = t66 * pkin(5);
t44 = t64 * t83 + t82 * t65;
t15 = t42 * pkin(4) - t44 * pkin(9) + t45;
t63 = sin(qJ(5));
t60 = sin(pkin(11));
t50 = t60 * pkin(1) + pkin(7);
t34 = (-pkin(8) - t50) * t65;
t71 = t83 * t50;
t35 = t83 * pkin(8) + t71;
t19 = t64 * t34 + t82 * t35;
t77 = t67 * t19;
t6 = t77 + (-pkin(10) * t44 + t15) * t63;
t87 = t66 * t6;
t85 = t67 * pkin(9);
t84 = pkin(4) - t55;
t18 = -t82 * t34 + t64 * t35;
t81 = t18 * t67;
t41 = t62 * t63 - t66 * t67;
t80 = t42 * t41;
t79 = t63 * t44;
t78 = t63 * t67;
t76 = t67 * t44;
t54 = pkin(9) + t89;
t75 = t67 * t54;
t74 = t46 + t56;
t73 = t44 * t97;
t7 = t67 * t15 - t63 * t19;
t4 = -pkin(10) * t76 + t7 + t91;
t1 = t66 * t4 - t62 * t6;
t70 = -pkin(4) * t44 - pkin(9) * t42;
t69 = -t42 * t54 + t44 * t55;
t43 = t62 * t67 + t66 * t63;
t59 = t67 ^ 2;
t58 = t63 ^ 2;
t57 = t67 * pkin(10);
t49 = 0.2e1 * t78;
t48 = t57 + t85;
t47 = (-pkin(9) - pkin(10)) * t63;
t40 = t44 ^ 2;
t39 = t43 ^ 2;
t38 = t42 ^ 2;
t37 = t57 + t75;
t36 = (-pkin(10) - t54) * t63;
t33 = t67 * t42;
t32 = t63 * t42;
t30 = t63 * t76;
t28 = t62 * t47 + t66 * t48;
t27 = t66 * t47 - t62 * t48;
t26 = t43 * t42;
t23 = -0.2e1 * t43 * t41;
t22 = (-t58 + t59) * t44;
t21 = t62 * t36 + t66 * t37;
t20 = t66 * t36 - t62 * t37;
t17 = -t62 * t79 + t66 * t76;
t16 = t43 * t44;
t14 = t18 * t63;
t12 = t17 * t43;
t11 = pkin(5) * t79 + t18;
t10 = t11 * t43;
t9 = t11 * t41;
t8 = t63 * t15 + t77;
t5 = -t43 * t16 - t17 * t41;
t2 = t62 * t4 + t87;
t3 = [1, 0, 0 (t60 ^ 2 + t61 ^ 2) * pkin(1) ^ 2, t65 ^ 2, t83 * t92, 0, 0, 0, -0.2e1 * t51 * t83, t51 * t92, t40, t73, 0, 0, 0, t42 * t95, t44 * t95, t59 * t40, -0.2e1 * t40 * t78, t76 * t96, t63 * t73, t38, 0.2e1 * t18 * t79 + 0.2e1 * t7 * t42, 0.2e1 * t18 * t76 - 0.2e1 * t8 * t42, t17 ^ 2, -0.2e1 * t17 * t16, t17 * t96, t16 * t97, t38, 0.2e1 * t1 * t42 + 0.2e1 * t11 * t16, 0.2e1 * t11 * t17 - 0.2e1 * t2 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t65, t83, 0, -t65 * t50, -t71, 0, 0, t44, -t42, 0, -t18, -t19, t30, t22, t32, t33, 0, t69 * t63 - t81, t69 * t67 + t14, t12, t5, t26, -t80, 0, t46 * t16 + t20 * t42 + t9, t46 * t17 - t21 * t42 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t65, 0, 0, 0, 0, 0, -t42, -t44, 0, 0, 0, 0, 0, -t33, t32, 0, 0, 0, 0, 0, t80, t26; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t72, -0.2e1 * t89, t58, t49, 0, 0, 0, -0.2e1 * t55 * t67, 0.2e1 * t55 * t63, t39, t23, 0, 0, 0, t41 * t94, t43 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t42, 0, -t18, -t19, t30, t22, t32, t33, 0, t70 * t63 - t81, t70 * t67 + t14, t12, t5, t26, -t80, 0, t56 * t16 + t27 * t42 + t9, t56 * t17 - t28 * t42 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t44, 0, 0, 0, 0, 0, -t33, t32, 0, 0, 0, 0, 0, t80, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t72, -t89, t58, t49, 0, 0, 0, t84 * t67, -t84 * t63, t39, t23, 0, 0, 0, t74 * t41, t74 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t58, t49, 0, 0, 0, 0.2e1 * pkin(4) * t67, -0.2e1 * pkin(4) * t63, t39, t23, 0, 0, 0, t41 * t93, t43 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t79, t42, t7, -t8, 0, 0, t17, -t16, t42, t42 * t88 + t1, -t87 + (-t4 - t91) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t76, 0, 0, 0, 0, 0, -t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t67, 0, -t63 * t54, -t75, 0, 0, t43, -t41, 0, t20, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t67, 0, -t63 * pkin(9), -t85, 0, 0, t43, -t41, 0, t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t88, -0.2e1 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t42, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t41, 0, t20, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t41, 0, t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t88, -t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
