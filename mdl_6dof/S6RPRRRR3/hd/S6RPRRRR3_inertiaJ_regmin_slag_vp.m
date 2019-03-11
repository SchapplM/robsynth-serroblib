% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRR3
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t57 = sin(qJ(5));
t58 = sin(qJ(4));
t61 = cos(qJ(5));
t62 = cos(qJ(4));
t34 = t57 * t58 - t61 * t62;
t47 = -t62 * pkin(4) - pkin(3);
t25 = t34 * pkin(5) + t47;
t90 = 0.2e1 * t25;
t89 = 0.2e1 * t47;
t63 = cos(qJ(3));
t88 = -0.2e1 * t63;
t87 = 0.2e1 * t63;
t86 = pkin(8) + pkin(9);
t85 = pkin(3) * t62;
t56 = sin(qJ(6));
t84 = t56 * pkin(5);
t83 = t57 * pkin(4);
t35 = t57 * t62 + t61 * t58;
t59 = sin(qJ(3));
t26 = t35 * t59;
t55 = cos(pkin(11));
t44 = -t55 * pkin(1) - pkin(2);
t33 = -t63 * pkin(3) - t59 * pkin(8) + t44;
t28 = t62 * t33;
t72 = t62 * t59;
t54 = sin(pkin(11));
t43 = t54 * pkin(1) + pkin(7);
t77 = t43 * t58;
t13 = -pkin(9) * t72 + t28 + (-pkin(4) - t77) * t63;
t68 = t63 * t43;
t65 = t62 * t68;
t16 = t65 + (-pkin(9) * t59 + t33) * t58;
t73 = t61 * t16;
t9 = t57 * t13 + t73;
t5 = -t26 * pkin(10) + t9;
t60 = cos(qJ(6));
t82 = t60 * t5;
t81 = t63 * pkin(4);
t80 = t63 * pkin(5);
t19 = -t56 * t34 + t60 * t35;
t79 = t19 * t63;
t78 = t35 * t63;
t76 = t58 * t59;
t75 = t58 * t62;
t74 = t58 * t63;
t71 = t62 * t63;
t18 = t60 * t34 + t56 * t35;
t70 = t63 * t18;
t69 = t63 * t34;
t37 = t59 * t43;
t29 = pkin(4) * t76 + t37;
t67 = t59 * t87;
t66 = t60 * t83;
t27 = -t57 * t76 + t61 * t72;
t8 = t61 * t13 - t57 * t16;
t4 = -t27 * pkin(10) + t8 - t80;
t1 = t60 * t4 - t56 * t5;
t38 = t86 * t58;
t39 = t86 * t62;
t22 = -t61 * t38 - t57 * t39;
t49 = t61 * pkin(4);
t46 = t49 + pkin(5);
t30 = t60 * t46 - t56 * t83;
t2 = t56 * t4 + t82;
t23 = -t57 * t38 + t61 * t39;
t53 = t63 ^ 2;
t52 = t62 ^ 2;
t51 = t59 ^ 2;
t50 = t58 ^ 2;
t48 = t60 * pkin(5);
t31 = t56 * t46 + t66;
t21 = t58 * t33 + t65;
t20 = -t58 * t68 + t28;
t17 = t26 * pkin(5) + t29;
t15 = -t34 * pkin(10) + t23;
t14 = -t35 * pkin(10) + t22;
t12 = -t56 * t26 + t60 * t27;
t11 = t60 * t26 + t56 * t27;
t7 = t56 * t14 + t60 * t15;
t6 = t60 * t14 - t56 * t15;
t3 = [1, 0, 0 (t54 ^ 2 + t55 ^ 2) * pkin(1) ^ 2, t51, t67, 0, 0, 0, t44 * t88, 0.2e1 * t44 * t59, t52 * t51, -0.2e1 * t51 * t75, -0.2e1 * t59 * t71, t58 * t67, t53, -0.2e1 * t20 * t63 + 0.2e1 * t51 * t77, 0.2e1 * t51 * t43 * t62 + 0.2e1 * t21 * t63, t27 ^ 2, -0.2e1 * t27 * t26, t27 * t88, t26 * t87, t53, 0.2e1 * t29 * t26 - 0.2e1 * t8 * t63, 0.2e1 * t29 * t27 + 0.2e1 * t9 * t63, t12 ^ 2, -0.2e1 * t12 * t11, t12 * t88, t11 * t87, t53, -0.2e1 * t1 * t63 + 0.2e1 * t17 * t11, 0.2e1 * t17 * t12 + 0.2e1 * t2 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t59, t63, 0, -t37, -t68, t58 * t72 (-t50 + t52) * t59, -t74, -t71, 0, -t43 * t72 + (-pkin(3) * t59 + pkin(8) * t63) * t58, pkin(8) * t71 + (t77 - t85) * t59, t27 * t35, -t35 * t26 - t27 * t34, -t78, t69, 0, -t22 * t63 + t47 * t26 + t29 * t34, t23 * t63 + t47 * t27 + t29 * t35, t12 * t19, -t19 * t11 - t12 * t18, -t79, t70, 0, t25 * t11 + t17 * t18 - t6 * t63, t25 * t12 + t17 * t19 + t7 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t59, 0, 0, 0, 0, 0, t71, -t74, 0, 0, 0, 0, 0, -t69, -t78, 0, 0, 0, 0, 0, -t70, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t50, 0.2e1 * t75, 0, 0, 0, 0.2e1 * t85, -0.2e1 * pkin(3) * t58, t35 ^ 2, -0.2e1 * t35 * t34, 0, 0, 0, t34 * t89, t35 * t89, t19 ^ 2, -0.2e1 * t19 * t18, 0, 0, 0, t18 * t90, t19 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t76, -t63, t20, -t21, 0, 0, t27, -t26, -t63, -t61 * t81 + t8, -t73 + (-t13 + t81) * t57, 0, 0, t12, -t11, -t63, -t30 * t63 + t1, t31 * t63 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t72, 0, 0, 0, 0, 0, -t26, -t27, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t62, 0, -t58 * pkin(8), -t62 * pkin(8), 0, 0, t35, -t34, 0, t22, -t23, 0, 0, t19, -t18, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t49, -0.2e1 * t83, 0, 0, 0, 0, 1, 0.2e1 * t30, -0.2e1 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, -t63, t8, -t9, 0, 0, t12, -t11, -t63, -t60 * t80 + t1, -t82 + (-t4 + t80) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34, 0, t22, -t23, 0, 0, t19, -t18, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t49, -t83, 0, 0, 0, 0, 1, t30 + t48, -t66 + (-pkin(5) - t46) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t48, -0.2e1 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, -t63, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t30, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t48, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
