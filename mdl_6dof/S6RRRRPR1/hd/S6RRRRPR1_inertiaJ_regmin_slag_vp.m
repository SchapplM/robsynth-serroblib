% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t64 = sin(qJ(3));
t65 = sin(qJ(2));
t68 = cos(qJ(2));
t81 = cos(qJ(3));
t40 = t64 * t65 - t81 * t68;
t55 = -t68 * pkin(2) - pkin(1);
t32 = t40 * pkin(3) + t55;
t89 = 0.2e1 * t32;
t88 = 0.2e1 * t55;
t66 = cos(qJ(6));
t87 = -0.2e1 * t66;
t86 = 0.2e1 * t68;
t85 = pkin(7) + pkin(8);
t63 = sin(qJ(4));
t67 = cos(qJ(4));
t45 = t85 * t65;
t46 = t85 * t68;
t28 = t64 * t45 - t81 * t46;
t71 = -t40 * pkin(9) - t28;
t27 = -t81 * t45 - t64 * t46;
t41 = t64 * t68 + t81 * t65;
t72 = -t41 * pkin(9) + t27;
t12 = -t63 * t72 - t67 * t71;
t73 = t67 * t40 + t63 * t41;
t10 = -qJ(5) * t73 - t12;
t60 = sin(pkin(11));
t61 = cos(pkin(11));
t11 = -t63 * t71 + t67 * t72;
t26 = -t63 * t40 + t67 * t41;
t70 = -t26 * qJ(5) + t11;
t4 = t60 * t10 - t61 * t70;
t84 = t4 * t66;
t83 = t63 * pkin(3);
t82 = t64 * pkin(2);
t17 = t60 * t26 + t61 * t73;
t62 = sin(qJ(6));
t14 = t62 * t17;
t18 = t61 * t26 - t60 * t73;
t80 = t62 * t18;
t79 = t62 * t66;
t78 = t66 * t18;
t57 = t81 * pkin(2);
t54 = t57 + pkin(3);
t38 = t67 * t54 - t63 * t82;
t35 = pkin(4) + t38;
t77 = t67 * t82;
t39 = t63 * t54 + t77;
t24 = t60 * t35 + t61 * t39;
t56 = t67 * pkin(3);
t53 = t56 + pkin(4);
t37 = t60 * t53 + t61 * t83;
t23 = t61 * t35 - t60 * t39;
t21 = -pkin(5) - t23;
t22 = pkin(10) + t24;
t76 = -t17 * t22 + t18 * t21;
t36 = t61 * t53 - t60 * t83;
t33 = -pkin(5) - t36;
t34 = pkin(10) + t37;
t75 = -t17 * t34 + t18 * t33;
t50 = t60 * pkin(4) + pkin(10);
t51 = -t61 * pkin(4) - pkin(5);
t74 = -t17 * t50 + t18 * t51;
t19 = pkin(4) * t73 + t32;
t59 = t66 ^ 2;
t58 = t62 ^ 2;
t49 = 0.2e1 * t79;
t44 = t51 * t62;
t29 = t33 * t62;
t20 = t21 * t62;
t16 = t18 ^ 2;
t15 = t66 * t17;
t13 = t62 * t78;
t8 = (-t58 + t59) * t18;
t7 = t17 * pkin(5) - t18 * pkin(10) + t19;
t6 = t61 * t10 + t60 * t70;
t3 = t4 * t62;
t2 = t66 * t6 + t62 * t7;
t1 = -t62 * t6 + t66 * t7;
t5 = [1, 0, 0, t65 ^ 2, t65 * t86, 0, 0, 0, pkin(1) * t86, -0.2e1 * pkin(1) * t65, t41 ^ 2, -0.2e1 * t41 * t40, 0, 0, 0, t40 * t88, t41 * t88, t26 ^ 2, -0.2e1 * t26 * t73, 0, 0, 0, t73 * t89, t26 * t89, -0.2e1 * t6 * t17 + 0.2e1 * t4 * t18, t19 ^ 2 + t4 ^ 2 + t6 ^ 2, t59 * t16, -0.2e1 * t16 * t79, 0.2e1 * t17 * t78, -0.2e1 * t17 * t80, t17 ^ 2, 0.2e1 * t1 * t17 + 0.2e1 * t4 * t80, -0.2e1 * t2 * t17 + 0.2e1 * t4 * t78; 0, 0, 0, 0, 0, t65, t68, 0, -t65 * pkin(7), -t68 * pkin(7), 0, 0, t41, -t40, 0, t27, t28, 0, 0, t26, -t73, 0, t11, t12, -t24 * t17 - t23 * t18, -t4 * t23 + t6 * t24, t13, t8, t14, t15, 0, t62 * t76 - t84, t66 * t76 + t3; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t82, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t39, 0, t23 ^ 2 + t24 ^ 2, t58, t49, 0, 0, 0, t21 * t87, 0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, 0, t27, t28, 0, 0, t26, -t73, 0, t11, t12, -t37 * t17 - t36 * t18, -t4 * t36 + t6 * t37, t13, t8, t14, t15, 0, t62 * t75 - t84, t66 * t75 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t57, -t82, 0, 0, 0, 0, 1, t38 + t56, -t77 + (-pkin(3) - t54) * t63, 0, t23 * t36 + t24 * t37, t58, t49, 0, 0, 0 (-t21 - t33) * t66, t29 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t83, 0, t36 ^ 2 + t37 ^ 2, t58, t49, 0, 0, 0, t33 * t87, 0.2e1 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t73, 0, t11, t12 (-t17 * t60 - t18 * t61) * pkin(4) (-t4 * t61 + t6 * t60) * pkin(4), t13, t8, t14, t15, 0, t62 * t74 - t84, t66 * t74 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, -t39, 0 (t23 * t61 + t24 * t60) * pkin(4), t58, t49, 0, 0, 0 (-t21 - t51) * t66, t44 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t56, -t83, 0 (t36 * t61 + t37 * t60) * pkin(4), t58, t49, 0, 0, 0 (-t33 - t51) * t66, t44 + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t60 ^ 2 + t61 ^ 2) * pkin(4) ^ 2, t58, t49, 0, 0, 0, t51 * t87, 0.2e1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, t15, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t80, t17, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t66, 0, -t62 * t22, -t66 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t66, 0, -t62 * t34, -t66 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t66, 0, -t62 * t50, -t66 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
