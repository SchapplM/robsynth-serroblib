% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t52 = cos(pkin(7));
t55 = sin(qJ(4));
t59 = cos(qJ(4));
t50 = sin(pkin(7));
t56 = sin(qJ(3));
t83 = t50 * t56;
t31 = t55 * t52 + t59 * t83;
t54 = sin(qJ(5));
t58 = cos(qJ(5));
t60 = cos(qJ(3));
t82 = t50 * t60;
t19 = t54 * t31 + t58 * t82;
t96 = -0.2e1 * t19;
t95 = -0.2e1 * t31;
t94 = -0.2e1 * t55;
t93 = 0.2e1 * t59;
t92 = pkin(2) * t56;
t91 = pkin(2) * t60;
t90 = pkin(4) * t58;
t89 = pkin(10) * t54;
t88 = t54 * pkin(5);
t64 = pkin(9) * t82;
t27 = t64 + (pkin(10) + t92) * t52;
t28 = (-pkin(3) * t60 - pkin(10) * t56 - pkin(2)) * t50;
t14 = -t55 * t27 + t59 * t28;
t12 = pkin(4) * t82 - t14;
t87 = t12 * t54;
t86 = t12 * t58;
t20 = t58 * t31 - t54 * t82;
t85 = t20 * t54;
t46 = t50 ^ 2;
t84 = t46 * t60;
t51 = sin(pkin(6));
t57 = sin(qJ(2));
t81 = t51 * t57;
t61 = cos(qJ(2));
t80 = t51 * t61;
t79 = t52 * t56;
t78 = t52 * t61;
t30 = -t59 * t52 + t55 * t83;
t77 = t54 * t30;
t76 = t54 * t55;
t75 = t54 * t58;
t74 = t54 * t59;
t73 = t55 * t30;
t72 = t58 * t30;
t71 = t58 * t55;
t70 = t58 * t59;
t69 = -qJ(6) - pkin(11);
t68 = qJ(6) * t55;
t67 = 0.2e1 * t82;
t66 = t55 * t93;
t65 = pkin(10) * t70;
t63 = t55 * t82;
t62 = t59 * t82;
t42 = pkin(9) * t83;
t26 = t42 + (-pkin(3) - t91) * t52;
t11 = t30 * pkin(4) - t31 * pkin(11) + t26;
t15 = t59 * t27 + t55 * t28;
t13 = -pkin(11) * t82 + t15;
t3 = t58 * t11 - t54 * t13;
t4 = t54 * t11 + t58 * t13;
t53 = cos(pkin(6));
t49 = t58 ^ 2;
t48 = t55 ^ 2;
t47 = t54 ^ 2;
t45 = -t58 * pkin(5) - pkin(4);
t39 = t69 * t58;
t38 = t69 * t54;
t37 = -t59 * pkin(4) - t55 * pkin(11) - pkin(3);
t36 = (pkin(10) + t88) * t55;
t34 = t58 * t37;
t33 = pkin(2) * t79 + t64;
t32 = t52 * t91 - t42;
t29 = -t50 * t80 + t53 * t52;
t24 = t54 * t37 + t65;
t23 = -pkin(10) * t74 + t34;
t21 = t65 + (t37 - t68) * t54;
t18 = t53 * t83 + (t56 * t78 + t57 * t60) * t51;
t17 = -t51 * t60 * t78 - t53 * t82 + t56 * t81;
t16 = -t58 * t68 + t34 + (-pkin(5) - t89) * t59;
t10 = t18 * t59 + t29 * t55;
t9 = t18 * t55 - t29 * t59;
t7 = t19 * pkin(5) + t12;
t6 = t10 * t58 + t17 * t54;
t5 = -t10 * t54 + t17 * t58;
t2 = -t19 * qJ(6) + t4;
t1 = t30 * pkin(5) - t20 * qJ(6) + t3;
t8 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t9 ^ 2; 0, 0, t80, -t81, 0, 0, 0, 0, 0, -t17 * t52 - t29 * t82, -t18 * t52 + t29 * t83, 0, 0, 0, 0, 0, t17 * t30 + t9 * t82, t10 * t82 + t17 * t31, 0, 0, 0, 0, 0, t9 * t19 + t5 * t30, t9 * t20 - t6 * t30, -t6 * t19 - t5 * t20, t5 * t1 + t6 * t2 + t9 * t7; 0, 1, 0, 0, t46 * t56 ^ 2, 0.2e1 * t56 * t84, 0.2e1 * t50 * t79, t52 * t67, t52 ^ 2, 0.2e1 * pkin(2) * t84 + 0.2e1 * t32 * t52, -0.2e1 * t33 * t52 - 0.2e1 * t46 * t92, t31 ^ 2, t30 * t95, t82 * t95, t30 * t67, t46 * t60 ^ 2, -0.2e1 * t14 * t82 + 0.2e1 * t26 * t30, 0.2e1 * t15 * t82 + 0.2e1 * t26 * t31, t20 ^ 2, t20 * t96, 0.2e1 * t20 * t30, t30 * t96, t30 ^ 2, 0.2e1 * t12 * t19 + 0.2e1 * t3 * t30, 0.2e1 * t12 * t20 - 0.2e1 * t4 * t30, -0.2e1 * t1 * t20 - 0.2e1 * t2 * t19, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, 0, 0, 0, 0, -t17 * t59, t17 * t55, 0, 0, 0, 0, 0, -t5 * t59 + t9 * t76, t6 * t59 + t9 * t71 (-t5 * t58 - t54 * t6) * t55, t5 * t16 + t6 * t21 + t9 * t36; 0, 0, 0, 0, 0, 0, t83, t82, t52, t32, -t33, t31 * t55, t31 * t59 - t73, -t63, -t62, 0, -pkin(3) * t30 + pkin(10) * t63 - t26 * t59, -pkin(3) * t31 + pkin(10) * t62 + t26 * t55, t20 * t71 (-t19 * t58 - t85) * t55, -t20 * t59 + t30 * t71, t19 * t59 - t54 * t73, -t30 * t59, t23 * t30 - t3 * t59 + (pkin(10) * t19 + t87) * t55, -t24 * t30 + t4 * t59 + (pkin(10) * t20 + t86) * t55, -t16 * t20 - t21 * t19 + (-t1 * t58 - t2 * t54) * t55, t1 * t16 + t2 * t21 + t7 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t48, t66, 0, 0, 0, pkin(3) * t93, pkin(3) * t94, t49 * t48, -0.2e1 * t48 * t75, t70 * t94, t54 * t66, t59 ^ 2, -0.2e1 * t23 * t59 + 0.2e1 * t48 * t89, 0.2e1 * t48 * pkin(10) * t58 + 0.2e1 * t24 * t59, 0.2e1 * (-t16 * t58 - t21 * t54) * t55, t16 ^ 2 + t21 ^ 2 + t36 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, 0, 0, 0, 0, -t9 * t58, t9 * t54, -t5 * t54 + t6 * t58, t5 * t38 - t6 * t39 + t9 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, -t82, t14, -t15, t85, -t54 * t19 + t20 * t58, t77, t72, 0, -pkin(4) * t19 - pkin(11) * t77 - t86, -pkin(4) * t20 - pkin(11) * t72 + t87, -t1 * t54 + t39 * t19 + t2 * t58 - t38 * t20, t1 * t38 - t2 * t39 + t7 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t59, 0, -t55 * pkin(10), -t59 * pkin(10), t54 * t71 (-t47 + t49) * t55, -t74, -t70, 0, -pkin(10) * t71 + (-pkin(4) * t55 + pkin(11) * t59) * t54, pkin(11) * t70 + (t89 - t90) * t55 (-t38 * t55 + t21) * t58 + (t39 * t55 - t16) * t54, t16 * t38 - t21 * t39 + t36 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t47, 0.2e1 * t75, 0, 0, 0, 0.2e1 * t90, -0.2e1 * pkin(4) * t54, -0.2e1 * t38 * t54 - 0.2e1 * t39 * t58, t38 ^ 2 + t39 ^ 2 + t45 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, t5 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t30, t3, -t4, -pkin(5) * t20, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t76, -t59, t23, -t24, -pkin(5) * t71, t16 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t58, 0, -t54 * pkin(11), -t58 * pkin(11), -t88, t38 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t8;
