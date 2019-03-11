% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t60 = sin(pkin(12));
t62 = cos(pkin(12));
t65 = sin(qJ(4));
t69 = cos(qJ(4));
t42 = t60 * t69 + t62 * t65;
t66 = sin(qJ(3));
t34 = t42 * t66;
t78 = t69 * t66;
t81 = t65 * t66;
t35 = -t60 * t81 + t62 * t78;
t64 = sin(qJ(6));
t68 = cos(qJ(6));
t15 = -t64 * t34 + t68 * t35;
t90 = -0.2e1 * t15;
t41 = -t60 * t65 + t62 * t69;
t54 = -t69 * pkin(4) - pkin(3);
t33 = -t41 * pkin(5) + t54;
t89 = 0.2e1 * t33;
t88 = -0.2e1 * t66;
t70 = cos(qJ(3));
t87 = 0.2e1 * t70;
t86 = pkin(3) * t69;
t85 = pkin(4) * t60;
t84 = pkin(8) * t65;
t61 = sin(pkin(6));
t83 = t61 * sin(qJ(2));
t82 = t61 * cos(qJ(2));
t80 = t65 * t69;
t79 = t65 * t70;
t77 = t69 * t70;
t76 = -qJ(5) - pkin(9);
t48 = -t70 * pkin(3) - t66 * pkin(9) - pkin(2);
t43 = t69 * t48;
t75 = qJ(5) * t66;
t22 = -t69 * t75 + t43 + (-pkin(4) - t84) * t70;
t73 = pkin(8) * t77;
t27 = t73 + (t48 - t75) * t65;
t11 = t60 * t22 + t62 * t27;
t49 = t76 * t65;
t50 = t76 * t69;
t29 = t60 * t49 - t62 * t50;
t55 = t66 * pkin(8);
t47 = pkin(4) * t81 + t55;
t74 = t66 * t87;
t10 = t62 * t22 - t60 * t27;
t6 = -t70 * pkin(5) - t35 * pkin(10) + t10;
t9 = -t34 * pkin(10) + t11;
t1 = t68 * t6 - t64 * t9;
t28 = t62 * t49 + t60 * t50;
t2 = t64 * t6 + t68 * t9;
t63 = cos(pkin(6));
t59 = t70 ^ 2;
t58 = t69 ^ 2;
t57 = t66 ^ 2;
t56 = t65 ^ 2;
t53 = t62 * pkin(4) + pkin(5);
t40 = t63 * t66 + t70 * t83;
t39 = -t63 * t70 + t66 * t83;
t37 = t64 * t53 + t68 * t85;
t36 = t68 * t53 - t64 * t85;
t32 = t65 * t48 + t73;
t31 = -pkin(8) * t79 + t43;
t26 = t40 * t69 - t65 * t82;
t25 = -t40 * t65 - t69 * t82;
t23 = t34 * pkin(5) + t47;
t21 = t64 * t41 + t68 * t42;
t20 = -t68 * t41 + t64 * t42;
t17 = t41 * pkin(10) + t29;
t16 = -t42 * pkin(10) + t28;
t14 = t68 * t34 + t64 * t35;
t13 = t60 * t25 + t62 * t26;
t12 = t62 * t25 - t60 * t26;
t8 = t64 * t16 + t68 * t17;
t7 = t68 * t16 - t64 * t17;
t4 = t64 * t12 + t68 * t13;
t3 = t68 * t12 - t64 * t13;
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 ^ 2 + t13 ^ 2 + t39 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t82, -t83, 0, 0, 0, 0, 0, t70 * t82, -t66 * t82, 0, 0, 0, 0, 0, -t25 * t70 + t39 * t81, t26 * t70 + t39 * t78, -t12 * t35 - t13 * t34, t12 * t10 + t13 * t11 + t39 * t47, 0, 0, 0, 0, 0, t39 * t14 - t3 * t70, t39 * t15 + t4 * t70; 0, 1, 0, 0, t57, t74, 0, 0, 0, pkin(2) * t87, pkin(2) * t88, t58 * t57, -0.2e1 * t57 * t80, t77 * t88, t65 * t74, t59, -0.2e1 * t31 * t70 + 0.2e1 * t57 * t84, 0.2e1 * t57 * pkin(8) * t69 + 0.2e1 * t32 * t70, -0.2e1 * t10 * t35 - 0.2e1 * t11 * t34, t10 ^ 2 + t11 ^ 2 + t47 ^ 2, t15 ^ 2, t14 * t90, t70 * t90, t14 * t87, t59, -0.2e1 * t1 * t70 + 0.2e1 * t23 * t14, 0.2e1 * t23 * t15 + 0.2e1 * t2 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t40, 0, 0, 0, 0, 0, -t39 * t69, t39 * t65, -t12 * t42 + t13 * t41, t12 * t28 + t13 * t29 + t39 * t54, 0, 0, 0, 0, 0, t39 * t20, t39 * t21; 0, 0, 0, 0, 0, 0, t66, t70, 0, -t55, -t70 * pkin(8), t65 * t78 (-t56 + t58) * t66, -t79, -t77, 0, -pkin(8) * t78 + (-pkin(3) * t66 + pkin(9) * t70) * t65, pkin(9) * t77 + (t84 - t86) * t66, -t10 * t42 + t11 * t41 - t28 * t35 - t29 * t34, t10 * t28 + t11 * t29 + t47 * t54, t15 * t21, -t21 * t14 - t15 * t20, -t21 * t70, t20 * t70, 0, t33 * t14 + t23 * t20 - t7 * t70, t33 * t15 + t23 * t21 + t8 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t56, 0.2e1 * t80, 0, 0, 0, 0.2e1 * t86, -0.2e1 * pkin(3) * t65, -0.2e1 * t28 * t42 + 0.2e1 * t29 * t41, t28 ^ 2 + t29 ^ 2 + t54 ^ 2, t21 ^ 2, -0.2e1 * t21 * t20, 0, 0, 0, t20 * t89, t21 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t26, 0 (t12 * t62 + t13 * t60) * pkin(4), 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t81, -t70, t31, -t32 (-t34 * t60 - t35 * t62) * pkin(4) (t10 * t62 + t11 * t60) * pkin(4), 0, 0, t15, -t14, -t70, -t36 * t70 + t1, t37 * t70 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t69, 0, -t65 * pkin(9), -t69 * pkin(9) (t41 * t60 - t42 * t62) * pkin(4) (t28 * t62 + t29 * t60) * pkin(4), 0, 0, t21, -t20, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t60 ^ 2 + t62 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t36, -0.2e1 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, t20, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, -t70, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t36, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
