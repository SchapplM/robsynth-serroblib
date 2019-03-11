% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t48 = sin(pkin(6));
t55 = sin(qJ(2));
t35 = t48 * t55;
t49 = cos(pkin(6));
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t21 = -t57 * t35 + t49 * t54;
t90 = -0.2e1 * t21;
t89 = 0.2e1 * t48;
t51 = (pkin(2) + qJ(4));
t88 = 2 * t51;
t56 = cos(qJ(6));
t87 = 0.2e1 * t56;
t58 = cos(qJ(2));
t36 = t48 * t58;
t50 = qJ(3) - pkin(9);
t61 = -t51 * t58 - pkin(1);
t10 = (-t50 * t55 + t61) * t48;
t81 = t49 * t55;
t24 = pkin(1) * t81 + pkin(8) * t36;
t41 = t49 * qJ(3);
t18 = -t41 - t24;
t31 = pkin(3) * t36;
t15 = t31 - t18;
t9 = pkin(4) * t36 - t49 * pkin(9) + t15;
t5 = -t54 * t10 + t57 * t9;
t3 = -pkin(5) * t36 - t5;
t53 = sin(qJ(6));
t86 = t3 * t53;
t85 = t3 * t56;
t22 = t54 * t35 + t49 * t57;
t13 = t22 * t56 + t53 * t36;
t84 = t13 * t53;
t43 = t48 ^ 2;
t83 = t43 * t58;
t47 = t57 ^ 2;
t82 = t47 * t50;
t80 = t53 * t21;
t37 = t53 * t54;
t79 = t53 * t56;
t78 = t53 * t57;
t77 = t54 * t50;
t76 = t56 * t21;
t75 = t56 * t54;
t38 = t56 * t57;
t74 = t57 * t13;
t73 = t57 * t21;
t72 = t57 * t50;
t71 = t57 * t54;
t70 = -t49 * t58 * pkin(1) + pkin(8) * t35;
t45 = t54 ^ 2;
t69 = -t45 - t47;
t68 = qJ(3) * t55;
t67 = 0.2e1 * t36;
t66 = -0.2e1 * t71;
t65 = t54 * t36;
t28 = t57 * t36;
t42 = t49 * pkin(2);
t20 = -t42 + t70;
t64 = 0.2e1 * t41 + t24;
t63 = t49 * qJ(4) - t20;
t62 = -pkin(5) * t57 - pkin(10) * t54;
t6 = t57 * t10 + t54 * t9;
t11 = pkin(3) * t35 - t63;
t8 = (-pkin(3) - pkin(4)) * t35 + t63;
t60 = qJ(3) ^ 2;
t59 = 0.2e1 * qJ(3);
t46 = t56 ^ 2;
t44 = t53 ^ 2;
t29 = qJ(3) * t36;
t25 = t54 * pkin(5) - t57 * pkin(10) + t51;
t19 = (-pkin(2) * t58 - pkin(1) - t68) * t48;
t17 = t53 * t25 + t50 * t75;
t16 = t56 * t25 - t53 * t77;
t14 = (t61 - t68) * t48;
t12 = t22 * t53 - t56 * t36;
t7 = t21 * pkin(5) - t22 * pkin(10) + t8;
t4 = pkin(10) * t36 + t6;
t2 = t56 * t4 + t53 * t7;
t1 = -t53 * t4 + t56 * t7;
t23 = [1, 0, 0, t43 * t55 ^ 2, 0.2e1 * t55 * t83, t81 * t89, t49 * t67, t49 ^ 2, 0.2e1 * pkin(1) * t83 - 0.2e1 * t49 * t70, -0.2e1 * t43 * pkin(1) * t55 - 0.2e1 * t24 * t49 (-t18 * t58 + t20 * t55) * t89, 0.2e1 * t19 * t36 + 0.2e1 * t20 * t49, -0.2e1 * t18 * t49 - 0.2e1 * t19 * t35, t18 ^ 2 + t19 ^ 2 + t20 ^ 2 (t11 * t55 + t15 * t58) * t89, -0.2e1 * t14 * t35 + 0.2e1 * t15 * t49, -0.2e1 * t11 * t49 - 0.2e1 * t14 * t36, t11 ^ 2 + t14 ^ 2 + t15 ^ 2, t22 ^ 2, t22 * t90, t22 * t67, t36 * t90, t43 * t58 ^ 2, 0.2e1 * t8 * t21 + 0.2e1 * t5 * t36, 0.2e1 * t8 * t22 - 0.2e1 * t6 * t36, t13 ^ 2, -0.2e1 * t13 * t12, 0.2e1 * t13 * t21, t12 * t90, t21 ^ 2, 0.2e1 * t1 * t21 + 0.2e1 * t3 * t12, 0.2e1 * t3 * t13 - 0.2e1 * t2 * t21; 0, 0, 0, 0, 0, t35, t36, t49, -t70, -t24, -pkin(2) * t35 + t29, -0.2e1 * t42 + t70, t64, -t20 * pkin(2) - t18 * qJ(3), -t51 * t35 + t29, t31 + t64, t51 * t49 - t11, t15 * qJ(3) - t11 * t51, t22 * t57, -t22 * t54 - t73, t28, -t65, 0, t51 * t21 + t28 * t50 + t8 * t54, t51 * t22 - t50 * t65 + t8 * t57, t56 * t74 (-t12 * t56 - t84) * t57, t13 * t54 + t56 * t73, -t12 * t54 - t53 * t73, t21 * t54, t1 * t54 + t16 * t21 + (-t12 * t50 + t86) * t57, -t17 * t21 - t2 * t54 + (-t13 * t50 + t85) * t57; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t59, pkin(2) ^ 2 + t60, 0, t59, t88 (t51 ^ 2) + t60, t47, t66, 0, 0, 0, t54 * t88, t57 * t88, t46 * t47, -0.2e1 * t47 * t79, t71 * t87, t53 * t66, t45, 0.2e1 * t16 * t54 - 0.2e1 * t53 * t82, -0.2e1 * t17 * t54 - 0.2e1 * t56 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t49, 0, t20, t35, 0, -t49, t11, 0, 0, 0, 0, 0, -t21, -t22, 0, 0, 0, 0, 0, -t76, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, -1, -t51, 0, 0, 0, 0, 0, -t54, -t57, 0, 0, 0, 0, 0, -t75, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t49, 0, t15, 0, 0, 0, 0, 0, t28, -t65, 0, 0, 0, 0, 0, -t57 * t12 - t21 * t37, -t21 * t75 - t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t53, t69 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, t36, t5, -t6, t84, -t53 * t12 + t13 * t56, t80, t76, 0, -pkin(5) * t12 - pkin(10) * t80 - t85, -pkin(5) * t13 - pkin(10) * t76 + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t54, 0, t72, -t77, t53 * t38 (-t44 + t46) * t57, t37, t75, 0, t62 * t53 + t56 * t72, -t53 * t72 + t56 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t54, 0, 0, 0, 0, 0, t38, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t44, 0.2e1 * t79, 0, 0, 0, pkin(5) * t87, -0.2e1 * pkin(5) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, t21, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t78, t54, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t56, 0, -t53 * pkin(10), -t56 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t23;
