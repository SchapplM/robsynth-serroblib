% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRR5
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t50 = sin(pkin(6));
t55 = sin(qJ(2));
t35 = t50 * t55;
t51 = cos(pkin(6));
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t21 = t54 * t35 + t51 * t57;
t92 = -0.2e1 * t21;
t91 = 0.2e1 * t50;
t90 = -0.2e1 * t54;
t89 = 0.2e1 * t57;
t59 = pkin(2) + pkin(3);
t56 = cos(qJ(6));
t88 = pkin(5) * t56;
t58 = cos(qJ(2));
t36 = t50 * t58;
t40 = t51 * qJ(3);
t82 = t51 * t55;
t24 = pkin(1) * t82 + pkin(8) * t36;
t70 = qJ(4) * t50;
t63 = -t58 * t70 + t24;
t14 = t40 + t63;
t10 = -t51 * pkin(9) + t14;
t19 = -t50 * pkin(1) - pkin(2) * t36 - qJ(3) * t35;
t15 = pkin(3) * t36 - t19;
t8 = (pkin(4) * t58 - pkin(9) * t55) * t50 + t15;
t5 = -t54 * t10 + t57 * t8;
t3 = -pkin(5) * t36 - t5;
t53 = sin(qJ(6));
t87 = t3 * t53;
t86 = t3 * t56;
t22 = t57 * t35 - t51 * t54;
t13 = t22 * t56 + t53 * t36;
t85 = t13 * t53;
t84 = t13 * t57;
t45 = t50 ^ 2;
t83 = t45 * t58;
t52 = qJ(3) - pkin(9);
t81 = t52 * t53;
t80 = t53 * t21;
t79 = t53 * t54;
t78 = t53 * t56;
t77 = t54 * t21;
t76 = t56 * t21;
t75 = t56 * t54;
t38 = t56 * t57;
t12 = t22 * t53 - t56 * t36;
t74 = t57 * t12;
t73 = t57 * t52;
t72 = -t51 * t58 * pkin(1) + pkin(8) * t35;
t71 = qJ(3) * t58;
t46 = pkin(4) + t59;
t69 = 0.2e1 * t36;
t68 = t54 * t89;
t67 = t54 * t36;
t66 = t57 * t36;
t65 = t53 * t77;
t64 = t21 * t75;
t44 = t51 * pkin(2);
t20 = -t44 + t72;
t11 = -t51 * pkin(3) - t55 * t70 + t20;
t6 = t57 * t10 + t54 * t8;
t9 = t51 * pkin(4) - t11;
t62 = qJ(3) ^ 2;
t61 = 0.2e1 * qJ(3);
t49 = t56 ^ 2;
t48 = t54 ^ 2;
t47 = t53 ^ 2;
t39 = 0.2e1 * t40;
t37 = t53 * t57;
t25 = t57 * pkin(5) + t54 * pkin(10) + t46;
t18 = t40 + t24;
t17 = t53 * t25 + t56 * t73;
t16 = t56 * t25 - t53 * t73;
t7 = t21 * pkin(5) - t22 * pkin(10) + t9;
t4 = pkin(10) * t36 + t6;
t2 = t56 * t4 + t53 * t7;
t1 = -t53 * t4 + t56 * t7;
t23 = [1, 0, 0, t45 * t55 ^ 2, 0.2e1 * t55 * t83, t82 * t91, t51 * t69, t51 ^ 2, 0.2e1 * pkin(1) * t83 - 0.2e1 * t51 * t72, -0.2e1 * t45 * pkin(1) * t55 - 0.2e1 * t24 * t51, -0.2e1 * t19 * t36 - 0.2e1 * t20 * t51 (t18 * t58 + t20 * t55) * t91, 0.2e1 * t18 * t51 - 0.2e1 * t19 * t35, t18 ^ 2 + t19 ^ 2 + t20 ^ 2, -0.2e1 * t11 * t51 + 0.2e1 * t15 * t36, 0.2e1 * t14 * t51 + 0.2e1 * t15 * t35 (-t11 * t55 - t14 * t58) * t91, t11 ^ 2 + t14 ^ 2 + t15 ^ 2, t22 ^ 2, t22 * t92, t22 * t69, t36 * t92, t45 * t58 ^ 2, 0.2e1 * t9 * t21 + 0.2e1 * t5 * t36, 0.2e1 * t9 * t22 - 0.2e1 * t6 * t36, t13 ^ 2, -0.2e1 * t13 * t12, 0.2e1 * t13 * t21, t12 * t92, t21 ^ 2, 0.2e1 * t1 * t21 + 0.2e1 * t3 * t12, 0.2e1 * t3 * t13 - 0.2e1 * t2 * t21; 0, 0, 0, 0, 0, t35, t36, t51, -t72, -t24, 0.2e1 * t44 - t72 (-pkin(2) * t55 + t71) * t50, t39 + t24, -t20 * pkin(2) + t18 * qJ(3), t59 * t51 - t11, t39 + t63 (t55 * t59 - t71) * t50, t14 * qJ(3) - t11 * t59, -t22 * t54, -t22 * t57 + t77, -t67, -t66, 0, t46 * t21 - t52 * t67 + t9 * t57, t46 * t22 - t52 * t66 - t9 * t54, -t13 * t75 (t12 * t56 + t85) * t54, -t64 + t84, t65 - t74, t21 * t57, t1 * t57 + t16 * t21 + (t12 * t52 - t87) * t54, -t17 * t21 - t2 * t57 + (t13 * t52 - t86) * t54; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, t61, pkin(2) ^ 2 + t62, 0.2e1 * t59, t61, 0, t59 ^ 2 + t62, t48, t68, 0, 0, 0, t46 * t89, t46 * t90, t49 * t48, -0.2e1 * t48 * t78, t38 * t90, t53 * t68, t57 ^ 2, 0.2e1 * t16 * t57 - 0.2e1 * t48 * t81, -0.2e1 * t48 * t52 * t56 - 0.2e1 * t17 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t35, 0, t20, -t51, 0, -t35, t11, 0, 0, 0, 0, 0, -t21, -t22, 0, 0, 0, 0, 0, -t76, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), -1, 0, 0, -t59, 0, 0, 0, 0, 0, -t57, t54, 0, 0, 0, 0, 0, -t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t35, 0, t15, 0, 0, 0, 0, 0, t66, -t67, 0, 0, 0, 0, 0, -t65 - t74, -t64 - t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, t36, t5, -t6, t85, -t53 * t12 + t13 * t56, t80, t76, 0, -pkin(5) * t12 - pkin(10) * t80 - t86, -pkin(5) * t13 - pkin(10) * t76 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t57, 0, -t54 * t52, -t73, -t53 * t75 (t47 - t49) * t54, t37, t38, 0, -t52 * t75 + (pkin(5) * t54 - pkin(10) * t57) * t53, -pkin(10) * t38 + (t81 + t88) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t54, 0, 0, 0, 0, 0, t38, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t47, 0.2e1 * t78, 0, 0, 0, 0.2e1 * t88, -0.2e1 * pkin(5) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, t21, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t79, t57, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t56, 0, -t53 * pkin(10), -t56 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t23;
