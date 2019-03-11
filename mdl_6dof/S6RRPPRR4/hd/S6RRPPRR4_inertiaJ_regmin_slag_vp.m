% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t51 = sin(pkin(11));
t53 = cos(pkin(11));
t52 = sin(pkin(6));
t60 = cos(qJ(2));
t82 = t52 * t60;
t57 = sin(qJ(2));
t83 = t52 * t57;
t31 = t51 * t83 - t53 * t82;
t54 = cos(pkin(6));
t56 = sin(qJ(5));
t59 = cos(qJ(5));
t22 = -t31 * t59 + t54 * t56;
t94 = -0.2e1 * t22;
t40 = t51 * pkin(2) + qJ(4);
t93 = 0.2e1 * t40;
t58 = cos(qJ(6));
t92 = 0.2e1 * t58;
t91 = pkin(3) + pkin(9);
t90 = pkin(1) * t57;
t32 = (t51 * t60 + t53 * t57) * t52;
t36 = (-pkin(2) * t60 - pkin(1)) * t52;
t62 = -t32 * qJ(4) + t36;
t10 = t91 * t31 + t62;
t38 = t54 * t60 * pkin(1);
t70 = pkin(8) + qJ(3);
t24 = t54 * pkin(2) - t70 * t83 + t38;
t66 = t54 * t90;
t27 = t70 * t82 + t66;
t14 = t53 * t24 - t51 * t27;
t8 = t32 * pkin(4) - t91 * t54 - t14;
t5 = -t56 * t10 + t59 * t8;
t3 = -t32 * pkin(5) - t5;
t55 = sin(qJ(6));
t89 = t3 * t55;
t88 = t3 * t58;
t23 = t31 * t56 + t54 * t59;
t17 = t23 * t58 + t32 * t55;
t87 = t17 * t55;
t86 = t32 * t56;
t30 = t32 * t59;
t42 = -t53 * pkin(2) - pkin(3);
t39 = -pkin(9) + t42;
t50 = t59 ^ 2;
t85 = t39 * t50;
t46 = t52 ^ 2;
t84 = t46 * t60;
t81 = t55 * t22;
t43 = t55 * t56;
t80 = t55 * t58;
t79 = t55 * t59;
t16 = t23 * t55 - t32 * t58;
t78 = t56 * t16;
t77 = t56 * t39;
t76 = t58 * t22;
t75 = t58 * t56;
t44 = t58 * t59;
t74 = t59 * t17;
t73 = t59 * t22;
t72 = t59 * t39;
t71 = t59 * t56;
t15 = t51 * t24 + t53 * t27;
t48 = t56 ^ 2;
t69 = -t48 - t50;
t68 = 0.2e1 * t52 * t54;
t67 = -0.2e1 * t71;
t65 = t55 * t73;
t64 = t58 * t73;
t11 = -t54 * qJ(4) - t15;
t63 = -pkin(5) * t59 - pkin(10) * t56;
t6 = t59 * t10 + t56 * t8;
t9 = -t31 * pkin(4) - t11;
t49 = t58 ^ 2;
t47 = t55 ^ 2;
t35 = t56 * pkin(5) - t59 * pkin(10) + t40;
t34 = pkin(8) * t82 + t66;
t33 = -pkin(8) * t83 + t38;
t20 = t55 * t35 + t39 * t75;
t19 = t58 * t35 - t55 * t77;
t18 = t31 * pkin(3) + t62;
t13 = t17 * t56;
t12 = -t54 * pkin(3) - t14;
t7 = t22 * pkin(5) - t23 * pkin(10) + t9;
t4 = t32 * pkin(10) + t6;
t2 = t58 * t4 + t55 * t7;
t1 = -t55 * t4 + t58 * t7;
t21 = [1, 0, 0, t46 * t57 ^ 2, 0.2e1 * t57 * t84, t57 * t68, t60 * t68, t54 ^ 2, 0.2e1 * pkin(1) * t84 + 0.2e1 * t33 * t54, -0.2e1 * t34 * t54 - 0.2e1 * t46 * t90, -0.2e1 * t14 * t32 - 0.2e1 * t15 * t31, t14 ^ 2 + t15 ^ 2 + t36 ^ 2, 0.2e1 * t11 * t31 + 0.2e1 * t12 * t32, 0.2e1 * t12 * t54 - 0.2e1 * t18 * t31, -0.2e1 * t11 * t54 - 0.2e1 * t18 * t32, t11 ^ 2 + t12 ^ 2 + t18 ^ 2, t23 ^ 2, t23 * t94, 0.2e1 * t23 * t32, t32 * t94, t32 ^ 2, 0.2e1 * t9 * t22 + 0.2e1 * t5 * t32, 0.2e1 * t9 * t23 - 0.2e1 * t6 * t32, t17 ^ 2, -0.2e1 * t17 * t16, 0.2e1 * t17 * t22, t16 * t94, t22 ^ 2, 0.2e1 * t1 * t22 + 0.2e1 * t3 * t16, 0.2e1 * t3 * t17 - 0.2e1 * t2 * t22; 0, 0, 0, 0, 0, t83, t82, t54, t33, -t34 (-t31 * t51 - t32 * t53) * pkin(2) (t14 * t53 + t15 * t51) * pkin(2), -t40 * t31 + t42 * t32 (-pkin(3) + t42) * t54 - t14, t40 * t54 - t11, -t11 * t40 + t12 * t42, t23 * t59, -t23 * t56 - t73, t30, -t86, 0, t40 * t22 + t32 * t72 + t9 * t56, t40 * t23 - t32 * t77 + t9 * t59, t58 * t74 (-t16 * t58 - t87) * t59, t13 + t64, -t65 - t78, t22 * t56, t1 * t56 + t19 * t22 + (-t16 * t39 + t89) * t59, -t2 * t56 - t20 * t22 + (-t17 * t39 + t88) * t59; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t51 ^ 2 + t53 ^ 2) * pkin(2) ^ 2, 0, 0.2e1 * t42, t93, t40 ^ 2 + t42 ^ 2, t50, t67, 0, 0, 0, t56 * t93, t59 * t93, t49 * t50, -0.2e1 * t50 * t80, t71 * t92, t55 * t67, t48, 0.2e1 * t19 * t56 - 0.2e1 * t55 * t85, -0.2e1 * t20 * t56 - 0.2e1 * t58 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, -t31, -t32, t18, 0, 0, 0, 0, 0, -t86, -t30, 0, 0, 0, 0, 0, -t65 + t78, t13 - t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t54, 0, t12, 0, 0, 0, 0, 0, t30, -t86, 0, 0, 0, 0, 0, -t59 * t16 - t22 * t43, -t22 * t75 - t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t55, t69 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, t32, t5, -t6, t87, -t55 * t16 + t17 * t58, t81, t76, 0, -pkin(5) * t16 - pkin(10) * t81 - t88, -pkin(5) * t17 - pkin(10) * t76 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t56, 0, t72, -t77, t55 * t44 (-t47 + t49) * t59, t43, t75, 0, t63 * t55 + t58 * t72, -t55 * t72 + t58 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t59, 0, 0, 0, 0, 0, -t75, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t56, 0, 0, 0, 0, 0, t44, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t47, 0.2e1 * t80, 0, 0, 0, pkin(5) * t92, -0.2e1 * pkin(5) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t79, t56, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t58, 0, -t55 * pkin(10), -t58 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t21;
