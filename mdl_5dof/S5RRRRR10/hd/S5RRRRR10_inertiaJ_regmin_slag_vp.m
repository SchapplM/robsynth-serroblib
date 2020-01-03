% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t58 = cos(pkin(5));
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t57 = sin(pkin(5));
t62 = sin(qJ(2));
t86 = t57 * t62;
t36 = -t58 * t64 + t61 * t86;
t37 = t58 * t61 + t64 * t86;
t60 = sin(qJ(4));
t89 = cos(qJ(4));
t22 = t89 * t36 + t60 * t37;
t98 = -0.2e1 * t22;
t53 = -t64 * pkin(3) - pkin(2);
t97 = 0.2e1 * t53;
t96 = 0.2e1 * t64;
t95 = pkin(8) + pkin(9);
t94 = pkin(1) * t62;
t65 = cos(qJ(2));
t93 = pkin(1) * t65;
t85 = t57 * t65;
t73 = pkin(7) * t85;
t32 = t73 + (pkin(8) + t94) * t58;
t33 = (-pkin(2) * t65 - pkin(8) * t62 - pkin(1)) * t57;
t18 = -t61 * t32 + t64 * t33;
t74 = pkin(3) * t85;
t11 = -t37 * pkin(9) + t18 - t74;
t19 = t64 * t32 + t61 * t33;
t14 = -t36 * pkin(9) + t19;
t6 = t89 * t11 - t60 * t14;
t4 = pkin(4) * t85 - t6;
t63 = cos(qJ(5));
t92 = t4 * t63;
t91 = t60 * pkin(3);
t70 = t89 * pkin(3);
t52 = -t70 - pkin(4);
t90 = pkin(4) - t52;
t23 = -t60 * t36 + t89 * t37;
t59 = sin(qJ(5));
t17 = t63 * t23 - t59 * t85;
t15 = t17 * t59;
t46 = t95 * t64;
t68 = t89 * t61;
t28 = t60 * t46 + t95 * t68;
t88 = t28 * t63;
t54 = t57 ^ 2;
t87 = t54 * t65;
t84 = t58 * t62;
t20 = t59 * t22;
t44 = t60 * t64 + t68;
t83 = t59 * t44;
t51 = pkin(10) + t91;
t82 = t59 * t51;
t81 = t59 * t63;
t80 = t60 * t61;
t21 = t63 * t22;
t79 = t63 * t44;
t78 = t63 * t51;
t43 = -t89 * t64 + t80;
t77 = -0.2e1 * t44 * t43;
t76 = -0.2e1 * t85;
t75 = 0.2e1 * t85;
t72 = t61 * t85;
t71 = t64 * t85;
t69 = t89 * t14;
t67 = -pkin(4) * t44 - pkin(10) * t43;
t66 = -t43 * t51 + t44 * t52;
t47 = pkin(7) * t86;
t31 = t47 + (-pkin(2) - t93) * t58;
t7 = t60 * t11 + t69;
t24 = t36 * pkin(3) + t31;
t56 = t63 ^ 2;
t55 = t59 ^ 2;
t49 = t54 * t65 ^ 2;
t48 = 0.2e1 * t81;
t42 = t44 ^ 2;
t41 = pkin(1) * t84 + t73;
t40 = t58 * t93 - t47;
t39 = t63 * t43;
t38 = t59 * t43;
t35 = t59 * t79;
t29 = t89 * t46 - t95 * t80;
t27 = t28 * t59;
t26 = (-t55 + t56) * t44;
t25 = t43 * pkin(4) - t44 * pkin(10) + t53;
t16 = t59 * t23 + t63 * t85;
t13 = t59 * t25 + t63 * t29;
t12 = t63 * t25 - t59 * t29;
t9 = -t59 * t16 + t17 * t63;
t8 = t22 * pkin(4) - t23 * pkin(10) + t24;
t5 = -pkin(10) * t85 + t7;
t3 = t4 * t59;
t2 = t63 * t5 + t59 * t8;
t1 = -t59 * t5 + t63 * t8;
t10 = [1, 0, 0, t54 * t62 ^ 2, 0.2e1 * t62 * t87, 0.2e1 * t57 * t84, t58 * t75, t58 ^ 2, 0.2e1 * pkin(1) * t87 + 0.2e1 * t40 * t58, -0.2e1 * t41 * t58 - 0.2e1 * t54 * t94, t37 ^ 2, -0.2e1 * t37 * t36, t37 * t76, t36 * t75, t49, -0.2e1 * t18 * t85 + 0.2e1 * t31 * t36, 0.2e1 * t19 * t85 + 0.2e1 * t31 * t37, t23 ^ 2, t23 * t98, t23 * t76, t22 * t75, t49, 0.2e1 * t24 * t22 - 0.2e1 * t6 * t85, 0.2e1 * t24 * t23 + 0.2e1 * t7 * t85, t17 ^ 2, -0.2e1 * t17 * t16, 0.2e1 * t17 * t22, t16 * t98, t22 ^ 2, 0.2e1 * t1 * t22 + 0.2e1 * t4 * t16, 0.2e1 * t4 * t17 - 0.2e1 * t2 * t22; 0, 0, 0, 0, 0, t86, t85, t58, t40, -t41, t37 * t61, -t61 * t36 + t37 * t64, -t72, -t71, 0, -pkin(2) * t36 + pkin(8) * t72 - t31 * t64, -pkin(2) * t37 + pkin(8) * t71 + t31 * t61, t23 * t44, -t44 * t22 - t23 * t43, -t44 * t85, t43 * t85, 0, t53 * t22 + t24 * t43 + t28 * t85, t53 * t23 + t24 * t44 + t29 * t85, t17 * t79, (-t16 * t63 - t15) * t44, t17 * t43 + t22 * t79, -t16 * t43 - t22 * t83, t22 * t43, t1 * t43 + t12 * t22 + t28 * t16 + t4 * t83, -t13 * t22 + t28 * t17 - t2 * t43 + t4 * t79; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t61 ^ 2, t61 * t96, 0, 0, 0, pkin(2) * t96, -0.2e1 * pkin(2) * t61, t42, t77, 0, 0, 0, t43 * t97, t44 * t97, t56 * t42, -0.2e1 * t42 * t81, 0.2e1 * t43 * t79, t59 * t77, t43 ^ 2, 0.2e1 * t12 * t43 + 0.2e1 * t28 * t83, -0.2e1 * t13 * t43 + 0.2e1 * t28 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, -t85, t18, -t19, 0, 0, t23, -t22, -t85, -t70 * t85 + t6, -t69 + (-t11 + t74) * t60, t15, t9, t20, t21, 0, t52 * t16 - t22 * t82 - t92, t52 * t17 - t22 * t78 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t64, 0, -t61 * pkin(8), -t64 * pkin(8), 0, 0, t44, -t43, 0, -t28, -t29, t35, t26, t38, t39, 0, t59 * t66 - t88, t63 * t66 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t70, -0.2e1 * t91, t55, t48, 0, 0, 0, -0.2e1 * t52 * t63, 0.2e1 * t52 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, -t85, t6, -t7, t15, t9, t20, t21, 0, -pkin(4) * t16 - pkin(10) * t20 - t92, -pkin(4) * t17 - pkin(10) * t21 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43, 0, -t28, -t29, t35, t26, t38, t39, 0, t59 * t67 - t88, t63 * t67 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t70, -t91, t55, t48, 0, 0, 0, t90 * t63, -t90 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t55, t48, 0, 0, 0, 0.2e1 * pkin(4) * t63, -0.2e1 * pkin(4) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t83, t43, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t63, 0, -t82, -t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t63, 0, -t59 * pkin(10), -t63 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t10;
