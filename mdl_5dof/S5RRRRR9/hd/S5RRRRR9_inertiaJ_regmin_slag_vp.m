% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t52 = sin(qJ(4));
t53 = sin(qJ(3));
t56 = cos(qJ(4));
t57 = cos(qJ(3));
t32 = t52 * t53 - t56 * t57;
t43 = -t57 * pkin(3) - pkin(2);
t25 = t32 * pkin(4) + t43;
t80 = 0.2e1 * t25;
t79 = 0.2e1 * t43;
t54 = sin(qJ(2));
t78 = -0.2e1 * t54;
t58 = cos(qJ(2));
t77 = -0.2e1 * t58;
t76 = 0.2e1 * t58;
t75 = pkin(7) + pkin(8);
t74 = pkin(2) * t57;
t73 = pkin(6) * t53;
t51 = sin(qJ(5));
t72 = t51 * pkin(4);
t71 = t52 * pkin(3);
t55 = cos(qJ(5));
t33 = t52 * t57 + t56 * t53;
t26 = t33 * t54;
t36 = -t58 * pkin(2) - t54 * pkin(7) - pkin(1);
t31 = t57 * t36;
t63 = t57 * t54;
t17 = -pkin(8) * t63 + t31 + (-pkin(3) - t73) * t58;
t62 = t57 * t58;
t59 = pkin(6) * t62;
t19 = t59 + (-pkin(8) * t54 + t36) * t53;
t64 = t56 * t19;
t9 = t52 * t17 + t64;
t7 = -t26 * pkin(9) + t9;
t70 = t55 * t7;
t69 = t58 * pkin(3);
t68 = t58 * pkin(4);
t67 = t53 * t54;
t66 = t53 * t57;
t65 = t53 * t58;
t44 = t54 * pkin(6);
t35 = pkin(3) * t67 + t44;
t61 = t54 * t76;
t60 = t55 * t71;
t27 = t32 * t54;
t8 = t56 * t17 - t52 * t19;
t4 = t27 * pkin(9) - t68 + t8;
t1 = t55 * t4 - t51 * t7;
t37 = t75 * t53;
t38 = t75 * t57;
t20 = -t56 * t37 - t52 * t38;
t46 = t56 * pkin(3);
t42 = t46 + pkin(4);
t28 = t55 * t42 - t51 * t71;
t2 = t51 * t4 + t70;
t21 = -t52 * t37 + t56 * t38;
t50 = t58 ^ 2;
t49 = t57 ^ 2;
t48 = t54 ^ 2;
t47 = t53 ^ 2;
t45 = t55 * pkin(4);
t29 = t51 * t42 + t60;
t24 = t53 * t36 + t59;
t23 = -pkin(6) * t65 + t31;
t18 = t26 * pkin(4) + t35;
t16 = -t51 * t32 + t55 * t33;
t15 = t55 * t32 + t51 * t33;
t13 = -t32 * pkin(9) + t21;
t12 = -t33 * pkin(9) + t20;
t11 = -t51 * t26 - t55 * t27;
t10 = t55 * t26 - t51 * t27;
t6 = t51 * t12 + t55 * t13;
t5 = t55 * t12 - t51 * t13;
t3 = [1, 0, 0, t48, t61, 0, 0, 0, pkin(1) * t76, pkin(1) * t78, t49 * t48, -0.2e1 * t48 * t66, t62 * t78, t53 * t61, t50, -0.2e1 * t23 * t58 + 0.2e1 * t48 * t73, 0.2e1 * t48 * pkin(6) * t57 + 0.2e1 * t24 * t58, t27 ^ 2, 0.2e1 * t27 * t26, -t27 * t77, t26 * t76, t50, 0.2e1 * t35 * t26 - 0.2e1 * t8 * t58, -0.2e1 * t35 * t27 + 0.2e1 * t9 * t58, t11 ^ 2, -0.2e1 * t11 * t10, t11 * t77, t10 * t76, t50, -0.2e1 * t1 * t58 + 0.2e1 * t18 * t10, 0.2e1 * t18 * t11 + 0.2e1 * t2 * t58; 0, 0, 0, 0, 0, t54, t58, 0, -t44, -t58 * pkin(6), t53 * t63, (-t47 + t49) * t54, -t65, -t62, 0, -pkin(6) * t63 + (-pkin(2) * t54 + pkin(7) * t58) * t53, pkin(7) * t62 + (t73 - t74) * t54, -t27 * t33, -t33 * t26 + t27 * t32, -t58 * t33, t58 * t32, 0, -t20 * t58 + t43 * t26 + t35 * t32, t21 * t58 - t43 * t27 + t35 * t33, t11 * t16, -t16 * t10 - t11 * t15, -t16 * t58, t15 * t58, 0, t25 * t10 + t18 * t15 - t5 * t58, t25 * t11 + t18 * t16 + t6 * t58; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t47, 0.2e1 * t66, 0, 0, 0, 0.2e1 * t74, -0.2e1 * pkin(2) * t53, t33 ^ 2, -0.2e1 * t33 * t32, 0, 0, 0, t32 * t79, t33 * t79, t16 ^ 2, -0.2e1 * t16 * t15, 0, 0, 0, t15 * t80, t16 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t67, -t58, t23, -t24, 0, 0, -t27, -t26, -t58, -t56 * t69 + t8, -t64 + (-t17 + t69) * t52, 0, 0, t11, -t10, -t58, -t28 * t58 + t1, t29 * t58 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t57, 0, -t53 * pkin(7), -t57 * pkin(7), 0, 0, t33, -t32, 0, t20, -t21, 0, 0, t16, -t15, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t46, -0.2e1 * t71, 0, 0, 0, 0, 1, 0.2e1 * t28, -0.2e1 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t26, -t58, t8, -t9, 0, 0, t11, -t10, -t58, -t55 * t68 + t1, -t70 + (-t4 + t68) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, 0, t20, -t21, 0, 0, t16, -t15, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t46, -t71, 0, 0, 0, 0, 1, t28 + t45, -t60 + (-pkin(4) - t42) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t45, -0.2e1 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, -t58, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t45, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
