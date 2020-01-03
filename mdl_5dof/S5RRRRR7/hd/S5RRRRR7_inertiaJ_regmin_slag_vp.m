% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR7
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t47 = cos(qJ(2));
t34 = -t47 * pkin(2) - pkin(1);
t42 = sin(qJ(3));
t43 = sin(qJ(2));
t46 = cos(qJ(3));
t49 = t42 * t43 - t46 * t47;
t19 = pkin(3) * t49 + t34;
t68 = 0.2e1 * t19;
t67 = 0.2e1 * t34;
t66 = 0.2e1 * t47;
t65 = pkin(6) + pkin(7);
t40 = sin(qJ(5));
t64 = pkin(4) * t40;
t41 = sin(qJ(4));
t63 = t41 * pkin(3);
t62 = t42 * pkin(2);
t44 = cos(qJ(5));
t45 = cos(qJ(4));
t26 = t65 * t43;
t27 = t65 * t47;
t16 = -t26 * t46 - t27 * t42;
t25 = t42 * t47 + t43 * t46;
t48 = -pkin(8) * t25 + t16;
t17 = t42 * t26 - t46 * t27;
t9 = -pkin(8) * t49 - t17;
t5 = t41 * t9 - t45 * t48;
t61 = t5 * t44;
t36 = t46 * pkin(2);
t33 = t36 + pkin(3);
t53 = -t45 * t33 + t41 * t62;
t20 = -pkin(4) + t53;
t60 = t20 * t44;
t35 = t45 * pkin(3);
t32 = -t35 - pkin(4);
t59 = t32 * t44;
t15 = t45 * t25 - t41 * t49;
t58 = t40 * t15;
t57 = t40 * t44;
t56 = t44 * t15;
t14 = t41 * t25 + t45 * t49;
t55 = -0.2e1 * t15 * t14;
t54 = t45 * t62;
t52 = -pkin(4) * t15 - pkin(9) * t14;
t23 = -t33 * t41 - t54;
t21 = pkin(9) - t23;
t51 = -t14 * t21 + t15 * t20;
t31 = pkin(9) + t63;
t50 = -t14 * t31 + t15 * t32;
t39 = t44 ^ 2;
t38 = t40 ^ 2;
t37 = pkin(4) * t44;
t30 = 0.2e1 * t57;
t29 = t32 * t40;
t18 = t20 * t40;
t13 = t15 ^ 2;
t12 = t44 * t14;
t11 = t40 * t14;
t10 = t40 * t56;
t7 = (-t38 + t39) * t15;
t6 = t41 * t48 + t45 * t9;
t4 = t14 * pkin(4) - t15 * pkin(9) + t19;
t3 = t5 * t40;
t2 = t4 * t40 + t44 * t6;
t1 = t4 * t44 - t40 * t6;
t8 = [1, 0, 0, t43 ^ 2, t43 * t66, 0, 0, 0, pkin(1) * t66, -0.2e1 * pkin(1) * t43, t25 ^ 2, -0.2e1 * t25 * t49, 0, 0, 0, t49 * t67, t25 * t67, t13, t55, 0, 0, 0, t14 * t68, t15 * t68, t39 * t13, -0.2e1 * t13 * t57, 0.2e1 * t14 * t56, t40 * t55, t14 ^ 2, 0.2e1 * t1 * t14 + 0.2e1 * t5 * t58, -0.2e1 * t14 * t2 + 0.2e1 * t5 * t56; 0, 0, 0, 0, 0, t43, t47, 0, -t43 * pkin(6), -t47 * pkin(6), 0, 0, t25, -t49, 0, t16, t17, 0, 0, t15, -t14, 0, -t5, -t6, t10, t7, t11, t12, 0, t40 * t51 - t61, t44 * t51 + t3; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t36, -0.2e1 * t62, 0, 0, 0, 0, 1, -0.2e1 * t53, 0.2e1 * t23, t38, t30, 0, 0, 0, -0.2e1 * t60, 0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t49, 0, t16, t17, 0, 0, t15, -t14, 0, -t5, -t6, t10, t7, t11, t12, 0, t40 * t50 - t61, t44 * t50 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t36, -t62, 0, 0, 0, 0, 1, t35 - t53, -t54 + (-pkin(3) - t33) * t41, t38, t30, 0, 0, 0, (-t20 - t32) * t44, t29 + t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t63, t38, t30, 0, 0, 0, -0.2e1 * t59, 0.2e1 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, -t5, -t6, t10, t7, t11, t12, 0, t40 * t52 - t61, t44 * t52 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t53, t23, t38, t30, 0, 0, 0, t37 - t60, t18 - t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t63, t38, t30, 0, 0, 0, t37 - t59, t29 - t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t38, t30, 0, 0, 0, 0.2e1 * t37, -0.2e1 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t58, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t44, 0, -t40 * t21, -t44 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t44, 0, -t40 * t31, -t44 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t44, 0, -t40 * pkin(9), -t44 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
