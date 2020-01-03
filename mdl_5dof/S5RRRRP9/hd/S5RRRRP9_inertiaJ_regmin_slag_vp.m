% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t40 = sin(qJ(2));
t42 = cos(qJ(3));
t52 = t42 * t40;
t39 = sin(qJ(3));
t55 = t39 * t40;
t17 = -t38 * t55 + t41 * t52;
t68 = -0.2e1 * t17;
t20 = t38 * t42 + t41 * t39;
t67 = -0.2e1 * t20;
t31 = -t42 * pkin(3) - pkin(2);
t66 = 0.2e1 * t31;
t65 = -0.2e1 * t40;
t43 = cos(qJ(2));
t64 = 0.2e1 * t43;
t63 = -pkin(8) - pkin(7);
t62 = pkin(2) * t42;
t61 = pkin(6) * t39;
t60 = t41 * pkin(3);
t59 = t43 * pkin(3);
t58 = t43 * pkin(4);
t23 = -t43 * pkin(2) - t40 * pkin(7) - pkin(1);
t51 = t42 * t43;
t48 = pkin(6) * t51;
t11 = t48 + (-pkin(8) * t40 + t23) * t39;
t18 = t42 * t23;
t9 = -pkin(8) * t52 + t18 + (-pkin(3) - t61) * t43;
t4 = t41 * t11 + t38 * t9;
t24 = t63 * t42;
t47 = t63 * t39;
t12 = -t38 * t24 - t41 * t47;
t57 = t12 * t43;
t13 = -t41 * t24 + t38 * t47;
t56 = t13 * t43;
t54 = t39 * t42;
t53 = t39 * t43;
t33 = t40 * pkin(6);
t22 = pkin(3) * t55 + t33;
t50 = t43 * qJ(5);
t49 = t40 * t64;
t46 = t38 * t11 - t41 * t9;
t45 = 0.2e1 * pkin(4);
t44 = 0.2e1 * qJ(5);
t37 = t43 ^ 2;
t36 = t42 ^ 2;
t35 = t40 ^ 2;
t34 = t39 ^ 2;
t32 = t38 * pkin(3);
t29 = pkin(4) + t60;
t27 = t32 + qJ(5);
t19 = t38 * t39 - t41 * t42;
t16 = t20 * t40;
t15 = t39 * t23 + t48;
t14 = -pkin(6) * t53 + t18;
t8 = t19 * pkin(4) - t20 * qJ(5) + t31;
t5 = t16 * pkin(4) - t17 * qJ(5) + t22;
t2 = t46 + t58;
t1 = -t50 + t4;
t3 = [1, 0, 0, t35, t49, 0, 0, 0, pkin(1) * t64, pkin(1) * t65, t36 * t35, -0.2e1 * t35 * t54, t51 * t65, t39 * t49, t37, -0.2e1 * t14 * t43 + 0.2e1 * t35 * t61, 0.2e1 * t35 * pkin(6) * t42 + 0.2e1 * t15 * t43, t17 ^ 2, t16 * t68, t43 * t68, t16 * t64, t37, 0.2e1 * t22 * t16 + 0.2e1 * t43 * t46, 0.2e1 * t22 * t17 + 0.2e1 * t4 * t43, 0.2e1 * t5 * t16 + 0.2e1 * t2 * t43, -0.2e1 * t1 * t16 + 0.2e1 * t2 * t17, -0.2e1 * t1 * t43 - 0.2e1 * t5 * t17, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t40, t43, 0, -t33, -t43 * pkin(6), t39 * t52, (-t34 + t36) * t40, -t53, -t51, 0, -pkin(6) * t52 + (-pkin(2) * t40 + pkin(7) * t43) * t39, pkin(7) * t51 + (t61 - t62) * t40, t17 * t20, -t20 * t16 - t17 * t19, -t20 * t43, t19 * t43, 0, t31 * t16 + t22 * t19 + t57, t31 * t17 + t22 * t20 + t56, t8 * t16 + t5 * t19 + t57, -t1 * t19 + t12 * t17 - t13 * t16 + t2 * t20, -t8 * t17 - t5 * t20 - t56, t1 * t13 + t2 * t12 + t5 * t8; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t34, 0.2e1 * t54, 0, 0, 0, 0.2e1 * t62, -0.2e1 * pkin(2) * t39, t20 ^ 2, t19 * t67, 0, 0, 0, t19 * t66, t20 * t66, 0.2e1 * t8 * t19, 0.2e1 * t12 * t20 - 0.2e1 * t13 * t19, t8 * t67, t12 ^ 2 + t13 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t55, -t43, t14, -t15, 0, 0, t17, -t16, -t43, -t41 * t59 - t46, t38 * t59 - t4, (-pkin(4) - t29) * t43 - t46, -t27 * t16 - t29 * t17, (-qJ(5) - t27) * t43 + t4, t1 * t27 - t2 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t42, 0, -t39 * pkin(7), -t42 * pkin(7), 0, 0, t20, -t19, 0, -t12, -t13, -t12, -t27 * t19 - t29 * t20, t13, -t12 * t29 + t13 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t60, -0.2e1 * t32, 0.2e1 * t29, 0, 0.2e1 * t27, t27 ^ 2 + t29 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t43, -t46, -t4, -t46 - 0.2e1 * t58, -pkin(4) * t17 - t16 * qJ(5), -0.2e1 * t50 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, -t12, -t13, -t12, -pkin(4) * t20 - t19 * qJ(5), t13, -t12 * pkin(4) + t13 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t60, -t32, t45 + t60, 0, t44 + t32, t29 * pkin(4) + t27 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t45, 0, t44, pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t17, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
