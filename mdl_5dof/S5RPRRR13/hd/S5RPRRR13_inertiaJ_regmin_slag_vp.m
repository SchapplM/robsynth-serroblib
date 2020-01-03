% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR13_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t29 = sin(qJ(5));
t30 = sin(qJ(4));
t32 = cos(qJ(5));
t33 = cos(qJ(4));
t17 = t29 * t33 + t32 * t30;
t34 = cos(qJ(3));
t11 = t34 * t17;
t58 = -0.2e1 * t11;
t24 = -t33 * pkin(4) - pkin(3);
t57 = 0.2e1 * t24;
t56 = 0.2e1 * t33;
t55 = 2 * qJ(2);
t54 = pkin(7) + pkin(8);
t53 = t29 * pkin(4);
t52 = t32 * pkin(4);
t31 = sin(qJ(3));
t18 = t31 * pkin(3) - t34 * pkin(7) + qJ(2);
t35 = -pkin(1) - pkin(6);
t43 = t33 * t35;
t37 = t31 * t43;
t5 = t37 + (-pkin(8) * t34 + t18) * t30;
t51 = t32 * t5;
t50 = t17 * t31;
t49 = t30 * t31;
t48 = t30 * t33;
t47 = t30 * t34;
t46 = t30 * t35;
t45 = t31 * t35;
t44 = t33 * t31;
t23 = t33 * t34;
t16 = t29 * t30 - t32 * t33;
t42 = t34 * t16;
t41 = t34 * t31;
t40 = t34 * t35;
t26 = t31 ^ 2;
t28 = t34 ^ 2;
t39 = -t26 - t28;
t38 = -0.2e1 * t41;
t14 = t33 * t18;
t4 = -pkin(8) * t23 + t14 + (pkin(4) - t46) * t31;
t1 = -t29 * t5 + t32 * t4;
t36 = -pkin(3) * t34 - pkin(7) * t31;
t27 = t33 ^ 2;
t25 = t30 ^ 2;
t20 = t54 * t33;
t19 = t54 * t30;
t15 = (pkin(4) * t30 - t35) * t34;
t12 = -t29 * t49 + t32 * t44;
t9 = t30 * t18 + t37;
t8 = -t30 * t45 + t14;
t7 = -t29 * t19 + t32 * t20;
t6 = -t32 * t19 - t29 * t20;
t2 = t29 * t4 + t51;
t3 = [1, 0, 0, -2 * pkin(1), t55, pkin(1) ^ 2 + qJ(2) ^ 2, t28, t38, 0, 0, 0, t31 * t55, t34 * t55, t27 * t28, -0.2e1 * t28 * t48, t41 * t56, t30 * t38, t26, -0.2e1 * t28 * t46 + 0.2e1 * t8 * t31, -0.2e1 * t28 * t43 - 0.2e1 * t9 * t31, t42 ^ 2, -t42 * t58, -0.2e1 * t42 * t31, t31 * t58, t26, 0.2e1 * t1 * t31 + 0.2e1 * t15 * t11, -0.2e1 * t15 * t42 - 0.2e1 * t2 * t31; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t30, t39 * t33, 0, 0, 0, 0, 0, -t34 * t11 - t31 * t50, -t12 * t31 + t34 * t42; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t34, -t31, 0, t40, -t45, t30 * t23, (-t25 + t27) * t34, t49, t44, 0, t30 * t36 + t33 * t40, -t30 * t40 + t36 * t33, -t42 * t17, -t17 * t11 + t16 * t42, t50, -t16 * t31, 0, t24 * t11 + t15 * t16 + t6 * t31, t15 * t17 - t24 * t42 - t7 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t31, 0, 0, 0, 0, 0, t23, -t47, 0, 0, 0, 0, 0, -t42, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t25, 0.2e1 * t48, 0, 0, 0, pkin(3) * t56, -0.2e1 * pkin(3) * t30, t17 ^ 2, -0.2e1 * t17 * t16, 0, 0, 0, t16 * t57, t17 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t47, t31, t8, -t9, 0, 0, -t42, -t11, t31, t31 * t52 + t1, -t51 + (-t31 * pkin(4) - t4) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t44, 0, 0, 0, 0, 0, -t50, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t33, 0, -t30 * pkin(7), -t33 * pkin(7), 0, 0, t17, -t16, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t52, -0.2e1 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t11, t31, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t52, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
