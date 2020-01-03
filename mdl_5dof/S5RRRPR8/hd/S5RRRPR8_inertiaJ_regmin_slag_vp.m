% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t32 = sin(qJ(3));
t33 = sin(qJ(2));
t35 = cos(qJ(3));
t36 = cos(qJ(2));
t17 = t32 * t33 - t35 * t36;
t18 = t32 * t36 + t35 * t33;
t28 = -t36 * pkin(2) - pkin(1);
t40 = -t18 * qJ(4) + t28;
t8 = t17 * pkin(3) + t40;
t59 = -0.2e1 * t8;
t53 = t32 * pkin(2);
t24 = qJ(4) + t53;
t58 = 0.2e1 * t24;
t57 = 0.2e1 * t28;
t56 = 0.2e1 * t36;
t38 = 0.2e1 * qJ(4);
t55 = pkin(3) + pkin(8);
t54 = -pkin(7) - pkin(6);
t52 = t35 * pkin(2);
t51 = t18 * t17;
t50 = t24 * t17;
t31 = sin(qJ(5));
t49 = t31 * t17;
t48 = t31 * t18;
t34 = cos(qJ(5));
t47 = t34 * t17;
t46 = t34 * t31;
t45 = qJ(4) * t17;
t44 = qJ(4) + t24;
t43 = 0.2e1 * t51;
t27 = -pkin(3) - t52;
t20 = t54 * t33;
t21 = t54 * t36;
t10 = -t35 * t20 - t32 * t21;
t22 = -pkin(8) + t27;
t42 = -t18 * t22 + t50;
t11 = t32 * t20 - t35 * t21;
t41 = t18 * t55 + t45;
t39 = -0.2e1 * pkin(3);
t30 = t34 ^ 2;
t29 = t31 ^ 2;
t23 = -0.2e1 * t46;
t16 = t18 ^ 2;
t15 = t17 ^ 2;
t14 = t34 * t18;
t13 = t17 * t46;
t9 = (-t29 + t30) * t17;
t7 = -t17 * pkin(4) + t11;
t6 = t18 * pkin(4) + t10;
t5 = t7 * t34;
t4 = t7 * t31;
t3 = t55 * t17 + t40;
t2 = t34 * t3 + t31 * t6;
t1 = -t31 * t3 + t34 * t6;
t12 = [1, 0, 0, t33 ^ 2, t33 * t56, 0, 0, 0, pkin(1) * t56, -0.2e1 * pkin(1) * t33, t16, -0.2e1 * t51, 0, 0, 0, t17 * t57, t18 * t57, 0.2e1 * t10 * t18 - 0.2e1 * t11 * t17, t17 * t59, t18 * t59, t10 ^ 2 + t11 ^ 2 + t8 ^ 2, t29 * t15, 0.2e1 * t15 * t46, t31 * t43, t34 * t43, t16, 0.2e1 * t1 * t18 - 0.2e1 * t7 * t47, -0.2e1 * t2 * t18 + 0.2e1 * t7 * t49; 0, 0, 0, 0, 0, t33, t36, 0, -t33 * pkin(6), -t36 * pkin(6), 0, 0, t18, -t17, 0, -t10, -t11, t27 * t18 - t50, t10, t11, t10 * t27 + t11 * t24, t13, t9, t14, -t48, 0, -t42 * t34 + t4, t42 * t31 + t5; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t52, -0.2e1 * t53, 0, 0.2e1 * t27, t58, t24 ^ 2 + t27 ^ 2, t30, t23, 0, 0, 0, t31 * t58, t34 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, -t10, -t11, -pkin(3) * t18 - t45, t10, t11, -t10 * pkin(3) + t11 * qJ(4), t13, t9, t14, -t48, 0, -t41 * t34 + t4, t41 * t31 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t52, -t53, 0, t39 - t52, t38 + t53, -t27 * pkin(3) + t24 * qJ(4), t30, t23, 0, 0, 0, t44 * t31, t44 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t39, t38, pkin(3) ^ 2 + qJ(4) ^ 2, t30, t23, 0, 0, 0, t31 * t38, t34 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, t10, 0, 0, 0, 0, 0, t14, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t47, t18, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t31, 0, t34 * t22, -t31 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t31, 0, -t34 * t55, t31 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t12;
