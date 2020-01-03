% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t31 = cos(pkin(9));
t24 = -t31 * pkin(4) - pkin(3);
t56 = 0.2e1 * t24;
t34 = sin(qJ(3));
t55 = 0.2e1 * t34;
t36 = cos(qJ(3));
t54 = -0.2e1 * t36;
t53 = 0.2e1 * t36;
t52 = t36 * pkin(3);
t29 = sin(pkin(9));
t33 = sin(qJ(5));
t35 = cos(qJ(5));
t16 = t35 * t29 + t33 * t31;
t51 = t16 * t36;
t30 = sin(pkin(8));
t22 = t30 * pkin(1) + pkin(6);
t50 = t22 * t29;
t49 = t29 * t34;
t48 = t31 * t34;
t20 = t34 * t22;
t15 = t33 * t29 - t35 * t31;
t47 = t36 * t15;
t46 = t36 * t22;
t45 = t36 * t29;
t44 = t36 * t31;
t43 = pkin(7) + qJ(4);
t32 = cos(pkin(8));
t23 = -t32 * pkin(1) - pkin(2);
t14 = -t34 * qJ(4) + t23 - t52;
t6 = t29 * t14 + t22 * t44;
t42 = t29 ^ 2 + t31 ^ 2;
t41 = t42 * qJ(4);
t12 = t31 * t14;
t5 = -t22 * t45 + t12;
t40 = -t5 * t29 + t6 * t31;
t39 = -pkin(3) * t34 + qJ(4) * t36;
t28 = t36 ^ 2;
t27 = t34 ^ 2;
t19 = t43 * t31;
t18 = t43 * t29;
t13 = pkin(4) * t49 + t20;
t10 = t15 * t34;
t9 = t16 * t34;
t8 = -t33 * t18 + t35 * t19;
t7 = -t35 * t18 - t33 * t19;
t4 = -pkin(7) * t49 + t6;
t3 = -pkin(7) * t48 + t12 + (-pkin(4) - t50) * t36;
t2 = t33 * t3 + t35 * t4;
t1 = t35 * t3 - t33 * t4;
t11 = [1, 0, 0, (t30 ^ 2 + t32 ^ 2) * pkin(1) ^ 2, t27, t34 * t53, 0, 0, 0, t23 * t54, t23 * t55, 0.2e1 * t27 * t50 - 0.2e1 * t5 * t36, 0.2e1 * t27 * t22 * t31 + 0.2e1 * t6 * t36, (-t29 * t6 - t31 * t5) * t55, t27 * t22 ^ 2 + t5 ^ 2 + t6 ^ 2, t10 ^ 2, 0.2e1 * t10 * t9, -t10 * t54, t9 * t53, t28, -0.2e1 * t1 * t36 + 0.2e1 * t13 * t9, -0.2e1 * t13 * t10 + 0.2e1 * t2 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t40 - t46) * t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t27 + t28, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t34, t36, 0, -t20, -t46, -t31 * t20 + t29 * t39, t29 * t20 + t31 * t39, t40, -pkin(3) * t20 + t40 * qJ(4), -t10 * t16, t10 * t15 - t16 * t9, -t51, t47, 0, t13 * t15 + t24 * t9 - t7 * t36, -t24 * t10 + t13 * t16 + t8 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t34, t44, -t45, t42 * t34, t34 * t41 + t52, 0, 0, 0, 0, 0, -t47, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t31, -0.2e1 * pkin(3) * t29, 0.2e1 * t41, t42 * qJ(4) ^ 2 + pkin(3) ^ 2, t16 ^ 2, -0.2e1 * t16 * t15, 0, 0, 0, t15 * t56, t16 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t48, 0, t20, 0, 0, 0, 0, 0, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t29, 0, -pkin(3), 0, 0, 0, 0, 0, t15, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t9, -t36, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t11;
