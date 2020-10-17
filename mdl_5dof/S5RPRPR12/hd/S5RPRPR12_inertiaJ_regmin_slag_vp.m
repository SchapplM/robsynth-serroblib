% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR12
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
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:15
% EndTime: 2019-12-31 18:30:16
% DurationCPUTime: 0.34s
% Computational Cost: add. (410->64), mult. (843->127), div. (0->0), fcn. (992->8), ass. (0->47)
t36 = sin(pkin(8));
t38 = cos(pkin(8));
t40 = sin(qJ(3));
t54 = cos(qJ(3));
t20 = t40 * t36 - t54 * t38;
t57 = -0.2e1 * t20;
t37 = cos(pkin(9));
t29 = -t37 * pkin(4) - pkin(3);
t56 = 0.2e1 * t29;
t30 = -t38 * pkin(2) - pkin(1);
t55 = 0.2e1 * t30;
t35 = sin(pkin(9));
t39 = sin(qJ(5));
t41 = cos(qJ(5));
t21 = t41 * t35 + t39 * t37;
t53 = t21 * t20;
t22 = t54 * t36 + t40 * t38;
t52 = t35 * t22;
t51 = t37 * t22;
t50 = pkin(6) + qJ(2);
t49 = pkin(7) + qJ(4);
t12 = t20 * pkin(3) - t22 * qJ(4) + t30;
t24 = t50 * t36;
t26 = t50 * t38;
t17 = -t40 * t24 + t54 * t26;
t6 = t35 * t12 + t37 * t17;
t48 = t35 ^ 2 + t37 ^ 2;
t47 = t36 ^ 2 + t38 ^ 2;
t5 = t37 * t12 - t35 * t17;
t46 = t6 * t35 + t5 * t37;
t45 = -t5 * t35 + t6 * t37;
t44 = -pkin(3) * t22 - qJ(4) * t20;
t19 = t39 * t35 - t41 * t37;
t15 = t54 * t24 + t40 * t26;
t25 = t49 * t37;
t23 = t49 * t35;
t18 = t19 * t20;
t16 = -t39 * t23 + t41 * t25;
t14 = -t41 * t23 - t39 * t25;
t9 = t19 * t22;
t8 = t21 * t22;
t7 = pkin(4) * t52 + t15;
t4 = -pkin(7) * t52 + t6;
t3 = t20 * pkin(4) - pkin(7) * t51 + t5;
t2 = t39 * t3 + t41 * t4;
t1 = t41 * t3 - t39 * t4;
t10 = [1, 0, 0, 0.2e1 * pkin(1) * t38, -0.2e1 * pkin(1) * t36, 0.2e1 * t47 * qJ(2), t47 * qJ(2) ^ 2 + pkin(1) ^ 2, t22 ^ 2, t22 * t57, 0, 0, 0, t20 * t55, t22 * t55, 0.2e1 * t15 * t52 + 0.2e1 * t5 * t20, 0.2e1 * t15 * t51 - 0.2e1 * t6 * t20, -0.2e1 * t46 * t22, t15 ^ 2 + t5 ^ 2 + t6 ^ 2, t9 ^ 2, 0.2e1 * t9 * t8, t9 * t57, t8 * t57, t20 ^ 2, 0.2e1 * t1 * t20 + 0.2e1 * t7 * t8, -0.2e1 * t2 * t20 - 0.2e1 * t7 * t9; 0, 0, 0, -t38, t36, 0, -pkin(1), 0, 0, 0, 0, 0, t20, t22, t37 * t20, -t35 * t20, -t48 * t22, t46, 0, 0, 0, 0, 0, -t18, -t53; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, 0, -t15, -t17, -t15 * t37 + t44 * t35, t15 * t35 + t44 * t37, t45, -t15 * pkin(3) + t45 * qJ(4), -t9 * t21, t9 * t19 - t21 * t8, t53, -t18, 0, t14 * t20 + t7 * t19 + t29 * t8, -t16 * t20 + t7 * t21 - t29 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t37, -0.2e1 * pkin(3) * t35, 0.2e1 * t48 * qJ(4), t48 * qJ(4) ^ 2 + pkin(3) ^ 2, t21 ^ 2, -0.2e1 * t21 * t19, 0, 0, 0, t19 * t56, t21 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t51, 0, t15, 0, 0, 0, 0, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t35, 0, -pkin(3), 0, 0, 0, 0, 0, t19, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t8, t20, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, 0, t14, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t10;
