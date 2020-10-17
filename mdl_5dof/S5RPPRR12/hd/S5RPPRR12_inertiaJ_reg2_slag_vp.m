% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR12_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:21
% EndTime: 2019-12-31 18:07:23
% DurationCPUTime: 0.51s
% Computational Cost: add. (386->54), mult. (671->101), div. (0->0), fcn. (754->6), ass. (0->49)
t26 = sin(pkin(8));
t27 = cos(pkin(8));
t49 = cos(qJ(4));
t39 = t49 * t27;
t48 = sin(qJ(4));
t14 = -t48 * t26 + t39;
t57 = -0.2e1 * t14;
t10 = t14 ^ 2;
t38 = t48 * t27;
t12 = t49 * t26 + t38;
t9 = t12 ^ 2;
t56 = t9 + t10;
t28 = -pkin(1) - qJ(3);
t50 = -pkin(6) + t28;
t16 = t50 * t26;
t4 = t48 * t16 - t50 * t39;
t55 = t4 ^ 2;
t18 = t26 * pkin(3) + qJ(2);
t54 = 0.2e1 * t18;
t53 = 0.2e1 * qJ(2);
t52 = t14 * pkin(4);
t51 = t4 * t14;
t29 = sin(qJ(5));
t47 = t29 * t12;
t46 = t29 * t14;
t30 = cos(qJ(5));
t45 = t29 * t30;
t44 = t30 * t14;
t21 = t26 ^ 2;
t22 = t27 ^ 2;
t17 = t21 + t22;
t24 = t29 ^ 2;
t25 = t30 ^ 2;
t43 = t24 + t25;
t42 = t12 * t57;
t41 = t29 * t44;
t37 = t43 * t12;
t36 = -pkin(7) * t12 - t52;
t3 = t12 * pkin(4) - t14 * pkin(7) + t18;
t6 = t49 * t16 + t50 * t38;
t1 = -t29 * t6 + t30 * t3;
t2 = t29 * t3 + t30 * t6;
t35 = t1 * t30 + t2 * t29;
t34 = -t1 * t29 + t2 * t30;
t33 = t6 * t12 - t51;
t31 = qJ(2) ^ 2;
t11 = t17 * t28;
t7 = t30 * t12;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t53, (pkin(1) ^ 2) + t31, t22, -0.2e1 * t27 * t26, 0, t21, 0, 0, t26 * t53, t27 * t53, -0.2e1 * t11, t17 * t28 ^ 2 + t31, t10, t42, 0, t9, 0, 0, t12 * t54, t14 * t54, -0.2e1 * t33, t18 ^ 2 + t6 ^ 2 + t55, t25 * t10, -0.2e1 * t10 * t45, 0.2e1 * t12 * t44, t24 * t10, t29 * t42, t9, 0.2e1 * t1 * t12 + 0.2e1 * t4 * t46, -0.2e1 * t2 * t12 + 0.2e1 * t4 * t44, t35 * t57, t1 ^ 2 + t2 ^ 2 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t17, t11, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t33, 0, 0, 0, 0, 0, 0, -t56 * t29, -t56 * t30, 0, t12 * t34 - t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 * t9 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t27, 0, qJ(2), 0, 0, 0, 0, 0, 0, t12, t14, 0, t18, 0, 0, 0, 0, 0, 0, t7, -t47, -t43 * t14, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, -t12, 0, -t4, -t6, 0, 0, t41, (-t24 + t25) * t14, t47, -t41, t7, 0, t36 * t29 - t4 * t30, t4 * t29 + t30 * t36, t34, -t4 * pkin(4) + pkin(7) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t12, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t46, t37, pkin(7) * t37 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t24, 0.2e1 * t45, 0, t25, 0, 0, 0.2e1 * pkin(4) * t30, -0.2e1 * pkin(4) * t29, 0.2e1 * t43 * pkin(7), pkin(7) ^ 2 * t43 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t46, t12, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t30, 0, -t29 * pkin(7), -t30 * pkin(7), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
