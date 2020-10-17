% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPRR2
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:56
% EndTime: 2019-12-05 17:39:58
% DurationCPUTime: 0.44s
% Computational Cost: add. (379->43), mult. (631->76), div. (0->0), fcn. (754->6), ass. (0->42)
t32 = sin(qJ(5));
t34 = cos(qJ(5));
t29 = sin(pkin(8));
t30 = cos(pkin(8));
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t17 = t35 * t29 + t33 * t30;
t19 = -t33 * t29 + t35 * t30;
t38 = -t32 * t17 + t34 * t19;
t39 = t34 * t17 + t32 * t19;
t54 = (t32 * t39 + t34 * t38) * pkin(4);
t53 = t38 ^ 2;
t52 = t39 ^ 2;
t15 = t19 ^ 2;
t47 = t17 ^ 2;
t51 = t15 + t47;
t50 = t52 + t53;
t31 = -pkin(1) - qJ(3);
t41 = -pkin(6) + t31;
t20 = t41 * t29;
t21 = t41 * t30;
t10 = -t33 * t20 + t35 * t21;
t3 = -t19 * pkin(7) + t10;
t11 = t35 * t20 + t33 * t21;
t4 = -t17 * pkin(7) + t11;
t1 = t34 * t3 - t32 * t4;
t2 = t32 * t3 + t34 * t4;
t49 = t1 * t38 + t2 * t39;
t23 = t29 * pkin(3) + qJ(2);
t12 = t17 * pkin(4) + t23;
t46 = 0.2e1 * t12;
t45 = 0.2e1 * t23;
t44 = 0.2e1 * qJ(2);
t43 = t32 * pkin(4);
t42 = t34 * pkin(4);
t26 = t29 ^ 2;
t27 = t30 ^ 2;
t22 = t26 + t27;
t40 = t10 * t19 + t11 * t17;
t36 = qJ(2) ^ 2;
t16 = t22 * t31;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t44, (pkin(1) ^ 2) + t36, t27, -0.2e1 * t30 * t29, 0, t26, 0, 0, t29 * t44, t30 * t44, -0.2e1 * t16, t22 * t31 ^ 2 + t36, t15, -0.2e1 * t19 * t17, 0, t47, 0, 0, t17 * t45, t19 * t45, -0.2e1 * t40, t10 ^ 2 + t11 ^ 2 + t23 ^ 2, t53, -0.2e1 * t38 * t39, 0, t52, 0, 0, t39 * t46, t38 * t46, -0.2e1 * t49, t1 ^ 2 + t12 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t22, t16, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t40, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t30, 0, qJ(2), 0, 0, 0, 0, 0, 0, t17, t19, 0, t23, 0, 0, 0, 0, 0, 0, t39, t38, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, 0, t10, -t11, 0, 0, 0, 0, t38, 0, -t39, 0, t1, -t2, -t54, (t1 * t34 + t2 * t32) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t17, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t42, -0.2e1 * t43, 0, (t32 ^ 2 + t34 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, -t39, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t42, -t43, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
