% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPPR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:51
% EndTime: 2019-12-31 17:47:53
% DurationCPUTime: 0.50s
% Computational Cost: add. (307->62), mult. (613->123), div. (0->0), fcn. (622->6), ass. (0->52)
t35 = sin(pkin(7));
t29 = t35 ^ 2;
t37 = cos(pkin(7));
t31 = t37 ^ 2;
t62 = t29 + t31;
t61 = -0.2e1 * t35;
t60 = 0.2e1 * t37;
t44 = -t35 * qJ(3) - pkin(1);
t12 = (-pkin(2) - qJ(4)) * t37 + t44;
t24 = t35 * qJ(2);
t16 = pkin(3) * t35 + t24;
t34 = sin(pkin(8));
t36 = cos(pkin(8));
t6 = t36 * t12 + t34 * t16;
t30 = t36 ^ 2;
t59 = t30 * t37;
t58 = t34 * t35;
t57 = t34 * t37;
t56 = t35 * t37;
t55 = t36 * t34;
t54 = t36 * t35;
t21 = t36 * t37;
t38 = sin(qJ(5));
t53 = t38 * t34;
t52 = t38 * t36;
t39 = cos(qJ(5));
t51 = t39 * t34;
t50 = t39 * t36;
t49 = t62 * qJ(2) ^ 2;
t17 = (pkin(3) + qJ(2)) * t37;
t28 = t34 ^ 2;
t18 = t30 + t28;
t48 = t38 ^ 2 + t39 ^ 2;
t47 = -0.2e1 * t21;
t46 = t37 * t53;
t45 = t37 * t51;
t4 = pkin(6) * t35 + t6;
t7 = (pkin(4) * t36 + pkin(6) * t34) * t37 + t17;
t1 = -t38 * t4 + t39 * t7;
t2 = t38 * t7 + t39 * t4;
t43 = -t1 * t38 + t2 * t39;
t5 = -t12 * t34 + t16 * t36;
t42 = t34 * t5 - t36 * t6;
t10 = -t35 * t38 + t45;
t9 = t35 * t39 + t46;
t41 = -t10 * t38 + t39 * t9;
t20 = t30 * t31;
t19 = 0.2e1 * t56;
t15 = -pkin(2) * t37 + t44;
t14 = 0.2e1 * t62 * qJ(2);
t3 = -pkin(4) * t35 - t5;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t29, t19, 0, t31, 0, 0, pkin(1) * t60, pkin(1) * t61, t14, pkin(1) ^ 2 + t49, 0, 0, 0, t29, t19, t31, t14, t15 * t60, t15 * t61, t15 ^ 2 + t49, t28 * t31, 0.2e1 * t31 * t55, -0.2e1 * t34 * t56, t20, t35 * t47, t29, 0.2e1 * t17 * t21 + 0.2e1 * t35 * t5, -0.2e1 * t17 * t57 - 0.2e1 * t35 * t6, t42 * t60, t17 ^ 2 + t5 ^ 2 + t6 ^ 2, t10 ^ 2, -0.2e1 * t10 * t9, t10 * t47, t9 ^ 2, 0.2e1 * t9 * t21, t20, 0.2e1 * t1 * t21 - 0.2e1 * t3 * t9, -0.2e1 * t10 * t3 - 0.2e1 * t2 * t21, 0.2e1 * t1 * t10 + 0.2e1 * t2 * t9, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t35, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, t37, -t35, t15, 0, 0, 0, 0, 0, 0, -t58, -t54, -t18 * t37, -t42, 0, 0, 0, 0, 0, 0, -t34 * t9 - t38 * t59, -t10 * t34 - t39 * t59, t41 * t36, t3 * t34 + t36 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t48 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, t24, 0, 0, 0, 0, 0, 0, t54, -t58, 0, t34 * t6 + t36 * t5, 0, 0, 0, 0, 0, 0, (t9 - t46) * t36, (t10 - t45) * t36, t41 * t34, -t3 * t36 + t34 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t48) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t48 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t57, 0, t17, 0, 0, 0, 0, 0, 0, t37 * t50, -t37 * t52, t10 * t39 + t38 * t9, t1 * t39 + t2 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, t9, t21, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t51, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
