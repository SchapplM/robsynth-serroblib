% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:57
% EndTime: 2019-12-31 19:54:58
% DurationCPUTime: 0.58s
% Computational Cost: add. (543->60), mult. (1050->116), div. (0->0), fcn. (1198->6), ass. (0->43)
t30 = sin(pkin(8));
t31 = cos(pkin(8));
t33 = sin(qJ(2));
t34 = cos(qJ(2));
t21 = t30 * t33 - t31 * t34;
t23 = t30 * t34 + t31 * t33;
t32 = sin(qJ(4));
t48 = cos(qJ(4));
t10 = t48 * t21 + t32 * t23;
t54 = t10 ^ 2;
t27 = -t34 * pkin(2) - pkin(1);
t15 = t21 * pkin(3) + t27;
t53 = 0.2e1 * t15;
t52 = 0.2e1 * t27;
t51 = 0.2e1 * t34;
t50 = t30 * pkin(2);
t49 = t31 * pkin(2);
t12 = -t32 * t21 + t48 * t23;
t47 = t12 * t10;
t46 = -qJ(3) - pkin(6);
t28 = t33 ^ 2;
t29 = t34 ^ 2;
t45 = t28 + t29;
t42 = t46 * t33;
t43 = t46 * t34;
t13 = t30 * t43 + t31 * t42;
t39 = -t23 * pkin(7) + t13;
t14 = t30 * t42 - t31 * t43;
t8 = -t21 * pkin(7) + t14;
t4 = t32 * t8 - t48 * t39;
t6 = t32 * t39 + t48 * t8;
t44 = t4 ^ 2 + t6 ^ 2;
t26 = pkin(3) + t49;
t41 = -t48 * t26 + t32 * t50;
t40 = -0.2e1 * t6 * t10 + 0.2e1 * t4 * t12;
t19 = t32 * t26 + t48 * t50;
t36 = 2 * pkin(4);
t35 = 2 * qJ(5);
t17 = -pkin(4) + t41;
t16 = qJ(5) + t19;
t9 = t12 ^ 2;
t2 = t10 * pkin(4) - t12 * qJ(5) + t15;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t28, t33 * t51, 0, t29, 0, 0, pkin(1) * t51, -0.2e1 * pkin(1) * t33, 0.2e1 * t45 * pkin(6), t45 * pkin(6) ^ 2 + pkin(1) ^ 2, t23 ^ 2, -0.2e1 * t23 * t21, 0, t21 ^ 2, 0, 0, t21 * t52, t23 * t52, -0.2e1 * t13 * t23 - 0.2e1 * t14 * t21, t13 ^ 2 + t14 ^ 2 + t27 ^ 2, t9, -0.2e1 * t47, 0, t54, 0, 0, t10 * t53, t12 * t53, t40, t15 ^ 2 + t44, t9, 0, 0.2e1 * t47, 0, 0, t54, 0.2e1 * t2 * t10, t40, -0.2e1 * t2 * t12, t2 ^ 2 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, t34, 0, -t33 * pkin(6), -t34 * pkin(6), 0, 0, 0, 0, t23, 0, -t21, 0, t13, -t14, (-t21 * t30 - t23 * t31) * pkin(2), (t13 * t31 + t14 * t30) * pkin(2), 0, 0, t12, 0, -t10, 0, -t4, -t6, -t19 * t10 + t12 * t41, t6 * t19 + t4 * t41, 0, t12, 0, 0, t10, 0, -t4, -t16 * t10 + t17 * t12, t6, t6 * t16 + t4 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t49, -0.2e1 * t50, 0, (t30 ^ 2 + t31 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t41, -0.2e1 * t19, 0, t19 ^ 2 + t41 ^ 2, 0, 0, 0, 1, 0, 0, -0.2e1 * t17, 0, 0.2e1 * t16, t16 ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t23, 0, t27, 0, 0, 0, 0, 0, 0, t10, t12, 0, t15, 0, 0, 0, 0, 0, 0, t10, 0, -t12, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, 0, -t4, -t6, 0, 0, 0, t12, 0, 0, t10, 0, -t4, -pkin(4) * t12 - t10 * qJ(5), t6, -t4 * pkin(4) + t6 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t41, -t19, 0, 0, 0, 0, 0, 1, 0, 0, t36 - t41, 0, t35 + t19, -t17 * pkin(4) + t16 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t36, 0, t35, pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
