% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRP5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:41
% EndTime: 2019-12-05 15:38:43
% DurationCPUTime: 0.43s
% Computational Cost: add. (182->47), mult. (403->74), div. (0->0), fcn. (463->6), ass. (0->37)
t25 = sin(pkin(8));
t26 = cos(pkin(8));
t27 = sin(qJ(4));
t44 = cos(qJ(4));
t17 = t25 * t44 + t27 * t26;
t20 = -t26 * pkin(3) - pkin(2);
t31 = -t27 * t25 + t26 * t44;
t3 = -pkin(4) * t31 - t17 * qJ(5) + t20;
t48 = -0.2e1 * t3;
t47 = t31 ^ 2;
t46 = 0.2e1 * t20;
t45 = 0.2e1 * t26;
t43 = t17 * t31;
t29 = cos(qJ(2));
t42 = t29 * t31;
t41 = t29 * t17;
t40 = pkin(6) + qJ(3);
t21 = t25 ^ 2;
t22 = t26 ^ 2;
t39 = t21 + t22;
t18 = t40 * t26;
t34 = t40 * t25;
t6 = t27 * t18 + t34 * t44;
t8 = t18 * t44 - t27 * t34;
t38 = t6 ^ 2 + t8 ^ 2;
t28 = sin(qJ(2));
t10 = t17 * t28;
t12 = t31 * t28;
t37 = t10 * t6 + t12 * t8;
t36 = t10 * t17 + t12 * t31;
t24 = t29 ^ 2;
t35 = t10 ^ 2 + t12 ^ 2 + t24;
t33 = t39 * qJ(3);
t32 = 0.2e1 * t6 * t17 + 0.2e1 * t31 * t8;
t23 = t28 ^ 2;
t13 = t17 ^ 2;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 + t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t39 + t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t26, -t29 * t25, t39 * t28, t29 * pkin(2) + t28 * t33, 0, 0, 0, 0, 0, 0, t42, -t41, t36, -t29 * t20 + t37, 0, 0, 0, 0, 0, 0, t42, t36, t41, -t29 * t3 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t21, t25 * t45, 0, t22, 0, 0, pkin(2) * t45, -0.2e1 * pkin(2) * t25, 0.2e1 * t33, qJ(3) ^ 2 * t39 + pkin(2) ^ 2, t13, 0.2e1 * t43, 0, t47, 0, 0, -t31 * t46, t17 * t46, t32, t20 ^ 2 + t38, t13, 0, -0.2e1 * t43, 0, 0, t47, t31 * t48, t32, t17 * t48, t3 ^ 2 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t25, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t31, t17, 0, t20, 0, 0, 0, 0, 0, 0, -t31, 0, -t17, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t12, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, t12, -t10 * pkin(4) + t12 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, t31, 0, -t6, -t8, 0, 0, 0, t17, 0, 0, -t31, 0, -t6, -pkin(4) * t17 + qJ(5) * t31, t8, -t6 * pkin(4) + t8 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
