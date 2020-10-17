% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPPR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:09
% EndTime: 2019-12-05 15:22:12
% DurationCPUTime: 0.42s
% Computational Cost: add. (235->53), mult. (546->102), div. (0->0), fcn. (578->6), ass. (0->43)
t28 = sin(pkin(9));
t30 = cos(pkin(9));
t32 = sin(qJ(5));
t33 = cos(qJ(5));
t14 = -t32 * t28 + t33 * t30;
t29 = sin(pkin(8));
t48 = -0.2e1 * t29;
t31 = cos(pkin(8));
t47 = -0.2e1 * t31;
t46 = 0.2e1 * t31;
t45 = t28 * t29;
t44 = t30 * t29;
t43 = t30 * t31;
t17 = -t31 * pkin(3) - t29 * qJ(4) - pkin(2);
t38 = qJ(3) * t31;
t6 = t28 * t17 + t30 * t38;
t24 = t28 ^ 2;
t26 = t30 ^ 2;
t40 = t24 + t26;
t25 = t29 ^ 2;
t27 = t31 ^ 2;
t39 = t25 + t27;
t37 = t25 * qJ(3);
t36 = t29 * t46;
t13 = t30 * t17;
t5 = -t28 * t38 + t13;
t35 = t6 * t28 + t5 * t30;
t15 = t33 * t28 + t32 * t30;
t34 = qJ(3) ^ 2;
t23 = t29 * qJ(3);
t22 = t25 * t34;
t21 = t26 * t25;
t20 = t24 * t25;
t16 = pkin(4) * t45 + t23;
t11 = t14 * t29;
t9 = t15 * t29;
t8 = t11 ^ 2;
t7 = t9 ^ 2;
t4 = -pkin(6) * t45 + t6;
t3 = -pkin(6) * t44 + t13 + (-qJ(3) * t28 - pkin(4)) * t31;
t2 = t32 * t3 + t33 * t4;
t1 = t33 * t3 - t32 * t4;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 + t20 + t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 + t7 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t28 * t5 + t30 * t6 - t38) * t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t1 + t11 * t2 - t31 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t25, t36, 0, t27, 0, 0, pkin(2) * t46, pkin(2) * t48, 0.2e1 * t39 * qJ(3), pkin(2) ^ 2 + t27 * t34 + t22, t21, -0.2e1 * t30 * t25 * t28, t43 * t48, t20, t28 * t36, t27, 0.2e1 * t28 * t37 - 0.2e1 * t5 * t31, 0.2e1 * t30 * t37 + 0.2e1 * t6 * t31, t35 * t48, t5 ^ 2 + t6 ^ 2 + t22, t8, -0.2e1 * t11 * t9, t11 * t47, t7, -t9 * t47, t27, -0.2e1 * t1 * t31 + 0.2e1 * t16 * t9, 0.2e1 * t16 * t11 + 0.2e1 * t2 * t31, -0.2e1 * t1 * t11 - 0.2e1 * t2 * t9, t1 ^ 2 + t16 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t15 - t9 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t29, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t43, t28 * t31, -t40 * t29, t35, 0, 0, 0, 0, 0, 0, -t14 * t31, t15 * t31, -t14 * t11 - t15 * t9, t1 * t14 + t2 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t44, 0, t23, 0, 0, 0, 0, 0, 0, t9, t11, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, -t31, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t10;
