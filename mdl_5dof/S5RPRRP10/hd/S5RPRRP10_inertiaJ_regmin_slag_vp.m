% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:14:48
% EndTime: 2021-01-15 19:14:52
% DurationCPUTime: 0.34s
% Computational Cost: add. (388->62), mult. (767->124), div. (0->0), fcn. (869->6), ass. (0->46)
t27 = sin(pkin(8));
t28 = cos(pkin(8));
t30 = sin(qJ(3));
t45 = cos(qJ(3));
t15 = t45 * t27 + t30 * t28;
t50 = -0.2e1 * t15;
t21 = -t28 * pkin(2) - pkin(1);
t49 = 0.2e1 * t21;
t14 = t30 * t27 - t45 * t28;
t48 = t14 * pkin(4);
t31 = cos(qJ(4));
t47 = t31 * pkin(4);
t42 = pkin(6) + qJ(2);
t16 = t42 * t27;
t17 = t42 * t28;
t9 = -t30 * t16 + t45 * t17;
t46 = t31 * t9;
t29 = sin(qJ(4));
t10 = t29 * t14;
t44 = t29 * t15;
t43 = t29 * t31;
t12 = t31 * t15;
t41 = -qJ(5) - pkin(7);
t40 = t27 ^ 2 + t28 ^ 2;
t25 = t29 ^ 2;
t26 = t31 ^ 2;
t39 = t25 + t26;
t38 = qJ(5) * t15;
t37 = t14 * t50;
t7 = t14 * pkin(3) - t15 * pkin(7) + t21;
t3 = -t29 * t9 + t31 * t7;
t36 = -pkin(3) * t15 - pkin(7) * t14;
t33 = -t31 * t38 + t3;
t1 = t33 + t48;
t2 = t46 + (t7 - t38) * t29;
t35 = t1 * t31 + t2 * t29;
t18 = t41 * t29;
t19 = t41 * t31;
t34 = t31 * t18 - t29 * t19;
t8 = t45 * t16 + t30 * t17;
t22 = -pkin(3) - t47;
t13 = t15 ^ 2;
t11 = t31 * t14;
t5 = pkin(4) * t44 + t8;
t4 = t29 * t7 + t46;
t6 = [1, 0, 0, 0.2e1 * pkin(1) * t28, 0.2e1 * t40 * qJ(2), t40 * qJ(2) ^ 2 + pkin(1) ^ 2, t13, t37, 0, 0, 0, t14 * t49, t15 * t49, t26 * t13, -0.2e1 * t13 * t43, 0.2e1 * t14 * t12, t29 * t37, t14 ^ 2, 0.2e1 * t3 * t14 + 0.2e1 * t8 * t44, 0.2e1 * t8 * t12 - 0.2e1 * t4 * t14, 0.2e1 * t1 * t14 + 0.2e1 * t5 * t44, 0.2e1 * t5 * t12 - 0.2e1 * t2 * t14, t35 * t50, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, -t28, 0, -pkin(1), 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t11, -t10, t11, -t10, -t39 * t15, t35; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, -t8, -t9, t29 * t12, (-t25 + t26) * t15, t10, t11, 0, t36 * t29 - t8 * t31, t8 * t29 + t36 * t31, t18 * t14 + t22 * t44 - t5 * t31, t22 * t12 + t19 * t14 + t5 * t29, -t1 * t29 - t34 * t15 + t2 * t31, t1 * t18 - t2 * t19 + t5 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t25, 0.2e1 * t43, 0, 0, 0, 0.2e1 * pkin(3) * t31, -0.2e1 * pkin(3) * t29, -0.2e1 * t22 * t31, 0.2e1 * t22 * t29, -0.2e1 * t18 * t29 - 0.2e1 * t19 * t31, t18 ^ 2 + t19 ^ 2 + t22 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t44, t14, t3, -t4, t33 + 0.2e1 * t48, -t2, -pkin(4) * t12, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, t31, -t29, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t31, 0, -t29 * pkin(7), -t31 * pkin(7), t18, t19, -t29 * pkin(4), t18 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t12, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t29, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
