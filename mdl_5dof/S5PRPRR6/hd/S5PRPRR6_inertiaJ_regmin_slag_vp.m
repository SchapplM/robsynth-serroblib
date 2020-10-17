% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:53
% EndTime: 2019-12-05 15:57:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (165->50), mult. (399->103), div. (0->0), fcn. (502->10), ass. (0->42)
t26 = cos(pkin(10));
t19 = -t26 * pkin(3) - pkin(2);
t46 = 0.2e1 * t19;
t45 = cos(qJ(4));
t25 = sin(pkin(5));
t44 = t25 * sin(qJ(2));
t32 = cos(qJ(2));
t43 = t25 * t32;
t24 = sin(pkin(10));
t29 = sin(qJ(4));
t14 = t29 * t24 - t45 * t26;
t28 = sin(qJ(5));
t42 = t28 * t14;
t15 = t45 * t24 + t29 * t26;
t41 = t28 * t15;
t31 = cos(qJ(5));
t40 = t28 * t31;
t39 = t31 * t15;
t38 = pkin(7) + qJ(3);
t37 = t24 ^ 2 + t26 ^ 2;
t36 = -0.2e1 * t15 * t14;
t35 = -pkin(4) * t15 - pkin(8) * t14;
t27 = cos(pkin(5));
t10 = -t24 * t44 + t27 * t26;
t11 = t27 * t24 + t26 * t44;
t34 = -t10 * t24 + t11 * t26;
t23 = t31 ^ 2;
t22 = t28 ^ 2;
t17 = t38 * t26;
t16 = t38 * t24;
t13 = t15 ^ 2;
t12 = t31 * t14;
t9 = -t29 * t16 + t45 * t17;
t8 = t45 * t16 + t29 * t17;
t7 = t14 * pkin(4) - t15 * pkin(8) + t19;
t6 = t29 * t10 + t45 * t11;
t5 = -t45 * t10 + t29 * t11;
t4 = -t28 * t43 + t31 * t6;
t3 = -t28 * t6 - t31 * t43;
t2 = t28 * t7 + t31 * t9;
t1 = -t28 * t9 + t31 * t7;
t18 = [1, 0, 0, 0, 0, 0, 0, t25 ^ 2 * t32 ^ 2 + t10 ^ 2 + t11 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t43, -t44, t26 * t43, -t24 * t43, t34, pkin(2) * t43 + t34 * qJ(3), 0, 0, 0, 0, 0, -t14 * t43, -t15 * t43, 0, 0, 0, 0, 0, t3 * t14 + t5 * t41, -t4 * t14 + t5 * t39; 0, 1, 0, 0, 0.2e1 * pkin(2) * t26, -0.2e1 * pkin(2) * t24, 0.2e1 * t37 * qJ(3), t37 * qJ(3) ^ 2 + pkin(2) ^ 2, t13, t36, 0, 0, 0, t14 * t46, t15 * t46, t23 * t13, -0.2e1 * t13 * t40, 0.2e1 * t14 * t39, t28 * t36, t14 ^ 2, 0.2e1 * t1 * t14 + 0.2e1 * t8 * t41, -0.2e1 * t2 * t14 + 0.2e1 * t8 * t39; 0, 0, 0, 0, 0, 0, 0, -t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t26, t24, 0, -pkin(2), 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t12, -t42; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t5 * t31, t5 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, -t8, -t9, t28 * t39, (-t22 + t23) * t15, t42, t12, 0, t35 * t28 - t8 * t31, t8 * t28 + t35 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t22, 0.2e1 * t40, 0, 0, 0, 0.2e1 * pkin(4) * t31, -0.2e1 * pkin(4) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t41, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t31, 0, -t28 * pkin(8), -t31 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t18;
