% Calculate inertial parameters regressor of gravitation load for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:54
% EndTime: 2019-12-05 15:57:55
% DurationCPUTime: 0.30s
% Computational Cost: add. (269->77), mult. (485->123), div. (0->0), fcn. (582->12), ass. (0->46)
t28 = sin(pkin(5));
t54 = g(3) * t28;
t27 = sin(pkin(9));
t32 = sin(qJ(2));
t34 = cos(qJ(2));
t42 = cos(pkin(9));
t43 = cos(pkin(5));
t38 = t43 * t42;
t14 = t27 * t32 - t34 * t38;
t15 = t27 * t34 + t32 * t38;
t29 = cos(pkin(10));
t22 = t29 * pkin(3) + pkin(2);
t30 = -pkin(7) - qJ(3);
t53 = -t14 * t22 - t15 * t30;
t41 = t27 * t43;
t16 = t42 * t32 + t34 * t41;
t17 = -t32 * t41 + t42 * t34;
t52 = -t16 * t22 - t17 * t30;
t25 = pkin(10) + qJ(4);
t24 = cos(t25);
t31 = sin(qJ(5));
t51 = t24 * t31;
t33 = cos(qJ(5));
t50 = t24 * t33;
t49 = t27 * t28;
t48 = t28 * t32;
t47 = t28 * t34;
t46 = t30 * t32;
t45 = t31 * t34;
t44 = t33 * t34;
t40 = t28 * t42;
t23 = sin(t25);
t39 = pkin(4) * t24 + pkin(8) * t23;
t10 = -t23 * t48 + t43 * t24;
t4 = -t15 * t23 - t24 * t40;
t6 = -t17 * t23 + t24 * t49;
t37 = g(1) * t6 + g(2) * t4 + g(3) * t10;
t11 = t43 * t23 + t24 * t48;
t5 = t15 * t24 - t23 * t40;
t7 = t17 * t24 + t23 * t49;
t36 = g(1) * t7 + g(2) * t5 + g(3) * t11;
t2 = -g(1) * t16 - g(2) * t14 + g(3) * t47;
t35 = g(1) * t17 + g(2) * t15 + g(3) * t48;
t18 = t22 * t47;
t1 = t2 * t23;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t35, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t29, t2 * sin(pkin(10)), -t35, -g(1) * (-t16 * pkin(2) + t17 * qJ(3)) - g(2) * (-t14 * pkin(2) + t15 * qJ(3)) - (pkin(2) * t34 + qJ(3) * t32) * t54, 0, 0, 0, 0, 0, 0, -t2 * t24, t1, -t35, -g(1) * t52 - g(2) * t53 - g(3) * (-t28 * t46 + t18), 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t50 + t17 * t31) - g(2) * (-t14 * t50 + t15 * t31) - (t24 * t44 + t31 * t32) * t54, -g(1) * (t16 * t51 + t17 * t33) - g(2) * (t14 * t51 + t15 * t33) - (-t24 * t45 + t32 * t33) * t54, -t1, -g(1) * (-t39 * t16 + t52) - g(2) * (-t39 * t14 + t53) - g(3) * t18 - (t39 * t34 - t46) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t36, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t33, t37 * t31, -t36, -g(1) * (t6 * pkin(4) + t7 * pkin(8)) - g(2) * (t4 * pkin(4) + t5 * pkin(8)) - g(3) * (t10 * pkin(4) + t11 * pkin(8)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t16 * t33 - t7 * t31) - g(2) * (t14 * t33 - t5 * t31) - g(3) * (-t11 * t31 - t28 * t44), -g(1) * (-t16 * t31 - t7 * t33) - g(2) * (-t14 * t31 - t5 * t33) - g(3) * (-t11 * t33 + t28 * t45), 0, 0;];
taug_reg = t3;
