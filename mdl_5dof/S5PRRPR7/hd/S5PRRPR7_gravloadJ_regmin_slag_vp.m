% Calculate minimal parameter regressor of gravitation load for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:34
% EndTime: 2019-12-05 16:37:35
% DurationCPUTime: 0.32s
% Computational Cost: add. (214->79), mult. (585->145), div. (0->0), fcn. (737->12), ass. (0->50)
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t43 = cos(pkin(9));
t44 = cos(pkin(5));
t37 = t44 * t43;
t42 = sin(pkin(9));
t12 = t42 * t28 - t31 * t37;
t36 = t44 * t42;
t14 = t43 * t28 + t31 * t36;
t57 = -g(1) * t14 - g(2) * t12;
t24 = sin(pkin(5));
t54 = g(3) * t24;
t23 = sin(pkin(10));
t30 = cos(qJ(3));
t53 = t23 * t30;
t52 = t24 * t28;
t51 = t24 * t31;
t25 = cos(pkin(10));
t26 = sin(qJ(5));
t50 = t25 * t26;
t29 = cos(qJ(5));
t49 = t25 * t29;
t48 = t25 * t30;
t27 = sin(qJ(3));
t47 = t26 * t27;
t46 = t27 * t29;
t45 = t30 * t31;
t41 = t27 * t51;
t40 = t24 * t43;
t39 = t24 * t42;
t13 = t28 * t37 + t42 * t31;
t15 = -t28 * t36 + t43 * t31;
t38 = -g(1) * t15 - g(2) * t13;
t16 = t27 * t52 - t44 * t30;
t7 = t13 * t27 + t30 * t40;
t9 = t15 * t27 - t30 * t39;
t34 = g(1) * t9 + g(2) * t7 + g(3) * t16;
t10 = t15 * t30 + t27 * t39;
t17 = t44 * t27 + t30 * t52;
t8 = t13 * t30 - t27 * t40;
t33 = g(1) * t10 + g(2) * t8 + g(3) * t17;
t32 = g(3) * t51 + t57;
t11 = (t23 * t28 + t25 * t45) * t24;
t6 = t17 * t25 - t23 * t51;
t5 = -t14 * t48 + t15 * t23;
t4 = -t12 * t48 + t13 * t23;
t3 = t32 * t27;
t2 = t10 * t25 + t14 * t23;
t1 = t12 * t23 + t8 * t25;
t18 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t32, g(3) * t52 - t38, 0, 0, 0, 0, 0, -t32 * t30, t3, -g(1) * t5 - g(2) * t4 - g(3) * t11, -g(1) * (t14 * t53 + t15 * t25) - g(2) * (t12 * t53 + t13 * t25) - (-t23 * t45 + t25 * t28) * t54, -t3, (-t28 * t54 + t38) * pkin(7) + (-t31 * t54 - t57) * (pkin(3) * t30 + qJ(4) * t27 + pkin(2)), 0, 0, 0, 0, 0, -g(1) * (-t14 * t47 + t5 * t29) - g(2) * (-t12 * t47 + t4 * t29) - g(3) * (t11 * t29 + t26 * t41), -g(1) * (-t14 * t46 - t5 * t26) - g(2) * (-t12 * t46 - t4 * t26) - g(3) * (-t11 * t26 + t29 * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t33, t34 * t25, -t34 * t23, -t33, -g(1) * (-t9 * pkin(3) + t10 * qJ(4)) - g(2) * (-t7 * pkin(3) + t8 * qJ(4)) - g(3) * (-t16 * pkin(3) + t17 * qJ(4)), 0, 0, 0, 0, 0, -g(1) * (t10 * t26 - t9 * t49) - g(2) * (t8 * t26 - t7 * t49) - g(3) * (-t16 * t49 + t17 * t26), -g(1) * (t10 * t29 + t9 * t50) - g(2) * (t8 * t29 + t7 * t50) - g(3) * (t16 * t50 + t17 * t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t2 * t26 + t9 * t29) - g(2) * (-t1 * t26 + t7 * t29) - g(3) * (t16 * t29 - t6 * t26), -g(1) * (-t2 * t29 - t9 * t26) - g(2) * (-t1 * t29 - t7 * t26) - g(3) * (-t16 * t26 - t6 * t29);];
taug_reg = t18;
