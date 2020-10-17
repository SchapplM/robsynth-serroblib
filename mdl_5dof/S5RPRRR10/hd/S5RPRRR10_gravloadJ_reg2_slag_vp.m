% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:44
% EndTime: 2019-12-31 19:10:46
% DurationCPUTime: 0.31s
% Computational Cost: add. (239->65), mult. (280->93), div. (0->0), fcn. (283->10), ass. (0->46)
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t15 = g(1) * t31 + g(2) * t29;
t23 = pkin(9) + qJ(3);
t19 = cos(t23);
t28 = sin(qJ(4));
t42 = t31 * t28;
t30 = cos(qJ(4));
t45 = t29 * t30;
t11 = -t19 * t42 + t45;
t18 = sin(t23);
t50 = g(3) * t18;
t41 = t31 * t30;
t46 = t29 * t28;
t9 = t19 * t46 + t41;
t56 = -g(1) * t11 + g(2) * t9 + t28 * t50;
t34 = -g(3) * t19 + t15 * t18;
t26 = cos(pkin(9));
t16 = t26 * pkin(2) + pkin(1);
t13 = t31 * t16;
t52 = g(2) * t13;
t24 = qJ(4) + qJ(5);
t20 = sin(t24);
t48 = t29 * t20;
t21 = cos(t24);
t47 = t29 * t21;
t44 = t31 * t20;
t43 = t31 * t21;
t27 = -pkin(6) - qJ(2);
t39 = pkin(4) * t28 - t27;
t38 = t19 * pkin(3) + t18 * pkin(7);
t14 = g(1) * t29 - g(2) * t31;
t17 = t30 * pkin(4) + pkin(3);
t32 = -pkin(8) - pkin(7);
t36 = t19 * t17 - t18 * t32;
t12 = t19 * t41 + t46;
t10 = -t19 * t45 + t42;
t8 = t14 * t18;
t7 = t19 * t43 + t48;
t6 = -t19 * t44 + t47;
t5 = -t19 * t47 + t44;
t4 = t19 * t48 + t43;
t3 = t15 * t19 + t50;
t2 = g(1) * t7 - g(2) * t5 + t21 * t50;
t1 = -g(1) * t6 + g(2) * t4 + t20 * t50;
t22 = [0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t26, -t14 * sin(pkin(9)), -t15, -g(1) * (-t29 * pkin(1) + t31 * qJ(2)) - g(2) * (t31 * pkin(1) + t29 * qJ(2)), 0, 0, 0, 0, 0, 0, t14 * t19, -t8, -t15, -g(1) * (-t29 * t16 - t31 * t27) - g(2) * (-t29 * t27 + t13), 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t8, -t52 + (g(1) * t27 - g(2) * t38) * t31 + (-g(1) * (-t16 - t38) + g(2) * t27) * t29, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, -g(1) * t4 - g(2) * t6, t8, -t52 + (-g(1) * t39 - g(2) * t36) * t31 + (-g(1) * (-t16 - t36) - g(2) * t39) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t3, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t30, -t34 * t28, -t3, -g(3) * t38 + t15 * (pkin(3) * t18 - pkin(7) * t19), 0, 0, 0, 0, 0, 0, t34 * t21, -t34 * t20, -t3, -g(3) * t36 + t15 * (t17 * t18 + t19 * t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, g(1) * t12 - g(2) * t10 + t30 * t50, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t56 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t22;
