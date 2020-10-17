% Calculate inertial parameters regressor of gravitation load for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:42
% EndTime: 2019-12-31 17:42:43
% DurationCPUTime: 0.16s
% Computational Cost: add. (195->43), mult. (182->56), div. (0->0), fcn. (171->10), ass. (0->38)
t20 = qJ(2) + qJ(3);
t16 = pkin(9) + t20;
t13 = sin(t16);
t14 = cos(t16);
t21 = sin(pkin(8));
t22 = cos(pkin(8));
t28 = g(1) * t22 + g(2) * t21;
t3 = -g(3) * t14 + t28 * t13;
t17 = sin(t20);
t40 = pkin(3) * t17;
t39 = pkin(4) * t13;
t38 = pkin(7) * t14;
t37 = g(3) * t13;
t23 = sin(qJ(5));
t35 = t21 * t23;
t25 = cos(qJ(5));
t34 = t21 * t25;
t33 = t22 * t23;
t32 = t22 * t25;
t18 = cos(t20);
t15 = pkin(3) * t18;
t31 = t14 * pkin(4) + t13 * pkin(7) + t15;
t24 = sin(qJ(2));
t7 = -t24 * pkin(2) - t40;
t30 = t7 - t39;
t29 = -t39 - t40;
t5 = -g(3) * t18 + t28 * t17;
t26 = cos(qJ(2));
t27 = -g(3) * t26 + t28 * t24;
t19 = t26 * pkin(2);
t10 = -g(1) * t21 + g(2) * t22;
t9 = t22 * t38;
t8 = t21 * t38;
t6 = g(3) * t17 + t28 * t18;
t4 = t28 * t14 + t37;
t2 = t3 * t25;
t1 = t3 * t23;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, g(3) * t24 + t28 * t26, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t27 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * (t15 + t19) - t28 * t7, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t30 * t22 + t9) - g(2) * (t30 * t21 + t8) - g(3) * (t19 + t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t29 * t22 + t9) - g(2) * (t29 * t21 + t8) - g(3) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t33 + t34) - g(2) * (-t14 * t35 - t32) + t23 * t37, -g(1) * (-t14 * t32 - t35) - g(2) * (-t14 * t34 + t33) + t25 * t37, 0, 0;];
taug_reg = t11;
