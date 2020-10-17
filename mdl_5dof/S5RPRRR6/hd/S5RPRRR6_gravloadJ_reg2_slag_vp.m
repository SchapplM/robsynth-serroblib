% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR6
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:40
% EndTime: 2019-12-31 19:01:41
% DurationCPUTime: 0.25s
% Computational Cost: add. (244->61), mult. (219->84), div. (0->0), fcn. (210->10), ass. (0->43)
t24 = qJ(3) + qJ(4);
t19 = sin(t24);
t20 = cos(t24);
t37 = t20 * pkin(4) + t19 * pkin(8);
t23 = qJ(1) + pkin(9);
t17 = sin(t23);
t18 = cos(t23);
t10 = g(1) * t18 + g(2) * t17;
t3 = -g(3) * t20 + t10 * t19;
t46 = pkin(4) * t19;
t45 = g(3) * t19;
t27 = sin(qJ(1));
t43 = t27 * pkin(1);
t42 = t18 * t19;
t41 = t18 * t20;
t25 = sin(qJ(5));
t40 = t20 * t25;
t28 = cos(qJ(5));
t39 = t20 * t28;
t29 = cos(qJ(3));
t21 = t29 * pkin(3);
t16 = t21 + pkin(2);
t30 = cos(qJ(1));
t22 = t30 * pkin(1);
t38 = t18 * t16 + t22;
t26 = sin(qJ(3));
t36 = -pkin(3) * t26 - t46;
t35 = g(1) * t17 - g(2) * t18;
t34 = g(1) * t27 - g(2) * t30;
t31 = -pkin(7) - pkin(6);
t33 = -t18 * t31 - t43;
t32 = -g(3) * t29 + t10 * t26;
t13 = pkin(8) * t41;
t12 = t17 * t20 * pkin(8);
t9 = t17 * t25 + t18 * t39;
t8 = t17 * t28 - t18 * t40;
t7 = -t17 * t39 + t18 * t25;
t6 = t17 * t40 + t18 * t28;
t5 = t35 * t19;
t4 = t10 * t20 + t45;
t2 = t3 * t28;
t1 = t3 * t25;
t11 = [0, 0, 0, 0, 0, 0, t34, g(1) * t30 + g(2) * t27, 0, 0, 0, 0, 0, 0, 0, 0, t35, t10, 0, t34 * pkin(1), 0, 0, 0, 0, 0, 0, t35 * t29, -t35 * t26, -t10, -g(1) * (-t17 * pkin(2) + t18 * pkin(6) - t43) - g(2) * (t18 * pkin(2) + t17 * pkin(6) + t22), 0, 0, 0, 0, 0, 0, t35 * t20, -t5, -t10, -g(1) * (-t17 * t16 + t33) - g(2) * (-t17 * t31 + t38), 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t5, -g(1) * t33 - g(2) * (pkin(4) * t41 + pkin(8) * t42 + t38) + (-g(1) * (-t16 - t37) + g(2) * t31) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, g(3) * t26 + t10 * t29, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t32 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t36 * t18 + t13) - g(2) * (t36 * t17 + t12) - g(3) * (t21 + t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-pkin(4) * t42 + t13) - g(2) * (-t17 * t46 + t12) - g(3) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t6 + t25 * t45, g(1) * t9 - g(2) * t7 + t28 * t45, 0, 0;];
taug_reg = t11;
