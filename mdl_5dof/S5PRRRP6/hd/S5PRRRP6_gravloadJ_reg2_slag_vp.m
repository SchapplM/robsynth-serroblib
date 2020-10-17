% Calculate inertial parameters regressor of gravitation load for
% S5PRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:18
% EndTime: 2019-12-05 16:52:19
% DurationCPUTime: 0.25s
% Computational Cost: add. (214->60), mult. (317->85), div. (0->0), fcn. (329->8), ass. (0->41)
t22 = sin(pkin(8));
t23 = cos(pkin(8));
t34 = g(1) * t23 + g(2) * t22;
t25 = sin(qJ(2));
t27 = cos(qJ(2));
t30 = -g(3) * t27 + t34 * t25;
t50 = g(3) * t25;
t21 = qJ(3) + qJ(4);
t19 = sin(t21);
t48 = t19 * t25;
t20 = cos(t21);
t47 = t20 * t25;
t26 = cos(qJ(3));
t46 = t22 * t26;
t45 = t22 * t27;
t44 = t23 * t26;
t43 = t23 * t27;
t24 = sin(qJ(3));
t42 = t24 * t27;
t41 = t26 * t27;
t28 = -pkin(7) - pkin(6);
t40 = t27 * t28;
t39 = t23 * t42;
t8 = t19 * t45 + t23 * t20;
t9 = -t23 * t19 + t20 * t45;
t38 = -t8 * pkin(4) + t9 * qJ(5);
t10 = t19 * t43 - t22 * t20;
t11 = t22 * t19 + t20 * t43;
t37 = -t10 * pkin(4) + t11 * qJ(5);
t18 = t26 * pkin(3) + pkin(2);
t36 = t27 * t18 - t25 * t28;
t33 = pkin(4) * t20 + qJ(5) * t19;
t31 = -t22 * t42 - t44;
t1 = g(1) * t10 + g(2) * t8 + g(3) * t48;
t3 = g(1) * t11 + g(2) * t9 + g(3) * t47;
t12 = t34 * t27 + t50;
t17 = pkin(3) * t46;
t13 = qJ(5) * t47;
t5 = t30 * t20;
t4 = t30 * t19;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t12, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t26, -t30 * t24, -t12, -g(3) * (t27 * pkin(2) + t25 * pkin(6)) + t34 * (pkin(2) * t25 - pkin(6) * t27), 0, 0, 0, 0, 0, 0, t5, -t4, -t12, -g(3) * t36 + t34 * (t18 * t25 + t40), 0, 0, 0, 0, 0, 0, t5, -t12, t4, -g(3) * (t33 * t27 + t36) + t34 * (t40 - (-t18 - t33) * t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t39 + t46) - g(2) * t31 + t24 * t50, -g(1) * (-t22 * t24 - t23 * t41) - g(2) * (-t22 * t41 + t23 * t24) + t26 * t50, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t17 + (g(2) * t44 + t12 * t24) * pkin(3), 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (-pkin(3) * t39 + t17 + t37) - g(2) * (t31 * pkin(3) + t38) - g(3) * (t13 + (-pkin(3) * t24 - pkin(4) * t19) * t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t37 - g(2) * t38 - g(3) * (-pkin(4) * t48 + t13); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t2;
