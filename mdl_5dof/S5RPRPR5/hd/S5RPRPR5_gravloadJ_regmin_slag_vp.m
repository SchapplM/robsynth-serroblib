% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:56
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:55:12
% EndTime: 2021-01-15 11:55:14
% DurationCPUTime: 0.21s
% Computational Cost: add. (153->50), mult. (188->80), div. (0->0), fcn. (209->10), ass. (0->45)
t25 = cos(pkin(8));
t29 = cos(qJ(3));
t30 = cos(qJ(1));
t34 = t30 * t29;
t27 = sin(qJ(3));
t28 = sin(qJ(1));
t41 = t28 * t27;
t11 = -t25 * t41 - t34;
t35 = t30 * t27;
t40 = t28 * t29;
t13 = t25 * t35 - t40;
t24 = sin(pkin(8));
t48 = g(1) * t24;
t50 = -g(2) * t11 - g(3) * t13 + t27 * t48;
t49 = (pkin(3) * t29 + pkin(2)) * t25 + (qJ(4) + pkin(6)) * t24 + pkin(1);
t23 = qJ(3) + pkin(9);
t22 = qJ(5) + t23;
t17 = sin(t22);
t45 = t28 * t17;
t18 = cos(t22);
t44 = t28 * t18;
t20 = sin(t23);
t43 = t28 * t20;
t21 = cos(t23);
t42 = t28 * t21;
t39 = t30 * t17;
t38 = t30 * t18;
t37 = t30 * t20;
t36 = t30 * t21;
t16 = g(2) * t30 + g(3) * t28;
t15 = g(2) * t28 - g(3) * t30;
t19 = -t27 * pkin(3) - qJ(2);
t14 = t25 * t34 + t41;
t12 = t25 * t40 - t35;
t10 = t25 * t36 + t43;
t9 = t25 * t37 - t42;
t8 = t25 * t42 - t37;
t7 = -t25 * t43 - t36;
t6 = t25 * t38 + t45;
t5 = t25 * t39 - t44;
t4 = t25 * t44 - t39;
t3 = -t25 * t45 - t38;
t2 = g(2) * t4 - g(3) * t6 + t18 * t48;
t1 = -g(2) * t3 - g(3) * t5 + t17 * t48;
t26 = [0, -t16, t15, -t16 * t25, -t15, -g(2) * (t30 * pkin(1) + t28 * qJ(2)) - g(3) * (t28 * pkin(1) - t30 * qJ(2)), 0, 0, 0, 0, 0, -g(2) * t14 - g(3) * t12, g(2) * t13 - g(3) * t11, -g(2) * t10 - g(3) * t8, g(2) * t9 - g(3) * t7, -t16 * t24, -g(2) * (-t19 * t28 + t49 * t30) - g(3) * (t19 * t30 + t49 * t28), 0, 0, 0, 0, 0, -g(2) * t6 - g(3) * t4, g(2) * t5 - g(3) * t3; 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, g(2) * t12 - g(3) * t14 + t29 * t48, -g(2) * t7 - g(3) * t9 + t20 * t48, g(2) * t8 - g(3) * t10 + t21 * t48, 0, t50 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t25 - t15 * t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t26;
