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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:26:00
% EndTime: 2022-01-23 09:26:01
% DurationCPUTime: 0.19s
% Computational Cost: add. (153->49), mult. (188->80), div. (0->0), fcn. (209->10), ass. (0->45)
t27 = cos(pkin(8));
t30 = cos(qJ(3));
t31 = cos(qJ(1));
t33 = t31 * t30;
t28 = sin(qJ(3));
t29 = sin(qJ(1));
t40 = t29 * t28;
t12 = t27 * t40 + t33;
t34 = t31 * t28;
t39 = t29 * t30;
t14 = -t27 * t34 + t39;
t26 = sin(pkin(8));
t45 = g(3) * t26;
t48 = -g(1) * t14 + g(2) * t12 + t28 * t45;
t25 = qJ(3) + pkin(9);
t23 = qJ(5) + t25;
t18 = sin(t23);
t44 = t29 * t18;
t19 = cos(t23);
t43 = t29 * t19;
t21 = sin(t25);
t42 = t29 * t21;
t22 = cos(t25);
t41 = t29 * t22;
t38 = t31 * t18;
t37 = t31 * t19;
t36 = t31 * t21;
t35 = t31 * t22;
t17 = g(1) * t31 + g(2) * t29;
t16 = g(1) * t29 - g(2) * t31;
t20 = t28 * pkin(3) + qJ(2);
t15 = t27 * t33 + t40;
t13 = -t27 * t39 + t34;
t11 = (qJ(4) + pkin(6)) * t26 + pkin(1) + (pkin(3) * t30 + pkin(2)) * t27;
t10 = t27 * t35 + t42;
t9 = -t27 * t36 + t41;
t8 = -t27 * t41 + t36;
t7 = t27 * t42 + t35;
t6 = t27 * t37 + t44;
t5 = -t27 * t38 + t43;
t4 = -t27 * t43 + t38;
t3 = t27 * t44 + t37;
t2 = g(1) * t6 - g(2) * t4 + t19 * t45;
t1 = -g(1) * t5 + g(2) * t3 + t18 * t45;
t24 = [0, t16, t17, t16 * t27, -t17, -g(1) * (-t29 * pkin(1) + t31 * qJ(2)) - g(2) * (t31 * pkin(1) + t29 * qJ(2)), 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t15, -g(1) * t12 - g(2) * t14, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t16 * t26, -g(1) * (-t11 * t29 + t20 * t31) - g(2) * (t11 * t31 + t20 * t29), 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, g(1) * t15 - g(2) * t13 + t30 * t45, -g(1) * t9 + g(2) * t7 + t21 * t45, g(1) * t10 - g(2) * t8 + t22 * t45, 0, t48 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t27 - t17 * t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t24;
