% Calculate minimal parameter regressor of gravitation load for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:25:21
% EndTime: 2019-05-05 16:25:21
% DurationCPUTime: 0.19s
% Computational Cost: add. (220->55), mult. (184->80), div. (0->0), fcn. (184->12), ass. (0->40)
t18 = qJ(1) + pkin(9);
t10 = sin(t18);
t13 = cos(t18);
t34 = g(1) * t13 + g(2) * t10;
t17 = qJ(3) + pkin(10);
t12 = cos(t17);
t9 = sin(t17);
t28 = -g(3) * t12 + t34 * t9;
t47 = g(3) * t9;
t25 = cos(qJ(1));
t24 = cos(qJ(3));
t14 = t24 * pkin(3);
t7 = t14 + pkin(2);
t43 = t25 * pkin(1) + t13 * t7;
t42 = t10 * t12;
t19 = sin(pkin(11));
t41 = t10 * t19;
t20 = cos(pkin(11));
t40 = t10 * t20;
t39 = t12 * t13;
t16 = pkin(11) + qJ(6);
t11 = cos(t16);
t38 = t13 * t11;
t37 = t13 * t19;
t36 = t13 * t20;
t35 = t9 * qJ(5);
t33 = g(1) * t10 - g(2) * t13;
t23 = sin(qJ(1));
t32 = g(1) * t23 - g(2) * t25;
t21 = -qJ(4) - pkin(7);
t31 = -t23 * pkin(1) - t13 * t21;
t30 = t12 * pkin(4) + t35;
t22 = sin(qJ(3));
t26 = -g(3) * t24 + t34 * t22;
t8 = sin(t16);
t4 = t10 * t8 + t12 * t38;
t3 = t10 * t11 - t8 * t39;
t2 = -t11 * t42 + t13 * t8;
t1 = t8 * t42 + t38;
t5 = [0, t32, g(1) * t25 + g(2) * t23, t32 * pkin(1), 0, 0, 0, 0, 0, t33 * t24, -t33 * t22, -t34, -g(1) * (-t10 * t7 + t31) - g(2) * (-t10 * t21 + t43) -g(1) * (-t12 * t40 + t37) - g(2) * (t12 * t36 + t41) -g(1) * (t12 * t41 + t36) - g(2) * (-t12 * t37 + t40) t33 * t9, -g(1) * t31 - g(2) * (pkin(4) * t39 + t13 * t35 + t43) + (-g(1) * (-t30 - t7) + g(2) * t21) * t10, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, g(3) * t22 + t34 * t24, 0, t26 * pkin(3), t28 * t20, -t28 * t19, -t34 * t12 - t47, -g(3) * (t14 + t30) + t34 * (pkin(3) * t22 + pkin(4) * t9 - qJ(5) * t12) 0, 0, 0, 0, 0, t28 * t11, -t28 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, 0, 0, -t33, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t8 * t47, g(1) * t4 - g(2) * t2 + t11 * t47;];
taug_reg  = t5;
