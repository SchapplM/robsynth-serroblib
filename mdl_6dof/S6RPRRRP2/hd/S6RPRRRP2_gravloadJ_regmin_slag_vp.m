% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:15:07
% EndTime: 2019-05-06 01:15:07
% DurationCPUTime: 0.23s
% Computational Cost: add. (246->52), mult. (244->84), div. (0->0), fcn. (251->10), ass. (0->43)
t24 = qJ(1) + pkin(10);
t18 = sin(t24);
t19 = cos(t24);
t39 = g(1) * t19 + g(2) * t18;
t25 = qJ(4) + qJ(5);
t20 = sin(t25);
t21 = cos(t25);
t30 = cos(qJ(3));
t44 = t20 * t30;
t3 = t18 * t44 + t19 * t21;
t27 = sin(qJ(3));
t48 = g(3) * t27;
t5 = t18 * t21 - t19 * t44;
t1 = -g(1) * t5 + g(2) * t3 + t20 * t48;
t7 = -g(3) * t30 + t39 * t27;
t26 = sin(qJ(4));
t15 = t26 * pkin(4) + pkin(5) * t20;
t46 = pkin(7) + t15;
t45 = t15 * t30;
t43 = t21 * t30;
t42 = t26 * t30;
t29 = cos(qJ(4));
t41 = t29 * t30;
t16 = t29 * pkin(4) + pkin(5) * t21;
t38 = g(1) * t18 - g(2) * t19;
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t37 = g(1) * t28 - g(2) * t31;
t14 = pkin(3) + t16;
t23 = -qJ(6) - pkin(9) - pkin(8);
t36 = t30 * t14 - t27 * t23;
t34 = pkin(2) + t36;
t33 = t37 * pkin(1);
t13 = t38 * t27;
t12 = t18 * t26 + t19 * t41;
t11 = t18 * t29 - t19 * t42;
t10 = -t18 * t41 + t19 * t26;
t9 = t18 * t42 + t19 * t29;
t8 = t39 * t30 + t48;
t6 = t18 * t20 + t19 * t43;
t4 = -t18 * t43 + t19 * t20;
t2 = g(1) * t6 - g(2) * t4 + t21 * t48;
t17 = [0, t37, g(1) * t31 + g(2) * t28, t33, 0, 0, 0, 0, 0, t38 * t30, -t13, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t13, t33 + (-g(1) * t46 - g(2) * t34) * t19 + (g(1) * t34 - g(2) * t46) * t18; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t7 * t29, -t7 * t26, 0, 0, 0, 0, 0, t7 * t21, -t7 * t20, -t8, -g(3) * t36 + t39 * (t14 * t27 + t23 * t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 + t26 * t48, g(1) * t12 - g(2) * t10 + t29 * t48, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t18 * t16 - t19 * t45) - g(2) * (-t19 * t16 - t18 * t45) + t15 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg  = t17;
