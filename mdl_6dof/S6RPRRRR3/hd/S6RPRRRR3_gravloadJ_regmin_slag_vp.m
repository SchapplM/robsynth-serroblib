% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:59:16
% EndTime: 2019-05-06 02:59:17
% DurationCPUTime: 0.20s
% Computational Cost: add. (302->44), mult. (242->81), div. (0->0), fcn. (276->12), ass. (0->43)
t27 = sin(qJ(3));
t30 = cos(qJ(3));
t24 = qJ(1) + pkin(11);
t19 = sin(t24);
t20 = cos(t24);
t36 = g(1) * t20 + g(2) * t19;
t33 = -g(3) * t30 + t36 * t27;
t44 = g(3) * t27;
t42 = t19 * t30;
t41 = t20 * t30;
t25 = qJ(4) + qJ(5);
t21 = sin(t25);
t40 = t21 * t30;
t22 = cos(t25);
t39 = t22 * t30;
t26 = sin(qJ(4));
t38 = t26 * t30;
t29 = cos(qJ(4));
t37 = t29 * t30;
t35 = g(1) * t19 - g(2) * t20;
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t34 = g(1) * t28 - g(2) * t31;
t23 = qJ(6) + t25;
t18 = cos(t23);
t17 = sin(t23);
t16 = t19 * t26 + t20 * t37;
t15 = t19 * t29 - t20 * t38;
t14 = -t19 * t37 + t20 * t26;
t13 = t19 * t38 + t20 * t29;
t12 = t19 * t21 + t20 * t39;
t11 = t19 * t22 - t20 * t40;
t10 = -t19 * t39 + t20 * t21;
t9 = t19 * t40 + t20 * t22;
t8 = t19 * t17 + t18 * t41;
t7 = -t17 * t41 + t19 * t18;
t6 = t20 * t17 - t18 * t42;
t5 = t17 * t42 + t20 * t18;
t4 = g(1) * t12 - g(2) * t10 + t22 * t44;
t3 = -g(1) * t11 + g(2) * t9 + t21 * t44;
t2 = g(1) * t8 - g(2) * t6 + t18 * t44;
t1 = -g(1) * t7 + g(2) * t5 + t17 * t44;
t32 = [0, t34, g(1) * t31 + g(2) * t28, t34 * pkin(1), 0, 0, 0, 0, 0, t35 * t30, -t35 * t27, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t36 * t30 + t44, 0, 0, 0, 0, 0, t33 * t29, -t33 * t26, 0, 0, 0, 0, 0, t33 * t22, -t33 * t21, 0, 0, 0, 0, 0, t33 * t18, -t33 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t26 * t44, g(1) * t16 - g(2) * t14 + t29 * t44, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t32;
