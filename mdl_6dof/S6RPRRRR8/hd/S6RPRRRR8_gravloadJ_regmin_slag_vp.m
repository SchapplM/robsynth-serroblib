% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x34]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:15:55
% EndTime: 2019-05-06 04:15:56
% DurationCPUTime: 0.18s
% Computational Cost: add. (182->44), mult. (218->62), div. (0->0), fcn. (236->10), ass. (0->40)
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t17 = g(1) * t28 - g(2) * t31;
t25 = qJ(3) + qJ(4);
t20 = sin(t25);
t22 = cos(t25);
t32 = -g(3) * t20 + t17 * t22;
t41 = g(3) * t22;
t24 = qJ(5) + qJ(6);
t19 = sin(t24);
t40 = t28 * t19;
t21 = cos(t24);
t39 = t28 * t21;
t26 = sin(qJ(5));
t38 = t28 * t26;
t29 = cos(qJ(5));
t37 = t28 * t29;
t36 = t31 * t19;
t35 = t31 * t21;
t34 = t31 * t26;
t33 = t31 * t29;
t18 = g(1) * t31 + g(2) * t28;
t30 = cos(qJ(3));
t27 = sin(qJ(3));
t16 = t20 * t33 - t38;
t15 = t20 * t34 + t37;
t14 = t20 * t37 + t34;
t13 = -t20 * t38 + t33;
t12 = t20 * t35 - t40;
t11 = t20 * t36 + t39;
t10 = t20 * t39 + t36;
t9 = -t20 * t40 + t35;
t7 = t17 * t20 + t41;
t6 = t32 * t29;
t5 = t32 * t26;
t4 = t32 * t21;
t3 = t32 * t19;
t2 = g(1) * t10 - g(2) * t12 + t21 * t41;
t1 = -g(1) * t9 - g(2) * t11 + t19 * t41;
t8 = [0, t17, t18, -t17, -t18, -g(1) * (-t28 * pkin(1) + t31 * qJ(2)) - g(2) * (t31 * pkin(1) + t28 * qJ(2)) 0, 0, 0, 0, 0, -t18 * t27, -t18 * t30, 0, 0, 0, 0, 0, -t18 * t20, -t18 * t22, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14, g(1) * t15 - g(2) * t13, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9; 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t27 - t17 * t30, g(3) * t30 + t17 * t27, 0, 0, 0, 0, 0, -t32, t7, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t7, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t15 + t26 * t41, g(1) * t14 - g(2) * t16 + t29 * t41, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t8;
