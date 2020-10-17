% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR5
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
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:32:07
% EndTime: 2019-05-06 03:32:08
% DurationCPUTime: 0.17s
% Computational Cost: add. (259->41), mult. (224->64), div. (0->0), fcn. (242->12), ass. (0->41)
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t18 = g(1) * t33 + g(2) * t31;
t26 = pkin(11) + qJ(3);
t23 = qJ(4) + t26;
t19 = sin(t23);
t20 = cos(t23);
t7 = -g(3) * t20 + t18 * t19;
t44 = g(3) * t19;
t27 = qJ(5) + qJ(6);
t24 = sin(t27);
t42 = t31 * t24;
t25 = cos(t27);
t41 = t31 * t25;
t30 = sin(qJ(5));
t40 = t31 * t30;
t32 = cos(qJ(5));
t39 = t31 * t32;
t38 = t33 * t24;
t37 = t33 * t25;
t36 = t33 * t30;
t35 = t33 * t32;
t17 = g(1) * t31 - g(2) * t33;
t22 = cos(t26);
t21 = sin(t26);
t16 = t20 * t35 + t40;
t15 = -t20 * t36 + t39;
t14 = -t20 * t39 + t36;
t13 = t20 * t40 + t35;
t12 = t20 * t37 + t42;
t11 = -t20 * t38 + t41;
t10 = -t20 * t41 + t38;
t9 = t20 * t42 + t37;
t8 = t18 * t20 + t44;
t6 = t7 * t32;
t5 = t7 * t30;
t4 = t7 * t25;
t3 = t7 * t24;
t2 = g(1) * t12 - g(2) * t10 + t25 * t44;
t1 = -g(1) * t11 + g(2) * t9 + t24 * t44;
t28 = [0, t17, t18, t17 * cos(pkin(11)) -t17 * sin(pkin(11)) -t18, -g(1) * (-t31 * pkin(1) + t33 * qJ(2)) - g(2) * (t33 * pkin(1) + t31 * qJ(2)) 0, 0, 0, 0, 0, t17 * t22, -t17 * t21, 0, 0, 0, 0, 0, t17 * t20, -t17 * t19, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11; 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t22 + t18 * t21, g(3) * t21 + t18 * t22, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t30 * t44, g(1) * t16 - g(2) * t14 + t32 * t44, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t28;
