% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR2
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:45:23
% EndTime: 2019-05-06 02:45:23
% DurationCPUTime: 0.16s
% Computational Cost: add. (254->38), mult. (210->65), div. (0->0), fcn. (228->12), ass. (0->40)
t25 = qJ(3) + qJ(4);
t20 = sin(t25);
t22 = cos(t25);
t23 = qJ(1) + pkin(11);
t17 = sin(t23);
t18 = cos(t23);
t35 = g(1) * t18 + g(2) * t17;
t7 = -g(3) * t22 + t35 * t20;
t41 = g(3) * t20;
t24 = qJ(5) + qJ(6);
t19 = sin(t24);
t39 = t19 * t22;
t21 = cos(t24);
t38 = t21 * t22;
t26 = sin(qJ(5));
t37 = t22 * t26;
t29 = cos(qJ(5));
t36 = t22 * t29;
t34 = g(1) * t17 - g(2) * t18;
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t33 = g(1) * t28 - g(2) * t31;
t30 = cos(qJ(3));
t27 = sin(qJ(3));
t16 = t17 * t26 + t18 * t36;
t15 = t17 * t29 - t18 * t37;
t14 = -t17 * t36 + t18 * t26;
t13 = t17 * t37 + t18 * t29;
t12 = t17 * t19 + t18 * t38;
t11 = t17 * t21 - t18 * t39;
t10 = -t17 * t38 + t18 * t19;
t9 = t17 * t39 + t18 * t21;
t8 = t35 * t22 + t41;
t6 = t7 * t29;
t5 = t7 * t26;
t4 = t7 * t21;
t3 = t7 * t19;
t2 = g(1) * t12 - g(2) * t10 + t21 * t41;
t1 = -g(1) * t11 + g(2) * t9 + t19 * t41;
t32 = [0, t33, g(1) * t31 + g(2) * t28, t33 * pkin(1), 0, 0, 0, 0, 0, t34 * t30, -t34 * t27, 0, 0, 0, 0, 0, t34 * t22, -t34 * t20, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t30 + t35 * t27, g(3) * t27 + t35 * t30, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t26 * t41, g(1) * t16 - g(2) * t14 + t29 * t41, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t32;
