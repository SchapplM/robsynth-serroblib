% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR1
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t19 = qJ(3) + qJ(4);
t17 = qJ(5) + t19;
t11 = sin(t17);
t12 = cos(t17);
t18 = qJ(1) + pkin(11);
t13 = sin(t18);
t14 = cos(t18);
t28 = g(1) * t14 + g(2) * t13;
t3 = -g(3) * t12 + t28 * t11;
t34 = g(3) * t11;
t20 = sin(qJ(6));
t32 = t13 * t20;
t23 = cos(qJ(6));
t31 = t13 * t23;
t30 = t14 * t20;
t29 = t14 * t23;
t27 = g(1) * t13 - g(2) * t14;
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t26 = g(1) * t22 - g(2) * t25;
t24 = cos(qJ(3));
t21 = sin(qJ(3));
t16 = cos(t19);
t15 = sin(t19);
t10 = t12 * t29 + t32;
t9 = -t12 * t30 + t31;
t8 = -t12 * t31 + t30;
t7 = t12 * t32 + t29;
t6 = g(3) * t15 + t28 * t16;
t5 = -g(3) * t16 + t28 * t15;
t4 = t28 * t12 + t34;
t2 = t3 * t23;
t1 = t3 * t20;
t33 = [0, t26, g(1) * t25 + g(2) * t22, t26 * pkin(1), 0, 0, 0, 0, 0, t27 * t24, -t27 * t21, 0, 0, 0, 0, 0, t27 * t16, -t27 * t15, 0, 0, 0, 0, 0, t27 * t12, -t27 * t11, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t24 + t28 * t21, g(3) * t21 + t28 * t24, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t20 * t34, g(1) * t10 - g(2) * t8 + t23 * t34;];
taug_reg  = t33;
