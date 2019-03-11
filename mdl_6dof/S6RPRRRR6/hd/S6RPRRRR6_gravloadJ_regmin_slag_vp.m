% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR6
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t18 = g(1) * t33 + g(2) * t31;
t26 = pkin(11) + qJ(3);
t21 = sin(t26);
t22 = cos(t26);
t35 = -g(3) * t22 + t18 * t21;
t49 = g(3) * t21;
t27 = qJ(4) + qJ(5);
t25 = qJ(6) + t27;
t19 = sin(t25);
t47 = t31 * t19;
t20 = cos(t25);
t46 = t31 * t20;
t23 = sin(t27);
t45 = t31 * t23;
t24 = cos(t27);
t44 = t31 * t24;
t30 = sin(qJ(4));
t43 = t31 * t30;
t32 = cos(qJ(4));
t42 = t31 * t32;
t41 = t33 * t19;
t40 = t33 * t20;
t39 = t33 * t23;
t38 = t33 * t24;
t37 = t33 * t30;
t36 = t33 * t32;
t17 = g(1) * t31 - g(2) * t33;
t16 = t22 * t36 + t43;
t15 = -t22 * t37 + t42;
t14 = -t22 * t42 + t37;
t13 = t22 * t43 + t36;
t12 = t22 * t38 + t45;
t11 = -t22 * t39 + t44;
t10 = -t22 * t44 + t39;
t9 = t22 * t45 + t38;
t8 = t22 * t40 + t47;
t7 = -t22 * t41 + t46;
t6 = -t22 * t46 + t41;
t5 = t22 * t47 + t40;
t4 = g(1) * t12 - g(2) * t10 + t24 * t49;
t3 = -g(1) * t11 + g(2) * t9 + t23 * t49;
t2 = g(1) * t8 - g(2) * t6 + t20 * t49;
t1 = -g(1) * t7 + g(2) * t5 + t19 * t49;
t28 = [0, t17, t18, t17 * cos(pkin(11)) -t17 * sin(pkin(11)) -t18, -g(1) * (-t31 * pkin(1) + t33 * qJ(2)) - g(2) * (t33 * pkin(1) + t31 * qJ(2)) 0, 0, 0, 0, 0, t17 * t22, -t17 * t21, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t18 * t22 + t49, 0, 0, 0, 0, 0, t35 * t32, -t35 * t30, 0, 0, 0, 0, 0, t35 * t24, -t35 * t23, 0, 0, 0, 0, 0, t35 * t20, -t35 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t30 * t49, g(1) * t16 - g(2) * t14 + t32 * t49, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t28;
