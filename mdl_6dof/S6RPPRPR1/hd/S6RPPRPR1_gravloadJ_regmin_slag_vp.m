% Calculate minimal parameter regressor of gravitation load for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t19 = qJ(1) + pkin(9);
t12 = sin(t19);
t15 = cos(t19);
t33 = g(1) * t15 + g(2) * t12;
t18 = pkin(10) + qJ(4);
t11 = sin(t18);
t14 = cos(t18);
t1 = -g(3) * t14 + t33 * t11;
t43 = g(3) * t11;
t25 = sin(qJ(1));
t41 = t25 * pkin(1);
t40 = t12 * t14;
t20 = sin(pkin(11));
t39 = t12 * t20;
t22 = cos(pkin(11));
t38 = t12 * t22;
t17 = pkin(11) + qJ(6);
t10 = sin(t17);
t37 = t15 * t10;
t13 = cos(t17);
t36 = t15 * t13;
t35 = t15 * t20;
t34 = t15 * t22;
t32 = g(1) * t12 - g(2) * t15;
t26 = cos(qJ(1));
t31 = g(1) * t25 - g(2) * t26;
t30 = t14 * pkin(4) + t11 * qJ(5);
t23 = cos(pkin(10));
t28 = t23 * pkin(3) + pkin(2) + t30;
t24 = -pkin(7) - qJ(3);
t16 = t26 * pkin(1);
t7 = t32 * t11;
t6 = t12 * t10 + t14 * t36;
t5 = t12 * t13 - t14 * t37;
t4 = -t13 * t40 + t37;
t3 = t10 * t40 + t36;
t2 = t33 * t14 + t43;
t8 = [0, t31, g(1) * t26 + g(2) * t25, t31 * pkin(1), t32 * t23, -t32 * sin(pkin(10)) -t33, -g(1) * (-t12 * pkin(2) + t15 * qJ(3) - t41) - g(2) * (t15 * pkin(2) + t12 * qJ(3) + t16) 0, 0, 0, 0, 0, t32 * t14, -t7, -g(1) * (-t14 * t38 + t35) - g(2) * (t14 * t34 + t39) -g(1) * (t14 * t39 + t34) - g(2) * (-t14 * t35 + t38) t7, g(1) * t41 - g(2) * t16 + (g(1) * t24 - g(2) * t28) * t15 + (g(1) * t28 + g(2) * t24) * t12, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1 * t22, -t1 * t20, -t2, -g(3) * t30 + t33 * (pkin(4) * t11 - qJ(5) * t14) 0, 0, 0, 0, 0, t1 * t13, -t1 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t10 * t43, g(1) * t6 - g(2) * t4 + t13 * t43;];
taug_reg  = t8;
