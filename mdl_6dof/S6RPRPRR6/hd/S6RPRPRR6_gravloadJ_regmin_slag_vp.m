% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t32 = sin(qJ(1));
t33 = cos(qJ(1));
t15 = g(1) * t33 + g(2) * t32;
t26 = pkin(10) + qJ(3);
t20 = sin(t26);
t22 = cos(t26);
t3 = -g(3) * t22 + t15 * t20;
t51 = g(3) * t20;
t25 = pkin(11) + qJ(5);
t23 = qJ(6) + t25;
t16 = sin(t23);
t49 = t32 * t16;
t17 = cos(t23);
t48 = t32 * t17;
t19 = sin(t25);
t47 = t32 * t19;
t21 = cos(t25);
t46 = t32 * t21;
t27 = sin(pkin(11));
t45 = t32 * t27;
t29 = cos(pkin(11));
t44 = t32 * t29;
t43 = t33 * t16;
t42 = t33 * t17;
t41 = t33 * t19;
t40 = t33 * t21;
t39 = t33 * t27;
t38 = t33 * t29;
t14 = g(1) * t32 - g(2) * t33;
t37 = t22 * pkin(3) + t20 * qJ(4);
t30 = cos(pkin(10));
t35 = t30 * pkin(2) + pkin(1) + t37;
t31 = -pkin(7) - qJ(2);
t13 = t14 * t20;
t12 = t22 * t40 + t47;
t11 = -t22 * t41 + t46;
t10 = -t22 * t46 + t41;
t9 = t22 * t47 + t40;
t8 = t22 * t42 + t49;
t7 = -t22 * t43 + t48;
t6 = -t22 * t48 + t43;
t5 = t22 * t49 + t42;
t4 = t15 * t22 + t51;
t2 = g(1) * t8 - g(2) * t6 + t17 * t51;
t1 = -g(1) * t7 + g(2) * t5 + t16 * t51;
t18 = [0, t14, t15, t14 * t30, -t14 * sin(pkin(10)) -t15, -g(1) * (-t32 * pkin(1) + t33 * qJ(2)) - g(2) * (t33 * pkin(1) + t32 * qJ(2)) 0, 0, 0, 0, 0, t14 * t22, -t13, -g(1) * (-t22 * t44 + t39) - g(2) * (t22 * t38 + t45) -g(1) * (t22 * t45 + t38) - g(2) * (-t22 * t39 + t44) t13 (g(1) * t31 - g(2) * t35) * t33 + (g(1) * t35 + g(2) * t31) * t32, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, t3 * t29, -t3 * t27, -t4, -g(3) * t37 + t15 * (pkin(3) * t20 - qJ(4) * t22) 0, 0, 0, 0, 0, t3 * t21, -t3 * t19, 0, 0, 0, 0, 0, t3 * t17, -t3 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 + t19 * t51, g(1) * t12 - g(2) * t10 + t21 * t51, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t18;
