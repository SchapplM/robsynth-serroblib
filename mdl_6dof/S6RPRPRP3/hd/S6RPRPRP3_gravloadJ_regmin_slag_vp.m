% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t27 = cos(pkin(10));
t17 = t27 * pkin(4) + pkin(3);
t31 = cos(qJ(3));
t28 = -pkin(8) - qJ(4);
t29 = sin(qJ(3));
t45 = t29 * t28;
t56 = t31 * t17 - t45;
t25 = qJ(1) + pkin(9);
t19 = sin(t25);
t21 = cos(t25);
t40 = g(1) * t21 + g(2) * t19;
t9 = -g(3) * t31 + t40 * t29;
t55 = g(1) * t19;
t52 = g(3) * t29;
t50 = t19 * t31;
t26 = sin(pkin(10));
t49 = t21 * t26;
t48 = t21 * t31;
t47 = t26 * t31;
t46 = t27 * t31;
t32 = cos(qJ(1));
t43 = t32 * pkin(1) + t21 * pkin(2) + t19 * pkin(7);
t30 = sin(qJ(1));
t42 = -t30 * pkin(1) + t21 * pkin(7);
t24 = pkin(10) + qJ(5);
t18 = sin(t24);
t20 = cos(t24);
t5 = t18 * t50 + t21 * t20;
t7 = t18 * t48 - t19 * t20;
t41 = g(1) * t5 - g(2) * t7;
t39 = -g(2) * t21 + t55;
t38 = g(1) * t30 - g(2) * t32;
t37 = t31 * pkin(3) + t29 * qJ(4);
t35 = pkin(5) * t20 + qJ(6) * t18 + t17;
t1 = g(1) * t7 + g(2) * t5 + t18 * t52;
t6 = -t21 * t18 + t20 * t50;
t8 = t19 * t18 + t20 * t48;
t34 = g(1) * t8 + g(2) * t6 + t20 * t52;
t11 = t39 * t29;
t10 = t40 * t31 + t52;
t4 = t9 * t20;
t3 = t9 * t18;
t2 = g(1) * t6 - g(2) * t8;
t12 = [0, t38, g(1) * t32 + g(2) * t30, t38 * pkin(1), 0, 0, 0, 0, 0, t39 * t31, -t11, -g(1) * (-t19 * t46 + t49) - g(2) * (t19 * t26 + t21 * t46) -g(1) * (t19 * t47 + t21 * t27) - g(2) * (t19 * t27 - t21 * t47) t11, -g(1) * t42 - g(2) * (t37 * t21 + t43) - (-pkin(2) - t37) * t55, 0, 0, 0, 0, 0, t2, -t41, t2, t11, t41, -g(1) * (pkin(4) * t49 - t6 * pkin(5) - t5 * qJ(6) + t42) - g(2) * (t8 * pkin(5) + t7 * qJ(6) + t56 * t21 + t43) + (-g(1) * (-pkin(2) - t56) - g(2) * pkin(4) * t26) * t19; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, t9 * t27, -t9 * t26, -t10, -g(3) * t37 + t40 * (pkin(3) * t29 - qJ(4) * t31) 0, 0, 0, 0, 0, t4, -t3, t4, -t10, t3, -g(3) * (t35 * t31 - t45) + t40 * (t28 * t31 + t35 * t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t34, t1, 0, -t34, -g(1) * (-t7 * pkin(5) + t8 * qJ(6)) - g(2) * (-t5 * pkin(5) + t6 * qJ(6)) - (-pkin(5) * t18 + qJ(6) * t20) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t12;
