% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t33 = t28 * pkin(3) - t30 * qJ(4);
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t54 = -g(1) * t29 + g(2) * t31;
t12 = -g(3) * t28 - t30 * t54;
t50 = g(3) * t30;
t25 = pkin(10) + qJ(5);
t20 = qJ(6) + t25;
t16 = sin(t20);
t48 = t29 * t16;
t17 = cos(t20);
t47 = t29 * t17;
t18 = sin(t25);
t46 = t29 * t18;
t19 = cos(t25);
t45 = t29 * t19;
t26 = sin(pkin(10));
t44 = t29 * t26;
t27 = cos(pkin(10));
t43 = t29 * t27;
t42 = t31 * t16;
t41 = t31 * t17;
t40 = t31 * t18;
t39 = t31 * t19;
t38 = t31 * t26;
t37 = t31 * t27;
t36 = t31 * pkin(1) + t29 * qJ(2);
t15 = g(1) * t31 + g(2) * t29;
t22 = t31 * qJ(2);
t13 = t15 * t30;
t11 = -t28 * t54 + t50;
t10 = t28 * t39 - t46;
t9 = t28 * t40 + t45;
t8 = t28 * t45 + t40;
t7 = -t28 * t46 + t39;
t6 = t28 * t41 - t48;
t5 = t28 * t42 + t47;
t4 = t28 * t47 + t42;
t3 = -t28 * t48 + t41;
t2 = g(1) * t4 - g(2) * t6 + t17 * t50;
t1 = -g(1) * t3 - g(2) * t5 + t16 * t50;
t14 = [0, -t54, t15, t54, -t15, -g(1) * (-t29 * pkin(1) + t22) - g(2) * t36, 0, 0, 0, 0, 0, -t15 * t28, -t13, -g(1) * (t28 * t37 - t44) - g(2) * (t28 * t43 + t38) -g(1) * (-t28 * t38 - t43) - g(2) * (-t28 * t44 + t37) t13, -g(1) * (t33 * t31 + t22) - g(2) * (t31 * pkin(7) + t36) + (-g(1) * (-pkin(1) - pkin(7)) - g(2) * t33) * t29, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, -t12 * t27, t12 * t26, -t11, g(3) * t33 + t54 * (pkin(3) * t30 + qJ(4) * t28) 0, 0, 0, 0, 0, -t12 * t19, t12 * t18, 0, 0, 0, 0, 0, -t12 * t17, t12 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t18 * t50, g(1) * t8 - g(2) * t10 + t19 * t50, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t14;
