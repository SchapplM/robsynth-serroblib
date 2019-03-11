% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t22 = qJ(5) + qJ(6);
t20 = sin(t22);
t29 = cos(qJ(1));
t24 = sin(qJ(4));
t25 = sin(qJ(2));
t28 = cos(qJ(2));
t42 = cos(qJ(4));
t46 = -t28 * t24 + t25 * t42;
t11 = t46 * t29;
t15 = t25 * t24 + t28 * t42;
t26 = sin(qJ(1));
t9 = t46 * t26;
t31 = -g(1) * t11 - g(2) * t9 + g(3) * t15;
t50 = t31 * t20;
t21 = cos(t22);
t49 = t31 * t21;
t23 = sin(qJ(5));
t48 = t31 * t23;
t27 = cos(qJ(5));
t47 = t31 * t27;
t17 = g(1) * t29 + g(2) * t26;
t43 = g(3) * t46;
t39 = g(1) * t26 - g(2) * t29;
t38 = t28 * pkin(2) + t25 * qJ(3);
t10 = t15 * t26;
t36 = t10 * t21 + t29 * t20;
t35 = t10 * t20 - t29 * t21;
t34 = t10 * t27 + t29 * t23;
t33 = t10 * t23 - t29 * t27;
t32 = pkin(1) + t38;
t12 = t15 * t29;
t30 = g(1) * t12 + g(2) * t10 + t43;
t14 = t39 * t28;
t13 = t39 * t25;
t8 = g(3) * t25 + t17 * t28;
t7 = -g(3) * t28 + t17 * t25;
t6 = t12 * t27 - t26 * t23;
t5 = -t12 * t23 - t26 * t27;
t4 = t12 * t21 - t26 * t20;
t3 = -t12 * t20 - t26 * t21;
t2 = g(1) * t4 + g(2) * t36 + t21 * t43;
t1 = -g(1) * t3 + g(2) * t35 + t20 * t43;
t16 = [0, t39, t17, 0, 0, 0, 0, 0, t14, -t13, t14, -t17, t13 (-g(1) * pkin(7) - g(2) * t32) * t29 + (-g(2) * pkin(7) + g(1) * t32) * t26, 0, 0, 0, 0, 0, g(1) * t10 - g(2) * t12, g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, g(1) * t34 - g(2) * t6, -g(1) * t33 - g(2) * t5, 0, 0, 0, 0, 0, g(1) * t36 - g(2) * t4, -g(1) * t35 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, t7, 0, -t8, -g(3) * t38 + t17 * (pkin(2) * t25 - qJ(3) * t28) 0, 0, 0, 0, 0, -t31, -t30, 0, 0, 0, 0, 0, -t47, t48, 0, 0, 0, 0, 0, -t49, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t30, 0, 0, 0, 0, 0, t47, -t48, 0, 0, 0, 0, 0, t49, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t33 + t23 * t43, g(1) * t6 + g(2) * t34 + t27 * t43, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t16;
