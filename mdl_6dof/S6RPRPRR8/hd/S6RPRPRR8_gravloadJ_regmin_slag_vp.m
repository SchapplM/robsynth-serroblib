% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR8
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
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t11 = g(1) * t26 - g(2) * t29;
t21 = qJ(3) + pkin(10);
t13 = sin(t21);
t14 = cos(t21);
t32 = -g(3) * t13 + t11 * t14;
t25 = sin(qJ(3));
t44 = pkin(3) * t25;
t42 = g(3) * t14;
t22 = qJ(5) + qJ(6);
t15 = sin(t22);
t41 = t26 * t15;
t16 = cos(t22);
t40 = t26 * t16;
t24 = sin(qJ(5));
t39 = t26 * t24;
t27 = cos(qJ(5));
t38 = t26 * t27;
t37 = t29 * t15;
t36 = t29 * t16;
t35 = t29 * t24;
t34 = t29 * t27;
t33 = t29 * pkin(1) + t26 * qJ(2);
t12 = g(1) * t29 + g(2) * t26;
t28 = cos(qJ(3));
t30 = g(3) * t25 - t11 * t28;
t23 = -qJ(4) - pkin(7);
t18 = t29 * qJ(2);
t10 = t13 * t34 - t39;
t9 = t13 * t35 + t38;
t8 = t13 * t38 + t35;
t7 = -t13 * t39 + t34;
t6 = t13 * t36 - t41;
t5 = t13 * t37 + t40;
t4 = t13 * t40 + t37;
t3 = -t13 * t41 + t36;
t2 = g(1) * t4 - g(2) * t6 + t16 * t42;
t1 = -g(1) * t3 - g(2) * t5 + t15 * t42;
t17 = [0, t11, t12, -t11, -t12, -g(1) * (-t26 * pkin(1) + t18) - g(2) * t33, 0, 0, 0, 0, 0, -t12 * t25, -t12 * t28, t11, -g(1) * (t29 * t44 + t18 + (-pkin(1) + t23) * t26) - g(2) * (-t29 * t23 + t26 * t44 + t33) 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, g(3) * t28 + t11 * t25, 0, t30 * pkin(3), 0, 0, 0, 0, 0, -t32 * t27, t32 * t24, 0, 0, 0, 0, 0, -t32 * t16, t32 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t24 * t42, g(1) * t8 - g(2) * t10 + t27 * t42, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t17;
