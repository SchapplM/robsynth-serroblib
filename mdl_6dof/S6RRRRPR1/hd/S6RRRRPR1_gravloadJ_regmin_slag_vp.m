% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t15 = g(1) * t34 + g(2) * t31;
t28 = qJ(2) + qJ(3);
t26 = qJ(4) + t28;
t20 = pkin(11) + t26;
t16 = sin(t20);
t17 = cos(t20);
t44 = -g(3) * t17 + t15 * t16;
t43 = g(3) * t16;
t29 = sin(qJ(6));
t41 = t31 * t29;
t32 = cos(qJ(6));
t40 = t31 * t32;
t39 = t34 * t29;
t38 = t34 * t32;
t22 = cos(t26);
t24 = cos(t28);
t37 = pkin(3) * t24 + pkin(4) * t22;
t33 = cos(qJ(2));
t36 = t33 * pkin(2) + t37;
t21 = sin(t26);
t23 = sin(t28);
t35 = -pkin(3) * t23 - pkin(4) * t21;
t14 = g(1) * t31 - g(2) * t34;
t3 = -g(3) * t22 + t15 * t21;
t30 = sin(qJ(2));
t25 = -qJ(5) - pkin(9) - pkin(8) - pkin(7);
t11 = pkin(1) + t36;
t10 = t17 * t38 + t41;
t9 = -t17 * t39 + t40;
t8 = -t17 * t40 + t39;
t7 = t17 * t41 + t38;
t6 = g(3) * t23 + t15 * t24;
t5 = -g(3) * t24 + t15 * t23;
t4 = g(3) * t21 + t15 * t22;
t2 = t44 * t32;
t1 = t44 * t29;
t12 = [0, t14, t15, 0, 0, 0, 0, 0, t14 * t33, -t14 * t30, 0, 0, 0, 0, 0, t14 * t24, -t14 * t23, 0, 0, 0, 0, 0, t14 * t22, -t14 * t21, -t15, -g(1) * (-t31 * t11 - t34 * t25) - g(2) * (t34 * t11 - t31 * t25) 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t33 + t15 * t30, g(3) * t30 + t15 * t33, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t36 - t15 * (-t30 * pkin(2) + t35) 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t37 - t15 * t35, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t3 * pkin(4), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t29 * t43, g(1) * t10 - g(2) * t8 + t32 * t43;];
taug_reg  = t12;
