% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t34 = cos(qJ(2));
t29 = qJ(2) + qJ(3);
t25 = cos(t29);
t26 = qJ(4) + t29;
t22 = sin(t26);
t23 = cos(t26);
t43 = t23 * pkin(4) + t22 * qJ(5);
t41 = pkin(3) * t25 + t43;
t52 = t34 * pkin(2) + t41;
t24 = sin(t29);
t50 = pkin(4) * t22;
t39 = -pkin(3) * t24 - t50;
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t17 = g(1) * t35 + g(2) * t32;
t4 = g(3) * t22 + t17 * t23;
t48 = g(3) * t23;
t30 = sin(qJ(6));
t47 = t32 * t30;
t33 = cos(qJ(6));
t46 = t32 * t33;
t45 = t35 * t30;
t44 = t35 * t33;
t42 = qJ(5) * t23;
t31 = sin(qJ(2));
t40 = -t31 * pkin(2) + t39;
t38 = g(1) * t32 - g(2) * t35;
t37 = pkin(1) + t52;
t28 = -pkin(9) - pkin(8) - pkin(7);
t16 = t35 * t42;
t15 = t32 * t42;
t12 = -t22 * t47 + t44;
t11 = t22 * t46 + t45;
t10 = t22 * t45 + t46;
t9 = t22 * t44 - t47;
t8 = t38 * t23;
t7 = t38 * t22;
t6 = g(3) * t24 + t17 * t25;
t5 = -g(3) * t25 + t17 * t24;
t3 = t17 * t22 - t48;
t2 = t4 * t33;
t1 = t4 * t30;
t13 = [0, t38, t17, 0, 0, 0, 0, 0, t38 * t34, -t38 * t31, 0, 0, 0, 0, 0, t38 * t25, -t38 * t24, 0, 0, 0, 0, 0, t8, -t7, -t17, -t8, t7 (g(1) * t28 - g(2) * t37) * t35 + (g(1) * t37 + g(2) * t28) * t32, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t34 + t17 * t31, g(3) * t31 + t17 * t34, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, -t3, -t4, -g(1) * (t40 * t35 + t16) - g(2) * (t40 * t32 + t15) - g(3) * t52, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, -t3, -t4, -g(1) * (t39 * t35 + t16) - g(2) * (t39 * t32 + t15) - g(3) * t41, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, -t3, -t4, -g(1) * (-t35 * t50 + t16) - g(2) * (-t32 * t50 + t15) - g(3) * t43, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11 + t33 * t48, g(1) * t10 - g(2) * t12 - t30 * t48;];
taug_reg  = t13;
