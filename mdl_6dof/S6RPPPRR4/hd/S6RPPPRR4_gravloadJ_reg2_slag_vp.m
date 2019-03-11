% Calculate inertial parameters regressor of gravitation load for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t39 = sin(pkin(9));
t40 = cos(pkin(9));
t44 = sin(qJ(1));
t45 = cos(qJ(1));
t8 = -t44 * t39 - t45 * t40;
t9 = t45 * t39 - t44 * t40;
t51 = -g(1) * t9 + g(2) * t8;
t21 = sin(qJ(5));
t23 = cos(qJ(5));
t50 = g(3) * t21 - t23 * t51;
t46 = g(3) * t23;
t20 = sin(qJ(6));
t43 = t20 * t21;
t22 = cos(qJ(6));
t42 = t21 * t22;
t41 = t45 * pkin(1) + t44 * qJ(2);
t38 = t45 * pkin(2) + t41;
t37 = -t8 * pkin(3) + t38;
t35 = g(1) * t8 + g(2) * t9;
t34 = -t44 * pkin(1) + t45 * qJ(2);
t32 = t21 * pkin(5) - t23 * pkin(8);
t31 = t9 * t20 + t8 * t42;
t30 = -t9 * t22 + t8 * t43;
t29 = qJ(4) + t32;
t28 = t9 * qJ(4) + t37;
t27 = -t44 * pkin(2) + t34;
t25 = t9 * pkin(3) + t27;
t24 = t8 * qJ(4) + t25;
t11 = g(1) * t45 + g(2) * t44;
t10 = g(1) * t44 - g(2) * t45;
t4 = t35 * t23;
t3 = -t8 * t20 + t9 * t42;
t2 = -t8 * t22 - t9 * t43;
t1 = -t21 * t51 - t46;
t5 = [0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t11, -g(1) * t34 - g(2) * t41, 0, 0, 0, 0, 0, 0, t51, t35, 0, -g(1) * t27 - g(2) * t38, 0, 0, 0, 0, 0, 0, 0, -t51, -t35, -g(1) * t24 - g(2) * t28, 0, 0, 0, 0, 0, 0, -t35 * t21, -t4, t51, -g(1) * (t9 * pkin(7) + t24) - g(2) * (-t8 * pkin(7) + t28) 0, 0, 0, 0, 0, 0, -g(1) * t31 - g(2) * t3, g(1) * t30 - g(2) * t2, t4, -g(1) * t25 - g(2) * t37 + (-g(1) * pkin(7) - g(2) * t29) * t9 + (g(2) * pkin(7) - g(1) * t29) * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t50 * t22, t50 * t20, -t1, -g(3) * t32 + t51 * (pkin(5) * t23 + pkin(8) * t21); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t30 - t20 * t46, g(1) * t3 - g(2) * t31 - t22 * t46, 0, 0;];
taug_reg  = t5;
