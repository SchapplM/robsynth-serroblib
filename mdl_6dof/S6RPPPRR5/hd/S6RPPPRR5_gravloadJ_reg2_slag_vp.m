% Calculate inertial parameters regressor of gravitation load for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPPRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = cos(pkin(9));
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t39 = sin(pkin(9));
t7 = t25 * t19 - t22 * t39;
t8 = t22 * t19 + t25 * t39;
t36 = g(1) * t8 - g(2) * t7;
t21 = sin(qJ(5));
t24 = cos(qJ(5));
t27 = -g(3) * t24 + t21 * t36;
t44 = g(3) * t21;
t20 = sin(qJ(6));
t42 = t20 * t24;
t23 = cos(qJ(6));
t41 = t23 * t24;
t40 = t25 * pkin(1) + t22 * qJ(2);
t38 = t25 * qJ(3) + t40;
t37 = t22 * pkin(3) + t38;
t35 = g(1) * t7 + g(2) * t8;
t34 = t24 * pkin(5) + t21 * pkin(8);
t14 = t25 * qJ(2);
t32 = t14 + (-pkin(1) - qJ(3)) * t22;
t31 = t8 * t20 + t41 * t7;
t30 = -t8 * t23 + t42 * t7;
t29 = t25 * pkin(3) + t32;
t28 = t8 * pkin(4) - t7 * pkin(7) + t37;
t26 = t7 * pkin(4) + t8 * pkin(7) + t29;
t10 = g(1) * t25 + g(2) * t22;
t9 = g(1) * t22 - g(2) * t25;
t4 = t35 * t21;
t3 = -t7 * t20 + t41 * t8;
t2 = -t7 * t23 - t42 * t8;
t1 = t24 * t36 + t44;
t5 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, -g(1) * (-t22 * pkin(1) + t14) - g(2) * t40, 0, 0, 0, 0, 0, 0, -t10, 0, t9, -g(1) * t32 - g(2) * t38, 0, 0, 0, 0, 0, 0, -t35, t36, 0, -g(1) * t29 - g(2) * t37, 0, 0, 0, 0, 0, 0, -t35 * t24, t4, -t36, -g(1) * t26 - g(2) * t28, 0, 0, 0, 0, 0, 0, -g(1) * t31 - g(2) * t3, g(1) * t30 - g(2) * t2, -t4, -g(1) * (t34 * t7 + t26) - g(2) * (t34 * t8 + t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t1, 0, 0, 0, 0, 0, 0, 0, 0, t27 * t23, -t27 * t20, -t1, -g(3) * t34 + t36 * (pkin(5) * t21 - pkin(8) * t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t30 + t20 * t44, g(1) * t3 - g(2) * t31 + t23 * t44, 0, 0;];
taug_reg  = t5;
