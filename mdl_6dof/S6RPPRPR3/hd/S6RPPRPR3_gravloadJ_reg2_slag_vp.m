% Calculate inertial parameters regressor of gravitation load for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t22 = qJ(1) + pkin(9);
t17 = sin(t22);
t19 = cos(t22);
t50 = -g(1) * t17 + g(2) * t19;
t21 = qJ(4) + pkin(10);
t16 = sin(t21);
t18 = cos(t21);
t49 = -g(3) * t16 - t18 * t50;
t45 = g(3) * t18;
t25 = sin(qJ(4));
t44 = t25 * pkin(4);
t24 = sin(qJ(6));
t43 = t17 * t24;
t27 = cos(qJ(6));
t42 = t17 * t27;
t41 = t19 * t24;
t40 = t19 * t27;
t29 = cos(qJ(1));
t39 = t29 * pkin(1) + t19 * pkin(2) + t17 * qJ(3);
t26 = sin(qJ(1));
t38 = -t26 * pkin(1) + t19 * qJ(3);
t37 = t16 * pkin(5) - t18 * pkin(8);
t8 = g(1) * t19 + g(2) * t17;
t36 = g(1) * t26 - g(2) * t29;
t35 = -t17 * pkin(2) + t38;
t23 = -qJ(5) - pkin(7);
t34 = t17 * t44 - t19 * t23 + t39;
t32 = t17 * t23 + t19 * t44 + t35;
t28 = cos(qJ(4));
t30 = g(3) * t25 + t28 * t50;
t6 = t16 * t40 - t43;
t5 = t16 * t41 + t42;
t4 = t16 * t42 + t41;
t3 = -t16 * t43 + t40;
t2 = t8 * t18;
t1 = -t16 * t50 + t45;
t7 = [0, 0, 0, 0, 0, 0, t36, g(1) * t29 + g(2) * t26, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t8, 0, t36 * pkin(1), 0, 0, 0, 0, 0, 0, 0, t50, -t8, -g(1) * t35 - g(2) * t39, 0, 0, 0, 0, 0, 0, -t8 * t25, -t8 * t28, -t50, -g(1) * ((-pkin(2) - pkin(7)) * t17 + t38) - g(2) * (t19 * pkin(7) + t39) 0, 0, 0, 0, 0, 0, -t8 * t16, -t2, -t50, -g(1) * t32 - g(2) * t34, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t2, -g(1) * (t19 * t37 + t32) - g(2) * (t17 * t37 + t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, g(3) * t28 - t25 * t50, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t1, 0, t30 * pkin(4), 0, 0, 0, 0, 0, 0, -t49 * t27, t49 * t24, -t1, -g(3) * (-t37 - t44) + t50 * (pkin(4) * t28 + pkin(5) * t18 + pkin(8) * t16); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t24 * t45, g(1) * t4 - g(2) * t6 + t27 * t45, 0, 0;];
taug_reg  = t7;
