% Calculate inertial parameters regressor of gravitation load for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPPRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t22 = qJ(1) + pkin(9);
t17 = sin(t22);
t19 = cos(t22);
t49 = -g(1) * t17 + g(2) * t19;
t21 = pkin(10) + qJ(5);
t16 = sin(t21);
t18 = cos(t21);
t48 = -g(3) * t16 - t18 * t49;
t23 = sin(pkin(10));
t47 = pkin(4) * t23;
t43 = g(3) * t18;
t26 = sin(qJ(6));
t42 = t17 * t26;
t28 = cos(qJ(6));
t41 = t17 * t28;
t40 = t19 * t26;
t39 = t19 * t28;
t29 = cos(qJ(1));
t38 = t29 * pkin(1) + t19 * pkin(2) + t17 * qJ(3);
t27 = sin(qJ(1));
t37 = -t27 * pkin(1) + t19 * qJ(3);
t35 = t16 * pkin(5) - t18 * pkin(8);
t8 = g(1) * t19 + g(2) * t17;
t34 = g(1) * t27 - g(2) * t29;
t33 = -t17 * pkin(2) + t37;
t25 = -pkin(7) - qJ(4);
t32 = t17 * t47 - t19 * t25 + t38;
t31 = t17 * t25 + t19 * t47 + t33;
t6 = t16 * t39 - t42;
t5 = t16 * t40 + t41;
t4 = t16 * t41 + t40;
t3 = -t16 * t42 + t39;
t2 = t8 * t18;
t1 = -t16 * t49 + t43;
t7 = [0, 0, 0, 0, 0, 0, t34, g(1) * t29 + g(2) * t27, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t8, 0, t34 * pkin(1), 0, 0, 0, 0, 0, 0, 0, t49, -t8, -g(1) * t33 - g(2) * t38, 0, 0, 0, 0, 0, 0, -t8 * t23, -t8 * cos(pkin(10)) -t49, -g(1) * ((-pkin(2) - qJ(4)) * t17 + t37) - g(2) * (t19 * qJ(4) + t38) 0, 0, 0, 0, 0, 0, -t8 * t16, -t2, -t49, -g(1) * t31 - g(2) * t32, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t2, -g(1) * (t35 * t19 + t31) - g(2) * (t35 * t17 + t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t28, t48 * t26, -t1, g(3) * t35 + t49 * (pkin(5) * t18 + pkin(8) * t16); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t26 * t43, g(1) * t4 - g(2) * t6 + t28 * t43, 0, 0;];
taug_reg  = t7;
