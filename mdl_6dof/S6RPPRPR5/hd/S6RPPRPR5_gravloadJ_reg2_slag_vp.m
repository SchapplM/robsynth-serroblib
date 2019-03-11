% Calculate inertial parameters regressor of gravitation load for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t30 = -t23 * pkin(4) + t25 * qJ(5);
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t8 = g(1) * t24 - g(2) * t26;
t9 = g(1) * t26 + g(2) * t24;
t6 = -g(3) * t23 + t9 * t25;
t49 = g(3) * t25;
t47 = t26 * pkin(7);
t19 = pkin(9) + qJ(6);
t12 = sin(t19);
t46 = t24 * t12;
t13 = cos(t19);
t45 = t24 * t13;
t20 = sin(pkin(9));
t44 = t24 * t20;
t21 = cos(pkin(9));
t43 = t24 * t21;
t42 = t26 * t12;
t41 = t26 * t13;
t40 = t26 * t20;
t39 = t26 * t21;
t38 = -pkin(1) - qJ(3);
t37 = t26 * pkin(1) + t24 * qJ(2);
t35 = t26 * qJ(3) + t37;
t34 = pkin(5) * t20 + pkin(7);
t33 = g(2) * t35;
t16 = t26 * qJ(2);
t32 = t38 * t24 + t16;
t11 = t21 * pkin(5) + pkin(4);
t22 = -pkin(8) - qJ(5);
t28 = t23 * t11 + t25 * t22;
t7 = t8 * t25;
t5 = t9 * t23 + t49;
t4 = t23 * t41 - t46;
t3 = -t23 * t42 - t45;
t2 = -t23 * t45 - t42;
t1 = t23 * t46 - t41;
t10 = [0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, -g(1) * (-t24 * pkin(1) + t16) - g(2) * t37, 0, 0, 0, 0, 0, 0, 0, -t9, t8, -g(1) * t32 - t33, 0, 0, 0, 0, 0, 0, t8 * t23, t7, t9, -g(1) * (t32 - t47) - g(2) * (-t24 * pkin(7) + t35) 0, 0, 0, 0, 0, 0, -g(1) * (-t23 * t43 - t40) - g(2) * (t23 * t39 - t44) -g(1) * (t23 * t44 - t39) - g(2) * (-t23 * t40 - t43) -t7, -g(1) * (t16 - t47) - g(2) * (-t30 * t26 + t35) + (-g(1) * (t30 + t38) + g(2) * pkin(7)) * t24, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3, -t7, -g(1) * t16 - t33 + (g(1) * t34 - g(2) * t28) * t26 + (-g(1) * (-t28 + t38) + g(2) * t34) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t21, t6 * t20, -t5, -g(3) * t30 - t9 * (pkin(4) * t25 + qJ(5) * t23) 0, 0, 0, 0, 0, 0, -t6 * t13, t6 * t12, -t5, g(3) * t28 - t9 * (t11 * t25 - t22 * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t12 * t49, g(1) * t4 - g(2) * t2 + t13 * t49, 0, 0;];
taug_reg  = t10;
