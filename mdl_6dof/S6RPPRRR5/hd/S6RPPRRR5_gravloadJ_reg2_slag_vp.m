% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = qJ(4) + qJ(5);
t20 = sin(t27);
t21 = cos(t27);
t60 = pkin(5) * t21 + pkin(9) * t20;
t40 = -t20 * pkin(5) + t21 * pkin(9);
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t11 = g(1) * t33 + g(2) * t30;
t59 = -g(3) * t20 + t11 * t21;
t32 = cos(qJ(4));
t58 = pkin(4) * t32;
t55 = g(1) * t30;
t53 = g(3) * t21;
t29 = sin(qJ(4));
t51 = t29 * pkin(4);
t28 = sin(qJ(6));
t50 = t30 * t28;
t31 = cos(qJ(6));
t49 = t30 * t31;
t48 = t33 * t28;
t47 = t33 * t31;
t46 = -pkin(1) - qJ(3);
t45 = t60 * t30;
t44 = t60 * t33;
t43 = t33 * pkin(1) + t30 * qJ(2);
t42 = t33 * qJ(3) + t43;
t24 = t33 * qJ(2);
t34 = -pkin(8) - pkin(7);
t41 = g(1) * (t33 * t34 + t24);
t39 = t46 - t51;
t38 = t30 * t34 + t33 * t51 + t42;
t10 = -g(2) * t33 + t55;
t36 = t46 * t30 + t24;
t35 = g(3) * t29 - t11 * t32;
t9 = t20 * t47 - t50;
t8 = -t20 * t48 - t49;
t7 = -t20 * t49 - t48;
t6 = t20 * t50 - t47;
t5 = t10 * t21;
t3 = t11 * t20 + t53;
t2 = t59 * t31;
t1 = t59 * t28;
t4 = [0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, -g(1) * (-t30 * pkin(1) + t24) - g(2) * t43, 0, 0, 0, 0, 0, 0, 0, -t11, t10, -g(1) * t36 - g(2) * t42, 0, 0, 0, 0, 0, 0, t10 * t29, t10 * t32, t11, -g(1) * (-t33 * pkin(7) + t36) - g(2) * (-t30 * pkin(7) + t42) 0, 0, 0, 0, 0, 0, t10 * t20, t5, t11, -g(2) * t38 - t39 * t55 - t41, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, -t5, -t41 - g(2) * (-t33 * t40 + t38) - (t40 + t39) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, g(3) * t32 + t11 * t29, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t3, 0, t35 * pkin(4), 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * (t33 * t58 + t44) - g(2) * (t30 * t58 + t45) - g(3) * (t40 - t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * t44 - g(2) * t45 - g(3) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t6 + t28 * t53, g(1) * t9 - g(2) * t7 + t31 * t53, 0, 0;];
taug_reg  = t4;
