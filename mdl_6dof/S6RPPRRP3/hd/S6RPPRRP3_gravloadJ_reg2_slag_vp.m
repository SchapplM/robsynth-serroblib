% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = qJ(1) + pkin(9);
t24 = sin(t28);
t25 = cos(t28);
t55 = g(2) * t25;
t11 = g(1) * t24 - t55;
t33 = cos(qJ(4));
t30 = sin(qJ(4));
t54 = g(3) * t30;
t57 = t11 * t33 - t54;
t56 = pkin(8) * t30;
t53 = g(3) * t33;
t26 = t33 * pkin(8);
t52 = t24 * t30;
t51 = t24 * t33;
t50 = t25 * t30;
t29 = sin(qJ(5));
t49 = t29 * t30;
t32 = cos(qJ(5));
t48 = t30 * t32;
t47 = pkin(4) * t51 + pkin(8) * t52;
t34 = cos(qJ(1));
t46 = t34 * pkin(1) + t25 * pkin(2) + t24 * qJ(3);
t31 = sin(qJ(1));
t45 = -t31 * pkin(1) + t25 * qJ(3);
t44 = t25 * pkin(7) + t46;
t6 = t24 * t49 - t25 * t32;
t8 = t24 * t32 + t25 * t49;
t43 = g(1) * t8 + g(2) * t6;
t12 = g(1) * t25 + g(2) * t24;
t42 = g(1) * t31 - g(2) * t34;
t41 = pkin(5) * t32 + qJ(6) * t29;
t40 = -pkin(4) - t41;
t39 = (-pkin(2) - pkin(7)) * t24 + t45;
t1 = g(1) * t6 - g(2) * t8 + t29 * t53;
t7 = t24 * t48 + t25 * t29;
t9 = -t24 * t29 + t25 * t48;
t38 = g(1) * t7 - g(2) * t9 + t32 * t53;
t37 = pkin(4) * t52 - pkin(8) * t51 + t44;
t35 = pkin(4) * t50 - t25 * t26 + t39;
t10 = t12 * t33;
t5 = g(1) * t52 - g(2) * t50 + t53;
t4 = t57 * t32;
t3 = t57 * t29;
t2 = -g(1) * t9 - g(2) * t7;
t13 = [0, 0, 0, 0, 0, 0, t42, g(1) * t34 + g(2) * t31, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, t42 * pkin(1), 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -g(1) * (-t24 * pkin(2) + t45) - g(2) * t46, 0, 0, 0, 0, 0, 0, -t12 * t30, -t10, t11, -g(1) * t39 - g(2) * t44, 0, 0, 0, 0, 0, 0, t2, t43, t10, -g(1) * t35 - g(2) * t37, 0, 0, 0, 0, 0, 0, t2, t10, -t43, -g(1) * (t9 * pkin(5) + t8 * qJ(6) + t35) - g(2) * (t7 * pkin(5) + t6 * qJ(6) + t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(1) * t47 - g(3) * (-t30 * pkin(4) + t26) - (-pkin(4) * t33 - t56) * t55, 0, 0, 0, 0, 0, 0, -t4, -t5, -t3, -g(1) * (t41 * t51 + t47) - g(3) * t26 - t40 * t54 - (t40 * t33 - t56) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t38, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t38, -g(1) * (-t6 * pkin(5) + t7 * qJ(6)) - g(2) * (t8 * pkin(5) - t9 * qJ(6)) - (-pkin(5) * t29 + qJ(6) * t32) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;
