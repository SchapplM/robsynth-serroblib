% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t12 = g(1) * t33 + g(2) * t30;
t32 = cos(qJ(4));
t29 = sin(qJ(4));
t56 = g(3) * t29;
t61 = -t12 * t32 + t56;
t55 = g(3) * t32;
t54 = t29 * pkin(4);
t53 = t33 * pkin(7);
t52 = t29 * t33;
t28 = sin(qJ(5));
t51 = t30 * t28;
t31 = cos(qJ(5));
t50 = t30 * t31;
t49 = t30 * t32;
t48 = t32 * t33;
t47 = t33 * t28;
t46 = t33 * t31;
t45 = -pkin(1) - qJ(3);
t44 = t33 * pkin(1) + t30 * qJ(2);
t43 = t33 * qJ(3) + t44;
t24 = t33 * qJ(2);
t42 = pkin(8) * t49 + t24 - t53;
t6 = t29 * t51 - t46;
t8 = t29 * t47 + t50;
t41 = g(1) * t6 - g(2) * t8;
t11 = g(1) * t30 - g(2) * t33;
t40 = t45 * t30 + t24;
t38 = pkin(4) * t52 - pkin(8) * t48 + t43;
t1 = g(1) * t8 + g(2) * t6 + t28 * t55;
t7 = t29 * t50 + t47;
t9 = t29 * t46 - t51;
t37 = g(1) * t9 + g(2) * t7 + t31 * t55;
t36 = -g(1) * (pkin(4) * t48 + pkin(8) * t52) - g(2) * (t30 * t29 * pkin(8) + pkin(4) * t49);
t34 = (-g(1) * (t45 - t54) + g(2) * pkin(7)) * t30;
t25 = t32 * pkin(8);
t10 = g(1) * t49 - g(2) * t48;
t5 = t12 * t29 + t55;
t4 = t61 * t31;
t3 = t61 * t28;
t2 = g(1) * t7 - g(2) * t9;
t13 = [0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -g(1) * (-t30 * pkin(1) + t24) - g(2) * t44, 0, 0, 0, 0, 0, 0, 0, -t12, t11, -g(1) * t40 - g(2) * t43, 0, 0, 0, 0, 0, 0, t11 * t29, t10, t12, -g(1) * (t40 - t53) - g(2) * (-t30 * pkin(7) + t43) 0, 0, 0, 0, 0, 0, t2, -t41, -t10, -g(1) * t42 - g(2) * t38 + t34, 0, 0, 0, 0, 0, 0, t2, -t10, t41, -g(1) * (-t7 * pkin(5) - t6 * qJ(6) + t42) - g(2) * (t9 * pkin(5) + t8 * qJ(6) + t38) + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t5, -g(3) * (t25 - t54) + t36, 0, 0, 0, 0, 0, 0, t4, -t5, t3, pkin(4) * t56 - g(3) * t25 + t36 + t61 * (pkin(5) * t31 + qJ(6) * t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t37, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t37, -g(1) * (-t8 * pkin(5) + t9 * qJ(6)) - g(2) * (-t6 * pkin(5) + t7 * qJ(6)) - (-pkin(5) * t28 + qJ(6) * t31) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;
