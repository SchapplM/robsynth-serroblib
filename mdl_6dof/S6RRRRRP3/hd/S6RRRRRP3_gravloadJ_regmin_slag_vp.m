% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t37 = cos(qJ(2));
t31 = qJ(4) + qJ(5);
t26 = cos(t31);
t36 = cos(qJ(4));
t21 = pkin(4) * t36 + pkin(5) * t26;
t19 = pkin(3) + t21;
t32 = qJ(2) + qJ(3);
t25 = sin(t32);
t27 = cos(t32);
t30 = -qJ(6) - pkin(10) - pkin(9);
t46 = t19 * t27 - t25 * t30;
t65 = pkin(2) * t37 + t46;
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t45 = g(1) * t38 + g(2) * t35;
t24 = sin(t31);
t52 = t38 * t24;
t55 = t35 * t26;
t11 = -t27 * t52 + t55;
t60 = g(3) * t25;
t51 = t38 * t26;
t56 = t35 * t24;
t9 = t27 * t56 + t51;
t1 = -g(1) * t11 + g(2) * t9 + t24 * t60;
t7 = -g(3) * t27 + t25 * t45;
t33 = sin(qJ(4));
t20 = pkin(4) * t33 + pkin(5) * t24;
t58 = t20 * t27;
t54 = t35 * t33;
t53 = t35 * t36;
t50 = t38 * t33;
t49 = t38 * t36;
t48 = t20 + pkin(8) + pkin(7);
t44 = g(1) * t35 - g(2) * t38;
t43 = t19 * t25 + t27 * t30;
t42 = pkin(1) + t65;
t34 = sin(qJ(2));
t17 = t27 * t49 + t54;
t16 = -t27 * t50 + t53;
t15 = -t27 * t53 + t50;
t14 = t27 * t54 + t49;
t13 = t44 * t25;
t12 = t27 * t51 + t56;
t10 = -t27 * t55 + t52;
t8 = t27 * t45 + t60;
t6 = t7 * t36;
t5 = t7 * t33;
t4 = t7 * t26;
t3 = t7 * t24;
t2 = g(1) * t12 - g(2) * t10 + t26 * t60;
t18 = [0, t44, t45, 0, 0, 0, 0, 0, t44 * t37, -t44 * t34, 0, 0, 0, 0, 0, t44 * t27, -t13, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13 (-g(1) * t48 - g(2) * t42) * t38 + (g(1) * t42 - g(2) * t48) * t35; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t37 + t34 * t45, g(3) * t34 + t37 * t45, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t65 + t45 * (pkin(2) * t34 + t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t46 + t43 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + g(2) * t14 + t33 * t60, g(1) * t17 - g(2) * t15 + t36 * t60, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t21 * t35 - t38 * t58) - g(2) * (-t21 * t38 - t35 * t58) + t20 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg  = t18;
