% Calculate minimal parameter regressor of gravitation load for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t29 = qJ(2) + qJ(3);
t23 = pkin(10) + t29;
t20 = sin(t23);
t21 = cos(t23);
t36 = t21 * pkin(4) + t20 * qJ(5);
t24 = sin(t29);
t48 = pkin(3) * t24;
t47 = pkin(4) * t20;
t46 = g(3) * t21;
t30 = sin(qJ(6));
t32 = sin(qJ(1));
t45 = t32 * t30;
t33 = cos(qJ(6));
t44 = t32 * t33;
t35 = cos(qJ(1));
t43 = t35 * t30;
t42 = t35 * t33;
t25 = cos(t29);
t22 = pkin(3) * t25;
t34 = cos(qJ(2));
t26 = t34 * pkin(2);
t41 = t22 + t26;
t40 = qJ(5) * t21;
t39 = t22 + t36;
t31 = sin(qJ(2));
t13 = -t31 * pkin(2) - t48;
t38 = t13 - t47;
t37 = -t47 - t48;
t17 = g(1) * t35 + g(2) * t32;
t16 = g(1) * t32 - g(2) * t35;
t4 = -g(3) * t20 - t17 * t21;
t5 = -g(3) * t25 + t17 * t24;
t28 = -qJ(4) - pkin(8) - pkin(7);
t15 = t35 * t40;
t14 = t32 * t40;
t12 = pkin(1) + t41;
t11 = t35 * t12;
t10 = -t20 * t45 + t42;
t9 = t20 * t44 + t43;
t8 = t20 * t43 + t44;
t7 = t20 * t42 - t45;
t6 = g(3) * t24 + t17 * t25;
t3 = -t17 * t20 + t46;
t2 = t4 * t33;
t1 = t4 * t30;
t18 = [0, t16, t17, 0, 0, 0, 0, 0, t16 * t34, -t16 * t31, 0, 0, 0, 0, 0, t16 * t25, -t16 * t24, -t17, -g(1) * (-t32 * t12 - t35 * t28) - g(2) * (-t32 * t28 + t11) -t17, -t16 * t21, t16 * t20, -g(2) * t11 + (g(1) * t28 - g(2) * t36) * t35 + (-g(1) * (-t12 - t36) + g(2) * t28) * t32, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t34 + t17 * t31, g(3) * t31 + t17 * t34, 0, 0, 0, 0, 0, t5, t6, 0, -g(3) * t41 - t17 * t13, 0, t3, t4, -g(1) * (t38 * t35 + t15) - g(2) * (t38 * t32 + t14) - g(3) * (t26 + t39) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t5 * pkin(3), 0, t3, t4, -g(1) * (t35 * t37 + t15) - g(2) * (t32 * t37 + t14) - g(3) * t39, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t33 * t46, g(1) * t8 - g(2) * t10 - t30 * t46;];
taug_reg  = t18;
