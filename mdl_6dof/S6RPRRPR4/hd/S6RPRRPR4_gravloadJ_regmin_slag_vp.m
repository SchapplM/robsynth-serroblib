% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t30 = pkin(10) + qJ(3);
t25 = cos(t30);
t26 = qJ(4) + t30;
t20 = sin(t26);
t21 = cos(t26);
t41 = t21 * pkin(4) + t20 * qJ(5);
t53 = pkin(3) * t25 + t41;
t35 = sin(qJ(1));
t36 = cos(qJ(1));
t16 = g(1) * t36 + g(2) * t35;
t5 = -g(3) * t21 + t16 * t20;
t52 = pkin(4) * t20;
t51 = g(3) * t20;
t29 = pkin(11) + qJ(6);
t22 = sin(t29);
t49 = t22 * t36;
t24 = cos(t29);
t48 = t24 * t36;
t31 = sin(pkin(11));
t47 = t31 * t36;
t33 = cos(pkin(11));
t46 = t33 * t36;
t45 = t35 * t22;
t44 = t35 * t24;
t43 = t35 * t31;
t42 = t35 * t33;
t40 = qJ(5) * t21;
t23 = sin(t30);
t39 = -pkin(3) * t23 - t52;
t15 = g(1) * t35 - g(2) * t36;
t34 = cos(pkin(10));
t38 = pkin(2) * t34 + pkin(1) + t53;
t28 = -pkin(8) - pkin(7) - qJ(2);
t14 = t36 * t40;
t13 = t35 * t40;
t11 = t15 * t20;
t10 = t21 * t48 + t45;
t9 = -t21 * t49 + t44;
t8 = -t21 * t44 + t49;
t7 = t21 * t45 + t48;
t6 = t16 * t21 + t51;
t4 = t5 * t33;
t3 = t5 * t31;
t2 = t5 * t24;
t1 = t5 * t22;
t12 = [0, t15, t16, t15 * t34, -t15 * sin(pkin(10)) -t16, -g(1) * (-t35 * pkin(1) + qJ(2) * t36) - g(2) * (pkin(1) * t36 + t35 * qJ(2)) 0, 0, 0, 0, 0, t15 * t25, -t15 * t23, 0, 0, 0, 0, 0, t15 * t21, -t11, -g(1) * (-t21 * t42 + t47) - g(2) * (t21 * t46 + t43) -g(1) * (t21 * t43 + t46) - g(2) * (-t21 * t47 + t42) t11 (g(1) * t28 - g(2) * t38) * t36 + (g(1) * t38 + g(2) * t28) * t35, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t25 + t16 * t23, g(3) * t23 + t16 * t25, 0, 0, 0, 0, 0, t5, t6, t4, -t3, -t6, -g(1) * (t39 * t36 + t14) - g(2) * (t39 * t35 + t13) - g(3) * t53, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t4, -t3, -t6, -g(1) * (-t36 * t52 + t14) - g(2) * (-t35 * t52 + t13) - g(3) * t41, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t22 * t51, g(1) * t10 - g(2) * t8 + t24 * t51;];
taug_reg  = t12;
