% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR6
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
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t30 = sin(qJ(1));
t32 = cos(qJ(1));
t15 = g(1) * t32 + g(2) * t30;
t24 = pkin(10) + qJ(3);
t21 = cos(t24);
t31 = cos(qJ(4));
t39 = t32 * t31;
t29 = sin(qJ(4));
t44 = t30 * t29;
t10 = t21 * t44 + t39;
t40 = t32 * t29;
t43 = t30 * t31;
t12 = -t21 * t40 + t43;
t20 = sin(t24);
t48 = g(3) * t20;
t53 = -g(1) * t12 + g(2) * t10 + t29 * t48;
t3 = -g(3) * t21 + t15 * t20;
t22 = qJ(4) + pkin(11) + qJ(6);
t16 = sin(t22);
t46 = t30 * t16;
t17 = cos(t22);
t45 = t30 * t17;
t42 = t32 * t16;
t41 = t32 * t17;
t37 = pkin(4) * t29 + pkin(7) + qJ(2);
t14 = g(1) * t30 - g(2) * t32;
t19 = t31 * pkin(4) + pkin(3);
t27 = -qJ(5) - pkin(8);
t36 = t21 * t19 - t20 * t27;
t26 = cos(pkin(10));
t34 = t26 * pkin(2) + pkin(1) + t36;
t13 = t21 * t39 + t44;
t11 = -t21 * t43 + t40;
t9 = t14 * t20;
t8 = t21 * t41 + t46;
t7 = -t21 * t42 + t45;
t6 = -t21 * t45 + t42;
t5 = t21 * t46 + t41;
t4 = t15 * t21 + t48;
t2 = g(1) * t8 - g(2) * t6 + t17 * t48;
t1 = -g(1) * t7 + g(2) * t5 + t16 * t48;
t18 = [0, t14, t15, t14 * t26, -t14 * sin(pkin(10)) -t15, -g(1) * (-t30 * pkin(1) + t32 * qJ(2)) - g(2) * (t32 * pkin(1) + t30 * qJ(2)) 0, 0, 0, 0, 0, t14 * t21, -t9, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13, -g(1) * t10 - g(2) * t12, t9 (-g(1) * t37 - g(2) * t34) * t32 + (g(1) * t34 - g(2) * t37) * t30, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t3 * t31, -t3 * t29, -t4, -g(3) * t36 + t15 * (t19 * t20 + t21 * t27) 0, 0, 0, 0, 0, t3 * t17, -t3 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, g(1) * t13 - g(2) * t11 + t31 * t48, 0, t53 * pkin(4), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t18;
