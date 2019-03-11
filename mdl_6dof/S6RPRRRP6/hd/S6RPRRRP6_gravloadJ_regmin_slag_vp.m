% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t18 = g(1) * t36 + g(2) * t34;
t29 = qJ(4) + qJ(5);
t23 = sin(t29);
t28 = pkin(10) + qJ(3);
t22 = cos(t28);
t47 = t34 * t23;
t24 = cos(t29);
t49 = t24 * t36;
t5 = t22 * t47 + t49;
t21 = sin(t28);
t53 = g(3) * t21;
t46 = t34 * t24;
t50 = t23 * t36;
t7 = -t22 * t50 + t46;
t1 = -g(1) * t7 + g(2) * t5 + t23 * t53;
t3 = -g(3) * t22 + t18 * t21;
t33 = sin(qJ(4));
t15 = pkin(4) * t33 + pkin(5) * t23;
t51 = t15 * t22;
t48 = t33 * t36;
t45 = t34 * t33;
t35 = cos(qJ(4));
t44 = t34 * t35;
t43 = t35 * t36;
t42 = t15 + pkin(7) + qJ(2);
t16 = t35 * pkin(4) + pkin(5) * t24;
t17 = g(1) * t34 - g(2) * t36;
t14 = pkin(3) + t16;
t27 = -qJ(6) - pkin(9) - pkin(8);
t40 = t14 * t22 - t21 * t27;
t31 = cos(pkin(10));
t38 = pkin(2) * t31 + pkin(1) + t40;
t13 = t22 * t43 + t45;
t12 = -t22 * t48 + t44;
t11 = -t22 * t44 + t48;
t10 = t22 * t45 + t43;
t9 = t17 * t21;
t8 = t22 * t49 + t47;
t6 = -t22 * t46 + t50;
t4 = t18 * t22 + t53;
t2 = g(1) * t8 - g(2) * t6 + t24 * t53;
t19 = [0, t17, t18, t17 * t31, -t17 * sin(pkin(10)) -t18, -g(1) * (-t34 * pkin(1) + qJ(2) * t36) - g(2) * (pkin(1) * t36 + t34 * qJ(2)) 0, 0, 0, 0, 0, t17 * t22, -t9, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13, -g(1) * t10 - g(2) * t12, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, t9 (-g(1) * t42 - g(2) * t38) * t36 + (g(1) * t38 - g(2) * t42) * t34; 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t3 * t35, -t3 * t33, 0, 0, 0, 0, 0, t3 * t24, -t3 * t23, -t4, -g(3) * t40 + t18 * (t14 * t21 + t22 * t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 + g(2) * t10 + t33 * t53, g(1) * t13 - g(2) * t11 + t35 * t53, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t34 * t16 - t36 * t51) - g(2) * (-t16 * t36 - t34 * t51) + t15 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3;];
taug_reg  = t19;
