% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t27 = sin(pkin(9));
t39 = g(1) * t27;
t26 = qJ(1) + qJ(2);
t22 = sin(t26);
t28 = cos(pkin(9));
t38 = t22 * t28;
t24 = cos(t26);
t37 = t24 * t28;
t29 = sin(qJ(4));
t36 = t28 * t29;
t31 = cos(qJ(4));
t35 = t28 * t31;
t34 = -t22 * pkin(2) + t24 * qJ(3);
t18 = g(2) * t24 + g(3) * t22;
t33 = -t24 * pkin(2) - t22 * qJ(3);
t32 = cos(qJ(1));
t30 = sin(qJ(1));
t25 = qJ(4) + qJ(5);
t23 = cos(t25);
t21 = sin(t25);
t17 = g(2) * t22 - g(3) * t24;
t16 = t18 * t28;
t15 = t18 * t27;
t14 = -t22 * t29 - t24 * t35;
t13 = -t22 * t31 + t24 * t36;
t12 = t22 * t35 - t24 * t29;
t11 = t22 * t36 + t24 * t31;
t10 = -t22 * t21 - t23 * t37;
t9 = t21 * t37 - t22 * t23;
t8 = -t24 * t21 + t23 * t38;
t7 = t21 * t38 + t24 * t23;
t6 = -g(2) * t14 + g(3) * t12;
t5 = -g(2) * t13 - g(3) * t11;
t4 = -g(2) * t10 + g(3) * t8;
t3 = -g(2) * t9 - g(3) * t7;
t2 = -g(2) * t8 - g(3) * t10 + t23 * t39;
t1 = -g(2) * t7 + g(3) * t9 + t21 * t39;
t19 = [0, g(2) * t32 + g(3) * t30, -g(2) * t30 + g(3) * t32, 0, t18, -t17, t16, -t15, t17, -g(2) * (-t32 * pkin(1) + t33) - g(3) * (-t30 * pkin(1) + t34), 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, t18, -t17, t16, -t15, t17, -g(2) * t33 - g(3) * t34, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t11 + g(3) * t13 + t29 * t39, -g(2) * t12 - g(3) * t14 + t31 * t39, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t19;
