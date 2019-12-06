% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t22 = sin(qJ(3));
t19 = sin(pkin(8));
t36 = g(1) * t19;
t20 = cos(pkin(8));
t24 = cos(qJ(3));
t25 = cos(qJ(1));
t28 = t25 * t24;
t23 = sin(qJ(1));
t33 = t23 * t22;
t7 = t20 * t33 + t28;
t29 = t25 * t22;
t32 = t23 * t24;
t9 = t20 * t29 - t32;
t40 = -g(2) * t7 + g(3) * t9 + t22 * t36;
t37 = pkin(3) * t22;
t17 = qJ(3) + pkin(9) + qJ(5);
t14 = sin(t17);
t35 = t23 * t14;
t15 = cos(t17);
t34 = t23 * t15;
t31 = t25 * t14;
t30 = t25 * t15;
t13 = g(2) * t25 + g(3) * t23;
t12 = g(2) * t23 - g(3) * t25;
t26 = (t24 * pkin(3) + pkin(2)) * t20 - t19 * (-qJ(4) - pkin(6)) + pkin(1);
t18 = t25 * qJ(2);
t11 = t13 * t19;
t10 = -t20 * t28 - t33;
t8 = t20 * t32 - t29;
t6 = -t20 * t30 - t35;
t5 = t20 * t31 - t34;
t4 = t20 * t34 - t31;
t3 = t20 * t35 + t30;
t2 = -g(2) * t4 - g(3) * t6 + t15 * t36;
t1 = -g(2) * t3 + g(3) * t5 + t14 * t36;
t16 = [0, t13, -t12, t13 * t20, -t11, t12, -g(2) * (-t25 * pkin(1) - t23 * qJ(2)) - g(3) * (-t23 * pkin(1) + t18), 0, 0, 0, 0, 0, -g(2) * t10 + g(3) * t8, -g(2) * t9 - g(3) * t7, t11, -g(3) * t18 + (g(2) * t26 - g(3) * t37) * t25 + (-g(2) * (-qJ(2) - t37) + g(3) * t26) * t23, 0, 0, 0, 0, 0, -g(2) * t6 + g(3) * t4, -g(2) * t5 - g(3) * t3; 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -g(2) * t8 - g(3) * t10 + t24 * t36, 0, t40 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t20 + t12 * t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t16;
