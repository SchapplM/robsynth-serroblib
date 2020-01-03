% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR16_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t21 = sin(qJ(2));
t22 = sin(qJ(1));
t25 = cos(qJ(2));
t26 = cos(qJ(1));
t35 = cos(pkin(5));
t33 = t26 * t35;
t13 = t21 * t33 + t22 * t25;
t19 = sin(qJ(5));
t23 = cos(qJ(5));
t12 = t22 * t21 - t25 * t33;
t20 = sin(qJ(4));
t24 = cos(qJ(4));
t18 = sin(pkin(5));
t40 = t18 * t26;
t29 = -t12 * t20 + t24 * t40;
t47 = t13 * t23 + t19 * t29;
t46 = -t13 * t19 + t23 * t29;
t45 = g(3) * t18;
t42 = t18 * t22;
t41 = t18 * t25;
t39 = t19 * t20;
t38 = t19 * t21;
t37 = t20 * t23;
t36 = t21 * t23;
t34 = t22 * t35;
t14 = t26 * t21 + t25 * t34;
t32 = g(1) * t12 - g(2) * t14;
t15 = -t21 * t34 + t26 * t25;
t31 = g(1) * t13 - g(2) * t15;
t30 = g(1) * t26 + g(2) * t22;
t6 = t12 * t24 + t20 * t40;
t4 = t14 * t24 - t20 * t42;
t28 = g(1) * t4 + g(2) * t6 + g(3) * (-t35 * t20 - t24 * t41);
t3 = -g(1) * t14 - g(2) * t12 + g(3) * t41;
t27 = g(1) * t15 + g(2) * t13 + t21 * t45;
t11 = -t20 * t41 + t35 * t24;
t5 = t14 * t20 + t24 * t42;
t2 = t15 * t19 + t5 * t23;
t1 = t15 * t23 - t5 * t19;
t7 = [0, g(1) * t22 - g(2) * t26, t30, 0, 0, 0, 0, 0, t31, -t32, -t30 * t18, -t31, t32, -g(1) * (-t22 * pkin(1) - t13 * pkin(2) + pkin(7) * t40 - t12 * qJ(3)) - g(2) * (t26 * pkin(1) + t15 * pkin(2) + pkin(7) * t42 + t14 * qJ(3)), 0, 0, 0, 0, 0, -g(1) * t29 - g(2) * t5, g(1) * t6 - g(2) * t4, 0, 0, 0, 0, 0, -g(1) * t46 - g(2) * t2, g(1) * t47 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t3, t27, 0, t3, -t27, -g(1) * (-t14 * pkin(2) + t15 * qJ(3)) - g(2) * (-t12 * pkin(2) + t13 * qJ(3)) - (pkin(2) * t25 + qJ(3) * t21) * t45, 0, 0, 0, 0, 0, -t27 * t20, -t27 * t24, 0, 0, 0, 0, 0, -g(1) * (-t14 * t19 + t15 * t37) - g(2) * (-t12 * t19 + t13 * t37) - (t19 * t25 + t20 * t36) * t45, -g(1) * (-t14 * t23 - t15 * t39) - g(2) * (-t12 * t23 - t13 * t39) - (-t20 * t38 + t23 * t25) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, g(1) * t5 - g(2) * t29 + g(3) * t11, 0, 0, 0, 0, 0, -t28 * t23, t28 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t47 - g(3) * (-t11 * t19 + t18 * t36), g(1) * t2 - g(2) * t46 - g(3) * (-t11 * t23 - t18 * t38);];
taug_reg = t7;
