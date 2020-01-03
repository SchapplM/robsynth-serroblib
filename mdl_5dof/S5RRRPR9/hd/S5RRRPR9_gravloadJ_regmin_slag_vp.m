% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t21 = sin(qJ(1));
t24 = cos(qJ(1));
t30 = g(1) * t24 + g(2) * t21;
t22 = cos(qJ(3));
t23 = cos(qJ(2));
t19 = sin(qJ(3));
t34 = t24 * t19;
t11 = t21 * t22 - t23 * t34;
t20 = sin(qJ(2));
t39 = g(3) * t20;
t33 = t24 * t22;
t37 = t21 * t23;
t9 = t19 * t37 + t33;
t44 = -g(1) * t11 + g(2) * t9 + t19 * t39;
t7 = -g(3) * t23 + t30 * t20;
t17 = qJ(3) + pkin(9) + qJ(5);
t14 = sin(t17);
t36 = t24 * t14;
t15 = cos(t17);
t35 = t24 * t15;
t31 = pkin(3) * t19 + pkin(6);
t29 = g(1) * t21 - g(2) * t24;
t16 = t22 * pkin(3) + pkin(2);
t18 = -qJ(4) - pkin(7);
t28 = t23 * t16 - t20 * t18;
t26 = pkin(1) + t28;
t13 = t29 * t20;
t12 = t21 * t19 + t23 * t33;
t10 = -t22 * t37 + t34;
t8 = t30 * t23 + t39;
t6 = t21 * t14 + t23 * t35;
t5 = t21 * t15 - t23 * t36;
t4 = -t15 * t37 + t36;
t3 = t14 * t37 + t35;
t2 = g(1) * t6 - g(2) * t4 + t15 * t39;
t1 = -g(1) * t5 + g(2) * t3 + t14 * t39;
t25 = [0, t29, t30, 0, 0, 0, 0, 0, t29 * t23, -t13, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, (-g(1) * t31 - g(2) * t26) * t24 + (g(1) * t26 - g(2) * t31) * t21, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t7 * t22, -t7 * t19, -t8, -g(3) * t28 + t30 * (t16 * t20 + t18 * t23), 0, 0, 0, 0, 0, t7 * t15, -t7 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, g(1) * t12 - g(2) * t10 + t22 * t39, 0, t44 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t25;
