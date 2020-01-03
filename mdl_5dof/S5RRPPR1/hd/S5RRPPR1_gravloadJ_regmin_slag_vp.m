% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t23 = qJ(1) + qJ(2);
t17 = pkin(8) + t23;
t11 = sin(t17);
t12 = cos(t17);
t19 = cos(t23);
t14 = pkin(2) * t19;
t31 = t12 * pkin(3) + t11 * qJ(4) + t14;
t18 = sin(t23);
t13 = pkin(2) * t18;
t30 = t11 * pkin(3) - t12 * qJ(4) + t13;
t29 = g(2) * t12 + g(3) * t11;
t28 = g(2) * t11 - g(3) * t12;
t7 = -g(2) * t19 - g(3) * t18;
t27 = cos(qJ(1));
t26 = sin(qJ(1));
t22 = pkin(9) + qJ(5);
t21 = t27 * pkin(1);
t20 = t26 * pkin(1);
t16 = cos(t22);
t15 = sin(t22);
t6 = g(2) * t18 - g(3) * t19;
t4 = t29 * cos(pkin(9));
t3 = t29 * sin(pkin(9));
t2 = t29 * t16;
t1 = t29 * t15;
t5 = [0, -g(2) * t27 - g(3) * t26, g(2) * t26 - g(3) * t27, 0, t7, t6, -g(2) * (t14 + t21) - g(3) * (t13 + t20), -t4, t3, -t28, -g(2) * (t21 + t31) - g(3) * (t20 + t30), 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, t7, t6, t7 * pkin(2), -t4, t3, -t28, -g(2) * t31 - g(3) * t30, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + t15 * t28, g(1) * t15 + t16 * t28;];
taug_reg = t5;
