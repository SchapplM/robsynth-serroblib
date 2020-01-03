% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t10 = qJ(1) + pkin(8);
t6 = sin(t10);
t8 = cos(t10);
t20 = g(1) * t8 + g(2) * t6;
t9 = pkin(9) + qJ(4);
t5 = sin(t9);
t7 = cos(t9);
t17 = -g(3) * t7 + t20 * t5;
t26 = g(3) * t5;
t13 = sin(qJ(5));
t24 = t6 * t13;
t15 = cos(qJ(5));
t23 = t6 * t15;
t22 = t8 * t13;
t21 = t8 * t15;
t19 = g(1) * t6 - g(2) * t8;
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t18 = g(1) * t14 - g(2) * t16;
t4 = t7 * t21 + t24;
t3 = -t7 * t22 + t23;
t2 = -t7 * t23 + t22;
t1 = t7 * t24 + t21;
t11 = [0, t18, g(1) * t16 + g(2) * t14, t18 * pkin(1), t19 * cos(pkin(9)), -t19 * sin(pkin(9)), -t20, -g(1) * (-t14 * pkin(1) - t6 * pkin(2) + t8 * qJ(3)) - g(2) * (t16 * pkin(1) + t8 * pkin(2) + t6 * qJ(3)), 0, 0, 0, 0, 0, t19 * t7, -t19 * t5, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3; 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t20 * t7 + t26, 0, 0, 0, 0, 0, t17 * t15, -t17 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t13 * t26, g(1) * t4 - g(2) * t2 + t15 * t26;];
taug_reg = t11;
