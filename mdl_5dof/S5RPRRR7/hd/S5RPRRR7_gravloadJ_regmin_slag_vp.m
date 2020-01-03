% Calculate minimal parameter regressor of gravitation load for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t18 = sin(qJ(3));
t21 = cos(qJ(3));
t15 = qJ(1) + pkin(9);
t11 = sin(t15);
t12 = cos(t15);
t27 = g(1) * t12 + g(2) * t11;
t24 = -g(3) * t21 + t27 * t18;
t33 = g(3) * t18;
t16 = qJ(4) + qJ(5);
t13 = sin(t16);
t31 = t13 * t21;
t14 = cos(t16);
t30 = t14 * t21;
t17 = sin(qJ(4));
t29 = t17 * t21;
t20 = cos(qJ(4));
t28 = t20 * t21;
t26 = g(1) * t11 - g(2) * t12;
t19 = sin(qJ(1));
t22 = cos(qJ(1));
t25 = g(1) * t19 - g(2) * t22;
t10 = t11 * t17 + t12 * t28;
t9 = t11 * t20 - t12 * t29;
t8 = -t11 * t28 + t12 * t17;
t7 = t11 * t29 + t12 * t20;
t6 = t11 * t13 + t12 * t30;
t5 = t11 * t14 - t12 * t31;
t4 = -t11 * t30 + t12 * t13;
t3 = t11 * t31 + t12 * t14;
t2 = g(1) * t6 - g(2) * t4 + t14 * t33;
t1 = -g(1) * t5 + g(2) * t3 + t13 * t33;
t23 = [0, t25, g(1) * t22 + g(2) * t19, t25 * pkin(1), 0, 0, 0, 0, 0, t26 * t21, -t26 * t18, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t27 * t21 + t33, 0, 0, 0, 0, 0, t24 * t20, -t24 * t17, 0, 0, 0, 0, 0, t24 * t14, -t24 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t17 * t33, g(1) * t10 - g(2) * t8 + t20 * t33, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t23;
