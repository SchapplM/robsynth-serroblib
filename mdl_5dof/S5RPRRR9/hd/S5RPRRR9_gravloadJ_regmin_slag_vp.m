% Calculate minimal parameter regressor of gravitation load for
% S5RPRRR9
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
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t10 = g(1) * t22 + g(2) * t20;
t16 = pkin(9) + qJ(3);
t15 = qJ(4) + t16;
t11 = sin(t15);
t12 = cos(t15);
t3 = -g(3) * t12 + t10 * t11;
t28 = g(3) * t11;
t19 = sin(qJ(5));
t26 = t20 * t19;
t21 = cos(qJ(5));
t25 = t20 * t21;
t24 = t22 * t19;
t23 = t22 * t21;
t9 = g(1) * t20 - g(2) * t22;
t14 = cos(t16);
t13 = sin(t16);
t8 = t12 * t23 + t26;
t7 = -t12 * t24 + t25;
t6 = -t12 * t25 + t24;
t5 = t12 * t26 + t23;
t4 = t10 * t12 + t28;
t2 = t3 * t21;
t1 = t3 * t19;
t17 = [0, t9, t10, t9 * cos(pkin(9)), -t9 * sin(pkin(9)), -t10, -g(1) * (-t20 * pkin(1) + t22 * qJ(2)) - g(2) * (t22 * pkin(1) + t20 * qJ(2)), 0, 0, 0, 0, 0, t9 * t14, -t9 * t13, 0, 0, 0, 0, 0, t9 * t12, -t9 * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t14 + t10 * t13, g(3) * t13 + t10 * t14, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 + t19 * t28, g(1) * t8 - g(2) * t6 + t21 * t28;];
taug_reg = t17;
