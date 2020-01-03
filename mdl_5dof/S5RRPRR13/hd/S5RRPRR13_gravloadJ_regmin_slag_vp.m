% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR13
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
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t23 = sin(qJ(1));
t25 = cos(qJ(1));
t31 = g(1) * t25 + g(2) * t23;
t22 = sin(qJ(2));
t24 = cos(qJ(2));
t11 = -g(3) * t24 + t31 * t22;
t40 = g(3) * t22;
t38 = t23 * t24;
t19 = pkin(9) + qJ(4);
t18 = qJ(5) + t19;
t14 = sin(t18);
t37 = t25 * t14;
t15 = cos(t18);
t36 = t25 * t15;
t16 = sin(t19);
t35 = t25 * t16;
t17 = cos(t19);
t34 = t25 * t17;
t20 = sin(pkin(9));
t33 = t25 * t20;
t21 = cos(pkin(9));
t32 = t25 * t21;
t30 = g(1) * t23 - g(2) * t25;
t29 = t24 * pkin(2) + t22 * qJ(3);
t27 = pkin(1) + t29;
t13 = t30 * t22;
t12 = t31 * t24 + t40;
t10 = t23 * t16 + t24 * t34;
t9 = t23 * t17 - t24 * t35;
t8 = -t17 * t38 + t35;
t7 = t16 * t38 + t34;
t6 = t23 * t14 + t24 * t36;
t5 = t23 * t15 - t24 * t37;
t4 = -t15 * t38 + t37;
t3 = t14 * t38 + t36;
t2 = g(1) * t6 - g(2) * t4 + t15 * t40;
t1 = -g(1) * t5 + g(2) * t3 + t14 * t40;
t26 = [0, t30, t31, 0, 0, 0, 0, 0, t30 * t24, -t13, -g(1) * (-t21 * t38 + t33) - g(2) * (t23 * t20 + t24 * t32), -g(1) * (t20 * t38 + t32) - g(2) * (t23 * t21 - t24 * t33), t13, (-g(1) * pkin(6) - g(2) * t27) * t25 + (-g(2) * pkin(6) + g(1) * t27) * t23, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, t11 * t21, -t11 * t20, -t12, -g(3) * t29 + t31 * (pkin(2) * t22 - qJ(3) * t24), 0, 0, 0, 0, 0, t11 * t17, -t11 * t16, 0, 0, 0, 0, 0, t11 * t15, -t11 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t16 * t40, g(1) * t10 - g(2) * t8 + t17 * t40, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t26;
