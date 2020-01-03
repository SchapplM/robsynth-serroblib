% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t24 = qJ(2) + pkin(8);
t20 = qJ(4) + t24;
t15 = sin(t20);
t16 = cos(t20);
t37 = t16 * pkin(4) + t15 * qJ(5);
t36 = pkin(4) * t15;
t25 = -qJ(3) - pkin(6);
t19 = cos(t24);
t28 = cos(qJ(2));
t21 = t28 * pkin(2);
t34 = pkin(3) * t19 + t21;
t33 = qJ(5) * t16;
t18 = sin(t24);
t26 = sin(qJ(2));
t7 = -t26 * pkin(2) - pkin(3) * t18;
t32 = t7 - t36;
t27 = sin(qJ(1));
t29 = cos(qJ(1));
t11 = g(1) * t29 + g(2) * t27;
t10 = g(1) * t27 - g(2) * t29;
t30 = -g(3) * t28 + t11 * t26;
t23 = -pkin(7) + t25;
t17 = t21 + pkin(1);
t9 = t29 * t33;
t8 = t27 * t33;
t6 = pkin(1) + t34;
t5 = t29 * t6;
t4 = t10 * t16;
t3 = t10 * t15;
t2 = g(3) * t15 + t11 * t16;
t1 = -g(3) * t16 + t11 * t15;
t12 = [0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t28, -t10 * t26, -t11, -g(1) * (-t27 * pkin(1) + t29 * pkin(6)) - g(2) * (t29 * pkin(1) + t27 * pkin(6)), 0, 0, 0, 0, 0, 0, t10 * t19, -t10 * t18, -t11, -g(1) * (-t27 * t17 - t29 * t25) - g(2) * (t29 * t17 - t27 * t25), 0, 0, 0, 0, 0, 0, t4, -t3, -t11, -g(1) * (-t29 * t23 - t27 * t6) - g(2) * (-t27 * t23 + t5), 0, 0, 0, 0, 0, 0, t4, -t11, t3, -g(2) * t5 + (g(1) * t23 - g(2) * t37) * t29 + (-g(1) * (-t37 - t6) + g(2) * t23) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, g(3) * t26 + t11 * t28, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t19 + t11 * t18, g(3) * t18 + t11 * t19, 0, t30 * pkin(2), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t34 - t11 * t7, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (t32 * t29 + t9) - g(2) * (t32 * t27 + t8) - g(3) * (t34 + t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-t29 * t36 + t9) - g(2) * (-t27 * t36 + t8) - g(3) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t12;
