% Calculate inertial parameters regressor of gravitation load for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t30 = t25 * pkin(4) + t23 * qJ(5);
t22 = qJ(1) + qJ(2);
t20 = qJ(3) + t22;
t17 = cos(t20);
t16 = sin(t20);
t44 = g(1) * t16;
t5 = -g(2) * t17 + t44;
t6 = g(1) * t17 + g(2) * t16;
t18 = sin(t22);
t45 = pkin(2) * t18;
t24 = sin(qJ(1));
t40 = t24 * pkin(1);
t38 = t17 * pkin(3) + t16 * pkin(8);
t19 = cos(t22);
t15 = pkin(2) * t19;
t36 = t15 + t38;
t13 = t17 * pkin(8);
t35 = -t16 * pkin(3) + t13;
t34 = t30 * t17 + t38;
t33 = t15 + t34;
t32 = -t40 - t45;
t7 = g(1) * t18 - g(2) * t19;
t26 = cos(qJ(1));
t31 = g(1) * t24 - g(2) * t26;
t28 = t35 - t45;
t27 = (-pkin(3) - t30) * t44;
t21 = t26 * pkin(1);
t8 = g(1) * t19 + g(2) * t18;
t4 = t5 * t25;
t3 = t5 * t23;
t2 = g(3) * t23 + t6 * t25;
t1 = -g(3) * t25 + t6 * t23;
t9 = [0, 0, 0, 0, 0, 0, t31, g(1) * t26 + g(2) * t24, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t31 * pkin(1), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(1) * t32 - g(2) * (t15 + t21), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t28 - t40) - g(2) * (t21 + t36), 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (t13 + t32) - g(2) * (t21 + t33) - t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t28 - g(2) * t36, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (t13 - t45) - g(2) * t33 - t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t35 - g(2) * t38, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t13 - g(2) * t34 - t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * t30 + t6 * (pkin(4) * t23 - qJ(5) * t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t9;
