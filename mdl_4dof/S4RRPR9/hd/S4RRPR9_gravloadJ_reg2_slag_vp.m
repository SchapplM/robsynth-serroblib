% Calculate inertial parameters regressor of gravitation load for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = -pkin(6) - qJ(3);
t20 = sin(qJ(2));
t22 = cos(qJ(2));
t18 = cos(pkin(7));
t9 = t18 * pkin(3) + pkin(2);
t27 = -t20 * t19 + t22 * t9;
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t8 = g(1) * t23 + g(2) * t21;
t5 = -g(3) * t22 + t8 * t20;
t42 = g(1) * t21;
t39 = g(3) * t20;
t35 = t21 * t22;
t16 = pkin(7) + qJ(4);
t10 = sin(t16);
t34 = t23 * t10;
t11 = cos(t16);
t33 = t23 * t11;
t17 = sin(pkin(7));
t32 = t23 * t17;
t31 = t23 * t18;
t30 = t23 * pkin(1) + t21 * pkin(5);
t29 = -g(2) * t23 + t42;
t26 = t22 * pkin(2) + t20 * qJ(3);
t13 = t23 * pkin(5);
t7 = t29 * t20;
t6 = t8 * t22 + t39;
t4 = t21 * t10 + t22 * t33;
t3 = t21 * t11 - t22 * t34;
t2 = -t11 * t35 + t34;
t1 = t10 * t35 + t33;
t12 = [0, 0, 0, 0, 0, 0, t29, t8, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t22, -t7, -t8, -g(1) * (-t21 * pkin(1) + t13) - g(2) * t30, 0, 0, 0, 0, 0, 0, -g(1) * (-t18 * t35 + t32) - g(2) * (t21 * t17 + t22 * t31), -g(1) * (t17 * t35 + t31) - g(2) * (t21 * t18 - t22 * t32), t7, -g(1) * t13 - g(2) * (t26 * t23 + t30) - (-pkin(1) - t26) * t42, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3, t7, -g(1) * (pkin(3) * t32 + t13) - g(2) * (t27 * t23 + t30) + (-g(1) * (-pkin(1) - t27) - g(2) * pkin(3) * t17) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t18, -t5 * t17, -t6, -g(3) * t26 + t8 * (pkin(2) * t20 - qJ(3) * t22), 0, 0, 0, 0, 0, 0, t5 * t11, -t5 * t10, -t6, -g(3) * t27 + t8 * (t19 * t22 + t20 * t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t10 * t39, g(1) * t4 - g(2) * t2 + t11 * t39, 0, 0;];
taug_reg = t12;
