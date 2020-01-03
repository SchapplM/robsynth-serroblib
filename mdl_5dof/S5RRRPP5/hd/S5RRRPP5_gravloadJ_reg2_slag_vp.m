% Calculate inertial parameters regressor of gravitation load for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t20 = qJ(2) + qJ(3);
t17 = sin(t20);
t18 = cos(t20);
t34 = t18 * pkin(3) + t17 * qJ(4);
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t6 = g(1) * t24 + g(2) * t22;
t41 = t6 * t17;
t21 = sin(qJ(2));
t39 = pkin(2) * t21;
t38 = pkin(3) * t17;
t13 = t18 * pkin(4);
t25 = -pkin(7) - pkin(6);
t35 = t24 * t25;
t33 = qJ(4) * t18;
t32 = qJ(5) + t25;
t23 = cos(qJ(2));
t19 = t23 * pkin(2);
t31 = t19 + t34;
t16 = t19 + pkin(1);
t10 = t24 * t16;
t30 = g(2) * (t34 * t24 + t10);
t29 = -t38 - t39;
t5 = g(1) * t22 - g(2) * t24;
t28 = -t16 - t34;
t27 = -g(3) * t23 + t6 * t21;
t26 = (pkin(3) + pkin(4)) * t41;
t9 = t24 * t33;
t7 = t22 * t33;
t4 = t5 * t18;
t3 = t5 * t17;
t2 = g(3) * t17 + t6 * t18;
t1 = -g(3) * t18 + t41;
t8 = [0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t23, -t5 * t21, -t6, -g(1) * (-t22 * pkin(1) + t24 * pkin(6)) - g(2) * (t24 * pkin(1) + t22 * pkin(6)), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t22 * t16 - t35) - g(2) * (-t22 * t25 + t10), 0, 0, 0, 0, 0, 0, t4, -t6, t3, g(1) * t35 - t30 + (-g(1) * t28 + g(2) * t25) * t22, 0, 0, 0, 0, 0, 0, t4, t3, t6, -t30 + (g(1) * t32 - g(2) * t13) * t24 + (-g(1) * (t28 - t13) + g(2) * t32) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, g(3) * t21 + t6 * t23, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t27 * pkin(2), 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (t29 * t24 + t9) - g(2) * (t29 * t22 + t7) - g(3) * t31, 0, 0, 0, 0, 0, 0, t1, -t2, 0, -g(1) * (-t24 * t39 + t9) - g(2) * (-t22 * t39 + t7) - g(3) * (t13 + t31) + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-t24 * t38 + t9) - g(2) * (-t22 * t38 + t7) - g(3) * t34, 0, 0, 0, 0, 0, 0, t1, -t2, 0, -g(1) * t9 - g(2) * t7 - g(3) * (t13 + t34) + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5;];
taug_reg = t8;
