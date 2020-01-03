% Calculate inertial parameters regressor of gravitation load for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = pkin(7) + qJ(3);
t16 = sin(t19);
t17 = cos(t19);
t29 = t17 * pkin(3) + t16 * qJ(4);
t23 = sin(qJ(1));
t24 = cos(qJ(1));
t10 = g(1) * t24 + g(2) * t23;
t36 = t10 * t16;
t35 = pkin(3) * t16;
t22 = -pkin(6) - qJ(2);
t32 = pkin(4) - t22;
t31 = t24 * t22;
t28 = qJ(4) * t17;
t27 = t17 * qJ(5);
t21 = cos(pkin(7));
t15 = t21 * pkin(2) + pkin(1);
t8 = t24 * t15;
t26 = g(2) * (t29 * t24 + t8);
t9 = g(1) * t23 - g(2) * t24;
t25 = -t15 - t29;
t7 = t24 * t28;
t5 = t23 * t28;
t4 = t9 * t17;
t3 = t9 * t16;
t2 = g(3) * t16 + t10 * t17;
t1 = -g(3) * t17 + t36;
t6 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t21, -t9 * sin(pkin(7)), -t10, -g(1) * (-t23 * pkin(1) + t24 * qJ(2)) - g(2) * (t24 * pkin(1) + t23 * qJ(2)), 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(1) * (-t23 * t15 - t31) - g(2) * (-t23 * t22 + t8), 0, 0, 0, 0, 0, 0, -t10, -t4, t3, g(1) * t31 - t26 + (-g(1) * t25 + g(2) * t22) * t23, 0, 0, 0, 0, 0, 0, -t10, t3, t4, -t26 + (-g(1) * t32 - g(2) * t27) * t24 + (-g(1) * (t25 - t27) - g(2) * t32) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (-t24 * t35 + t7) - g(2) * (-t23 * t35 + t5) - g(3) * t29, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -g(1) * t7 - g(2) * t5 - g(3) * (t27 + t29) + (pkin(3) + qJ(5)) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2;];
taug_reg = t6;
