% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t22 = sin(qJ(5));
t34 = pkin(8) + qJ(4);
t29 = sin(t34);
t30 = cos(t34);
t36 = sin(qJ(1));
t37 = cos(qJ(1));
t3 = -t36 * t29 - t37 * t30;
t4 = t37 * t29 - t36 * t30;
t26 = g(1) * t4 - g(2) * t3;
t39 = t26 * t22;
t23 = cos(qJ(5));
t38 = t26 * t23;
t35 = t37 * pkin(1) + t36 * qJ(2);
t20 = sin(pkin(8));
t33 = t37 * t20;
t32 = t36 * t20;
t21 = cos(pkin(8));
t15 = t21 * pkin(3) + pkin(2);
t31 = pkin(3) * t32 + t37 * t15 + t35;
t28 = -t4 * pkin(4) - t3 * pkin(7);
t27 = t3 * pkin(4) - t4 * pkin(7);
t2 = g(1) * t3 + g(2) * t4;
t25 = -t36 * pkin(1) + t37 * qJ(2);
t24 = pkin(3) * t33 - t36 * t15 + t25;
t9 = g(1) * t37 + g(2) * t36;
t8 = g(1) * t36 - g(2) * t37;
t6 = t37 * t21 + t32;
t5 = -t36 * t21 + t33;
t1 = [0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, -t9, -g(1) * t25 - g(2) * t35, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t6, -g(1) * t6 + g(2) * t5, 0, -g(1) * (-t36 * pkin(2) + t25) - g(2) * (t37 * pkin(2) + t35), 0, 0, 0, 0, 0, 0, -t26, t2, 0, -g(1) * t24 - g(2) * t31, 0, 0, 0, 0, 0, 0, -t38, t39, -t2, -g(1) * (t24 - t28) - g(2) * (-t27 + t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t2, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, t2, -g(1) * t28 - g(2) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t23 - t2 * t22, -g(3) * t22 - t2 * t23, 0, 0;];
taug_reg = t1;
