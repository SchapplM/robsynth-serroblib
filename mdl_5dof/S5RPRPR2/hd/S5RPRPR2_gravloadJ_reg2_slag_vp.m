% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = qJ(1) + pkin(8);
t23 = qJ(3) + t27;
t16 = sin(t23);
t17 = cos(t23);
t29 = cos(pkin(9));
t18 = t29 * pkin(4) + pkin(3);
t30 = -pkin(7) - qJ(4);
t39 = t16 * t18 + t17 * t30;
t38 = t17 * pkin(3) + t16 * qJ(4);
t20 = sin(t27);
t31 = sin(qJ(1));
t37 = t31 * pkin(1) + pkin(2) * t20;
t22 = cos(t27);
t32 = cos(qJ(1));
t36 = t32 * pkin(1) + pkin(2) * t22;
t35 = -t16 * t30 + t17 * t18;
t34 = t16 * pkin(3) - t17 * qJ(4);
t6 = g(2) * t17 + g(3) * t16;
t5 = g(2) * t16 - g(3) * t17;
t33 = -g(2) * t32 - g(3) * t31;
t26 = pkin(9) + qJ(5);
t21 = cos(t26);
t19 = sin(t26);
t4 = t6 * t29;
t3 = t6 * sin(pkin(9));
t2 = t6 * t21;
t1 = t6 * t19;
t7 = [0, 0, 0, 0, 0, 0, t33, g(2) * t31 - g(3) * t32, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t22 - g(3) * t20, g(2) * t20 - g(3) * t22, 0, t33 * pkin(1), 0, 0, 0, 0, 0, 0, -t6, t5, 0, -g(2) * t36 - g(3) * t37, 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * (t36 + t38) - g(3) * (t34 + t37), 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(2) * (t35 + t36) - g(3) * (t37 + t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * t38 - g(3) * t34, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(2) * t35 - g(3) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t21 + t5 * t19, g(1) * t19 + t5 * t21, 0, 0;];
taug_reg = t7;
