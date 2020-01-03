% Calculate inertial parameters regressor of gravitation load for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = qJ(1) + pkin(7);
t16 = sin(t19);
t35 = g(1) * t16;
t17 = cos(t19);
t21 = cos(pkin(8));
t34 = t17 * t21;
t20 = sin(pkin(8));
t33 = qJ(4) * t20;
t25 = cos(qJ(1));
t32 = t25 * pkin(1) + t17 * pkin(2) + t16 * qJ(3);
t23 = sin(qJ(1));
t31 = -t23 * pkin(1) + t17 * qJ(3);
t30 = pkin(3) * t34 + t17 * t33 + t32;
t9 = g(1) * t17 + g(2) * t16;
t8 = -g(2) * t17 + t35;
t29 = g(1) * t23 - g(2) * t25;
t22 = sin(qJ(5));
t24 = cos(qJ(5));
t28 = t20 * t24 - t21 * t22;
t27 = t20 * t22 + t21 * t24;
t26 = -pkin(3) * t21 - pkin(2) - t33;
t7 = t8 * t21;
t6 = t8 * t20;
t5 = t27 * t17;
t4 = t28 * t17;
t3 = t27 * t16;
t2 = t28 * t16;
t1 = g(3) * t21 - t9 * t20;
t10 = [0, 0, 0, 0, 0, 0, t29, g(1) * t25 + g(2) * t23, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t29 * pkin(1), 0, 0, 0, 0, 0, 0, t7, -t6, -t9, -g(1) * (-t16 * pkin(2) + t31) - g(2) * t32, 0, 0, 0, 0, 0, 0, t7, -t9, t6, -g(1) * t31 - g(2) * t30 - t26 * t35, 0, 0, 0, 0, 0, 0, g(1) * t3 - g(2) * t5, g(1) * t2 - g(2) * t4, t9, -g(1) * (-t17 * pkin(6) + t31) - g(2) * (pkin(4) * t34 + t30) + (-g(1) * (-pkin(4) * t21 + t26) + g(2) * pkin(6)) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2 + g(3) * t27, g(1) * t5 + g(2) * t3 + g(3) * t28, 0, 0;];
taug_reg = t10;
