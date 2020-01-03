% Calculate inertial parameters regressor of gravitation load for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t15 = sin(qJ(5));
t25 = pkin(8) + qJ(2);
t22 = sin(t25);
t23 = cos(t25);
t27 = sin(qJ(4));
t28 = cos(qJ(4));
t3 = -t22 * t27 - t23 * t28;
t4 = -t22 * t28 + t23 * t27;
t19 = g(1) * t4 - g(2) * t3;
t30 = t19 * t15;
t16 = cos(qJ(5));
t29 = t19 * t16;
t26 = t23 * pkin(2) + t22 * qJ(3);
t24 = t23 * pkin(3) + t26;
t21 = -t4 * pkin(4) - t3 * pkin(7);
t20 = t3 * pkin(4) - t4 * pkin(7);
t2 = g(1) * t3 + g(2) * t4;
t18 = -t22 * pkin(2) + t23 * qJ(3);
t17 = -t22 * pkin(3) + t18;
t6 = g(1) * t23 + g(2) * t22;
t5 = g(1) * t22 - g(2) * t23;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * t18 - g(2) * t26, 0, 0, 0, 0, 0, 0, -t19, t2, 0, -g(1) * t17 - g(2) * t24, 0, 0, 0, 0, 0, 0, -t29, t30, -t2, -g(1) * (t17 - t21) - g(2) * (-t20 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t2, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t30, t2, -g(1) * t21 - g(2) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t16 - t2 * t15, -g(3) * t15 - t2 * t16, 0, 0;];
taug_reg = t1;
