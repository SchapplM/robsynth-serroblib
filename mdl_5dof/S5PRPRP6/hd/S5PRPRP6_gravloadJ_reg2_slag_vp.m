% Calculate inertial parameters regressor of gravitation load for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t20 = sin(qJ(2));
t40 = (pkin(2) + pkin(6)) * t20;
t17 = sin(pkin(7));
t18 = cos(pkin(7));
t39 = -g(1) * t18 - g(2) * t17;
t22 = cos(qJ(2));
t5 = g(3) * t20 - t39 * t22;
t37 = pkin(2) * t20;
t33 = g(3) * t22;
t19 = sin(qJ(4));
t32 = t19 * t20;
t21 = cos(qJ(4));
t31 = t20 * t21;
t30 = t22 * pkin(2) + t20 * qJ(3);
t29 = qJ(3) * t22;
t28 = t22 * pkin(6) + t30;
t11 = t17 * t29;
t12 = t18 * t29;
t27 = -g(1) * t12 - g(2) * t11;
t25 = pkin(4) * t19 - qJ(5) * t21;
t7 = t17 * t21 + t18 * t32;
t9 = -t17 * t32 + t18 * t21;
t24 = -g(1) * t7 + g(2) * t9 + t19 * t33;
t6 = t17 * t19 - t18 * t31;
t8 = t17 * t31 + t18 * t19;
t1 = g(1) * t6 - g(2) * t8 + t21 * t33;
t4 = -t20 * t39 - t33;
t3 = t5 * t21;
t2 = t5 * t19;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t5, -g(1) * (-t18 * t37 + t12) - g(2) * (-t17 * t37 + t11) - g(3) * t30, 0, 0, 0, 0, 0, 0, -t2, -t3, t4, -g(3) * t28 - t39 * t40 + t27, 0, 0, 0, 0, 0, 0, -t2, t4, t3, -g(3) * (t25 * t20 + t28) + t27 + t39 * (t25 * t22 - t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t24, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, t24, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (t8 * pkin(4) - t9 * qJ(5)) - (-pkin(4) * t21 - qJ(5) * t19) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t10;
