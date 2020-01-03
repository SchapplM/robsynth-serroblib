% Calculate inertial parameters regressor of gravitation load for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t13 = sin(qJ(1));
t20 = t13 * pkin(1);
t10 = qJ(1) + pkin(7);
t5 = sin(t10);
t6 = cos(t10);
t3 = g(1) * t6 + g(2) * t5;
t19 = g(1) * t5 - g(2) * t6;
t15 = cos(qJ(1));
t18 = g(1) * t13 - g(2) * t15;
t12 = sin(qJ(3));
t14 = cos(qJ(3));
t17 = -g(3) * t14 + t3 * t12;
t16 = -pkin(6) - pkin(5);
t11 = qJ(3) + qJ(4);
t9 = t15 * pkin(1);
t8 = cos(t11);
t7 = sin(t11);
t4 = t14 * pkin(3) + pkin(2);
t2 = g(3) * t7 + t3 * t8;
t1 = -g(3) * t8 + t3 * t7;
t21 = [0, 0, 0, 0, 0, 0, t18, g(1) * t15 + g(2) * t13, 0, 0, 0, 0, 0, 0, 0, 0, t19, t3, 0, t18 * pkin(1), 0, 0, 0, 0, 0, 0, t19 * t14, -t19 * t12, -t3, -g(1) * (-t5 * pkin(2) + t6 * pkin(5) - t20) - g(2) * (t6 * pkin(2) + t5 * pkin(5) + t9), 0, 0, 0, 0, 0, 0, t19 * t8, -t19 * t7, -t3, -g(1) * (-t6 * t16 - t5 * t4 - t20) - g(2) * (-t5 * t16 + t6 * t4 + t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, g(3) * t12 + t3 * t14, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t17 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t21;
