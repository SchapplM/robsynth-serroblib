% Calculate inertial parameters regressor of gravitation load for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t20 = -pkin(6) - pkin(5);
t15 = qJ(2) + qJ(3);
t12 = cos(t15);
t18 = cos(qJ(2));
t13 = t18 * pkin(2);
t22 = pkin(3) * t12 + t13;
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t8 = g(1) * t19 + g(2) * t17;
t7 = g(1) * t17 - g(2) * t19;
t11 = sin(t15);
t1 = -g(3) * t12 + t8 * t11;
t16 = sin(qJ(2));
t21 = -g(3) * t18 + t8 * t16;
t14 = -qJ(4) + t20;
t10 = t13 + pkin(1);
t5 = pkin(1) + t22;
t4 = t7 * t12;
t3 = t7 * t11;
t2 = g(3) * t11 + t8 * t12;
t6 = [0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t18, -t7 * t16, -t8, -g(1) * (-t17 * pkin(1) + t19 * pkin(5)) - g(2) * (t19 * pkin(1) + t17 * pkin(5)), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (-t17 * t10 - t19 * t20) - g(2) * (t19 * t10 - t17 * t20), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (-t19 * t14 - t17 * t5) - g(2) * (-t17 * t14 + t19 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, g(3) * t16 + t8 * t18, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t21 * pkin(2), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t22 - t8 * (-t16 * pkin(2) - pkin(3) * t11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t6;
