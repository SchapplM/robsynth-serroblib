% Calculate inertial parameters regressor of gravitation load for
% S4RRRP2
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t18 = sin(qJ(1));
t26 = t18 * pkin(1);
t15 = qJ(1) + qJ(2);
t12 = sin(t15);
t13 = cos(t15);
t25 = t13 * pkin(2) + t12 * pkin(6);
t24 = -t12 * pkin(2) + t13 * pkin(6);
t19 = cos(qJ(3));
t11 = t19 * pkin(3) + pkin(2);
t16 = -qJ(4) - pkin(6);
t23 = t13 * t11 - t12 * t16;
t6 = g(1) * t13 + g(2) * t12;
t5 = g(1) * t12 - g(2) * t13;
t20 = cos(qJ(1));
t22 = g(1) * t18 - g(2) * t20;
t21 = -t12 * t11 - t13 * t16;
t17 = sin(qJ(3));
t1 = -g(3) * t19 + t6 * t17;
t14 = t20 * pkin(1);
t4 = t5 * t19;
t3 = t5 * t17;
t2 = g(3) * t17 + t6 * t19;
t7 = [0, 0, 0, 0, 0, 0, t22, g(1) * t20 + g(2) * t18, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t22 * pkin(1), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t24 - t26) - g(2) * (t14 + t25), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t21 - t26) - g(2) * (t14 + t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t24 - g(2) * t25, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t21 - g(2) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t7;
