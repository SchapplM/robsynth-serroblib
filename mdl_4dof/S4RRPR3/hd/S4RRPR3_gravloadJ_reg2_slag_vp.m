% Calculate inertial parameters regressor of gravitation load for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:34
% DurationCPUTime: 0.13s
% Computational Cost: add. (119->27), mult. (82->34), div. (0->0), fcn. (72->8), ass. (0->24)
t17 = qJ(1) + qJ(2);
t14 = sin(t17);
t26 = pkin(2) * t14;
t19 = sin(qJ(1));
t25 = t19 * pkin(1);
t13 = pkin(7) + t17;
t10 = sin(t13);
t11 = cos(t13);
t15 = cos(t17);
t12 = pkin(2) * t15;
t24 = t11 * pkin(3) + t10 * pkin(6) + t12;
t4 = g(1) * t11 + g(2) * t10;
t3 = g(1) * t10 - g(2) * t11;
t5 = g(1) * t14 - g(2) * t15;
t21 = cos(qJ(1));
t23 = g(1) * t19 - g(2) * t21;
t22 = -t10 * pkin(3) + t11 * pkin(6) - t26;
t20 = cos(qJ(4));
t18 = sin(qJ(4));
t16 = t21 * pkin(1);
t6 = g(1) * t15 + g(2) * t14;
t2 = t3 * t20;
t1 = t3 * t18;
t7 = [0, 0, 0, 0, 0, 0, t23, g(1) * t21 + g(2) * t19, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t23 * pkin(1), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (-t25 - t26) - g(2) * (t12 + t16), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t22 - t25) - g(2) * (t16 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(2), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t22 - g(2) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t20 + t4 * t18, g(3) * t18 + t4 * t20, 0, 0;];
taug_reg = t7;
