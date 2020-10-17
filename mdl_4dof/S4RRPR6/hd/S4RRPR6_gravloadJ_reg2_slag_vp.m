% Calculate inertial parameters regressor of gravitation load for
% S4RRPR6
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:54
% EndTime: 2019-12-31 17:04:55
% DurationCPUTime: 0.11s
% Computational Cost: add. (104->34), mult. (121->47), div. (0->0), fcn. (109->8), ass. (0->22)
t17 = qJ(2) + pkin(7);
t12 = cos(t17);
t21 = cos(qJ(2));
t14 = t21 * pkin(2);
t24 = pkin(3) * t12 + t14;
t18 = -qJ(3) - pkin(5);
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t6 = g(1) * t22 + g(2) * t20;
t5 = g(1) * t20 - g(2) * t22;
t19 = sin(qJ(2));
t23 = -g(3) * t21 + t6 * t19;
t16 = -pkin(6) + t18;
t13 = qJ(4) + t17;
t11 = sin(t17);
t10 = t14 + pkin(1);
t9 = cos(t13);
t8 = sin(t13);
t3 = pkin(1) + t24;
t2 = g(3) * t8 + t6 * t9;
t1 = -g(3) * t9 + t6 * t8;
t4 = [0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t21, -t5 * t19, -t6, -g(1) * (-t20 * pkin(1) + t22 * pkin(5)) - g(2) * (t22 * pkin(1) + t20 * pkin(5)), 0, 0, 0, 0, 0, 0, t5 * t12, -t5 * t11, -t6, -g(1) * (-t20 * t10 - t22 * t18) - g(2) * (t22 * t10 - t20 * t18), 0, 0, 0, 0, 0, 0, t5 * t9, -t5 * t8, -t6, -g(1) * (-t22 * t16 - t20 * t3) - g(2) * (-t20 * t16 + t22 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, g(3) * t19 + t6 * t21, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t12 + t6 * t11, g(3) * t11 + t6 * t12, 0, t23 * pkin(2), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t24 - t6 * (-t19 * pkin(2) - pkin(3) * t11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t4;
