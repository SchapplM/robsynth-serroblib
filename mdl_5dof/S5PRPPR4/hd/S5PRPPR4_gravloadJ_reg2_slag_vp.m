% Calculate inertial parameters regressor of gravitation load for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:55
% EndTime: 2019-12-31 17:36:55
% DurationCPUTime: 0.13s
% Computational Cost: add. (126->38), mult. (128->48), div. (0->0), fcn. (132->6), ass. (0->26)
t18 = pkin(7) + qJ(2);
t16 = sin(t18);
t30 = g(1) * t16;
t17 = cos(t18);
t20 = cos(pkin(8));
t29 = t17 * t20;
t28 = t17 * pkin(2) + t16 * qJ(3);
t19 = sin(pkin(8));
t27 = qJ(4) * t19;
t26 = pkin(3) * t29 + t17 * t27 + t28;
t9 = g(1) * t17 + g(2) * t16;
t8 = -g(2) * t17 + t30;
t21 = sin(qJ(5));
t22 = cos(qJ(5));
t25 = t19 * t22 - t20 * t21;
t24 = t19 * t21 + t20 * t22;
t23 = -pkin(3) * t20 - pkin(2) - t27;
t13 = t17 * qJ(3);
t7 = t8 * t20;
t6 = t8 * t19;
t5 = t24 * t17;
t4 = t25 * t17;
t3 = t24 * t16;
t2 = t25 * t16;
t1 = g(3) * t20 - t9 * t19;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, -t9, -g(1) * (-t16 * pkin(2) + t13) - g(2) * t28, 0, 0, 0, 0, 0, 0, t7, -t9, t6, -g(1) * t13 - g(2) * t26 - t23 * t30, 0, 0, 0, 0, 0, 0, g(1) * t3 - g(2) * t5, g(1) * t2 - g(2) * t4, t9, -g(1) * (-t17 * pkin(6) + t13) - g(2) * (pkin(4) * t29 + t26) + (-g(1) * (-pkin(4) * t20 + t23) + g(2) * pkin(6)) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2 + g(3) * t24, g(1) * t5 + g(2) * t3 + g(3) * t25, 0, 0;];
taug_reg = t10;
