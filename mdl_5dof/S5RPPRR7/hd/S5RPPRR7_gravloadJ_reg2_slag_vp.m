% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:53
% EndTime: 2019-12-31 17:59:54
% DurationCPUTime: 0.21s
% Computational Cost: add. (142->45), mult. (152->60), div. (0->0), fcn. (149->8), ass. (0->28)
t19 = sin(qJ(4));
t22 = cos(qJ(4));
t27 = t19 * pkin(4) - t22 * pkin(7);
t17 = qJ(1) + pkin(8);
t14 = sin(t17);
t15 = cos(t17);
t41 = -g(1) * t14 + g(2) * t15;
t40 = -g(3) * t19 - t22 * t41;
t39 = -pkin(2) - pkin(6);
t35 = g(3) * t22;
t18 = sin(qJ(5));
t32 = t18 * t19;
t21 = cos(qJ(5));
t31 = t19 * t21;
t23 = cos(qJ(1));
t30 = t23 * pkin(1) + t15 * pkin(2) + t14 * qJ(3);
t20 = sin(qJ(1));
t29 = -t20 * pkin(1) + t15 * qJ(3);
t8 = g(1) * t15 + g(2) * t14;
t26 = g(1) * t20 - g(2) * t23;
t25 = g(2) * (t15 * pkin(6) + t30);
t6 = t8 * t22;
t5 = -t14 * t18 + t15 * t31;
t4 = t14 * t21 + t15 * t32;
t3 = t14 * t31 + t15 * t18;
t2 = -t14 * t32 + t15 * t21;
t1 = -t19 * t41 + t35;
t7 = [0, 0, 0, 0, 0, 0, t26, g(1) * t23 + g(2) * t20, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t8, 0, t26 * pkin(1), 0, 0, 0, 0, 0, 0, 0, t41, -t8, -g(1) * (-t14 * pkin(2) + t29) - g(2) * t30, 0, 0, 0, 0, 0, 0, -t8 * t19, -t6, -t41, -g(1) * (t39 * t14 + t29) - t25, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t3, g(1) * t4 - g(2) * t2, t6, -g(1) * (t27 * t15 + t29) - t25 + (-g(1) * t39 - g(2) * t27) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t40 * t21, t40 * t18, -t1, g(3) * t27 + t41 * (pkin(4) * t22 + pkin(7) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4 + t18 * t35, g(1) * t3 - g(2) * t5 + t21 * t35, 0, 0;];
taug_reg = t7;
