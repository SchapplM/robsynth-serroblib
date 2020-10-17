% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:40
% EndTime: 2019-12-31 20:15:41
% DurationCPUTime: 0.19s
% Computational Cost: add. (196->47), mult. (166->47), div. (0->0), fcn. (149->8), ass. (0->34)
t23 = qJ(1) + qJ(2);
t18 = sin(t23);
t20 = cos(t23);
t28 = -pkin(8) - pkin(7);
t24 = sin(qJ(4));
t38 = pkin(4) * t24;
t39 = t18 * t28 + t20 * t38;
t25 = sin(qJ(1));
t37 = t25 * pkin(1);
t36 = t20 * pkin(2) + t18 * qJ(3);
t27 = cos(qJ(1));
t21 = t27 * pkin(1);
t35 = t21 + t36;
t13 = t20 * qJ(3);
t34 = -t18 * pkin(2) + t13;
t33 = t13 + (-pkin(2) - pkin(7)) * t18;
t8 = g(1) * t20 + g(2) * t18;
t7 = g(1) * t18 - g(2) * t20;
t32 = g(1) * t25 - g(2) * t27;
t31 = t18 * t38 - t20 * t28 + t36;
t30 = t34 - t37;
t26 = cos(qJ(4));
t29 = g(3) * t24 - t7 * t26;
t22 = qJ(4) + qJ(5);
t19 = cos(t22);
t17 = sin(t22);
t14 = t20 * pkin(7);
t6 = t8 * t26;
t5 = t8 * t24;
t4 = t8 * t19;
t3 = t8 * t17;
t2 = g(3) * t17 - t7 * t19;
t1 = g(3) * t19 + t7 * t17;
t9 = [0, 0, 0, 0, 0, 0, t32, g(1) * t27 + g(2) * t25, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t32 * pkin(1), 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * t30 - g(2) * t35, 0, 0, 0, 0, 0, 0, -t5, -t6, t7, -g(1) * (t33 - t37) - g(2) * (t14 + t35), 0, 0, 0, 0, 0, 0, -t3, -t4, t7, -g(1) * (t30 + t39) - g(2) * (t21 + t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * t34 - g(2) * t36, 0, 0, 0, 0, 0, 0, -t5, -t6, t7, -g(1) * t33 - g(2) * (t14 + t36), 0, 0, 0, 0, 0, 0, -t3, -t4, t7, -g(1) * (t34 + t39) - g(2) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, g(3) * t26 + t7 * t24, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t29 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t9;
