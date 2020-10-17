% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:43
% EndTime: 2019-12-31 19:49:43
% DurationCPUTime: 0.18s
% Computational Cost: add. (237->43), mult. (164->46), div. (0->0), fcn. (147->8), ass. (0->33)
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t30 = t25 * pkin(4) + t23 * qJ(5);
t22 = qJ(1) + qJ(2);
t18 = pkin(8) + t22;
t16 = cos(t18);
t15 = sin(t18);
t42 = g(1) * t15;
t5 = -g(2) * t16 + t42;
t6 = g(1) * t16 + g(2) * t15;
t19 = sin(t22);
t43 = pkin(2) * t19;
t24 = sin(qJ(1));
t38 = t24 * pkin(1);
t20 = cos(t22);
t17 = pkin(2) * t20;
t35 = t16 * pkin(3) + t15 * pkin(7) + t17;
t13 = t16 * pkin(7);
t34 = t13 - t43;
t33 = t30 * t16 + t35;
t32 = -t38 - t43;
t7 = g(1) * t19 - g(2) * t20;
t26 = cos(qJ(1));
t31 = g(1) * t24 - g(2) * t26;
t28 = -t15 * pkin(3) + t34;
t27 = (-pkin(3) - t30) * t42;
t21 = t26 * pkin(1);
t8 = g(1) * t20 + g(2) * t19;
t4 = t5 * t25;
t3 = t5 * t23;
t2 = g(3) * t23 + t6 * t25;
t1 = -g(3) * t25 + t6 * t23;
t9 = [0, 0, 0, 0, 0, 0, t31, g(1) * t26 + g(2) * t24, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t31 * pkin(1), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(1) * t32 - g(2) * (t17 + t21), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t28 - t38) - g(2) * (t21 + t35), 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (t13 + t32) - g(2) * (t21 + t33) - t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t28 - g(2) * t35, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t34 - g(2) * t33 - t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * t30 + t6 * (pkin(4) * t23 - qJ(5) * t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t9;
