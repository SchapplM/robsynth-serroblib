% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:54
% EndTime: 2020-01-03 12:05:55
% DurationCPUTime: 0.12s
% Computational Cost: add. (191->38), mult. (172->61), div. (0->0), fcn. (196->10), ass. (0->38)
t29 = sin(pkin(9));
t41 = g(1) * t29;
t28 = qJ(1) + qJ(2);
t24 = sin(t28);
t30 = cos(pkin(9));
t40 = t24 * t30;
t26 = cos(t28);
t39 = t26 * t30;
t31 = sin(qJ(4));
t38 = t30 * t31;
t33 = cos(qJ(4));
t37 = t30 * t33;
t36 = t26 * pkin(2) + t24 * qJ(3);
t35 = t24 * pkin(2) - t26 * qJ(3);
t18 = g(2) * t26 + g(3) * t24;
t34 = cos(qJ(1));
t32 = sin(qJ(1));
t27 = qJ(4) + qJ(5);
t25 = cos(t27);
t23 = sin(t27);
t17 = g(2) * t24 - g(3) * t26;
t16 = t18 * t30;
t15 = t18 * t29;
t14 = t24 * t31 + t26 * t37;
t13 = -t24 * t33 + t26 * t38;
t12 = t24 * t37 - t26 * t31;
t11 = -t24 * t38 - t26 * t33;
t10 = t24 * t23 + t25 * t39;
t9 = t23 * t39 - t24 * t25;
t8 = -t26 * t23 + t25 * t40;
t7 = -t23 * t40 - t26 * t25;
t6 = -g(2) * t14 - g(3) * t12;
t5 = g(2) * t13 - g(3) * t11;
t4 = -g(2) * t10 - g(3) * t8;
t3 = g(2) * t9 - g(3) * t7;
t2 = g(2) * t8 - g(3) * t10 + t25 * t41;
t1 = -g(2) * t7 - g(3) * t9 + t23 * t41;
t19 = [0, -g(2) * t34 - g(3) * t32, g(2) * t32 - g(3) * t34, 0, -t18, t17, -t16, t15, -t17, -g(2) * (t34 * pkin(1) + t36) - g(3) * (t32 * pkin(1) + t35), 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, -t18, t17, -t16, t15, -t17, -g(2) * t36 - g(3) * t35, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t11 - g(3) * t13 + t31 * t41, g(2) * t12 - g(3) * t14 + t33 * t41, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t19;
