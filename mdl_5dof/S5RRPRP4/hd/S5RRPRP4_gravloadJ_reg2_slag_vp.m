% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:56
% EndTime: 2019-12-31 19:52:56
% DurationCPUTime: 0.16s
% Computational Cost: add. (181->46), mult. (174->43), div. (0->0), fcn. (157->6), ass. (0->27)
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t39 = -t20 * pkin(4) + t22 * qJ(5);
t19 = qJ(1) + qJ(2);
t15 = sin(t19);
t16 = cos(t19);
t38 = -g(1) * t15 + g(2) * t16;
t21 = sin(qJ(1));
t34 = t21 * pkin(1);
t33 = t16 * pkin(2) + t15 * qJ(3);
t12 = t16 * pkin(7);
t31 = t12 + t33;
t23 = cos(qJ(1));
t17 = t23 * pkin(1);
t30 = t17 + t33;
t11 = t16 * qJ(3);
t29 = -t15 * pkin(2) + t11;
t28 = t11 + (-pkin(2) - pkin(7)) * t15;
t6 = g(1) * t16 + g(2) * t15;
t27 = g(1) * t21 - g(2) * t23;
t25 = -t39 * t15 + t31;
t24 = -t39 * t16 + t28;
t4 = t6 * t22;
t3 = t6 * t20;
t2 = -g(3) * t20 - t38 * t22;
t1 = g(3) * t22 - t20 * t38;
t5 = [0, 0, 0, 0, 0, 0, t27, g(1) * t23 + g(2) * t21, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t6, 0, t27 * pkin(1), 0, 0, 0, 0, 0, 0, 0, t38, -t6, -g(1) * (t29 - t34) - g(2) * t30, 0, 0, 0, 0, 0, 0, -t3, -t4, -t38, -g(1) * (t28 - t34) - g(2) * (t12 + t30), 0, 0, 0, 0, 0, 0, -t3, -t38, t4, -g(1) * (t24 - t34) - g(2) * (t17 + t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t6, -g(1) * t29 - g(2) * t33, 0, 0, 0, 0, 0, 0, -t3, -t4, -t38, -g(1) * t28 - g(2) * t31, 0, 0, 0, 0, 0, 0, -t3, -t38, t4, -g(1) * t24 - g(2) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, -t1, -g(3) * t39 + t38 * (pkin(4) * t22 + qJ(5) * t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2;];
taug_reg = t5;
