% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:57
% EndTime: 2019-12-31 18:40:58
% DurationCPUTime: 0.17s
% Computational Cost: add. (221->42), mult. (152->45), div. (0->0), fcn. (137->8), ass. (0->29)
t21 = sin(qJ(4));
t23 = cos(qJ(4));
t28 = t23 * pkin(4) + t21 * qJ(5);
t20 = qJ(1) + pkin(8);
t18 = qJ(3) + t20;
t15 = cos(t18);
t41 = t28 * t15;
t14 = sin(t18);
t40 = g(1) * t14;
t5 = -g(2) * t15 + t40;
t6 = g(1) * t15 + g(2) * t14;
t36 = t14 * pkin(3);
t34 = t15 * pkin(3) + t14 * pkin(7);
t17 = cos(t20);
t24 = cos(qJ(1));
t33 = t24 * pkin(1) + pkin(2) * t17;
t31 = t33 + t34;
t16 = sin(t20);
t22 = sin(qJ(1));
t30 = -t22 * pkin(1) - pkin(2) * t16;
t29 = g(1) * t22 - g(2) * t24;
t11 = t15 * pkin(7);
t26 = t11 + t30;
t25 = (-pkin(3) - t28) * t40;
t4 = t5 * t23;
t3 = t5 * t21;
t2 = g(3) * t21 + t6 * t23;
t1 = -g(3) * t23 + t6 * t21;
t7 = [0, 0, 0, 0, 0, 0, t29, g(1) * t24 + g(2) * t22, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t16 - g(2) * t17, g(1) * t17 + g(2) * t16, 0, t29 * pkin(1), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(1) * t30 - g(2) * t33, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t26 - t36) - g(2) * t31, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t26 - g(2) * (t31 + t41) - t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t11 - t36) - g(2) * t34, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t11 - g(2) * (t34 + t41) - t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * t28 + t6 * (pkin(4) * t21 - qJ(5) * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t7;
