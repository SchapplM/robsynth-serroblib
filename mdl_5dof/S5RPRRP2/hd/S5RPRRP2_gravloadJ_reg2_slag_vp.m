% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP2
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:29
% EndTime: 2020-01-03 11:45:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (197->41), mult. (132->43), div. (0->0), fcn. (117->8), ass. (0->26)
t24 = qJ(1) + pkin(8);
t21 = qJ(3) + t24;
t16 = sin(t21);
t17 = cos(t21);
t28 = cos(qJ(4));
t18 = t28 * pkin(4) + pkin(3);
t25 = -qJ(5) - pkin(7);
t36 = t16 * t18 + t17 * t25;
t35 = t17 * pkin(3) + t16 * pkin(7);
t19 = sin(t24);
t27 = sin(qJ(1));
t34 = t27 * pkin(1) + pkin(2) * t19;
t20 = cos(t24);
t29 = cos(qJ(1));
t33 = t29 * pkin(1) + pkin(2) * t20;
t32 = t16 * pkin(3) - t17 * pkin(7);
t31 = -t16 * t25 + t17 * t18;
t6 = g(2) * t17 + g(3) * t16;
t5 = g(2) * t16 - g(3) * t17;
t30 = -g(2) * t29 - g(3) * t27;
t26 = sin(qJ(4));
t2 = -g(1) * t28 + t5 * t26;
t4 = t6 * t28;
t3 = t6 * t26;
t1 = g(1) * t26 + t5 * t28;
t7 = [0, 0, 0, 0, 0, 0, t30, g(2) * t27 - g(3) * t29, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t20 - g(3) * t19, g(2) * t19 - g(3) * t20, 0, t30 * pkin(1), 0, 0, 0, 0, 0, 0, -t6, t5, 0, -g(2) * t33 - g(3) * t34, 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * (t33 + t35) - g(3) * (t32 + t34), 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * (t31 + t33) - g(3) * (t34 + t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * t35 - g(3) * t32, 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * t31 - g(3) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6;];
taug_reg = t7;
