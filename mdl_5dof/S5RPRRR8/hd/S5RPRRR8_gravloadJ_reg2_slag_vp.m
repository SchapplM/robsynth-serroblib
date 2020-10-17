% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:06
% EndTime: 2019-12-31 19:06:06
% DurationCPUTime: 0.20s
% Computational Cost: add. (160->44), mult. (280->52), div. (0->0), fcn. (327->8), ass. (0->33)
t20 = qJ(4) + qJ(5);
t13 = sin(t20);
t34 = sin(qJ(3));
t35 = sin(qJ(1));
t36 = cos(qJ(3));
t37 = cos(qJ(1));
t5 = -t35 * t34 - t37 * t36;
t6 = t37 * t34 - t35 * t36;
t29 = g(1) * t6 - g(2) * t5;
t41 = t29 * t13;
t14 = cos(t20);
t40 = t29 * t14;
t21 = sin(qJ(4));
t39 = t29 * t21;
t22 = cos(qJ(4));
t38 = t29 * t22;
t33 = t37 * pkin(1) + t35 * qJ(2);
t32 = t37 * pkin(2) + t33;
t31 = -t6 * pkin(3) - t5 * pkin(7);
t30 = t5 * pkin(3) - t6 * pkin(7);
t4 = g(1) * t5 + g(2) * t6;
t28 = -t35 * pkin(1) + t37 * qJ(2);
t12 = t22 * pkin(4) + pkin(3);
t23 = -pkin(8) - pkin(7);
t27 = -t6 * t12 + t5 * t23;
t26 = t5 * t12 + t6 * t23;
t25 = -t35 * pkin(2) + t28;
t24 = g(3) * t22 - t4 * t21;
t8 = g(1) * t37 + g(2) * t35;
t7 = g(1) * t35 - g(2) * t37;
t2 = -g(3) * t13 - t4 * t14;
t1 = g(3) * t14 - t4 * t13;
t3 = [0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t8, -g(1) * t28 - g(2) * t33, 0, 0, 0, 0, 0, 0, -t29, t4, 0, -g(1) * t25 - g(2) * t32, 0, 0, 0, 0, 0, 0, -t38, t39, -t4, -g(1) * (t25 - t31) - g(2) * (-t30 + t32), 0, 0, 0, 0, 0, 0, -t40, t41, -t4, -g(1) * (t25 - t27) - g(2) * (-t26 + t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t4, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, t4, -g(1) * t31 - g(2) * t30, 0, 0, 0, 0, 0, 0, t40, -t41, t4, -g(1) * t27 - g(2) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -g(3) * t21 - t4 * t22, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t24 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t3;
