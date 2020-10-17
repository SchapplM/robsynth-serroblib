% Calculate inertial parameters regressor of gravitation load for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:15
% EndTime: 2019-12-31 17:26:15
% DurationCPUTime: 0.19s
% Computational Cost: add. (145->49), mult. (205->74), div. (0->0), fcn. (198->8), ass. (0->37)
t20 = qJ(2) + qJ(3);
t17 = sin(t20);
t18 = cos(t20);
t41 = t18 * pkin(3) + t17 * pkin(7);
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t10 = g(1) * t26 + g(2) * t23;
t3 = -g(3) * t18 + t10 * t17;
t40 = pkin(3) * t17;
t39 = pkin(7) * t18;
t38 = g(3) * t17;
t21 = sin(qJ(4));
t36 = t23 * t21;
t24 = cos(qJ(4));
t35 = t23 * t24;
t34 = t26 * t21;
t33 = t26 * t24;
t22 = sin(qJ(2));
t31 = -pkin(2) * t22 - t40;
t29 = g(1) * t23 - g(2) * t26;
t25 = cos(qJ(2));
t28 = -g(3) * t25 + t10 * t22;
t27 = -pkin(6) - pkin(5);
t19 = t25 * pkin(2);
t16 = t19 + pkin(1);
t13 = t26 * t39;
t12 = t23 * t39;
t11 = t26 * t16;
t9 = t18 * t33 + t36;
t8 = -t18 * t34 + t35;
t7 = -t18 * t35 + t34;
t6 = t18 * t36 + t33;
t5 = t29 * t17;
t4 = t10 * t18 + t38;
t2 = t3 * t24;
t1 = t3 * t21;
t14 = [0, 0, 0, 0, 0, 0, t29, t10, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t25, -t29 * t22, -t10, -g(1) * (-t23 * pkin(1) + t26 * pkin(5)) - g(2) * (t26 * pkin(1) + t23 * pkin(5)), 0, 0, 0, 0, 0, 0, t29 * t18, -t5, -t10, -g(1) * (-t23 * t16 - t26 * t27) - g(2) * (-t23 * t27 + t11), 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t5, -g(2) * t11 + (g(1) * t27 - g(2) * t41) * t26 + (-g(1) * (-t16 - t41) + g(2) * t27) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, g(3) * t22 + t10 * t25, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t28 * pkin(2), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t26 * t31 + t13) - g(2) * (t23 * t31 + t12) - g(3) * (t19 + t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t26 * t40 + t13) - g(2) * (-t23 * t40 + t12) - g(3) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t6 + t21 * t38, g(1) * t9 - g(2) * t7 + t24 * t38, 0, 0;];
taug_reg = t14;
