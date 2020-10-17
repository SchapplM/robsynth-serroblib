% Calculate inertial parameters regressor of gravitation load for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:46
% EndTime: 2019-12-31 20:55:46
% DurationCPUTime: 0.22s
% Computational Cost: add. (244->56), mult. (231->71), div. (0->0), fcn. (208->8), ass. (0->38)
t26 = qJ(2) + qJ(3);
t20 = pkin(8) + t26;
t16 = sin(t20);
t17 = cos(t20);
t33 = t17 * pkin(4) + t16 * qJ(5);
t31 = -pkin(7) - pkin(6);
t21 = sin(t26);
t40 = pkin(3) * t21;
t39 = pkin(4) * t16;
t22 = cos(t26);
t18 = pkin(3) * t22;
t29 = cos(qJ(2));
t23 = t29 * pkin(2);
t38 = t18 + t23;
t37 = qJ(5) * t17;
t36 = t18 + t33;
t27 = sin(qJ(2));
t9 = -t27 * pkin(2) - t40;
t35 = t9 - t39;
t34 = -t39 - t40;
t28 = sin(qJ(1));
t30 = cos(qJ(1));
t13 = g(1) * t30 + g(2) * t28;
t12 = g(1) * t28 - g(2) * t30;
t3 = -g(3) * t22 + t13 * t21;
t32 = -g(3) * t29 + t13 * t27;
t25 = -qJ(4) + t31;
t19 = t23 + pkin(1);
t11 = t30 * t37;
t10 = t28 * t37;
t8 = pkin(1) + t38;
t7 = t30 * t8;
t6 = t12 * t17;
t5 = t12 * t16;
t4 = g(3) * t21 + t13 * t22;
t2 = g(3) * t16 + t13 * t17;
t1 = -g(3) * t17 + t13 * t16;
t14 = [0, 0, 0, 0, 0, 0, t12, t13, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t29, -t12 * t27, -t13, -g(1) * (-t28 * pkin(1) + t30 * pkin(6)) - g(2) * (t30 * pkin(1) + t28 * pkin(6)), 0, 0, 0, 0, 0, 0, t12 * t22, -t12 * t21, -t13, -g(1) * (-t28 * t19 - t30 * t31) - g(2) * (t30 * t19 - t28 * t31), 0, 0, 0, 0, 0, 0, t6, -t5, -t13, -g(1) * (-t30 * t25 - t28 * t8) - g(2) * (-t28 * t25 + t7), 0, 0, 0, 0, 0, 0, t6, -t13, t5, -g(2) * t7 + (g(1) * t25 - g(2) * t33) * t30 + (-g(1) * (-t33 - t8) + g(2) * t25) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, g(3) * t27 + t13 * t29, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t32 * pkin(2), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t38 - t13 * t9, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (t35 * t30 + t11) - g(2) * (t35 * t28 + t10) - g(3) * (t23 + t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(3), 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (t34 * t30 + t11) - g(2) * (t34 * t28 + t10) - g(3) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t14;
