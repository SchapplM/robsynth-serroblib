% Calculate inertial parameters regressor of gravitation load for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:08
% EndTime: 2019-12-05 16:10:10
% DurationCPUTime: 0.24s
% Computational Cost: add. (163->58), mult. (264->78), div. (0->0), fcn. (267->8), ass. (0->36)
t19 = sin(pkin(7));
t20 = cos(pkin(7));
t31 = g(1) * t20 + g(2) * t19;
t23 = sin(qJ(2));
t25 = cos(qJ(2));
t8 = -g(3) * t25 + t31 * t23;
t43 = g(3) * t23;
t24 = cos(qJ(3));
t41 = t19 * t24;
t40 = t19 * t25;
t39 = t20 * t24;
t38 = t20 * t25;
t21 = -qJ(4) - pkin(6);
t37 = t21 * t25;
t22 = sin(qJ(3));
t36 = t22 * t25;
t35 = t24 * t25;
t34 = t20 * t36;
t14 = t24 * pkin(3) + pkin(2);
t33 = t25 * t14 - t23 * t21;
t18 = qJ(3) + pkin(8);
t15 = sin(t18);
t16 = cos(t18);
t30 = pkin(4) * t16 + qJ(5) * t15;
t28 = -t19 * t36 - t39;
t4 = t15 * t40 + t20 * t16;
t6 = t15 * t38 - t19 * t16;
t1 = g(1) * t6 + g(2) * t4 + t15 * t43;
t5 = -t20 * t15 + t16 * t40;
t7 = t19 * t15 + t16 * t38;
t27 = g(1) * t7 + g(2) * t5 + t16 * t43;
t9 = t31 * t25 + t43;
t13 = pkin(3) * t41;
t3 = t8 * t16;
t2 = t8 * t15;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t24, -t8 * t22, -t9, -g(3) * (t25 * pkin(2) + t23 * pkin(6)) + t31 * (pkin(2) * t23 - pkin(6) * t25), 0, 0, 0, 0, 0, 0, t3, -t2, -t9, -g(3) * t33 + t31 * (t14 * t23 + t37), 0, 0, 0, 0, 0, 0, t3, -t9, t2, -g(3) * (t30 * t25 + t33) + t31 * (t37 - (-t14 - t30) * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t34 + t41) - g(2) * t28 + t22 * t43, -g(1) * (-t19 * t22 - t20 * t35) - g(2) * (-t19 * t35 + t20 * t22) + t24 * t43, 0, 0, 0, 0, 0, 0, 0, 0, t1, t27, 0, -g(1) * t13 + (g(2) * t39 + t22 * t9) * pkin(3), 0, 0, 0, 0, 0, 0, t1, 0, -t27, -g(1) * (-pkin(3) * t34 - t6 * pkin(4) + t7 * qJ(5) + t13) - g(2) * (pkin(3) * t28 - t4 * pkin(4) + t5 * qJ(5)) - (-pkin(3) * t22 - pkin(4) * t15 + qJ(5) * t16) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t10;
