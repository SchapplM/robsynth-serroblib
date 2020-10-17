% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR14_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:16
% EndTime: 2019-12-31 18:35:17
% DurationCPUTime: 0.22s
% Computational Cost: add. (126->51), mult. (193->65), div. (0->0), fcn. (187->8), ass. (0->33)
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t44 = -g(1) * t22 + g(2) * t25;
t18 = qJ(3) + pkin(8);
t12 = sin(t18);
t13 = cos(t18);
t43 = -g(3) * t12 - t13 * t44;
t39 = g(3) * t13;
t21 = sin(qJ(3));
t38 = t21 * pkin(3);
t20 = sin(qJ(5));
t37 = t22 * t20;
t23 = cos(qJ(5));
t36 = t22 * t23;
t35 = t25 * t20;
t34 = t25 * t23;
t33 = t25 * pkin(1) + t22 * qJ(2);
t15 = t25 * qJ(2);
t32 = -t22 * pkin(1) + t15;
t31 = t12 * pkin(4) - t13 * pkin(7);
t8 = g(1) * t25 + g(2) * t22;
t19 = -qJ(4) - pkin(6);
t30 = t22 * t19 + t25 * t38 + t32;
t29 = -t25 * t19 + t22 * t38 + t33;
t24 = cos(qJ(3));
t26 = g(3) * t21 + t24 * t44;
t6 = t12 * t34 - t37;
t5 = t12 * t35 + t36;
t4 = t12 * t36 + t35;
t3 = -t12 * t37 + t34;
t2 = t8 * t13;
t1 = -t12 * t44 + t39;
t7 = [0, 0, 0, 0, 0, 0, -t44, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t8, -g(1) * t32 - g(2) * t33, 0, 0, 0, 0, 0, 0, -t8 * t21, -t8 * t24, -t44, -g(1) * (t15 + (-pkin(1) - pkin(6)) * t22) - g(2) * (t25 * pkin(6) + t33), 0, 0, 0, 0, 0, 0, -t8 * t12, -t2, -t44, -g(1) * t30 - g(2) * t29, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t2, -g(1) * (t31 * t25 + t30) - g(2) * (t31 * t22 + t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, g(3) * t24 - t21 * t44, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t1, 0, t26 * pkin(3), 0, 0, 0, 0, 0, 0, -t43 * t23, t43 * t20, -t1, -g(3) * (-t31 - t38) + t44 * (pkin(3) * t24 + pkin(4) * t13 + pkin(7) * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t20 * t39, g(1) * t4 - g(2) * t6 + t23 * t39, 0, 0;];
taug_reg = t7;
