% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR10
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:05
% EndTime: 2019-12-31 18:26:06
% DurationCPUTime: 0.21s
% Computational Cost: add. (170->51), mult. (212->55), div. (0->0), fcn. (240->8), ass. (0->37)
t28 = cos(qJ(3));
t19 = t28 * pkin(3) + pkin(2);
t29 = cos(qJ(1));
t26 = sin(qJ(1));
t39 = t29 * pkin(1) + t26 * qJ(2);
t50 = t29 * t19 + t39;
t24 = sin(qJ(5));
t38 = qJ(3) + pkin(8);
t35 = sin(t38);
t36 = cos(t38);
t3 = -t26 * t35 - t29 * t36;
t4 = -t26 * t36 + t29 * t35;
t34 = g(1) * t4 - g(2) * t3;
t49 = t34 * t24;
t27 = cos(qJ(5));
t48 = t34 * t27;
t47 = t4 * pkin(4) + t3 * pkin(7);
t25 = sin(qJ(3));
t44 = t26 * t25;
t43 = t26 * t28;
t42 = t29 * t25;
t41 = t29 * t28;
t15 = pkin(3) * t42;
t40 = pkin(3) * t43 - t15;
t37 = pkin(3) * t41;
t2 = g(1) * t3 + g(2) * t4;
t5 = -t41 - t44;
t6 = t42 - t43;
t33 = g(1) * t6 - g(2) * t5;
t32 = g(1) * t5 + g(2) * t6;
t13 = pkin(3) * t44;
t31 = t3 * pkin(4) - t4 * pkin(7) - t13;
t21 = t29 * qJ(2);
t30 = t15 + t21 + (-pkin(1) - t19) * t26;
t8 = g(1) * t29 + g(2) * t26;
t7 = g(1) * t26 - g(2) * t29;
t1 = [0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t8, -g(1) * (-t26 * pkin(1) + t21) - g(2) * t39, 0, 0, 0, 0, 0, 0, -t33, t32, 0, -g(1) * (t21 + (-pkin(1) - pkin(2)) * t26) - g(2) * (t29 * pkin(2) + t39), 0, 0, 0, 0, 0, 0, -t34, t2, 0, -g(1) * t30 - g(2) * (t13 + t50), 0, 0, 0, 0, 0, 0, -t48, t49, -t2, -g(1) * (t30 + t47) - g(2) * (-t31 + t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t2, 0, -g(1) * t40 - g(2) * (-t13 - t37), 0, 0, 0, 0, 0, 0, t48, -t49, t2, -g(1) * (t40 - t47) - g(2) * (t31 - t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t27 - t2 * t24, -g(3) * t24 - t2 * t27, 0, 0;];
taug_reg = t1;
