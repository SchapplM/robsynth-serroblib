% Calculate inertial parameters regressor of gravitation load for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:13
% EndTime: 2019-12-31 21:51:14
% DurationCPUTime: 0.23s
% Computational Cost: add. (283->58), mult. (236->69), div. (0->0), fcn. (213->8), ass. (0->41)
t27 = qJ(3) + qJ(4);
t21 = sin(t27);
t23 = cos(t27);
t44 = t23 * pkin(4) + t21 * qJ(5);
t30 = sin(qJ(1));
t48 = t30 * pkin(1);
t28 = qJ(1) + qJ(2);
t22 = sin(t28);
t47 = t21 * t22;
t24 = cos(t28);
t46 = t21 * t24;
t33 = -pkin(8) - pkin(7);
t45 = t24 * t33;
t43 = t24 * pkin(2) + t22 * pkin(7);
t42 = qJ(5) * t23;
t31 = cos(qJ(3));
t25 = t31 * pkin(3);
t20 = t25 + pkin(2);
t12 = t24 * t20;
t41 = t44 * t24 + t12;
t40 = -t22 * pkin(2) + t24 * pkin(7);
t39 = -t22 * t33 + t12;
t29 = sin(qJ(3));
t38 = -pkin(3) * t29 - pkin(4) * t21;
t8 = g(1) * t24 + g(2) * t22;
t7 = g(1) * t22 - g(2) * t24;
t32 = cos(qJ(1));
t37 = g(1) * t30 - g(2) * t32;
t36 = -t22 * t20 - t45;
t35 = -g(3) * t31 + t8 * t29;
t34 = (-g(1) * (-t20 - t44) + g(2) * t33) * t22;
t26 = t32 * pkin(1);
t11 = t24 * t42;
t9 = t22 * t42;
t6 = t7 * t31;
t5 = t7 * t29;
t4 = t7 * t23;
t3 = g(1) * t47 - g(2) * t46;
t2 = g(3) * t21 + t8 * t23;
t1 = -g(3) * t23 + t8 * t21;
t10 = [0, 0, 0, 0, 0, 0, t37, g(1) * t32 + g(2) * t30, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t37 * pkin(1), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t40 - t48) - g(2) * (t26 + t43), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t36 - t48) - g(2) * (t26 + t39), 0, 0, 0, 0, 0, 0, t4, -t8, t3, -g(1) * (-t45 - t48) - g(2) * (t26 + t41) + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t40 - g(2) * t43, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * t36 - g(2) * t39, 0, 0, 0, 0, 0, 0, t4, -t8, t3, g(1) * t45 - g(2) * t41 + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, g(3) * t29 + t8 * t31, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t35 * pkin(3), 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (t38 * t24 + t11) - g(2) * (t38 * t22 + t9) - g(3) * (t25 + t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-pkin(4) * t46 + t11) - g(2) * (-pkin(4) * t47 + t9) - g(3) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t10;
