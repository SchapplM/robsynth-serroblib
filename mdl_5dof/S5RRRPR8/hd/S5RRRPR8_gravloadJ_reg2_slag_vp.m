% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:26
% EndTime: 2019-12-31 21:20:27
% DurationCPUTime: 0.32s
% Computational Cost: add. (240->74), mult. (308->89), div. (0->0), fcn. (293->8), ass. (0->47)
t25 = qJ(2) + qJ(3);
t22 = sin(t25);
t23 = cos(t25);
t42 = t23 * pkin(3) + t22 * qJ(4);
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t11 = g(1) * t31 + g(2) * t28;
t56 = t11 * t22;
t4 = g(3) * t22 + t11 * t23;
t27 = sin(qJ(2));
t54 = pkin(2) * t27;
t53 = pkin(3) * t22;
t49 = g(3) * t23;
t18 = t23 * pkin(8);
t32 = -pkin(7) - pkin(6);
t48 = pkin(4) - t32;
t26 = sin(qJ(5));
t47 = t28 * t26;
t29 = cos(qJ(5));
t46 = t28 * t29;
t45 = t31 * t26;
t44 = t31 * t29;
t43 = t31 * t32;
t41 = qJ(4) * t23;
t30 = cos(qJ(2));
t24 = t30 * pkin(2);
t40 = t24 + t42;
t21 = t24 + pkin(1);
t15 = t31 * t21;
t39 = g(2) * (t42 * t31 + t15);
t38 = -t53 - t54;
t37 = g(1) * t28 - g(2) * t31;
t36 = -t21 - t42;
t34 = -g(3) * t30 + t11 * t27;
t33 = (pkin(3) + pkin(8)) * t56;
t14 = t31 * t41;
t12 = t28 * t41;
t10 = -t22 * t47 + t44;
t9 = t22 * t46 + t45;
t8 = t22 * t45 + t46;
t7 = t22 * t44 - t47;
t6 = t37 * t23;
t5 = t37 * t22;
t3 = -t49 + t56;
t2 = t4 * t29;
t1 = t4 * t26;
t13 = [0, 0, 0, 0, 0, 0, t37, t11, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t30, -t37 * t27, -t11, -g(1) * (-t28 * pkin(1) + t31 * pkin(6)) - g(2) * (t31 * pkin(1) + t28 * pkin(6)), 0, 0, 0, 0, 0, 0, t6, -t5, -t11, -g(1) * (-t28 * t21 - t43) - g(2) * (-t28 * t32 + t15), 0, 0, 0, 0, 0, 0, -t11, -t6, t5, g(1) * t43 - t39 + (-g(1) * t36 + g(2) * t32) * t28, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, t6, -t39 + (-g(1) * t48 - g(2) * t18) * t31 + (-g(1) * (t36 - t18) - g(2) * t48) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, g(3) * t27 + t11 * t30, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t34 * pkin(2), 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (t38 * t31 + t14) - g(2) * (t38 * t28 + t12) - g(3) * t40, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * (-t31 * t54 + t14) - g(2) * (-t28 * t54 + t12) - g(3) * (t18 + t40) + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (-t31 * t53 + t14) - g(2) * (-t28 * t53 + t12) - g(3) * t42, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * t14 - g(2) * t12 - g(3) * (t18 + t42) + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t29 * t49, g(1) * t8 - g(2) * t10 - t26 * t49, 0, 0;];
taug_reg = t13;
