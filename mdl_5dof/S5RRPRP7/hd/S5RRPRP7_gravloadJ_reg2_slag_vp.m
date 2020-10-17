% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:31
% EndTime: 2019-12-31 20:01:32
% DurationCPUTime: 0.28s
% Computational Cost: add. (229->70), mult. (327->93), div. (0->0), fcn. (333->8), ass. (0->48)
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t12 = g(1) * t34 + g(2) * t31;
t27 = qJ(2) + pkin(8);
t23 = sin(t27);
t56 = t12 * t23;
t18 = t23 * pkin(7);
t24 = cos(t27);
t55 = -t24 * pkin(3) - t18;
t30 = sin(qJ(2));
t54 = pkin(2) * t30;
t51 = g(3) * t23;
t50 = t24 * t34;
t29 = sin(qJ(4));
t49 = t31 * t29;
t32 = cos(qJ(4));
t48 = t31 * t32;
t28 = -qJ(3) - pkin(6);
t47 = t34 * t28;
t46 = t34 * t29;
t45 = t34 * t32;
t33 = cos(qJ(2));
t25 = t33 * pkin(2);
t22 = t25 + pkin(1);
t17 = t34 * t22;
t44 = pkin(3) * t50 + t34 * t18 + t17;
t43 = t25 - t55;
t7 = t24 * t49 + t45;
t9 = t24 * t46 - t48;
t42 = g(1) * t7 - g(2) * t9;
t41 = -pkin(3) * t23 - t54;
t11 = g(1) * t31 - g(2) * t34;
t40 = pkin(4) * t32 + qJ(5) * t29;
t1 = g(1) * t9 + g(2) * t7 + t29 * t51;
t10 = t24 * t45 + t49;
t8 = t24 * t48 - t46;
t38 = g(1) * t10 + g(2) * t8 + t32 * t51;
t37 = -g(3) * t24 + t56;
t36 = -g(3) * t33 + t12 * t30;
t35 = (-g(1) * (-t22 + t55) + g(2) * t28) * t31;
t15 = pkin(7) * t50;
t13 = t31 * t24 * pkin(7);
t6 = t11 * t23;
t5 = t12 * t24 + t51;
t4 = t37 * t32;
t3 = t37 * t29;
t2 = g(1) * t8 - g(2) * t10;
t14 = [0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t33, -t11 * t30, -t12, -g(1) * (-t31 * pkin(1) + t34 * pkin(6)) - g(2) * (t34 * pkin(1) + t31 * pkin(6)), 0, 0, 0, 0, 0, 0, t11 * t24, -t6, -t12, -g(1) * (-t31 * t22 - t47) - g(2) * (-t31 * t28 + t17), 0, 0, 0, 0, 0, 0, t2, -t42, t6, g(1) * t47 - g(2) * t44 + t35, 0, 0, 0, 0, 0, 0, t2, t6, t42, -g(1) * (-t8 * pkin(4) - t7 * qJ(5) - t47) - g(2) * (t10 * pkin(4) + t9 * qJ(5) + t44) + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, g(3) * t30 + t12 * t33, 0, 0, 0, 0, 0, 0, 0, 0, t37, t5, 0, t36 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, -t5, -g(1) * (t41 * t34 + t15) - g(2) * (t41 * t31 + t13) - g(3) * t43, 0, 0, 0, 0, 0, 0, t4, -t5, t3, -g(1) * (-t34 * t54 + t15) - g(2) * (-t31 * t54 + t13) - g(3) * (t24 * t40 + t43) + (pkin(3) + t40) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t38, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t38, -g(1) * (-t9 * pkin(4) + t10 * qJ(5)) - g(2) * (-t7 * pkin(4) + t8 * qJ(5)) - (-pkin(4) * t29 + qJ(5) * t32) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t14;
