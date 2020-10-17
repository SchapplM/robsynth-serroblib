% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:18
% EndTime: 2020-01-03 12:15:19
% DurationCPUTime: 0.21s
% Computational Cost: add. (286->52), mult. (217->61), div. (0->0), fcn. (194->10), ass. (0->42)
t44 = -pkin(8) - pkin(7);
t38 = qJ(3) + qJ(4);
t31 = cos(t38);
t42 = cos(qJ(3));
t35 = t42 * pkin(3);
t51 = pkin(4) * t31 + t35;
t15 = pkin(2) + t51;
t39 = qJ(1) + qJ(2);
t30 = sin(t39);
t32 = cos(t39);
t37 = -pkin(9) + t44;
t54 = t30 * t15 + t32 * t37;
t28 = t35 + pkin(2);
t53 = t30 * t28 + t32 * t44;
t52 = t32 * pkin(2) + t30 * pkin(7);
t50 = t30 * pkin(2) - t32 * pkin(7);
t49 = t32 * t15 - t30 * t37;
t48 = t32 * t28 - t30 * t44;
t47 = g(2) * t32 + g(3) * t30;
t13 = g(2) * t30 - g(3) * t32;
t41 = sin(qJ(1));
t43 = cos(qJ(1));
t46 = -g(2) * t43 - g(3) * t41;
t29 = sin(t38);
t4 = -g(1) * t31 + t13 * t29;
t40 = sin(qJ(3));
t45 = -g(1) * t42 + t13 * t40;
t36 = t43 * pkin(1);
t34 = t41 * pkin(1);
t33 = qJ(5) + t38;
t27 = cos(t33);
t26 = sin(t33);
t10 = t47 * t42;
t9 = t47 * t40;
t8 = t47 * t31;
t7 = t47 * t29;
t6 = t47 * t27;
t5 = t47 * t26;
t3 = g(1) * t29 + t13 * t31;
t2 = -g(1) * t27 + t13 * t26;
t1 = g(1) * t26 + t13 * t27;
t11 = [0, 0, 0, 0, 0, 0, t46, g(2) * t41 - g(3) * t43, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t13, 0, t46 * pkin(1), 0, 0, 0, 0, 0, 0, -t10, t9, -t13, -g(2) * (t36 + t52) - g(3) * (t34 + t50), 0, 0, 0, 0, 0, 0, -t8, t7, -t13, -g(2) * (t36 + t48) - g(3) * (t34 + t53), 0, 0, 0, 0, 0, 0, -t6, t5, -t13, -g(2) * (t36 + t49) - g(3) * (t34 + t54); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t13, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, -t13, -g(2) * t52 - g(3) * t50, 0, 0, 0, 0, 0, 0, -t8, t7, -t13, -g(2) * t48 - g(3) * t53, 0, 0, 0, 0, 0, 0, -t6, t5, -t13, -g(2) * t49 - g(3) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, g(1) * t40 + t13 * t42, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, t45 * pkin(3), 0, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t51 - t13 * (-t40 * pkin(3) - pkin(4) * t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t11;
