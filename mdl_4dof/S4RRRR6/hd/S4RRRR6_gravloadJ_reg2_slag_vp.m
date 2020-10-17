% Calculate inertial parameters regressor of gravitation load for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:51
% EndTime: 2019-12-31 17:30:52
% DurationCPUTime: 0.37s
% Computational Cost: add. (226->84), mult. (596->144), div. (0->0), fcn. (724->10), ass. (0->52)
t33 = sin(qJ(2));
t34 = sin(qJ(1));
t37 = cos(qJ(2));
t55 = cos(pkin(4));
t67 = cos(qJ(1));
t45 = t55 * t67;
t16 = t34 * t33 - t37 * t45;
t31 = sin(qJ(4));
t35 = cos(qJ(4));
t17 = t33 * t45 + t34 * t37;
t32 = sin(qJ(3));
t36 = cos(qJ(3));
t30 = sin(pkin(4));
t54 = t30 * t67;
t5 = t17 * t36 - t32 * t54;
t70 = -t16 * t35 + t5 * t31;
t69 = t16 * t31 + t5 * t35;
t68 = g(3) * t30;
t64 = t30 * t33;
t63 = t30 * t34;
t62 = t30 * t36;
t61 = t30 * t37;
t60 = t31 * t36;
t59 = t35 * t36;
t58 = t36 * t37;
t57 = pkin(2) * t61 + pkin(7) * t64;
t56 = t67 * pkin(1) + pkin(6) * t63;
t53 = -t34 * pkin(1) + pkin(6) * t54;
t52 = -t16 * pkin(2) + t17 * pkin(7);
t49 = t34 * t55;
t18 = t67 * t33 + t37 * t49;
t19 = -t33 * t49 + t67 * t37;
t51 = -t18 * pkin(2) + t19 * pkin(7);
t50 = -t17 * t32 - t36 * t54;
t8 = t19 * t32 - t34 * t62;
t48 = g(1) * t50 + g(2) * t8;
t47 = pkin(3) * t36 + pkin(8) * t32;
t46 = g(1) * t16 - g(2) * t18;
t44 = t19 * pkin(2) + t18 * pkin(7) + t56;
t43 = g(1) * t67 + g(2) * t34;
t42 = -t17 * pkin(2) - t16 * pkin(7) + t53;
t14 = -t32 * t64 + t55 * t36;
t41 = g(1) * t8 - g(2) * t50 - g(3) * t14;
t15 = t55 * t32 + t33 * t62;
t9 = t19 * t36 + t32 * t63;
t40 = g(1) * t9 + g(2) * t5 + g(3) * t15;
t39 = -g(1) * t18 - g(2) * t16 + g(3) * t61;
t38 = g(1) * t19 + g(2) * t17 + g(3) * t64;
t3 = t39 * t32;
t2 = t18 * t31 + t9 * t35;
t1 = t18 * t35 - t9 * t31;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t34 - g(2) * t67, t43, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t19, -t46, -t43 * t30, -g(1) * t53 - g(2) * t56, 0, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t48, t46, -g(1) * t42 - g(2) * t44, 0, 0, 0, 0, 0, 0, g(1) * t69 - g(2) * t2, -g(1) * t70 - g(2) * t1, -t48, -g(1) * (-pkin(3) * t5 + pkin(8) * t50 + t42) - g(2) * (t9 * pkin(3) + t8 * pkin(8) + t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, 0, 0, 0, 0, 0, 0, 0, 0, -t39 * t36, t3, -t38, -g(1) * t51 - g(2) * t52 - g(3) * t57, 0, 0, 0, 0, 0, 0, -g(1) * (-t18 * t59 + t19 * t31) - g(2) * (-t16 * t59 + t17 * t31) - (t31 * t33 + t35 * t58) * t68, -g(1) * (t18 * t60 + t19 * t35) - g(2) * (t16 * t60 + t17 * t35) - (-t31 * t58 + t33 * t35) * t68, -t3, -g(1) * (-t47 * t18 + t51) - g(2) * (-t47 * t16 + t52) - g(3) * (t47 * t61 + t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t35, -t41 * t31, -t40, -g(1) * (-t8 * pkin(3) + t9 * pkin(8)) - g(2) * (pkin(3) * t50 + t5 * pkin(8)) - g(3) * (t14 * pkin(3) + t15 * pkin(8)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t70 - g(3) * (-t15 * t31 - t35 * t61), g(1) * t2 + g(2) * t69 - g(3) * (-t15 * t35 + t31 * t61), 0, 0;];
taug_reg = t4;
