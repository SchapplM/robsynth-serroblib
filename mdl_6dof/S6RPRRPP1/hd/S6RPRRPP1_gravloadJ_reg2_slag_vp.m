% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:17:25
% EndTime: 2019-05-05 21:17:26
% DurationCPUTime: 0.41s
% Computational Cost: add. (434->89), mult. (430->119), div. (0->0), fcn. (439->10), ass. (0->55)
t41 = cos(qJ(4));
t28 = t41 * pkin(4) + pkin(3);
t37 = -qJ(5) - pkin(8);
t39 = sin(qJ(3));
t42 = cos(qJ(3));
t55 = t42 * t28 - t39 * t37;
t36 = qJ(1) + pkin(9);
t30 = sin(t36);
t32 = cos(t36);
t17 = g(1) * t32 + g(2) * t30;
t9 = -g(3) * t42 + t17 * t39;
t73 = g(1) * t30;
t70 = g(3) * t39;
t38 = sin(qJ(4));
t68 = t30 * t38;
t67 = t30 * t41;
t66 = t30 * t42;
t65 = t32 * t38;
t64 = t32 * t41;
t63 = t32 * t42;
t62 = t37 * t42;
t61 = t38 * t42;
t59 = t41 * t42;
t58 = t32 * t61;
t43 = cos(qJ(1));
t57 = t43 * pkin(1) + t32 * pkin(2) + t30 * pkin(7);
t40 = sin(qJ(1));
t56 = -t40 * pkin(1) + t32 * pkin(7);
t35 = qJ(4) + pkin(10);
t29 = sin(t35);
t31 = cos(t35);
t5 = t29 * t66 + t32 * t31;
t7 = t29 * t63 - t30 * t31;
t54 = g(1) * t5 - g(2) * t7;
t53 = t42 * pkin(3) + t39 * pkin(8);
t51 = -g(2) * t32 + t73;
t50 = g(1) * t40 - g(2) * t43;
t49 = pkin(5) * t31 + qJ(6) * t29;
t11 = t30 * t61 + t64;
t1 = g(1) * t7 + g(2) * t5 + t29 * t70;
t6 = -t32 * t29 + t31 * t66;
t8 = t30 * t29 + t31 * t63;
t47 = g(1) * t8 + g(2) * t6 + t31 * t70;
t46 = pkin(4) * t68 + t55 * t32 + t57;
t10 = t17 * t42 + t70;
t44 = pkin(4) * t65 + t56 + (-pkin(2) - t55) * t30;
t22 = pkin(4) * t67;
t15 = t51 * t39;
t14 = t32 * t59 + t68;
t13 = -t58 + t67;
t12 = -t30 * t59 + t65;
t4 = t9 * t31;
t3 = t9 * t29;
t2 = g(1) * t6 - g(2) * t8;
t16 = [0, 0, 0, 0, 0, 0, t50, g(1) * t43 + g(2) * t40, 0, 0, 0, 0, 0, 0, 0, 0, t51, t17, 0, t50 * pkin(1), 0, 0, 0, 0, 0, 0, t51 * t42, -t15, -t17, -g(1) * (-t30 * pkin(2) + t56) - g(2) * t57, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t15, -g(1) * t56 - g(2) * (t53 * t32 + t57) - (-pkin(2) - t53) * t73, 0, 0, 0, 0, 0, 0, t2, -t54, t15, -g(1) * t44 - g(2) * t46, 0, 0, 0, 0, 0, 0, t2, t15, t54, -g(1) * (-t6 * pkin(5) - t5 * qJ(6) + t44) - g(2) * (t8 * pkin(5) + t7 * qJ(6) + t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t41, -t9 * t38, -t10, -g(3) * t53 + t17 * (pkin(3) * t39 - pkin(8) * t42) 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(3) * t55 + t17 * (t28 * t39 + t62) 0, 0, 0, 0, 0, 0, t4, -t10, t3, -g(3) * (t49 * t42 + t55) + t17 * (t62 - (-t28 - t49) * t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 + g(2) * t11 + t38 * t70, g(1) * t14 - g(2) * t12 + t41 * t70, 0, 0, 0, 0, 0, 0, 0, 0, t1, t47, 0, -g(1) * t22 + (g(2) * t64 + t10 * t38) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t47, -g(1) * (-pkin(4) * t58 - t7 * pkin(5) + t8 * qJ(6) + t22) - g(2) * (-t11 * pkin(4) - t5 * pkin(5) + t6 * qJ(6)) - (-pkin(4) * t38 - pkin(5) * t29 + qJ(6) * t31) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t16;
