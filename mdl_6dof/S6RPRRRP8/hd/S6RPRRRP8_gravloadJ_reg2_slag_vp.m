% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:48:05
% EndTime: 2019-05-06 01:48:06
% DurationCPUTime: 0.31s
% Computational Cost: add. (318->92), mult. (449->111), div. (0->0), fcn. (447->8), ass. (0->58)
t40 = cos(qJ(3));
t73 = pkin(3) * t40;
t35 = qJ(3) + qJ(4);
t29 = sin(t35);
t72 = pkin(9) * t29;
t41 = cos(qJ(1));
t71 = g(2) * t41;
t70 = g(3) * t29;
t30 = cos(t35);
t69 = g(3) * t30;
t36 = sin(qJ(5));
t68 = g(3) * t36;
t26 = t30 * pkin(9);
t37 = sin(qJ(3));
t67 = t37 * pkin(3);
t38 = sin(qJ(1));
t66 = t29 * t38;
t65 = t29 * t41;
t64 = t30 * t38;
t63 = t38 * t36;
t39 = cos(qJ(5));
t62 = t38 * t39;
t61 = t41 * t36;
t60 = t41 * t39;
t59 = pkin(4) * t64 + pkin(9) * t66;
t58 = pkin(1) * t41 + qJ(2) * t38;
t57 = t30 * t63;
t32 = t41 * qJ(2);
t56 = -pkin(1) * t38 + t32;
t55 = -pkin(4) * t29 + t26;
t54 = pkin(5) * t30 * t62 + qJ(6) * t57 + t59;
t10 = t29 * t61 + t62;
t8 = t29 * t63 - t60;
t53 = g(1) * t10 + g(2) * t8;
t52 = -pkin(4) * t30 - t72;
t16 = g(1) * t41 + g(2) * t38;
t15 = g(1) * t38 - t71;
t42 = -pkin(8) - pkin(7);
t51 = t38 * t42 + t41 * t67 + t56;
t50 = t38 * t67 - t41 * t42 + t58;
t49 = -pkin(5) * t39 - qJ(6) * t36 - pkin(4);
t1 = g(1) * t8 - g(2) * t10 + t30 * t68;
t11 = t29 * t60 - t63;
t9 = t29 * t62 + t61;
t48 = g(1) * t9 - g(2) * t11 + t39 * t69;
t47 = t49 * t70;
t6 = -t15 * t30 + t70;
t46 = g(3) * t37 - t15 * t40;
t45 = pkin(4) * t65 - t26 * t41 + t51;
t44 = pkin(4) * t66 - pkin(9) * t64 + t50;
t43 = t30 * t49 - t72;
t24 = t38 * t73;
t7 = t16 * t30;
t5 = g(1) * t66 - g(2) * t65 + t69;
t4 = t6 * t39;
t3 = -g(2) * t30 * t61 + g(1) * t57 - t29 * t68;
t2 = -g(1) * t11 - g(2) * t9;
t12 = [0, 0, 0, 0, 0, 0, t15, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, -g(1) * t56 - g(2) * t58, 0, 0, 0, 0, 0, 0, -t16 * t37, -t16 * t40, t15, -g(1) * (t32 + (-pkin(1) - pkin(7)) * t38) - g(2) * (pkin(7) * t41 + t58) 0, 0, 0, 0, 0, 0, -t16 * t29, -t7, t15, -g(1) * t51 - g(2) * t50, 0, 0, 0, 0, 0, 0, t2, t53, t7, -g(1) * t45 - g(2) * t44, 0, 0, 0, 0, 0, 0, t2, t7, -t53, -g(1) * (pkin(5) * t11 + qJ(6) * t10 + t45) - g(2) * (pkin(5) * t9 + qJ(6) * t8 + t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, g(3) * t40 + t15 * t37, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, t46 * pkin(3), 0, 0, 0, 0, 0, 0, t4, t3, -t5, -g(1) * (t24 + t59) - g(3) * (t55 - t67) - (t52 - t73) * t71, 0, 0, 0, 0, 0, 0, t4, -t5, -t3, -g(1) * (t24 + t54) - g(3) * (t26 - t67) - t47 - (t43 - t73) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, -t5, -g(1) * t59 - g(3) * t55 - t52 * t71, 0, 0, 0, 0, 0, 0, t4, -t5, -t3, -g(1) * t54 - g(3) * t26 - t43 * t71 - t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t48, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t48, -g(1) * (-pkin(5) * t8 + qJ(6) * t9) - g(2) * (pkin(5) * t10 - qJ(6) * t11) - (-pkin(5) * t36 + qJ(6) * t39) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t12;
