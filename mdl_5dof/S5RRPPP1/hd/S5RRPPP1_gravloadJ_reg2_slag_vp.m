% Calculate inertial parameters regressor of gravitation load for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:31
% EndTime: 2019-12-31 19:24:32
% DurationCPUTime: 0.33s
% Computational Cost: add. (229->84), mult. (644->122), div. (0->0), fcn. (733->8), ass. (0->67)
t48 = sin(qJ(2));
t45 = sin(pkin(5));
t68 = qJ(3) * t45;
t35 = t48 * t68;
t50 = cos(qJ(2));
t70 = t50 * pkin(2) + t35;
t47 = cos(pkin(5));
t44 = sin(pkin(8));
t79 = t48 * t44;
t64 = t47 * t79;
t46 = cos(pkin(8));
t73 = t50 * t46;
t23 = -t64 + t73;
t49 = sin(qJ(1));
t51 = cos(qJ(1));
t34 = g(1) * t51 + g(2) * t49;
t85 = g(3) * t48 + t34 * t50;
t84 = g(1) * t49;
t82 = t45 * t49;
t81 = t45 * t50;
t80 = t45 * t51;
t78 = t48 * t46;
t77 = t48 * t49;
t76 = t48 * t51;
t75 = t49 * t47;
t74 = t50 * t44;
t72 = t50 * t51;
t42 = t51 * pkin(7);
t67 = qJ(3) * t47;
t71 = t51 * t67 + t42;
t69 = t51 * pkin(1) + t49 * pkin(7);
t66 = pkin(2) * t77;
t65 = pkin(2) * t76;
t63 = t50 * t68;
t55 = -t47 * t73 + t79;
t17 = t55 * t49;
t54 = t47 * t74 + t78;
t18 = t54 * t49;
t28 = t49 * t63;
t62 = -t18 * pkin(3) - t17 * qJ(4) + t28;
t19 = t55 * t51;
t20 = t54 * t51;
t30 = t51 * t63;
t61 = -t20 * pkin(3) - t19 * qJ(4) + t30;
t60 = pkin(2) * t72 + t51 * t35 + t49 * t67 + t69;
t10 = t49 * t74 + (t48 * t75 + t80) * t46;
t22 = t47 * t78 + t74;
t12 = t22 * t51 - t46 * t82;
t2 = g(1) * t10 - g(2) * t12;
t11 = t23 * t49 - t44 * t80;
t13 = t44 * t82 + t46 * t72 - t51 * t64;
t3 = g(1) * t11 - g(2) * t13;
t59 = -g(2) * t51 + t84;
t58 = -t11 * pkin(3) - t10 * qJ(4) + t71;
t57 = t23 * pkin(3) + t22 * qJ(4) + t70;
t56 = -pkin(2) * t48 + pkin(4) * t81;
t4 = g(1) * t19 + g(2) * t17 - g(3) * t22;
t5 = g(1) * t20 + g(2) * t18 - g(3) * t23;
t53 = t13 * pkin(3) + t12 * qJ(4) + t60;
t52 = (-pkin(1) - t70) * t84;
t25 = t45 * t76 + t75;
t24 = t45 * t77 - t51 * t47;
t14 = t85 * t45;
t7 = g(1) * t24 - g(2) * t25;
t6 = -g(1) * t25 - g(2) * t24 + g(3) * t81;
t1 = -g(1) * t12 - g(2) * t10 - g(3) * t55;
t8 = [0, 0, 0, 0, 0, 0, t59, t34, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t50, -t59 * t48, -t34, -g(1) * (-t49 * pkin(1) + t42) - g(2) * t69, 0, 0, 0, 0, 0, 0, t3, -t2, t7, -g(1) * t71 - g(2) * t60 - t52, 0, 0, 0, 0, 0, 0, t7, -t3, t2, -g(1) * t58 - g(2) * t53 - t52, 0, 0, 0, 0, 0, 0, t7, t2, t3, -g(1) * (-t24 * pkin(4) - t11 * qJ(5) + t58) - g(2) * (t25 * pkin(4) + t13 * qJ(5) + t53) - t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t50 + t34 * t48, t85, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, -t14, -g(1) * (t30 - t65) - g(2) * (t28 - t66) - g(3) * t70, 0, 0, 0, 0, 0, 0, -t14, -t5, t4, -g(1) * (t61 - t65) - g(2) * (t62 - t66) - g(3) * t57, 0, 0, 0, 0, 0, 0, -t14, t4, t5, -g(1) * (-t20 * qJ(5) + t56 * t51 + t61) - g(2) * (-t18 * qJ(5) + t56 * t49 + t62) - g(3) * (t48 * t45 * pkin(4) + t23 * qJ(5) + t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t54;];
taug_reg = t8;
