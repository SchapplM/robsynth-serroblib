% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:03:21
% EndTime: 2019-05-06 16:03:23
% DurationCPUTime: 0.56s
% Computational Cost: add. (324->120), mult. (498->156), div. (0->0), fcn. (498->10), ass. (0->77)
t46 = sin(qJ(1));
t49 = cos(qJ(1));
t22 = g(1) * t49 + g(2) * t46;
t45 = sin(qJ(2));
t93 = t22 * t45;
t34 = t45 * qJ(3);
t48 = cos(qJ(2));
t59 = t48 * pkin(2) + t34;
t47 = cos(qJ(4));
t60 = t49 * t47;
t44 = sin(qJ(4));
t71 = t46 * t44;
t13 = t45 * t60 - t71;
t61 = t49 * t44;
t70 = t46 * t47;
t15 = t45 * t70 + t61;
t84 = g(3) * t48;
t92 = -g(1) * t13 - g(2) * t15 + t47 * t84;
t12 = g(3) * t45 + t22 * t48;
t89 = g(1) * t46;
t83 = t44 * pkin(4);
t82 = t48 * pkin(8);
t79 = t44 * t48;
t78 = t45 * t46;
t77 = t45 * t49;
t42 = qJ(4) + pkin(10);
t31 = sin(t42);
t20 = pkin(5) * t31 + t83;
t76 = t46 * t20;
t33 = qJ(6) + t42;
t28 = sin(t33);
t75 = t46 * t28;
t29 = cos(t33);
t74 = t46 * t29;
t73 = t46 * t31;
t32 = cos(t42);
t72 = t46 * t32;
t43 = -qJ(5) - pkin(8);
t41 = -pkin(9) + t43;
t69 = t48 * t41;
t68 = t48 * t43;
t67 = t48 * t49;
t66 = t49 * t20;
t65 = t49 * t28;
t64 = t49 * t29;
t63 = t49 * t31;
t62 = t49 * t32;
t36 = t47 * pkin(4);
t21 = pkin(5) * t32 + t36;
t58 = t49 * pkin(1) + t46 * pkin(7);
t57 = qJ(3) * t48;
t56 = pkin(4) * t79;
t54 = t45 * t61;
t53 = pkin(2) * t67 + t49 * t34 + t58;
t52 = -g(2) * t49 + t89;
t51 = -pkin(1) - t59;
t38 = t49 * pkin(7);
t30 = t36 + pkin(3);
t25 = t49 * t57;
t23 = t46 * t57;
t19 = pkin(3) + t21;
t18 = t52 * t48;
t17 = t52 * t45;
t16 = -t45 * t71 + t60;
t14 = t54 + t70;
t11 = -t84 + t93;
t10 = -t45 * t73 + t62;
t9 = t45 * t72 + t63;
t8 = t45 * t63 + t72;
t7 = t45 * t62 - t73;
t6 = -t45 * t75 + t64;
t5 = t45 * t74 + t65;
t4 = t45 * t65 + t74;
t3 = t45 * t64 - t75;
t2 = g(1) * t4 - g(2) * t6 - t28 * t84;
t1 = -g(1) * t3 - g(2) * t5 + t29 * t84;
t24 = [0, 0, 0, 0, 0, 0, t52, t22, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, -t22, -g(1) * (-t46 * pkin(1) + t38) - g(2) * t58, 0, 0, 0, 0, 0, 0, -t22, -t18, t17, -g(1) * t38 - g(2) * t53 - t51 * t89, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14, g(1) * t15 - g(2) * t13, t18, -g(1) * (t49 * pkin(3) + t38) - g(2) * (pkin(8) * t67 + t53) + (-g(1) * (t51 - t82) - g(2) * pkin(3)) * t46, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, t18, -g(1) * (t49 * t30 + t38) - g(2) * (pkin(4) * t54 - t43 * t67 + t53) + (-g(1) * (-t45 * t83 + t51 + t68) - g(2) * t30) * t46, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t18, -g(1) * (t49 * t19 + t38) - g(2) * (-t41 * t67 + t45 * t66 + t53) + (-g(1) * (-t45 * t20 + t51 + t69) - g(2) * t19) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -g(1) * (-pkin(2) * t77 + t25) - g(2) * (-pkin(2) * t78 + t23) - g(3) * t59, 0, 0, 0, 0, 0, 0, -t12 * t44, -t12 * t47, t11, -g(1) * t25 - g(2) * t23 - g(3) * (t59 + t82) + (pkin(2) + pkin(8)) * t93, 0, 0, 0, 0, 0, 0, -t12 * t31, -t12 * t32, t11, -g(1) * (t49 * t56 + t25) - g(2) * (t46 * t56 + t23) - g(3) * (t59 - t68) + (-g(3) * t83 + t22 * (pkin(2) - t43)) * t45, 0, 0, 0, 0, 0, 0, -t12 * t28, -t12 * t29, t11, -g(1) * (t48 * t66 + t25) - g(2) * (t48 * t76 + t23) - g(3) * (t59 - t69) + (-g(3) * t20 + t22 * (pkin(2) - t41)) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, g(1) * t14 - g(2) * t16 - g(3) * t79, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t32 * t84, g(1) * t8 - g(2) * t10 - t31 * t84, 0, t92 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t21 * t77 - t76) - g(2) * (t21 * t78 + t66) + t21 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t24;
