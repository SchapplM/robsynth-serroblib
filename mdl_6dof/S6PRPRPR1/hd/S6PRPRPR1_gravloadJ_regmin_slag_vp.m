% Calculate minimal parameter regressor of gravitation load for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% taug_reg [6x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:06:19
% EndTime: 2021-01-16 01:06:22
% DurationCPUTime: 0.68s
% Computational Cost: add. (429->116), mult. (462->192), div. (0->0), fcn. (512->18), ass. (0->90)
t37 = qJ(2) + pkin(11);
t35 = cos(t37);
t40 = cos(pkin(10));
t26 = t40 * t35;
t33 = sin(t37);
t38 = sin(pkin(10));
t41 = cos(pkin(6));
t96 = t38 * t41;
t105 = t33 * t96 - t26;
t25 = t38 * t35;
t87 = t40 * t41;
t106 = t33 * t87 + t25;
t36 = qJ(4) + pkin(12);
t34 = cos(t36);
t27 = t41 * t34;
t32 = sin(t36);
t39 = sin(pkin(6));
t107 = (g(1) * t105 - g(2) * t106) * t32 + g(3) * (-t39 * t33 * t32 + t27);
t31 = pkin(6) - t37;
t102 = sin(t31);
t101 = g(3) * t39;
t100 = t34 * t39;
t99 = t35 * t41;
t98 = t38 * t33;
t97 = t38 * t39;
t45 = sin(qJ(2));
t95 = t38 * t45;
t94 = t39 * t40;
t43 = sin(qJ(6));
t93 = t39 * t43;
t44 = sin(qJ(4));
t92 = t39 * t44;
t46 = cos(qJ(6));
t91 = t39 * t46;
t47 = cos(qJ(4));
t90 = t39 * t47;
t48 = cos(qJ(2));
t89 = t39 * t48;
t88 = t40 * t33;
t86 = t41 * t32;
t85 = t41 * t45;
t84 = t41 * t47;
t83 = t41 * t48;
t82 = t43 * t38;
t81 = t43 * t40;
t80 = t46 * t38;
t79 = t46 * t40;
t78 = g(3) * t93;
t77 = g(3) * t91;
t76 = t34 * t82;
t75 = t34 * t80;
t74 = t34 * t81;
t73 = t34 * t79;
t72 = t38 * t90;
t70 = t41 * t82;
t69 = t41 * t80;
t68 = t40 * t90;
t66 = t41 * t81;
t65 = t41 * t79;
t64 = t40 * t83;
t63 = cos(t31) / 0.2e1;
t60 = pkin(6) + t37;
t53 = sin(t60) / 0.2e1;
t21 = t53 - t102 / 0.2e1;
t62 = t40 * t21 + t25;
t61 = -t38 * t21 + t26;
t59 = g(1) * t38 - g(2) * t40;
t58 = cos(t60);
t55 = -t38 * t83 - t40 * t45;
t54 = t59 * t39;
t50 = t58 / 0.2e1 + t63;
t1 = -t40 * t50 + t98;
t22 = t102 / 0.2e1 + t53;
t4 = t38 * t50 + t88;
t52 = g(1) * t4 + g(2) * t1 - g(3) * t22;
t51 = -g(1) * (t35 * t96 + t88) + g(2) * (t35 * t87 - t98) + t35 * t101;
t49 = -g(1) * t55 - g(3) * t89;
t42 = -qJ(5) - pkin(8);
t30 = t47 * pkin(4) + pkin(3);
t24 = pkin(2) * t64;
t23 = t63 - t58 / 0.2e1;
t19 = -g(3) * t41 - t54;
t14 = t34 * t69 - t81;
t13 = t34 * t65 + t82;
t12 = t34 * t70 + t79;
t11 = -t34 * t66 + t80;
t10 = t33 * t100 + t86;
t8 = t32 * t91 + t43 * t99;
t7 = -t32 * t93 + t46 * t99;
t2 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(2) * (t64 - t95) + t49, -g(1) * (t38 * t85 - t40 * t48) - g(2) * (-t38 * t48 - t40 * t85) + t45 * t101, -g(2) * t24 + (g(2) * t95 + t49) * pkin(2), 0, 0, 0, 0, 0, -t51 * t47, t51 * t44, t52 * t34, -t52 * t32, -g(1) * t61 - g(2) * t62 - g(3) * t23, -g(1) * (pkin(2) * t55 - t4 * t30 - t42 * t61) - g(2) * (-pkin(2) * t95 - t1 * t30 - t42 * t62 + t24) - g(3) * (pkin(2) * t89 + t22 * t30 - t23 * t42), 0, 0, 0, 0, 0, (g(1) * t14 - g(2) * t13 - t34 * t77) * t35 + (-g(1) * (-t70 - t73) - g(2) * (t66 - t75) - t78) * t33, (-g(1) * t12 - g(2) * t11 + t34 * t78) * t35 + (-g(1) * (-t69 + t74) - g(2) * (t65 + t76) - t77) * t33; 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t44 + t72) - g(2) * (-t106 * t44 - t68) - g(3) * (-t33 * t92 + t84), -g(1) * (t105 * t47 - t38 * t92) - g(2) * (-t106 * t47 + t40 * t92) - g(3) * (-t33 * t90 - t41 * t44), -g(1) * (-t32 * t61 + t34 * t97) - g(2) * (-t32 * t62 - t34 * t94) - g(3) * (-t23 * t32 + t27), -g(1) * (-t32 * t97 - t34 * t61) - g(2) * (t32 * t94 - t34 * t62) - g(3) * (-t23 * t34 - t86), 0, (-g(1) * (-t44 * t61 + t72) - g(2) * (-t44 * t62 - t68) - g(3) * (-t23 * t44 + t84)) * pkin(4), 0, 0, 0, 0, 0, (-t34 * t54 - t107) * t46, (t59 * t100 + t107) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t33 - t35 * t74 + t38 * t7) - g(2) * (t11 * t33 - t35 * t76 - t40 * t7) - g(3) * (-t10 * t43 - t35 * t91), -g(1) * (t14 * t33 - t35 * t73 - t38 * t8) - g(2) * (-t13 * t33 - t35 * t75 + t40 * t8) - g(3) * (-t10 * t46 + t35 * t93);];
taug_reg = t2;
