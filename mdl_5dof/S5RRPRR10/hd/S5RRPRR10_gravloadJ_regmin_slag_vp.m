% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:02:42
% EndTime: 2021-01-15 22:02:46
% DurationCPUTime: 0.50s
% Computational Cost: add. (280->94), mult. (442->167), div. (0->0), fcn. (533->16), ass. (0->84)
t43 = cos(pkin(5));
t46 = sin(qJ(2));
t51 = cos(qJ(1));
t72 = t51 * t46;
t47 = sin(qJ(1));
t50 = cos(qJ(2));
t78 = t47 * t50;
t23 = -t43 * t78 - t72;
t71 = t51 * t50;
t80 = t47 * t46;
t56 = t43 * t71 - t80;
t42 = sin(pkin(5));
t89 = g(3) * t42;
t98 = -g(1) * t23 - g(2) * t56 - t50 * t89;
t41 = qJ(2) + pkin(10);
t39 = sin(t41);
t73 = t51 * t39;
t40 = cos(t41);
t81 = t47 * t40;
t10 = t43 * t73 + t81;
t45 = sin(qJ(4));
t49 = cos(qJ(4));
t86 = t42 * t45;
t90 = g(3) * (-t39 * t86 + t43 * t49);
t29 = t51 * t40;
t82 = t47 * t39;
t95 = t43 * t82 - t29;
t97 = -(g(1) * t95 - g(2) * t10) * t45 - t90;
t38 = pkin(5) - t41;
t96 = sin(t38) / 0.2e1;
t37 = pkin(5) + t41;
t30 = sin(t37);
t92 = -t30 / 0.2e1 + t96;
t88 = t40 * t42;
t44 = sin(qJ(5));
t87 = t40 * t44;
t85 = t42 * t49;
t84 = t44 * t47;
t83 = t44 * t51;
t79 = t47 * t49;
t48 = cos(qJ(5));
t77 = t48 * t47;
t76 = t48 * t49;
t75 = t48 * t51;
t74 = t49 * t51;
t69 = t47 * t86;
t67 = t47 * t76;
t66 = t44 * t79;
t28 = t51 * t86;
t65 = t44 * t74;
t64 = t48 * t74;
t63 = t10 * t49 - t28;
t62 = g(1) * t51 + g(2) * t47;
t61 = g(1) * t47 - g(2) * t51;
t32 = cos(t37);
t33 = cos(t38);
t27 = t32 + t33;
t60 = t82 - t51 * t27 / 0.2e1;
t59 = t73 + t47 * t27 / 0.2e1;
t58 = -t49 * t95 + t69;
t55 = t61 * t42;
t54 = -g(1) * (t43 * t81 + t73) + g(2) * (t43 * t29 - t82) + g(3) * t88;
t12 = t43 * t67 - t83;
t20 = t43 * t84 + t64;
t53 = -t12 * t39 + t20 * t40 + t48 * t69;
t13 = t43 * t65 - t77;
t18 = -t43 * t75 - t66;
t52 = -t13 * t39 + t18 * t40 + t44 * t28;
t36 = t50 * pkin(2) + pkin(1);
t25 = t43 * t46 * pkin(2) - t42 * (pkin(7) + qJ(3));
t24 = -t43 * t80 + t71;
t22 = -t43 * t72 - t78;
t19 = t43 * t77 - t65;
t17 = t43 * t83 - t67;
t16 = t39 * t85 + t43 * t45;
t14 = t43 * t64 + t84;
t11 = t43 * t66 + t75;
t6 = t47 * t92 + t29;
t5 = t51 * t92 - t81;
t4 = t42 * t79 + t45 * t95;
t3 = t10 * t45 + t42 * t74;
t2 = t11 * t39 + t19 * t40 - t44 * t69;
t1 = -t14 * t39 + t17 * t40 + t48 * t28;
t7 = [0, t61, t62, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t24, g(1) * t56 - g(2) * t23, -g(1) * t5 - g(2) * t6, -g(1) * t60 + g(2) * t59, -t62 * t42, -g(1) * (-t25 * t51 - t47 * t36) - g(2) * (-t47 * t25 + t51 * t36), 0, 0, 0, 0, 0, g(1) * t63 - g(2) * t58, -g(1) * t3 - g(2) * t4, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t53, g(1) * t52 - g(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, t98, g(1) * t24 - g(2) * t22 + t46 * t89, g(1) * t59 + g(2) * t60 - g(3) * (t96 + t30 / 0.2e1), g(1) * t6 - g(2) * t5 - g(3) * (-t33 / 0.2e1 + t32 / 0.2e1), 0, t98 * pkin(2), 0, 0, 0, 0, 0, -t54 * t49, t54 * t45, 0, 0, 0, 0, 0, -g(1) * (-t12 * t40 - t20 * t39) - g(2) * (t14 * t40 + t17 * t39) - (t39 * t44 + t40 * t76) * t89, -g(1) * (t11 * t40 - t19 * t39) - g(2) * (-t13 * t40 - t18 * t39) - (t39 * t48 - t49 * t87) * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t43 - t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 + g(2) * t3 - t90, g(1) * t58 + g(2) * t63 + g(3) * t16, 0, 0, 0, 0, 0, (-t49 * t55 + t97) * t48, (t61 * t85 - t97) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t52 - g(3) * (-t16 * t44 - t48 * t88), g(1) * t53 - g(2) * t1 - g(3) * (-t16 * t48 + t42 * t87);];
taug_reg = t7;
