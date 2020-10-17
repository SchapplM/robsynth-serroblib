% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:05:18
% EndTime: 2019-05-07 11:05:22
% DurationCPUTime: 0.92s
% Computational Cost: add. (1717->138), mult. (3581->244), div. (0->0), fcn. (4203->10), ass. (0->90)
t76 = sin(pkin(11));
t77 = cos(pkin(11));
t80 = sin(qJ(3));
t84 = cos(qJ(3));
t55 = -t76 * t80 + t77 * t84;
t56 = t76 * t84 + t77 * t80;
t79 = sin(qJ(5));
t83 = cos(qJ(5));
t36 = -t83 * t55 + t79 * t56;
t69 = -t84 * pkin(3) - pkin(2);
t48 = -t55 * pkin(4) + t69;
t24 = t36 * pkin(5) + t48;
t107 = 0.2e1 * t24;
t106 = 0.2e1 * t48;
t81 = sin(qJ(2));
t105 = -0.2e1 * t81;
t85 = cos(qJ(2));
t104 = -0.2e1 * t85;
t103 = 0.2e1 * t85;
t102 = pkin(2) * t84;
t101 = pkin(3) * t76;
t100 = pkin(7) * t80;
t78 = sin(qJ(6));
t99 = t78 * pkin(5);
t49 = t56 * t81;
t92 = t84 * t81;
t96 = t80 * t81;
t50 = -t76 * t96 + t77 * t92;
t27 = t83 * t49 + t79 * t50;
t62 = -t85 * pkin(2) - t81 * pkin(8) - pkin(1);
t57 = t84 * t62;
t89 = qJ(4) * t81;
t38 = -t84 * t89 + t57 + (-pkin(3) - t100) * t85;
t91 = t84 * t85;
t87 = pkin(7) * t91;
t41 = t87 + (t62 - t89) * t80;
t22 = t77 * t38 - t76 * t41;
t15 = -t85 * pkin(4) - t50 * pkin(9) + t22;
t23 = t76 * t38 + t77 * t41;
t18 = -t49 * pkin(9) + t23;
t9 = t79 * t15 + t83 * t18;
t7 = -t27 * pkin(10) + t9;
t82 = cos(qJ(6));
t98 = t82 * t7;
t97 = t85 * pkin(5);
t95 = t80 * t84;
t94 = t80 * t85;
t68 = t77 * pkin(3) + pkin(4);
t53 = t83 * t101 + t79 * t68;
t93 = t82 * t53;
t90 = -qJ(4) - pkin(8);
t63 = t90 * t80;
t64 = t90 * t84;
t43 = t76 * t63 - t77 * t64;
t70 = t81 * pkin(7);
t61 = pkin(3) * t96 + t70;
t88 = t81 * t103;
t28 = -t79 * t49 + t83 * t50;
t8 = t83 * t15 - t79 * t18;
t6 = -t28 * pkin(10) + t8 - t97;
t1 = t82 * t6 - t78 * t7;
t42 = t77 * t63 + t76 * t64;
t31 = -t56 * pkin(9) + t42;
t32 = t55 * pkin(9) + t43;
t16 = t83 * t31 - t79 * t32;
t52 = -t79 * t101 + t83 * t68;
t51 = pkin(5) + t52;
t29 = t82 * t51 - t78 * t53;
t39 = t49 * pkin(4) + t61;
t2 = t78 * t6 + t98;
t17 = t79 * t31 + t83 * t32;
t75 = t85 ^ 2;
t74 = t84 ^ 2;
t73 = t81 ^ 2;
t72 = t80 ^ 2;
t71 = t82 * pkin(5);
t47 = t80 * t62 + t87;
t46 = -pkin(7) * t94 + t57;
t37 = t79 * t55 + t83 * t56;
t30 = t78 * t51 + t93;
t21 = -t78 * t36 + t82 * t37;
t20 = t82 * t36 + t78 * t37;
t19 = t27 * pkin(5) + t39;
t13 = -t78 * t27 + t82 * t28;
t12 = t82 * t27 + t78 * t28;
t11 = -t36 * pkin(10) + t17;
t10 = -t37 * pkin(10) + t16;
t4 = t78 * t10 + t82 * t11;
t3 = t82 * t10 - t78 * t11;
t5 = [1, 0, 0, t73, t88, 0, 0, 0, pkin(1) * t103, pkin(1) * t105, t74 * t73, -0.2e1 * t73 * t95, t91 * t105, t80 * t88, t75, 0.2e1 * t73 * t100 - 0.2e1 * t46 * t85, 0.2e1 * t73 * pkin(7) * t84 + 0.2e1 * t47 * t85, -0.2e1 * t22 * t50 - 0.2e1 * t23 * t49, t22 ^ 2 + t23 ^ 2 + t61 ^ 2, t28 ^ 2, -0.2e1 * t28 * t27, t28 * t104, t27 * t103, t75, 0.2e1 * t39 * t27 - 0.2e1 * t8 * t85, 0.2e1 * t39 * t28 + 0.2e1 * t9 * t85, t13 ^ 2, -0.2e1 * t13 * t12, t13 * t104, t12 * t103, t75, -0.2e1 * t1 * t85 + 0.2e1 * t19 * t12, 0.2e1 * t19 * t13 + 0.2e1 * t2 * t85; 0, 0, 0, 0, 0, t81, t85, 0, -t70, -t85 * pkin(7), t80 * t92 (-t72 + t74) * t81, -t94, -t91, 0, -pkin(7) * t92 + (-pkin(2) * t81 + pkin(8) * t85) * t80, pkin(8) * t91 + (t100 - t102) * t81, -t22 * t56 + t23 * t55 - t42 * t50 - t43 * t49, t22 * t42 + t23 * t43 + t61 * t69, t28 * t37, -t37 * t27 - t28 * t36, -t37 * t85, t36 * t85, 0, -t16 * t85 + t48 * t27 + t39 * t36, t17 * t85 + t48 * t28 + t39 * t37, t13 * t21, -t21 * t12 - t13 * t20, -t21 * t85, t20 * t85, 0, t24 * t12 + t19 * t20 - t3 * t85, t24 * t13 + t19 * t21 + t4 * t85; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t72, 0.2e1 * t95, 0, 0, 0, 0.2e1 * t102, -0.2e1 * pkin(2) * t80, -0.2e1 * t42 * t56 + 0.2e1 * t43 * t55, t42 ^ 2 + t43 ^ 2 + t69 ^ 2, t37 ^ 2, -0.2e1 * t37 * t36, 0, 0, 0, t36 * t106, t37 * t106, t21 ^ 2, -0.2e1 * t21 * t20, 0, 0, 0, t20 * t107, t21 * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, -t96, -t85, t46, -t47 (-t49 * t76 - t50 * t77) * pkin(3) (t22 * t77 + t23 * t76) * pkin(3), 0, 0, t28, -t27, -t85, -t52 * t85 + t8, t53 * t85 - t9, 0, 0, t13, -t12, -t85, -t29 * t85 + t1, t30 * t85 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t84, 0, -t80 * pkin(8), -t84 * pkin(8) (t55 * t76 - t56 * t77) * pkin(3) (t42 * t77 + t43 * t76) * pkin(3), 0, 0, t37, -t36, 0, t16, -t17, 0, 0, t21, -t20, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t76 ^ 2 + t77 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t52, -0.2e1 * t53, 0, 0, 0, 0, 1, 0.2e1 * t29, -0.2e1 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, t27, t28, 0, 0, 0, 0, 0, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, t36, t37, 0, 0, 0, 0, 0, t20, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, -t85, t8, -t9, 0, 0, t13, -t12, -t85, -t82 * t97 + t1, -t98 + (-t6 + t97) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, 0, t16, -t17, 0, 0, t21, -t20, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t52, -t53, 0, 0, 0, 0, 1, t29 + t71, -t93 + (-pkin(5) - t51) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t71, -0.2e1 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, -t85, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t71, -t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
