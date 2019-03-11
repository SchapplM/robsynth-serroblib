% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t68 = cos(pkin(6));
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t66 = sin(pkin(6));
t71 = sin(qJ(2));
t99 = t66 * t71;
t42 = -t68 * t73 + t70 * t99;
t43 = t68 * t70 + t73 * t99;
t65 = sin(pkin(11));
t67 = cos(pkin(11));
t29 = -t65 * t42 + t67 * t43;
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t74 = cos(qJ(2));
t98 = t66 * t74;
t22 = t69 * t29 + t72 * t98;
t108 = -0.2e1 * t22;
t107 = -0.2e1 * t43;
t106 = 0.2e1 * t73;
t105 = pkin(1) * t71;
t104 = pkin(1) * t74;
t103 = t72 * pkin(5);
t23 = t72 * t29 - t69 * t98;
t102 = t23 * t69;
t101 = t23 * t72;
t62 = t66 ^ 2;
t100 = t62 * t74;
t97 = t68 * t71;
t28 = t67 * t42 + t65 * t43;
t96 = t69 * t28;
t49 = t65 * t70 - t67 * t73;
t95 = t69 * t49;
t50 = t65 * t73 + t67 * t70;
t94 = t69 * t50;
t59 = t65 * pkin(3) + pkin(10);
t93 = t69 * t59;
t92 = t69 * t72;
t88 = -qJ(4) - pkin(9);
t53 = t88 * t73;
t80 = t88 * t70;
t35 = -t67 * t53 + t65 * t80;
t91 = t72 * t35;
t90 = t72 * t50;
t89 = t72 * t59;
t83 = pkin(8) * t98;
t38 = t83 + (pkin(9) + t105) * t68;
t39 = (-pkin(2) * t74 - pkin(9) * t71 - pkin(1)) * t66;
t24 = -t70 * t38 + t73 * t39;
t15 = -pkin(3) * t98 - t43 * qJ(4) + t24;
t25 = t73 * t38 + t70 * t39;
t20 = -t42 * qJ(4) + t25;
t9 = t65 * t15 + t67 * t20;
t63 = t69 ^ 2;
t64 = t72 ^ 2;
t87 = t63 + t64;
t86 = qJ(6) * t50;
t85 = qJ(6) + t59;
t84 = 0.2e1 * t98;
t82 = t70 * t98;
t81 = t73 * t98;
t61 = -t73 * pkin(3) - pkin(2);
t60 = -t67 * pkin(3) - pkin(4);
t55 = pkin(8) * t99;
t37 = t55 + (-pkin(2) - t104) * t68;
t31 = t42 * pkin(3) + t37;
t11 = t28 * pkin(4) - t29 * pkin(10) + t31;
t7 = -pkin(10) * t98 + t9;
t3 = t72 * t11 - t69 * t7;
t8 = t67 * t15 - t65 * t20;
t32 = t49 * pkin(4) - t50 * pkin(10) + t61;
t16 = t72 * t32 - t69 * t35;
t33 = -t65 * t53 - t67 * t80;
t6 = pkin(4) * t98 - t8;
t1 = t28 * pkin(5) - t23 * qJ(6) + t3;
t4 = t69 * t11 + t72 * t7;
t2 = -t22 * qJ(6) + t4;
t79 = t1 * t72 + t2 * t69;
t12 = t49 * pkin(5) - t72 * t86 + t16;
t13 = t91 + (t32 - t86) * t69;
t78 = t12 * t72 + t13 * t69;
t46 = t85 * t69;
t47 = t85 * t72;
t77 = -t46 * t72 + t47 * t69;
t76 = -t49 * t59 + t50 * t60;
t52 = t60 - t103;
t48 = t50 ^ 2;
t45 = pkin(1) * t97 + t83;
t44 = t68 * t104 - t55;
t41 = t72 * t49;
t27 = t72 * t28;
t26 = pkin(5) * t94 + t33;
t21 = t69 * t22;
t17 = t69 * t32 + t91;
t5 = t22 * pkin(5) + t6;
t10 = [1, 0, 0, t62 * t71 ^ 2, 0.2e1 * t71 * t100, 0.2e1 * t66 * t97, t68 * t84, t68 ^ 2, 0.2e1 * pkin(1) * t100 + 0.2e1 * t44 * t68, -0.2e1 * t62 * t105 - 0.2e1 * t45 * t68, t43 ^ 2, t42 * t107, t98 * t107, t42 * t84, t62 * t74 ^ 2, -0.2e1 * t24 * t98 + 0.2e1 * t37 * t42, 0.2e1 * t25 * t98 + 0.2e1 * t37 * t43, -0.2e1 * t9 * t28 - 0.2e1 * t8 * t29, t31 ^ 2 + t8 ^ 2 + t9 ^ 2, t23 ^ 2, t23 * t108, 0.2e1 * t23 * t28, t28 * t108, t28 ^ 2, 0.2e1 * t6 * t22 + 0.2e1 * t3 * t28, 0.2e1 * t6 * t23 - 0.2e1 * t4 * t28, -0.2e1 * t1 * t23 - 0.2e1 * t2 * t22, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t99, t98, t68, t44, -t45, t43 * t70, -t70 * t42 + t43 * t73, -t82, -t81, 0, -pkin(2) * t42 + pkin(9) * t82 - t37 * t73, -pkin(2) * t43 + pkin(9) * t81 + t37 * t70, -t35 * t28 + t33 * t29 - t9 * t49 - t8 * t50, t31 * t61 - t8 * t33 + t9 * t35, t23 * t90 (-t22 * t72 - t102) * t50, t23 * t49 + t28 * t90, -t22 * t49 - t28 * t94, t28 * t49, t16 * t28 + t33 * t22 + t3 * t49 + t6 * t94, -t17 * t28 + t33 * t23 - t4 * t49 + t6 * t90, -t12 * t23 - t13 * t22 - t79 * t50, t1 * t12 + t2 * t13 + t5 * t26; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t70 ^ 2, t70 * t106, 0, 0, 0, pkin(2) * t106, -0.2e1 * pkin(2) * t70, 0.2e1 * t33 * t50 - 0.2e1 * t35 * t49, t33 ^ 2 + t35 ^ 2 + t61 ^ 2, t64 * t48, -0.2e1 * t48 * t92, 0.2e1 * t49 * t90, -0.2e1 * t49 * t94, t49 ^ 2, 0.2e1 * t16 * t49 + 0.2e1 * t33 * t94, -0.2e1 * t17 * t49 + 0.2e1 * t33 * t90, -0.2e1 * t78 * t50, t12 ^ 2 + t13 ^ 2 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t42, -t98, t24, -t25 (-t28 * t65 - t29 * t67) * pkin(3) (t65 * t9 + t67 * t8) * pkin(3), t102, -t21 + t101, t96, t27, 0, t60 * t22 - t28 * t93 - t6 * t72, t60 * t23 - t28 * t89 + t6 * t69, -t1 * t69 + t2 * t72 - t47 * t22 + t46 * t23, -t1 * t46 + t2 * t47 + t5 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t73, 0, -t70 * pkin(9), -t73 * pkin(9) (-t49 * t65 - t50 * t67) * pkin(3) (-t33 * t67 + t35 * t65) * pkin(3), t69 * t90 (-t63 + t64) * t50, t95, t41, 0, -t33 * t72 + t76 * t69, t33 * t69 + t76 * t72, -t12 * t69 + t13 * t72 - t77 * t50, -t12 * t46 + t13 * t47 + t26 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t65 ^ 2 + t67 ^ 2) * pkin(3) ^ 2, t63, 0.2e1 * t92, 0, 0, 0, -0.2e1 * t60 * t72, 0.2e1 * t60 * t69, 0.2e1 * t46 * t69 + 0.2e1 * t47 * t72, t46 ^ 2 + t47 ^ 2 + t52 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, t27, -t96, -t21 - t101, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, t41, -t95, -t87 * t50, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, t28, t3, -t4, -pkin(5) * t23, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t94, t49, t16, -t17, -pkin(5) * t90, t12 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t72, 0, -t93, -t89, -t69 * pkin(5), -t46 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t69, 0, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t10;
