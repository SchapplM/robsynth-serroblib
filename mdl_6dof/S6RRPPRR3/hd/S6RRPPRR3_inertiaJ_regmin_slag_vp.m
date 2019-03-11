% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t102 = cos(qJ(5));
t68 = sin(pkin(11));
t69 = sin(pkin(6));
t71 = cos(pkin(11));
t75 = sin(qJ(2));
t77 = cos(qJ(2));
t44 = (t68 * t77 + t71 * t75) * t69;
t67 = sin(pkin(12));
t70 = cos(pkin(12));
t72 = cos(pkin(6));
t35 = t44 * t67 - t72 * t70;
t36 = t44 * t70 + t72 * t67;
t74 = sin(qJ(5));
t19 = t102 * t35 + t74 * t36;
t106 = -0.2e1 * t19;
t61 = -t71 * pkin(2) - pkin(3);
t55 = -t70 * pkin(4) + t61;
t105 = 0.2e1 * t55;
t104 = pkin(1) * t75;
t59 = t68 * pkin(2) + qJ(4);
t103 = pkin(9) + t59;
t20 = t102 * t36 - t74 * t35;
t95 = t69 * t77;
t96 = t69 * t75;
t43 = t68 * t96 - t71 * t95;
t73 = sin(qJ(6));
t76 = cos(qJ(6));
t14 = t76 * t20 + t43 * t73;
t101 = t14 * t73;
t100 = t43 * t67;
t13 = t73 * t20 - t43 * t76;
t90 = t74 * t67;
t52 = -t102 * t70 + t90;
t99 = t52 * t13;
t81 = t102 * t67;
t53 = t74 * t70 + t81;
t98 = t53 * t43;
t63 = t69 ^ 2;
t97 = t63 * t77;
t94 = t70 * t43;
t93 = t73 * t19;
t46 = t73 * t52;
t92 = t73 * t53;
t91 = t73 * t76;
t18 = t76 * t19;
t89 = t76 * t53;
t88 = pkin(8) + qJ(3);
t57 = t72 * t77 * pkin(1);
t37 = t72 * pkin(2) - t88 * t96 + t57;
t84 = t72 * t104;
t41 = t88 * t95 + t84;
t27 = t68 * t37 + t71 * t41;
t22 = t72 * qJ(4) + t27;
t54 = (-pkin(2) * t77 - pkin(1)) * t69;
t28 = t43 * pkin(3) - t44 * qJ(4) + t54;
t11 = t70 * t22 + t67 * t28;
t87 = t67 ^ 2 + t70 ^ 2;
t86 = -0.2e1 * t53 * t52;
t85 = 0.2e1 * t69 * t72;
t83 = t19 * t92;
t82 = t19 * t89;
t10 = -t67 * t22 + t70 * t28;
t26 = t71 * t37 - t68 * t41;
t80 = -pkin(5) * t53 - pkin(10) * t52;
t79 = -t10 * t67 + t11 * t70;
t23 = -t72 * pkin(3) - t26;
t8 = t43 * pkin(4) - t36 * pkin(9) + t10;
t9 = -t35 * pkin(9) + t11;
t5 = t102 * t8 - t74 * t9;
t6 = t102 * t9 + t74 * t8;
t15 = t35 * pkin(4) + t23;
t66 = t76 ^ 2;
t65 = t73 ^ 2;
t51 = t53 ^ 2;
t50 = pkin(8) * t95 + t84;
t49 = -pkin(8) * t96 + t57;
t48 = t103 * t70;
t47 = t76 * t52;
t33 = t52 * t43;
t31 = t102 * t48 - t103 * t90;
t30 = t103 * t81 + t74 * t48;
t29 = t52 * pkin(5) - t53 * pkin(10) + t55;
t17 = t73 * t29 + t76 * t31;
t16 = t76 * t29 - t73 * t31;
t12 = t14 * t52;
t7 = t19 * pkin(5) - t20 * pkin(10) + t15;
t4 = t43 * pkin(10) + t6;
t3 = -t43 * pkin(5) - t5;
t2 = t76 * t4 + t73 * t7;
t1 = -t73 * t4 + t76 * t7;
t21 = [1, 0, 0, t63 * t75 ^ 2, 0.2e1 * t75 * t97, t75 * t85, t77 * t85, t72 ^ 2, 0.2e1 * pkin(1) * t97 + 0.2e1 * t49 * t72, -0.2e1 * t63 * t104 - 0.2e1 * t50 * t72, -0.2e1 * t26 * t44 - 0.2e1 * t27 * t43, t26 ^ 2 + t27 ^ 2 + t54 ^ 2, 0.2e1 * t10 * t43 + 0.2e1 * t23 * t35, -0.2e1 * t11 * t43 + 0.2e1 * t23 * t36, -0.2e1 * t10 * t36 - 0.2e1 * t11 * t35, t10 ^ 2 + t11 ^ 2 + t23 ^ 2, t20 ^ 2, t20 * t106, 0.2e1 * t20 * t43, t43 * t106, t43 ^ 2, 0.2e1 * t15 * t19 + 0.2e1 * t5 * t43, 0.2e1 * t15 * t20 - 0.2e1 * t6 * t43, t14 ^ 2, -0.2e1 * t14 * t13, 0.2e1 * t14 * t19, t13 * t106, t19 ^ 2, 0.2e1 * t1 * t19 + 0.2e1 * t3 * t13, 0.2e1 * t3 * t14 - 0.2e1 * t2 * t19; 0, 0, 0, 0, 0, t96, t95, t72, t49, -t50 (-t43 * t68 - t44 * t71) * pkin(2) (t26 * t71 + t27 * t68) * pkin(2), -t59 * t100 - t23 * t70 + t61 * t35, t23 * t67 + t61 * t36 - t59 * t94 (-t35 * t70 + t36 * t67) * t59 + t79, t23 * t61 + t79 * t59, t20 * t53, -t53 * t19 - t20 * t52, t98, -t33, 0, t15 * t52 + t55 * t19 - t30 * t43, t15 * t53 + t55 * t20 - t31 * t43, t14 * t89 (-t13 * t76 - t101) * t53, t12 + t82, -t83 - t99, t19 * t52, t1 * t52 + t30 * t13 + t16 * t19 + t3 * t92, t30 * t14 - t17 * t19 - t2 * t52 + t3 * t89; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t68 ^ 2 + t71 ^ 2) * pkin(2) ^ 2, -0.2e1 * t61 * t70, 0.2e1 * t61 * t67, 0.2e1 * t87 * t59, t87 * t59 ^ 2 + t61 ^ 2, t51, t86, 0, 0, 0, t52 * t105, t53 * t105, t66 * t51, -0.2e1 * t51 * t91, 0.2e1 * t52 * t89, t73 * t86, t52 ^ 2, 0.2e1 * t16 * t52 + 0.2e1 * t30 * t92, -0.2e1 * t17 * t52 + 0.2e1 * t30 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t94, -t100, -t67 * t35 - t70 * t36, t10 * t70 + t11 * t67, 0, 0, 0, 0, 0, -t33, -t98, 0, 0, 0, 0, 0, -t83 + t99, t12 - t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t36, 0, t23, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, t18, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t67, 0, t61, 0, 0, 0, 0, 0, t52, t53, 0, 0, 0, 0, 0, t47, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t43, t5, -t6, t101, -t73 * t13 + t14 * t76, t93, t18, 0, -pkin(5) * t13 - pkin(10) * t93 - t3 * t76, -pkin(5) * t14 - pkin(10) * t18 + t3 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t52, 0, -t30, -t31, t73 * t89 (-t65 + t66) * t53, t46, t47, 0, -t30 * t76 + t80 * t73, t30 * t73 + t80 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t53, 0, 0, 0, 0, 0, -t47, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t65, 0.2e1 * t91, 0, 0, 0, 0.2e1 * pkin(5) * t76, -0.2e1 * pkin(5) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t19, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -t92, t52, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t76, 0, -t73 * pkin(10), -t76 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t21;
