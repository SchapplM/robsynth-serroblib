% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x34]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t67 = sin(qJ(4));
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t84 = cos(qJ(4));
t45 = t67 * t71 + t84 * t68;
t40 = t45 ^ 2;
t44 = t67 * t68 - t84 * t71;
t42 = t44 ^ 2;
t99 = -t40 - t42;
t65 = sin(qJ(6));
t66 = sin(qJ(5));
t69 = cos(qJ(6));
t70 = cos(qJ(5));
t46 = t65 * t70 + t69 * t66;
t28 = t44 * t46;
t98 = 0.2e1 * t28;
t97 = -0.2e1 * t44;
t96 = 0.2e1 * t45;
t75 = t84 * pkin(3);
t58 = -t75 - pkin(4);
t87 = t70 * pkin(5);
t50 = t58 - t87;
t95 = 0.2e1 * t50;
t59 = -pkin(4) - t87;
t94 = 0.2e1 * t59;
t93 = 2 * qJ(2);
t92 = t45 * pkin(5);
t91 = t65 * pkin(5);
t90 = t67 * pkin(3);
t89 = t69 * pkin(5);
t54 = t68 * pkin(3) + qJ(2);
t20 = t45 * pkin(4) + t44 * pkin(9) + t54;
t72 = -pkin(1) - pkin(7);
t48 = (-pkin(8) + t72) * t68;
t60 = t71 * t72;
t49 = -t71 * pkin(8) + t60;
t24 = t84 * t48 + t67 * t49;
t79 = t70 * t24;
t6 = t79 + (pkin(10) * t44 + t20) * t66;
t88 = t69 * t6;
t86 = t70 * pkin(9);
t85 = pkin(4) - t58;
t23 = t67 * t48 - t84 * t49;
t83 = t23 * t70;
t43 = t65 * t66 - t69 * t70;
t82 = t44 * t43;
t36 = t44 * t66;
t81 = t44 * t70;
t29 = t46 * t45;
t34 = t66 * t45;
t80 = t66 * t70;
t35 = t70 * t45;
t57 = pkin(9) + t90;
t78 = t70 * t57;
t77 = t50 + t59;
t76 = t44 * t96;
t7 = t70 * t20 - t66 * t24;
t5 = pkin(10) * t81 + t7 + t92;
t1 = t69 * t5 - t65 * t6;
t74 = pkin(4) * t44 - pkin(9) * t45;
t73 = -t44 * t58 - t45 * t57;
t64 = t70 ^ 2;
t63 = t66 ^ 2;
t62 = t70 * pkin(10);
t53 = 0.2e1 * t80;
t52 = t62 + t86;
t51 = (-pkin(9) - pkin(10)) * t66;
t41 = t46 ^ 2;
t38 = t62 + t78;
t37 = (-pkin(10) - t57) * t66;
t33 = t44 * t80;
t31 = t65 * t51 + t69 * t52;
t30 = t69 * t51 - t65 * t52;
t26 = t43 * t45;
t25 = -0.2e1 * t46 * t43;
t22 = t23 * t66;
t21 = (t63 - t64) * t44;
t19 = t65 * t37 + t69 * t38;
t18 = t69 * t37 - t65 * t38;
t15 = -t65 * t34 + t69 * t35;
t12 = -pkin(5) * t36 + t23;
t11 = t82 * t46;
t10 = t12 * t46;
t9 = t12 * t43;
t8 = t66 * t20 + t79;
t4 = t28 * t46 - t43 * t82;
t2 = t65 * t5 + t88;
t3 = [1, 0, 0, -2 * pkin(1), t93, pkin(1) ^ 2 + qJ(2) ^ 2, t71 ^ 2, -0.2e1 * t71 * t68, 0, 0, 0, t68 * t93, t71 * t93, t42, t76, 0, 0, 0, t54 * t96, t54 * t97, t64 * t42, -0.2e1 * t42 * t80, t35 * t97, t66 * t76, t40, -0.2e1 * t23 * t36 + 0.2e1 * t7 * t45, -0.2e1 * t23 * t81 - 0.2e1 * t8 * t45, t82 ^ 2, t82 * t98, t82 * t96, t45 * t98, t40, 0.2e1 * t1 * t45 - 0.2e1 * t12 * t28, 0.2e1 * t12 * t82 - 0.2e1 * t2 * t45; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99 * t66, t99 * t70, 0, 0, 0, 0, 0, -t28 * t44 - t29 * t45, -t15 * t45 + t44 * t82; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t71, -t68, 0, t60, -t68 * t72, 0, 0, -t44, -t45, 0, -t23, -t24, -t33, t21, t34, t35, 0, t73 * t66 - t83, t73 * t70 + t22, t11, t4, t29, -t26, 0, t18 * t45 - t28 * t50 + t9, -t19 * t45 + t50 * t82 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t68, 0, 0, 0, 0, 0, -t44, -t45, 0, 0, 0, 0, 0, -t81, t36, 0, 0, 0, 0, 0, t82, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t75, -0.2e1 * t90, t63, t53, 0, 0, 0, -0.2e1 * t58 * t70, 0.2e1 * t58 * t66, t41, t25, 0, 0, 0, t43 * t95, t46 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t45, 0, -t23, -t24, -t33, t21, t34, t35, 0, t74 * t66 - t83, t74 * t70 + t22, t11, t4, t29, -t26, 0, -t28 * t59 + t30 * t45 + t9, -t31 * t45 + t59 * t82 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t45, 0, 0, 0, 0, 0, -t81, t36, 0, 0, 0, 0, 0, t82, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t75, -t90, t63, t53, 0, 0, 0, t85 * t70, -t85 * t66, t41, t25, 0, 0, 0, t77 * t43, t77 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t63, t53, 0, 0, 0, 0.2e1 * pkin(4) * t70, -0.2e1 * pkin(4) * t66, t41, t25, 0, 0, 0, t43 * t94, t46 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t36, t45, t7, -t8, 0, 0, t82, t28, t45, t45 * t89 + t1, -t88 + (-t5 - t92) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t35, 0, 0, 0, 0, 0, -t29, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t70, 0, -t66 * t57, -t78, 0, 0, t46, -t43, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t70, 0, -t66 * pkin(9), -t86, 0, 0, t46, -t43, 0, t30, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t89, -0.2e1 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t28, t45, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t43, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t43, 0, t30, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t89, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
