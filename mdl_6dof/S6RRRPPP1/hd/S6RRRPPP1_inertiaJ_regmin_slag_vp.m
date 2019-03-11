% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t70 = sin(qJ(3));
t72 = cos(qJ(3));
t67 = sin(pkin(6));
t86 = qJ(4) * t67;
t47 = -pkin(3) * t72 - t70 * t86 - pkin(2);
t69 = cos(pkin(6));
t80 = qJ(4) * t69 + pkin(9);
t49 = t80 * t70;
t68 = cos(pkin(10));
t106 = (t47 * t67 - t49 * t69) * t68;
t73 = cos(qJ(2));
t71 = sin(qJ(2));
t93 = t70 * t71;
t77 = t67 * t73 + t69 * t93;
t52 = -pkin(2) * t73 - pkin(9) * t71 - pkin(1);
t89 = t72 * t73;
t34 = pkin(8) * t89 + t70 * t52;
t21 = -qJ(4) * t77 + t34;
t100 = pkin(8) * t70;
t48 = t72 * t52;
t85 = qJ(4) * t72;
t22 = -t71 * t69 * t85 + t48 + (-pkin(3) - t100) * t73;
t31 = (pkin(3) * t70 - t67 * t85 + pkin(8)) * t71;
t66 = sin(pkin(10));
t7 = -t66 * t21 + (t22 * t69 + t31 * t67) * t68;
t105 = 0.2e1 * t67;
t104 = -0.2e1 * t71;
t103 = 0.2e1 * t73;
t102 = pkin(2) * t72;
t101 = pkin(3) * t68;
t99 = t67 * pkin(3);
t98 = t66 * t69;
t59 = t67 * t66;
t60 = t67 * t68;
t97 = t67 * t72;
t95 = t68 * t69;
t94 = t69 * t72;
t92 = t70 * t72;
t91 = t70 * t73;
t90 = t72 * t71;
t88 = pkin(4) + qJ(6);
t50 = t80 * t72;
t44 = t66 * t50;
t87 = pkin(4) * t97 + t44;
t43 = pkin(3) * t98 + t68 * t86;
t84 = t71 * t103;
t8 = t21 * t68 + t22 * t98 + t31 * t59;
t17 = t47 * t59 - t49 * t98 + t50 * t68;
t82 = -pkin(4) - t101;
t81 = -qJ(5) * t66 - pkin(3);
t13 = -t22 * t67 + t31 * t69;
t23 = t47 * t69 + t49 * t67;
t32 = -qJ(5) * t69 - t43;
t46 = t67 * t93 - t69 * t73;
t4 = -qJ(5) * t46 - t8;
t26 = -t66 * t77 + t68 * t90;
t76 = -qJ(5) * t26 + t13;
t41 = t66 * t94 + t68 * t70;
t75 = -qJ(5) * t41 + t23;
t14 = qJ(5) * t97 - t17;
t65 = t72 ^ 2;
t64 = t71 ^ 2;
t63 = t70 ^ 2;
t62 = t67 ^ 2;
t55 = t66 * t86;
t42 = pkin(3) * t95 - t55;
t40 = t66 * t70 - t68 * t94;
t36 = (-pkin(4) * t68 + t81) * t67;
t35 = t69 * t82 + t55;
t33 = -pkin(8) * t91 + t48;
t28 = (-t68 * t88 + t81) * t67;
t27 = pkin(5) * t60 - t32;
t25 = t66 * t90 + t68 * t77;
t24 = pkin(5) * t59 + t55 + (-qJ(6) + t82) * t69;
t16 = -t44 + t106;
t15 = t87 - t106;
t12 = pkin(4) * t40 + t75;
t11 = -pkin(5) * t40 - t14;
t10 = t40 * t88 + t75;
t9 = t49 * t95 + pkin(5) * t41 + (qJ(6) * t72 - t47 * t68) * t67 + t87;
t6 = pkin(4) * t25 + t76;
t5 = -t46 * pkin(4) - t7;
t3 = t25 * t88 + t76;
t2 = -pkin(5) * t25 - t4;
t1 = pkin(5) * t26 - t46 * t88 - t7;
t18 = [1, 0, 0, t64, t84, 0, 0, 0, pkin(1) * t103, pkin(1) * t104, t65 * t64, -0.2e1 * t64 * t92, t89 * t104, t70 * t84, t73 ^ 2, 0.2e1 * t100 * t64 - 0.2e1 * t33 * t73, 0.2e1 * pkin(8) * t64 * t72 + 0.2e1 * t34 * t73, 0.2e1 * t13 * t25 + 0.2e1 * t46 * t7, 0.2e1 * t13 * t26 - 0.2e1 * t46 * t8, -0.2e1 * t25 * t8 - 0.2e1 * t26 * t7, t13 ^ 2 + t7 ^ 2 + t8 ^ 2, 0.2e1 * t25 * t4 + 0.2e1 * t26 * t5, -0.2e1 * t25 * t6 + 0.2e1 * t46 * t5, -0.2e1 * t26 * t6 - 0.2e1 * t4 * t46, t4 ^ 2 + t5 ^ 2 + t6 ^ 2, 0.2e1 * t1 * t26 - 0.2e1 * t2 * t25, 0.2e1 * t2 * t46 - 0.2e1 * t26 * t3, -0.2e1 * t1 * t46 + 0.2e1 * t25 * t3, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t71, t73, 0, -t71 * pkin(8), -t73 * pkin(8), t70 * t90 (-t63 + t65) * t71, -t91, -t89, 0, -pkin(8) * t90 + (-pkin(2) * t71 + pkin(9) * t73) * t70, pkin(9) * t89 + (t100 - t102) * t71, t13 * t40 + t16 * t46 + t23 * t25 - t7 * t97, t13 * t41 - t17 * t46 + t23 * t26 + t8 * t97, -t16 * t26 - t17 * t25 - t40 * t8 - t41 * t7, t13 * t23 + t16 * t7 + t17 * t8, t14 * t25 + t15 * t26 + t4 * t40 + t41 * t5, -t12 * t25 + t15 * t46 - t40 * t6 - t5 * t97, -t12 * t26 - t14 * t46 + t4 * t97 - t41 * t6, t12 * t6 + t14 * t4 + t15 * t5, t1 * t41 - t11 * t25 - t2 * t40 + t26 * t9, -t10 * t26 + t11 * t46 - t2 * t97 - t3 * t41, t1 * t97 + t10 * t25 + t3 * t40 - t46 * t9, t1 * t9 + t10 * t3 + t11 * t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t63, 0.2e1 * t92, 0, 0, 0, 0.2e1 * t102, -0.2e1 * pkin(2) * t70, -0.2e1 * t16 * t97 + 0.2e1 * t23 * t40, 0.2e1 * t17 * t97 + 0.2e1 * t23 * t41, -0.2e1 * t16 * t41 - 0.2e1 * t17 * t40, t16 ^ 2 + t17 ^ 2 + t23 ^ 2, 0.2e1 * t14 * t40 + 0.2e1 * t15 * t41, -0.2e1 * t12 * t40 - 0.2e1 * t15 * t97, -0.2e1 * t12 * t41 + 0.2e1 * t14 * t97, t12 ^ 2 + t14 ^ 2 + t15 ^ 2, -0.2e1 * t11 * t40 + 0.2e1 * t41 * t9, -0.2e1 * t10 * t41 - 0.2e1 * t11 * t97, 0.2e1 * t10 * t40 + 0.2e1 * t9 * t97, t10 ^ 2 + t11 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t93, -t73, t33, -t34, t42 * t46 + t7 * t69 + (-pkin(3) * t25 - t13 * t68) * t67, -t43 * t46 - t69 * t8 + (-pkin(3) * t26 + t13 * t66) * t67, -t25 * t43 - t26 * t42 + (-t66 * t7 + t68 * t8) * t67, -t13 * t99 + t42 * t7 + t43 * t8, t25 * t32 + t26 * t35 + (-t4 * t68 + t5 * t66) * t67, -t25 * t36 + t35 * t46 + t5 * t69 + t6 * t60, -t26 * t36 - t32 * t46 - t4 * t69 - t59 * t6, t32 * t4 + t35 * t5 + t36 * t6, t24 * t26 - t25 * t27 + (t1 * t66 + t2 * t68) * t67, t2 * t69 - t26 * t28 + t27 * t46 - t3 * t59, -t1 * t69 - t24 * t46 + t25 * t28 - t3 * t60, t1 * t24 + t2 * t27 + t28 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t72, 0, -t70 * pkin(9), -t72 * pkin(9), t16 * t69 + (-pkin(3) * t40 - t23 * t68 - t42 * t72) * t67, -t17 * t69 + (-pkin(3) * t41 + t23 * t66 + t43 * t72) * t67, -t40 * t43 - t41 * t42 + (-t16 * t66 + t17 * t68) * t67, t16 * t42 + t17 * t43 - t23 * t99, t32 * t40 + t35 * t41 + (-t14 * t68 + t15 * t66) * t67, t15 * t69 - t36 * t40 + (t12 * t68 - t35 * t72) * t67, -t14 * t69 - t36 * t41 + (-t12 * t66 + t32 * t72) * t67, t12 * t36 + t14 * t32 + t15 * t35, t24 * t41 - t27 * t40 + (t11 * t68 + t66 * t9) * t67, t11 * t69 - t28 * t41 + (-t10 * t66 - t27 * t72) * t67, t28 * t40 - t9 * t69 + (-t10 * t68 + t24 * t72) * t67, t10 * t28 + t11 * t27 + t24 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t101 * t62 + 0.2e1 * t42 * t69, -0.2e1 * pkin(3) * t62 * t66 - 0.2e1 * t43 * t69 (-t42 * t66 + t43 * t68) * t105, pkin(3) ^ 2 * t62 + t42 ^ 2 + t43 ^ 2 (-t32 * t68 + t35 * t66) * t105, 0.2e1 * t35 * t69 + 0.2e1 * t36 * t60, -0.2e1 * t32 * t69 - 0.2e1 * t36 * t59, t32 ^ 2 + t35 ^ 2 + t36 ^ 2 (t24 * t66 + t27 * t68) * t105, 0.2e1 * t27 * t69 - 0.2e1 * t28 * t59, -0.2e1 * t24 * t69 - 0.2e1 * t28 * t60, t24 ^ 2 + t27 ^ 2 + t28 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t26, 0, t13, 0, -t25, -t26, t6, 0, -t26, t25, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t41, 0, t23, 0, -t40, -t41, t12, 0, -t41, t40, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t59, 0, -t99, 0, t60, -t59, t36, 0, -t59, -t60, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t46, 0, t5, t26, 0, -t46, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t97, 0, t15, t41, 0, t97, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t69, 0, t35, t59, 0, -t69, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t46, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t97, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t69, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t18;
