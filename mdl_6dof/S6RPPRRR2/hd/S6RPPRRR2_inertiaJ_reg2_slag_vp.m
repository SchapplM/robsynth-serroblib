% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t64 = sin(qJ(6));
t67 = cos(qJ(6));
t62 = cos(pkin(11));
t66 = sin(qJ(4));
t60 = sin(pkin(11));
t89 = cos(qJ(4));
t75 = t89 * t60;
t41 = t66 * t62 + t75;
t68 = cos(qJ(5));
t81 = t68 * t41;
t65 = sin(qJ(5));
t85 = t65 * t41;
t19 = -t64 * t85 + t67 * t81;
t44 = t64 * t65 - t67 * t68;
t10 = t19 * t44;
t46 = t64 * t68 + t67 * t65;
t17 = t46 * t41;
t87 = t46 * t17;
t106 = t10 - t87;
t105 = -0.2e1 * t41;
t61 = sin(pkin(10));
t95 = t61 * pkin(1);
t51 = qJ(3) + t95;
t90 = pkin(7) + t51;
t35 = t90 * t62;
t20 = t66 * t35 + t90 * t75;
t104 = t20 ^ 2;
t83 = t66 * t60;
t39 = -t89 * t62 + t83;
t103 = t39 ^ 2;
t102 = 0.2e1 * t39;
t63 = cos(pkin(10));
t94 = t63 * pkin(1);
t53 = -pkin(2) - t94;
t47 = -t62 * pkin(3) + t53;
t101 = 0.2e1 * t47;
t55 = -t68 * pkin(5) - pkin(4);
t100 = 0.2e1 * t55;
t99 = 0.2e1 * t60;
t98 = -pkin(9) - pkin(8);
t97 = t39 * pkin(4);
t96 = t39 * pkin(5);
t93 = t64 * pkin(5);
t92 = t67 * pkin(5);
t14 = -t41 * pkin(8) + t47 + t97;
t22 = t89 * t35 - t90 * t83;
t82 = t68 * t22;
t5 = t82 + (-pkin(9) * t41 + t14) * t65;
t91 = t67 * t5;
t88 = t20 * t39;
t24 = t39 * t44;
t25 = t46 * t39;
t58 = t65 ^ 2;
t86 = t58 * t41;
t31 = t65 * t39;
t84 = t65 * t68;
t56 = t60 ^ 2;
t57 = t62 ^ 2;
t80 = t56 + t57;
t59 = t68 ^ 2;
t79 = t58 + t59;
t78 = t39 * t105;
t77 = t65 * t81;
t6 = t68 * t14 - t65 * t22;
t4 = -pkin(9) * t81 + t6 + t96;
t1 = t67 * t4 - t64 * t5;
t76 = t79 * pkin(8);
t74 = -pkin(4) * t41 - pkin(8) * t39;
t7 = t65 * t14 + t82;
t73 = t6 * t68 + t7 * t65;
t72 = -t6 * t65 + t7 * t68;
t49 = t98 * t68;
t48 = t98 * t65;
t43 = t46 ^ 2;
t42 = t44 ^ 2;
t38 = t41 ^ 2;
t34 = t68 * t39;
t33 = t59 * t41;
t32 = t59 * t38;
t30 = t58 * t38;
t27 = t64 * t48 - t67 * t49;
t26 = t67 * t48 + t64 * t49;
t23 = -t33 - t86;
t16 = t19 ^ 2;
t15 = t17 ^ 2;
t11 = t19 * t46;
t9 = t17 * t44;
t8 = pkin(5) * t85 + t20;
t2 = t64 * t4 + t91;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t94, -0.2e1 * t95, 0 (t61 ^ 2 + t63 ^ 2) * pkin(1) ^ 2, t56, t62 * t99, 0, t57, 0, 0, -0.2e1 * t53 * t62, t53 * t99, 0.2e1 * t80 * t51, t80 * t51 ^ 2 + t53 ^ 2, t38, t78, 0, t103, 0, 0, t39 * t101, t41 * t101, 0.2e1 * t20 * t41 - 0.2e1 * t22 * t39, t22 ^ 2 + t47 ^ 2 + t104, t32, -0.2e1 * t38 * t84, t81 * t102, t30, t65 * t78, t103, 0.2e1 * t20 * t85 + 0.2e1 * t6 * t39, 0.2e1 * t20 * t81 - 0.2e1 * t7 * t39, t73 * t105, t6 ^ 2 + t7 ^ 2 + t104, t16, -0.2e1 * t19 * t17, t19 * t102, t15, -t17 * t102, t103, 0.2e1 * t1 * t39 + 0.2e1 * t8 * t17, 0.2e1 * t8 * t19 - 0.2e1 * t2 * t39, -0.2e1 * t1 * t19 - 0.2e1 * t2 * t17, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t41 + t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t41 + t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t17 + t2 * t19 + t8 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 + t103, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 + t30 + t103, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 + t15 + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t60, 0, t53, 0, 0, 0, 0, 0, 0, t39, t41, 0, t47, 0, 0, 0, 0, 0, 0, t34, -t31, t23, t73, 0, 0, 0, 0, 0, 0, -t24, -t25, t106, -t1 * t44 + t2 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 + t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, -t39, 0, -t20, -t22, 0, 0, t77, t33 - t86, t31, -t77, t34, 0, -t20 * t68 + t74 * t65, t20 * t65 + t74 * t68, t72, -t20 * pkin(4) + t72 * pkin(8), t11, -t10 - t87, t25, t9, -t24, 0, t55 * t17 + t26 * t39 + t8 * t44, t55 * t19 - t27 * t39 + t8 * t46, -t1 * t46 - t27 * t17 - t26 * t19 - t2 * t44, t1 * t26 + t2 * t27 + t8 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t41, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t31, -t23, t41 * t76 - t97, 0, 0, 0, 0, 0, 0, t24, t25, -t106, -t17 * t26 + t19 * t27 + t39 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t26 + t46 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t58, 0.2e1 * t84, 0, t59, 0, 0, 0.2e1 * pkin(4) * t68, -0.2e1 * pkin(4) * t65, 0.2e1 * t76, t79 * pkin(8) ^ 2 + pkin(4) ^ 2, t43, -0.2e1 * t46 * t44, 0, t42, 0, 0, t44 * t100, t46 * t100, -0.2e1 * t26 * t46 - 0.2e1 * t27 * t44, t26 ^ 2 + t27 ^ 2 + t55 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, -t85, t39, t6, -t7, 0, 0, 0, 0, t19, 0, -t17, t39, t39 * t92 + t1, -t91 + (-t4 - t96) * t64 (-t17 * t64 - t19 * t67) * pkin(5) (t1 * t67 + t2 * t64) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t81, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t19, 0 (-t17 * t67 + t19 * t64) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t65, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t46, 0 (-t44 * t67 + t46 * t64) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t68, 0, -t65 * pkin(8), -t68 * pkin(8), 0, 0, 0, 0, t46, 0, -t44, 0, t26, -t27 (-t44 * t64 - t46 * t67) * pkin(5) (t26 * t67 + t27 * t64) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t92, -0.2e1 * t93, 0 (t64 ^ 2 + t67 ^ 2) * pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, t39, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, -t44, 0, t26, -t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t92, -t93, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t3;
