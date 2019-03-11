% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t67 = sin(qJ(4));
t68 = sin(qJ(3));
t71 = cos(qJ(4));
t72 = cos(qJ(3));
t41 = t67 * t68 - t71 * t72;
t74 = -pkin(4) - pkin(5);
t42 = t67 * t72 + t71 * t68;
t57 = -t72 * pkin(3) - pkin(2);
t78 = t42 * qJ(5) - t57;
t16 = t74 * t41 + t78;
t103 = 0.2e1 * t16;
t102 = -0.2e1 * t42;
t101 = 0.2e1 * t57;
t69 = sin(qJ(2));
t100 = -0.2e1 * t69;
t73 = cos(qJ(2));
t99 = -0.2e1 * t73;
t98 = 0.2e1 * t73;
t97 = -pkin(9) - pkin(8);
t61 = t73 * pkin(4);
t47 = -t73 * pkin(2) - t69 * pkin(8) - pkin(1);
t40 = t72 * t47;
t87 = t72 * t69;
t95 = pkin(7) * t68;
t26 = -pkin(9) * t87 + t40 + (-pkin(3) - t95) * t73;
t86 = t72 * t73;
t81 = pkin(7) * t86;
t29 = t81 + (-pkin(9) * t69 + t47) * t68;
t85 = -t71 * t26 + t67 * t29;
t11 = t61 + t85;
t90 = t68 * t69;
t38 = -t67 * t90 + t71 * t87;
t4 = t73 * pkin(5) - t38 * pkin(10) + t11;
t66 = sin(qJ(6));
t14 = t67 * t26 + t71 * t29;
t84 = t73 * qJ(5);
t10 = -t84 + t14;
t37 = t42 * t69;
t7 = t37 * pkin(10) + t10;
t70 = cos(qJ(6));
t2 = t66 * t4 + t70 * t7;
t96 = pkin(2) * t72;
t94 = t71 * pkin(3);
t93 = t73 * pkin(3);
t49 = t97 * t72;
t80 = t97 * t68;
t30 = -t67 * t49 - t71 * t80;
t92 = t30 * t73;
t31 = -t71 * t49 + t67 * t80;
t91 = t31 * t73;
t89 = t68 * t72;
t88 = t68 * t73;
t60 = t69 * pkin(7);
t44 = pkin(3) * t90 + t60;
t59 = t67 * pkin(3);
t53 = t59 + qJ(5);
t83 = qJ(5) + t53;
t82 = t69 * t98;
t55 = pkin(4) + t94;
t1 = t70 * t4 - t66 * t7;
t79 = t38 * qJ(5) - t44;
t77 = -t42 * pkin(10) + t30;
t76 = 0.2e1 * pkin(4);
t75 = 0.2e1 * qJ(5);
t65 = t73 ^ 2;
t64 = t72 ^ 2;
t63 = t69 ^ 2;
t62 = t68 ^ 2;
t58 = t70 * t74;
t52 = -pkin(5) - t55;
t48 = t70 * t52;
t46 = t70 * qJ(5) + t66 * t74;
t45 = t66 * qJ(5) - t58;
t36 = t66 * t52 + t70 * t53;
t35 = t66 * t53 - t48;
t34 = t68 * t47 + t81;
t33 = -pkin(7) * t88 + t40;
t25 = t66 * t41 + t70 * t42;
t24 = -t70 * t41 + t66 * t42;
t23 = t41 * pkin(4) - t78;
t20 = t41 * pkin(10) + t31;
t18 = t66 * t37 + t70 * t38;
t17 = -t70 * t37 + t66 * t38;
t15 = t37 * pkin(4) - t79;
t12 = t74 * t37 + t79;
t9 = t70 * t20 + t66 * t77;
t8 = t66 * t20 - t70 * t77;
t3 = [1, 0, 0, t63, t82, 0, 0, 0, pkin(1) * t98, pkin(1) * t100, t64 * t63, -0.2e1 * t63 * t89, t86 * t100, t68 * t82, t65, -0.2e1 * t33 * t73 + 0.2e1 * t63 * t95, 0.2e1 * t63 * pkin(7) * t72 + 0.2e1 * t34 * t73, t38 ^ 2, -0.2e1 * t38 * t37, t38 * t99, t37 * t98, t65, 0.2e1 * t44 * t37 + 0.2e1 * t73 * t85, 0.2e1 * t14 * t73 + 0.2e1 * t44 * t38, 0.2e1 * t11 * t73 + 0.2e1 * t15 * t37, -0.2e1 * t10 * t37 + 0.2e1 * t11 * t38, -0.2e1 * t10 * t73 - 0.2e1 * t15 * t38, t10 ^ 2 + t11 ^ 2 + t15 ^ 2, t18 ^ 2, -0.2e1 * t18 * t17, t18 * t98, t17 * t99, t65, 0.2e1 * t1 * t73 + 0.2e1 * t12 * t17, 0.2e1 * t12 * t18 - 0.2e1 * t2 * t73; 0, 0, 0, 0, 0, t69, t73, 0, -t60, -t73 * pkin(7), t68 * t87 (-t62 + t64) * t69, -t88, -t86, 0, -pkin(7) * t87 + (-pkin(2) * t69 + pkin(8) * t73) * t68, pkin(8) * t86 + (t95 - t96) * t69, t38 * t42, -t42 * t37 - t38 * t41, -t42 * t73, t41 * t73, 0, t57 * t37 + t44 * t41 + t92, t57 * t38 + t44 * t42 + t91, t15 * t41 + t23 * t37 + t92, -t10 * t41 + t11 * t42 + t30 * t38 - t31 * t37, -t15 * t42 - t23 * t38 - t91, t10 * t31 + t11 * t30 + t15 * t23, t18 * t25, -t25 * t17 - t18 * t24, t25 * t73, -t24 * t73, 0, t12 * t24 + t16 * t17 - t8 * t73, t12 * t25 + t16 * t18 - t9 * t73; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t62, 0.2e1 * t89, 0, 0, 0, 0.2e1 * t96, -0.2e1 * pkin(2) * t68, t42 ^ 2, t41 * t102, 0, 0, 0, t41 * t101, t42 * t101, 0.2e1 * t23 * t41, 0.2e1 * t30 * t42 - 0.2e1 * t31 * t41, t23 * t102, t23 ^ 2 + t30 ^ 2 + t31 ^ 2, t25 ^ 2, -0.2e1 * t25 * t24, 0, 0, 0, t24 * t103, t25 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t90, -t73, t33, -t34, 0, 0, t38, -t37, -t73, -t71 * t93 - t85, t67 * t93 - t14, -t55 * t73 - t11, -t53 * t37 - t55 * t38, -t83 * t73 + t14, t10 * t53 - t11 * t55, 0, 0, -t18, t17, -t73, -t35 * t73 - t1, -t36 * t73 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t72, 0, -t68 * pkin(8), -t72 * pkin(8), 0, 0, t42, -t41, 0, -t30, -t31, -t30, -t53 * t41 - t55 * t42, t31, -t30 * t55 + t31 * t53, 0, 0, -t25, t24, 0, t8, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t94, -0.2e1 * t59, 0.2e1 * t55, 0, 0.2e1 * t53, t53 ^ 2 + t55 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t35, 0.2e1 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, -t73, -t85, -t14, -0.2e1 * t61 - t85, -pkin(4) * t38 - t37 * qJ(5), -0.2e1 * t84 + t14, -t11 * pkin(4) + t10 * qJ(5), 0, 0, -t18, t17, -t73, -t45 * t73 - t1, -t46 * t73 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t41, 0, -t30, -t31, -t30, -pkin(4) * t42 - t41 * qJ(5), t31, -t30 * pkin(4) + t31 * qJ(5), 0, 0, -t25, t24, 0, t8, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t94, -t59, t76 + t94, 0, t75 + t59, t55 * pkin(4) + t53 * qJ(5), 0, 0, 0, 0, 1, t83 * t66 - t48 - t58, t83 * t70 + (t52 + t74) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t76, 0, t75, pkin(4) ^ 2 + qJ(5) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t45, 0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t38, 0, t11, 0, 0, 0, 0, 0, t70 * t73, -t66 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t55, 0, 0, 0, 0, 0, -t70, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), 0, 0, 0, 0, 0, -t70, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, t73, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, -t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t35, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t45, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
