% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t67 = sin(qJ(4));
t68 = sin(qJ(3));
t70 = cos(qJ(4));
t71 = cos(qJ(3));
t44 = t67 * t68 - t70 * t71;
t69 = sin(qJ(2));
t39 = t44 * t69;
t98 = 0.2e1 * t39;
t57 = -t71 * pkin(3) - pkin(2);
t97 = 0.2e1 * t57;
t96 = -0.2e1 * t69;
t72 = cos(qJ(2));
t95 = 0.2e1 * t72;
t94 = pkin(8) + pkin(9);
t93 = pkin(2) * t71;
t92 = pkin(7) * t68;
t66 = cos(pkin(10));
t91 = t66 * pkin(4);
t90 = t67 * pkin(3);
t89 = t72 * pkin(3);
t47 = -t72 * pkin(2) - t69 * pkin(8) - pkin(1);
t43 = t71 * t47;
t84 = t71 * t69;
t27 = -pkin(9) * t84 + t43 + (-pkin(3) - t92) * t72;
t83 = t71 * t72;
t80 = pkin(7) * t83;
t29 = t80 + (-pkin(9) * t69 + t47) * t68;
t18 = t70 * t27 - t67 * t29;
t10 = -t72 * pkin(4) + t39 * qJ(5) + t18;
t85 = t70 * t29;
t19 = t67 * t27 + t85;
t45 = t67 * t71 + t70 * t68;
t38 = t45 * t69;
t17 = -t38 * qJ(5) + t19;
t65 = sin(pkin(10));
t4 = t65 * t10 + t66 * t17;
t88 = t68 * t69;
t87 = t68 * t71;
t86 = t68 * t72;
t60 = t70 * pkin(3);
t56 = t60 + pkin(4);
t41 = t65 * t56 + t66 * t90;
t59 = t69 * pkin(7);
t46 = pkin(3) * t88 + t59;
t82 = t69 * t95;
t78 = t94 * t71;
t79 = t94 * t68;
t31 = -t67 * t79 + t70 * t78;
t23 = -t44 * qJ(5) + t31;
t30 = -t67 * t78 - t70 * t79;
t74 = -t45 * qJ(5) + t30;
t14 = t65 * t23 - t66 * t74;
t16 = t66 * t23 + t65 * t74;
t81 = t14 ^ 2 + t16 ^ 2;
t20 = t66 * t38 - t65 * t39;
t21 = -t65 * t38 - t66 * t39;
t77 = t14 * t21 - t16 * t20;
t3 = t66 * t10 - t65 * t17;
t76 = -t66 * t56 + t65 * t90;
t28 = t38 * pkin(4) + t46;
t35 = t44 * pkin(4) + t57;
t25 = t66 * t44 + t65 * t45;
t26 = -t65 * t44 + t66 * t45;
t75 = 0.2e1 * t14 * t26 - 0.2e1 * t16 * t25;
t64 = t72 ^ 2;
t63 = t71 ^ 2;
t62 = t69 ^ 2;
t61 = t68 ^ 2;
t58 = t65 * pkin(4);
t53 = pkin(5) + t91;
t52 = t58 + qJ(6);
t37 = -pkin(5) + t76;
t36 = qJ(6) + t41;
t34 = t68 * t47 + t80;
t33 = -pkin(7) * t86 + t43;
t12 = t25 * pkin(5) - t26 * qJ(6) + t35;
t6 = t20 * pkin(5) - t21 * qJ(6) + t28;
t2 = t72 * pkin(5) - t3;
t1 = -t72 * qJ(6) + t4;
t5 = [1, 0, 0, t62, t82, 0, 0, 0, pkin(1) * t95, pkin(1) * t96, t63 * t62, -0.2e1 * t62 * t87, t83 * t96, t68 * t82, t64, -0.2e1 * t33 * t72 + 0.2e1 * t62 * t92, 0.2e1 * t62 * pkin(7) * t71 + 0.2e1 * t34 * t72, t39 ^ 2, t38 * t98, t72 * t98, t38 * t95, t64, -0.2e1 * t18 * t72 + 0.2e1 * t46 * t38, 0.2e1 * t19 * t72 - 0.2e1 * t46 * t39, -0.2e1 * t4 * t20 - 0.2e1 * t3 * t21, t28 ^ 2 + t3 ^ 2 + t4 ^ 2, 0.2e1 * t2 * t72 + 0.2e1 * t6 * t20, -0.2e1 * t1 * t20 + 0.2e1 * t2 * t21, -0.2e1 * t1 * t72 - 0.2e1 * t6 * t21, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t69, t72, 0, -t59, -t72 * pkin(7), t68 * t84 (-t61 + t63) * t69, -t86, -t83, 0, -pkin(7) * t84 + (-pkin(2) * t69 + pkin(8) * t72) * t68, pkin(8) * t83 + (t92 - t93) * t69, -t39 * t45, -t45 * t38 + t39 * t44, -t45 * t72, t44 * t72, 0, -t30 * t72 + t57 * t38 + t46 * t44, t31 * t72 - t57 * t39 + t46 * t45, -t4 * t25 - t3 * t26 + t77, -t3 * t14 + t4 * t16 + t28 * t35, t12 * t20 + t14 * t72 + t6 * t25, -t1 * t25 + t2 * t26 + t77, -t12 * t21 - t16 * t72 - t6 * t26, t1 * t16 + t6 * t12 + t2 * t14; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t61, 0.2e1 * t87, 0, 0, 0, 0.2e1 * t93, -0.2e1 * pkin(2) * t68, t45 ^ 2, -0.2e1 * t45 * t44, 0, 0, 0, t44 * t97, t45 * t97, t75, t35 ^ 2 + t81, 0.2e1 * t12 * t25, t75, -0.2e1 * t12 * t26, t12 ^ 2 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t88, -t72, t33, -t34, 0, 0, -t39, -t38, -t72, -t70 * t89 + t18, -t85 + (-t27 + t89) * t67, -t41 * t20 + t21 * t76, -t3 * t76 + t4 * t41 (-pkin(5) + t37) * t72 + t3, -t36 * t20 + t37 * t21 (-qJ(6) - t36) * t72 + t4, t1 * t36 + t2 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t71, 0, -t68 * pkin(8), -t71 * pkin(8), 0, 0, t45, -t44, 0, t30, -t31, -t41 * t25 + t26 * t76, t14 * t76 + t16 * t41, -t14, -t36 * t25 + t37 * t26, t16, t14 * t37 + t16 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t60, -0.2e1 * t90, 0, t41 ^ 2 + t76 ^ 2, -0.2e1 * t37, 0, 0.2e1 * t36, t36 ^ 2 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t38, -t72, t18, -t19 (-t20 * t65 - t21 * t66) * pkin(4) (t3 * t66 + t4 * t65) * pkin(4) (-pkin(5) - t53) * t72 + t3, -t52 * t20 - t53 * t21 (-qJ(6) - t52) * t72 + t4, t1 * t52 - t2 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t44, 0, t30, -t31 (-t25 * t65 - t26 * t66) * pkin(4) (-t14 * t66 + t16 * t65) * pkin(4), -t14, -t52 * t25 - t53 * t26, t16, -t14 * t53 + t16 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t60, -t90, 0 (t41 * t65 - t66 * t76) * pkin(4), 0.2e1 * pkin(5) - t76 + t91, 0, t58 + 0.2e1 * qJ(6) + t41, t36 * t52 - t37 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t65 ^ 2 + t66 ^ 2) * pkin(4) ^ 2, 0.2e1 * t53, 0, 0.2e1 * t52, t52 ^ 2 + t53 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t20, 0, -t21, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t25, 0, -t26, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t21, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t5;
