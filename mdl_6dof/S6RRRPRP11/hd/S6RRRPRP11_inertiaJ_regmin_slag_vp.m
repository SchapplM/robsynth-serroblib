% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t61 = cos(pkin(6));
t63 = sin(qJ(3));
t66 = cos(qJ(3));
t60 = sin(pkin(6));
t64 = sin(qJ(2));
t98 = t60 * t64;
t35 = t61 * t63 + t66 * t98;
t108 = -0.2e1 * t35;
t107 = -0.2e1 * t63;
t106 = 0.2e1 * t66;
t105 = 2 * qJ(4);
t104 = pkin(3) + pkin(10);
t103 = pkin(1) * t64;
t67 = cos(qJ(2));
t102 = pkin(1) * t67;
t65 = cos(qJ(5));
t101 = t65 * pkin(5);
t62 = sin(qJ(5));
t100 = t35 * t62;
t30 = t35 * t63;
t55 = t60 ^ 2;
t99 = t55 * t67;
t97 = t60 * t67;
t96 = t61 * t64;
t34 = -t61 * t66 + t63 * t98;
t20 = t34 * t65 + t62 * t97;
t95 = t62 * t20;
t52 = t63 * pkin(9);
t43 = t63 * pkin(4) + t52;
t94 = t62 * t43;
t93 = t62 * t63;
t92 = t62 * t66;
t91 = t62 * t104;
t90 = t63 * t66;
t21 = -t34 * t62 + t65 * t97;
t89 = t65 * t21;
t88 = t65 * t62;
t87 = t65 * t66;
t51 = t65 * t104;
t82 = pkin(8) * t97;
t27 = t82 + (pkin(9) + t103) * t61;
t28 = (-pkin(2) * t67 - pkin(9) * t64 - pkin(1)) * t60;
t14 = t66 * t27 + t63 * t28;
t53 = t66 * pkin(9);
t44 = t66 * pkin(4) + t53;
t57 = t63 ^ 2;
t59 = t66 ^ 2;
t86 = t57 + t59;
t85 = qJ(4) * t66;
t84 = 0.2e1 * t97;
t83 = -0.2e1 * t90;
t81 = t63 * t97;
t80 = t66 * t97;
t79 = qJ(4) * t97;
t13 = -t63 * t27 + t66 * t28;
t47 = pkin(3) * t97;
t12 = -t13 + t47;
t7 = t35 * pkin(4) + pkin(10) * t97 + t12;
t46 = pkin(8) * t98;
t26 = t46 + (-pkin(2) - t102) * t61;
t70 = -t35 * qJ(4) + t26;
t8 = t104 * t34 + t70;
t3 = -t62 * t8 + t65 * t7;
t78 = -t63 * qJ(4) - pkin(2);
t77 = pkin(9) * t81;
t76 = pkin(9) * t80;
t38 = -t104 * t66 + t78;
t75 = qJ(6) * t66 - t38;
t1 = t35 * pkin(5) + t21 * qJ(6) + t3;
t4 = t62 * t7 + t65 * t8;
t2 = t20 * qJ(6) + t4;
t74 = t1 * t65 + t2 * t62;
t73 = -pkin(3) * t63 + t85;
t11 = t79 - t14;
t72 = -t11 * t66 + t12 * t63;
t71 = -t104 * t63 + t85;
t9 = -t34 * pkin(4) - t11;
t58 = t65 ^ 2;
t56 = t62 ^ 2;
t50 = t65 * t63;
t49 = t62 * pkin(5) + qJ(4);
t45 = t56 + t58;
t42 = -t66 * pkin(3) + t78;
t41 = -t65 * qJ(6) - t51;
t40 = (-qJ(6) - t104) * t62;
t39 = t65 * t43;
t37 = pkin(1) * t96 + t82;
t36 = t61 * t102 - t46;
t33 = pkin(5) * t87 + t44;
t32 = t35 ^ 2;
t31 = t35 * t65;
t19 = t40 * t62 + t41 * t65;
t18 = t65 * t38 + t94;
t17 = -t62 * t38 + t39;
t16 = -t75 * t65 + t94;
t15 = t63 * pkin(5) + t75 * t62 + t39;
t10 = t34 * pkin(3) + t70;
t5 = -t20 * pkin(5) + t9;
t6 = [1, 0, 0, t55 * t64 ^ 2, 0.2e1 * t64 * t99, 0.2e1 * t60 * t96, t61 * t84, t61 ^ 2, 0.2e1 * pkin(1) * t99 + 0.2e1 * t36 * t61, -0.2e1 * t55 * t103 - 0.2e1 * t37 * t61, t32, t34 * t108, t97 * t108, t34 * t84, t55 * t67 ^ 2, -0.2e1 * t13 * t97 + 0.2e1 * t26 * t34, 0.2e1 * t14 * t97 + 0.2e1 * t26 * t35, 0.2e1 * t11 * t34 + 0.2e1 * t12 * t35, -0.2e1 * t10 * t34 - 0.2e1 * t12 * t97, -0.2e1 * t10 * t35 + 0.2e1 * t11 * t97, t10 ^ 2 + t11 ^ 2 + t12 ^ 2, t21 ^ 2, -0.2e1 * t21 * t20, t21 * t108, 0.2e1 * t20 * t35, t32, -0.2e1 * t9 * t20 + 0.2e1 * t3 * t35, -0.2e1 * t9 * t21 - 0.2e1 * t4 * t35, 0.2e1 * t1 * t21 + 0.2e1 * t2 * t20, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t98, t97, t61, t36, -t37, t30, -t63 * t34 + t35 * t66, -t81, -t80, 0, -pkin(2) * t34 - t26 * t66 + t77, -pkin(2) * t35 + t26 * t63 + t76 (-t34 * t66 + t30) * pkin(9) + t72, t10 * t66 - t42 * t34 - t77, -t10 * t63 - t42 * t35 - t76, t72 * pkin(9) + t10 * t42, t21 * t92 (t89 - t95) * t66, -t21 * t63 - t35 * t92, t20 * t63 - t35 * t87, t30, t17 * t35 - t44 * t20 + t3 * t63 + t9 * t87, -t18 * t35 - t44 * t21 - t4 * t63 - t9 * t92, t15 * t21 + t16 * t20 + (t1 * t62 - t2 * t65) * t66, t1 * t15 + t2 * t16 + t5 * t33; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t57, 0.2e1 * t90, 0, 0, 0, pkin(2) * t106, pkin(2) * t107, 0.2e1 * t86 * pkin(9), t42 * t106, t42 * t107, t86 * pkin(9) ^ 2 + t42 ^ 2, t56 * t59, 0.2e1 * t59 * t88, t62 * t83, t65 * t83, t57, 0.2e1 * t17 * t63 + 0.2e1 * t44 * t87, -0.2e1 * t18 * t63 - 0.2e1 * t44 * t92 (t15 * t62 - t16 * t65) * t106, t15 ^ 2 + t16 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34, -t97, t13, -t14, -t35 * pkin(3) - qJ(4) * t34, -t13 + 0.2e1 * t47, -0.2e1 * t79 + t14, -t12 * pkin(3) - t11 * qJ(4), -t89, t65 * t20 + t21 * t62, t31, -t100, 0, -qJ(4) * t20 - t35 * t51 + t9 * t62, -qJ(4) * t21 + t35 * t91 + t9 * t65, t40 * t20 + t41 * t21 - t74, t1 * t41 + t2 * t40 + t5 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t66, 0, -t52, -t53, t73, t52, t53, t73 * pkin(9), -t62 * t87 (t56 - t58) * t66, t50, -t93, 0, t44 * t62 + t71 * t65, t44 * t65 - t71 * t62 (-t40 * t66 - t15) * t65 + (t41 * t66 - t16) * t62, t15 * t41 + t16 * t40 + t33 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t105, pkin(3) ^ 2 + (qJ(4) ^ 2) t58, -0.2e1 * t88, 0, 0, 0, t62 * t105, t65 * t105, -0.2e1 * t19, t40 ^ 2 + t41 ^ 2 + t49 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t97, 0, t12, 0, 0, 0, 0, 0, t31, -t100, t89 + t95, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, t52, 0, 0, 0, 0, 0, t50, -t93, 0, t15 * t65 + t16 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, -t45, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t20, t35, t3, -t4, pkin(5) * t21, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, -t87, t63, t17, -t18, pkin(5) * t92, t15 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t62, 0, -t51, t91, -t101, t41 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t62, 0, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
