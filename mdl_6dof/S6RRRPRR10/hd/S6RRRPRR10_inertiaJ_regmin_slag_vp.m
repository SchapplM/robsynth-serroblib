% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t73 = sin(qJ(3));
t77 = cos(qJ(3));
t112 = -pkin(3) * t77 - qJ(4) * t73;
t49 = -pkin(2) + t112;
t38 = pkin(4) * t77 - t49;
t72 = sin(qJ(5));
t76 = cos(qJ(5));
t42 = t72 * t73 + t76 * t77;
t24 = pkin(5) * t42 + t38;
t111 = 0.2e1 * t24;
t110 = 0.2e1 * t38;
t109 = -0.2e1 * t49;
t74 = sin(qJ(2));
t108 = -0.2e1 * t74;
t78 = cos(qJ(2));
t107 = -0.2e1 * t78;
t106 = 0.2e1 * t78;
t105 = -pkin(3) - pkin(4);
t104 = pkin(2) * t73;
t103 = pkin(2) * t77;
t102 = pkin(3) * t73;
t101 = pkin(7) * t73;
t100 = pkin(7) * t77;
t71 = sin(qJ(6));
t99 = t71 * pkin(5);
t75 = cos(qJ(6));
t98 = t75 * pkin(5);
t59 = t77 * t74;
t91 = t76 * t73;
t34 = t59 * t72 - t74 * t91;
t66 = t78 * pkin(3);
t50 = -pkin(2) * t78 - pkin(8) * t74 - pkin(1);
t93 = t73 * t78;
t89 = pkin(7) * t93 - t50 * t77;
t29 = t66 + t89;
t17 = pkin(4) * t78 - pkin(9) * t59 + t29;
t90 = t77 * t78;
t32 = pkin(7) * t90 + t50 * t73;
t86 = t78 * qJ(4);
t28 = -t86 + t32;
t95 = t73 * t74;
t21 = pkin(9) * t95 + t28;
t9 = t17 * t72 + t21 * t76;
t5 = -pkin(10) * t34 + t9;
t97 = t75 * t5;
t96 = t78 * pkin(5);
t94 = t73 * t77;
t48 = qJ(4) * t76 + t105 * t72;
t92 = t75 * t48;
t67 = t73 ^ 2;
t69 = t77 ^ 2;
t88 = t67 + t69;
t87 = t77 * qJ(4);
t85 = t74 * t106;
t35 = t42 * t74;
t83 = -t17 * t76 + t21 * t72;
t4 = -pkin(10) * t35 - t83 + t96;
t84 = -t4 * t75 + t5 * t71;
t47 = qJ(4) * t72 - t105 * t76;
t46 = -pkin(5) - t47;
t22 = -t46 * t75 + t48 * t71;
t62 = t73 * pkin(8);
t51 = -pkin(9) * t73 + t62;
t63 = t77 * pkin(8);
t52 = -pkin(9) * t77 + t63;
t26 = -t51 * t76 + t52 * t72;
t2 = t4 * t71 + t97;
t82 = t87 - t102;
t81 = t28 * t77 + t29 * t73;
t27 = t51 * t72 + t52 * t76;
t44 = -t72 * t77 + t91;
t80 = -pkin(10) * t44 - t26;
t54 = t74 * t87;
t25 = t54 + (t105 * t73 - pkin(7)) * t74;
t70 = t78 ^ 2;
t68 = t74 ^ 2;
t56 = pkin(8) * t93;
t43 = t71 * t76 + t72 * t75;
t41 = t71 * t72 - t75 * t76;
t33 = -t54 + (pkin(7) + t102) * t74;
t23 = t46 * t71 + t92;
t20 = -t42 * t71 + t44 * t75;
t19 = t42 * t75 + t44 * t71;
t14 = -pkin(10) * t42 + t27;
t12 = -t34 * t71 + t35 * t75;
t11 = t34 * t75 + t35 * t71;
t10 = t34 * pkin(5) + t25;
t7 = t14 * t75 + t71 * t80;
t6 = t14 * t71 - t75 * t80;
t1 = [1, 0, 0, t68, t85, 0, 0, 0, pkin(1) * t106, pkin(1) * t108, t69 * t68, -0.2e1 * t68 * t94, t90 * t108, t73 * t85, t70, 0.2e1 * t101 * t68 + 0.2e1 * t78 * t89, 0.2e1 * t100 * t68 + 0.2e1 * t32 * t78, 0.2e1 * t29 * t78 + 0.2e1 * t33 * t95, 0.2e1 * (-t28 * t73 + t29 * t77) * t74, -0.2e1 * t28 * t78 - 0.2e1 * t33 * t59, t28 ^ 2 + t29 ^ 2 + t33 ^ 2, t35 ^ 2, -0.2e1 * t35 * t34, t35 * t106, t34 * t107, t70, 0.2e1 * t25 * t34 - 0.2e1 * t78 * t83, 0.2e1 * t25 * t35 - 0.2e1 * t78 * t9, t12 ^ 2, -0.2e1 * t12 * t11, t12 * t106, t11 * t107, t70, 0.2e1 * t10 * t11 - 0.2e1 * t78 * t84, 0.2e1 * t10 * t12 - 0.2e1 * t2 * t78; 0, 0, 0, 0, 0, t74, t78, 0, -t74 * pkin(7), -t78 * pkin(7), t73 * t59 (-t67 + t69) * t74, -t93, -t90, 0, t56 + (-t100 - t104) * t74, pkin(8) * t90 + (t101 - t103) * t74, -t33 * t77 + t49 * t95 + t56, t81, -t33 * t73 + (-pkin(8) * t78 - t49 * t74) * t77, pkin(8) * t81 + t33 * t49, t35 * t44, -t34 * t44 - t35 * t42, t44 * t78, -t42 * t78, 0, t25 * t42 - t26 * t78 + t34 * t38, t25 * t44 - t27 * t78 + t35 * t38, t12 * t20, -t11 * t20 - t12 * t19, t20 * t78, -t19 * t78, 0, t10 * t19 + t11 * t24 - t6 * t78, t10 * t20 + t12 * t24 - t7 * t78; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t67, 0.2e1 * t94, 0, 0, 0, 0.2e1 * t103, -0.2e1 * t104, t77 * t109, 0.2e1 * t88 * pkin(8), t73 * t109, pkin(8) ^ 2 * t88 + t49 ^ 2, t44 ^ 2, -0.2e1 * t44 * t42, 0, 0, 0, t42 * t110, t44 * t110, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t111, t20 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t95, -t78, -t89, -t32, -0.2e1 * t66 - t89, t112 * t74, -0.2e1 * t86 + t32, -pkin(3) * t29 + qJ(4) * t28, 0, 0, -t35, t34, -t78, -t47 * t78 + t83, -t48 * t78 + t9, 0, 0, -t12, t11, -t78, -t22 * t78 + t84, -t23 * t78 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t77, 0, -t62, -t63, -t62, t82, t63, t82 * pkin(8), 0, 0, -t44, t42, 0, t26, t27, 0, 0, -t20, t19, 0, t6, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t47, 0.2e1 * t48, 0, 0, 0, 0, 1, 0.2e1 * t22, 0.2e1 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t59, 0, t29, 0, 0, 0, 0, 0, t76 * t78, -t72 * t78, 0, 0, 0, 0, 0, -t41 * t78, -t43 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, -t76, t72, 0, 0, 0, 0, 0, t41, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34, t78, -t83, -t9, 0, 0, t12, -t11, t78, t75 * t96 - t84, -t97 + (-t4 - t96) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t42, 0, -t26, -t27, 0, 0, t20, -t19, 0, -t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t47, -t48, 0, 0, 0, 0, -1, -t22 - t98, -t92 + (pkin(5) - t46) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t72, 0, 0, 0, 0, 0, -t41, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t98, -0.2e1 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, t78, -t84, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, -t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t98, -t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t1;
