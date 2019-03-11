% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t72 = sin(pkin(7));
t83 = cos(qJ(3));
t102 = t72 * t83;
t79 = sin(qJ(3));
t103 = t72 * t79;
t75 = cos(pkin(7));
t78 = sin(qJ(4));
t82 = cos(qJ(4));
t46 = t82 * t103 + t78 * t75;
t71 = sin(pkin(13));
t74 = cos(pkin(13));
t30 = t74 * t102 + t71 * t46;
t31 = -t71 * t102 + t74 * t46;
t77 = sin(qJ(6));
t81 = cos(qJ(6));
t14 = t81 * t30 + t77 * t31;
t116 = -0.2e1 * t14;
t51 = t77 * t71 - t81 * t74;
t43 = t51 * t78;
t115 = 0.2e1 * t43;
t114 = -0.2e1 * t46;
t65 = -t74 * pkin(5) - pkin(4);
t113 = 0.2e1 * t65;
t112 = 0.2e1 * t82;
t111 = pkin(2) * t79;
t110 = pkin(2) * t83;
t109 = pkin(10) * t71;
t66 = t78 * pkin(10);
t108 = t82 * pkin(10);
t92 = pkin(9) * t102;
t40 = t92 + (pkin(10) + t111) * t75;
t41 = (-pkin(3) * t83 - pkin(10) * t79 - pkin(2)) * t72;
t24 = -t78 * t40 + t82 * t41;
t23 = pkin(4) * t102 - t24;
t107 = t23 * t71;
t106 = t23 * t74;
t68 = t72 ^ 2;
t105 = t68 * t83;
t104 = t71 * t78;
t73 = sin(pkin(6));
t80 = sin(qJ(2));
t101 = t73 * t80;
t84 = cos(qJ(2));
t100 = t73 * t84;
t99 = t74 * t78;
t98 = t75 * t79;
t97 = t75 * t84;
t96 = pkin(11) + qJ(5);
t60 = pkin(9) * t103;
t39 = t60 + (-pkin(3) - t110) * t75;
t45 = t78 * t103 - t82 * t75;
t19 = t45 * pkin(4) - t46 * qJ(5) + t39;
t25 = t82 * t40 + t78 * t41;
t22 = -qJ(5) * t102 + t25;
t8 = t71 * t19 + t74 * t22;
t55 = -t82 * pkin(4) - t78 * qJ(5) - pkin(3);
t37 = t74 * t108 + t71 * t55;
t95 = t71 ^ 2 + t74 ^ 2;
t94 = qJ(5) * t45;
t93 = 0.2e1 * t102;
t91 = t78 * t102;
t90 = t82 * t102;
t7 = t74 * t19 - t71 * t22;
t89 = -t7 * t71 + t8 * t74;
t76 = cos(pkin(6));
t29 = t76 * t103 + (t79 * t97 + t80 * t83) * t73;
t44 = -t72 * t100 + t76 * t75;
t21 = t29 * t82 + t44 * t78;
t28 = -t73 * t83 * t97 + t79 * t101 - t76 * t102;
t10 = t21 * t74 + t28 * t71;
t9 = -t21 * t71 + t28 * t74;
t88 = t10 * t74 - t9 * t71;
t87 = -pkin(4) * t78 + qJ(5) * t82;
t50 = t74 * t55;
t36 = -t71 * t108 + t50;
t86 = -t36 * t71 + t37 * t74;
t52 = t81 * t71 + t77 * t74;
t70 = t78 ^ 2;
t57 = t96 * t74;
t56 = t96 * t71;
t53 = pkin(5) * t104 + t66;
t48 = pkin(2) * t98 + t92;
t47 = t75 * t110 - t60;
t42 = t52 * t78;
t34 = -t77 * t56 + t81 * t57;
t33 = -t81 * t56 - t77 * t57;
t32 = -pkin(11) * t104 + t37;
t26 = -pkin(11) * t99 + t50 + (-pkin(5) - t109) * t82;
t20 = t29 * t78 - t44 * t82;
t15 = -t77 * t30 + t81 * t31;
t13 = t77 * t26 + t81 * t32;
t12 = t81 * t26 - t77 * t32;
t11 = t30 * pkin(5) + t23;
t6 = -t30 * pkin(11) + t8;
t5 = t45 * pkin(5) - t31 * pkin(11) + t7;
t4 = t81 * t10 + t77 * t9;
t3 = -t77 * t10 + t81 * t9;
t2 = t77 * t5 + t81 * t6;
t1 = t81 * t5 - t77 * t6;
t16 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t20 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t100, -t101, 0, 0, 0, 0, 0, -t44 * t102 - t28 * t75, t44 * t103 - t29 * t75, 0, 0, 0, 0, 0, t20 * t102 + t28 * t45, t21 * t102 + t28 * t46, t20 * t30 + t9 * t45, -t10 * t45 + t20 * t31, -t10 * t30 - t9 * t31, t10 * t8 + t20 * t23 + t9 * t7, 0, 0, 0, 0, 0, t20 * t14 + t3 * t45, t20 * t15 - t4 * t45; 0, 1, 0, 0, t68 * t79 ^ 2, 0.2e1 * t79 * t105, 0.2e1 * t72 * t98, t75 * t93, t75 ^ 2, 0.2e1 * pkin(2) * t105 + 0.2e1 * t47 * t75, -0.2e1 * t68 * t111 - 0.2e1 * t48 * t75, t46 ^ 2, t45 * t114, t102 * t114, t45 * t93, t68 * t83 ^ 2, -0.2e1 * t24 * t102 + 0.2e1 * t39 * t45, 0.2e1 * t25 * t102 + 0.2e1 * t39 * t46, 0.2e1 * t23 * t30 + 0.2e1 * t7 * t45, 0.2e1 * t23 * t31 - 0.2e1 * t8 * t45, -0.2e1 * t8 * t30 - 0.2e1 * t7 * t31, t23 ^ 2 + t7 ^ 2 + t8 ^ 2, t15 ^ 2, t15 * t116, 0.2e1 * t15 * t45, t45 * t116, t45 ^ 2, 0.2e1 * t1 * t45 + 0.2e1 * t11 * t14, 0.2e1 * t11 * t15 - 0.2e1 * t2 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29, 0, 0, 0, 0, 0, -t28 * t82, t28 * t78, t20 * t104 - t9 * t82, t10 * t82 + t20 * t99 (-t10 * t71 - t74 * t9) * t78, t10 * t37 + t20 * t66 + t9 * t36, 0, 0, 0, 0, 0, t20 * t42 - t3 * t82, -t20 * t43 + t4 * t82; 0, 0, 0, 0, 0, 0, t103, t102, t75, t47, -t48, t46 * t78, -t78 * t45 + t46 * t82, -t91, -t90, 0, -pkin(3) * t45 + pkin(10) * t91 - t39 * t82, -pkin(3) * t46 + pkin(10) * t90 + t39 * t78, t36 * t45 - t7 * t82 + (pkin(10) * t30 + t107) * t78, -t37 * t45 + t8 * t82 + (pkin(10) * t31 + t106) * t78, -t37 * t30 - t36 * t31 + (-t7 * t74 - t71 * t8) * t78, t23 * t66 + t7 * t36 + t8 * t37, -t15 * t43, t43 * t14 - t15 * t42, -t15 * t82 - t43 * t45, t14 * t82 - t42 * t45, -t45 * t82, -t1 * t82 + t11 * t42 + t12 * t45 + t53 * t14, -t11 * t43 - t13 * t45 + t53 * t15 + t2 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t70, t78 * t112, 0, 0, 0, pkin(3) * t112, -0.2e1 * pkin(3) * t78, 0.2e1 * t70 * t109 - 0.2e1 * t36 * t82, 0.2e1 * t70 * pkin(10) * t74 + 0.2e1 * t37 * t82, 0.2e1 * (-t36 * t74 - t37 * t71) * t78, t70 * pkin(10) ^ 2 + t36 ^ 2 + t37 ^ 2, t43 ^ 2, t42 * t115, t82 * t115, t42 * t112, t82 ^ 2, -0.2e1 * t12 * t82 + 0.2e1 * t53 * t42, 0.2e1 * t13 * t82 - 0.2e1 * t53 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, -t20 * t74, t20 * t71, t88, -t20 * pkin(4) + t88 * qJ(5), 0, 0, 0, 0, 0, t20 * t51, t20 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, -t102, t24, -t25, -pkin(4) * t30 - t71 * t94 - t106, -pkin(4) * t31 - t74 * t94 + t107 (-t30 * t74 + t31 * t71) * qJ(5) + t89, -t23 * pkin(4) + t89 * qJ(5), t15 * t52, -t52 * t14 - t15 * t51, t52 * t45, -t51 * t45, 0, t11 * t51 + t65 * t14 + t33 * t45, t11 * t52 + t65 * t15 - t34 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t82, 0, -t66, -t108, -pkin(10) * t99 + t87 * t71, pkin(10) * t104 + t87 * t74, t86, -pkin(4) * t66 + t86 * qJ(5), -t43 * t52, -t52 * t42 + t43 * t51, -t52 * t82, t51 * t82, 0, -t33 * t82 + t65 * t42 + t53 * t51, t34 * t82 - t65 * t43 + t53 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t74, -0.2e1 * pkin(4) * t71, 0.2e1 * t95 * qJ(5), t95 * qJ(5) ^ 2 + pkin(4) ^ 2, t52 ^ 2, -0.2e1 * t52 * t51, 0, 0, 0, t51 * t113, t52 * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t31, 0, t23, 0, 0, 0, 0, 0, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t99, 0, t66, 0, 0, 0, 0, 0, t42, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t71, 0, -pkin(4), 0, 0, 0, 0, 0, t51, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, t45, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t42, -t82, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t51, 0, t33, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t16;
