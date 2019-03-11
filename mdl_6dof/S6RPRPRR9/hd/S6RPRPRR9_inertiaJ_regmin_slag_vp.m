% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t67 = sin(pkin(12));
t69 = sin(pkin(6));
t71 = cos(pkin(12));
t73 = cos(pkin(6));
t76 = sin(qJ(3));
t72 = cos(pkin(7));
t79 = cos(qJ(3));
t95 = t72 * t79;
t68 = sin(pkin(7));
t98 = t68 * t79;
t35 = -t73 * t98 + (t67 * t76 - t71 * t95) * t69;
t96 = t71 * t72;
t99 = t68 * t76;
t36 = t73 * t99 + (t67 * t79 + t76 * t96) * t69;
t66 = sin(pkin(13));
t70 = cos(pkin(13));
t28 = -t35 * t66 + t36 * t70;
t97 = t69 * t71;
t49 = t68 * t97 - t73 * t72;
t75 = sin(qJ(5));
t78 = cos(qJ(5));
t23 = t28 * t75 + t49 * t78;
t112 = -0.2e1 * t23;
t111 = -0.2e1 * t36;
t110 = 0.2e1 * t75;
t62 = t69 ^ 2;
t109 = pkin(1) * t62;
t108 = pkin(1) * t73;
t77 = cos(qJ(6));
t107 = pkin(5) * t77;
t27 = t70 * t35 + t36 * t66;
t100 = t67 * t69;
t57 = t71 * t108;
t37 = pkin(2) * t73 + t57 + (-pkin(9) * t72 - qJ(2)) * t100;
t44 = (-pkin(9) * t67 * t68 - pkin(2) * t71 - pkin(1)) * t69;
t29 = -t37 * t68 + t72 * t44;
t22 = pkin(3) * t35 + t29;
t12 = pkin(4) * t27 - pkin(10) * t28 + t22;
t84 = qJ(2) * t69;
t52 = t67 * t108 + t71 * t84;
t34 = (t68 * t73 + t69 * t96) * pkin(9) + t52;
t20 = -t34 * t76 + t37 * t95 + t44 * t98;
t16 = -pkin(3) * t49 - qJ(4) * t36 + t20;
t21 = t34 * t79 + (t37 * t72 + t44 * t68) * t76;
t19 = -qJ(4) * t35 + t21;
t11 = t66 * t16 + t70 * t19;
t9 = -pkin(10) * t49 + t11;
t5 = t12 * t78 - t75 * t9;
t3 = -pkin(5) * t27 - t5;
t74 = sin(qJ(6));
t106 = t3 * t74;
t105 = t3 * t77;
t24 = t28 * t78 - t49 * t75;
t13 = t24 * t74 - t77 * t27;
t104 = t13 * t78;
t14 = t24 * t77 + t27 * t74;
t103 = t14 * t74;
t102 = t14 * t78;
t60 = pkin(3) * t66 + pkin(10);
t101 = t60 * t74;
t94 = t74 * t23;
t93 = t74 * t75;
t92 = t74 * t77;
t91 = t74 * t78;
t90 = t75 * t27;
t89 = t75 * t60;
t88 = t77 * t23;
t87 = t77 * t75;
t86 = t77 * t78;
t85 = t78 * t60;
t83 = t78 * t110;
t82 = t23 * t93;
t81 = t23 * t87;
t61 = -pkin(3) * t70 - pkin(4);
t10 = t16 * t70 - t66 * t19;
t6 = t12 * t75 + t78 * t9;
t8 = pkin(4) * t49 - t10;
t65 = t77 ^ 2;
t64 = t75 ^ 2;
t63 = t74 ^ 2;
t53 = -pkin(5) * t78 - pkin(11) * t75 + t61;
t51 = -t67 * t84 + t57;
t48 = (t66 * t79 + t70 * t76) * t68;
t46 = t66 * t99 - t70 * t98;
t41 = t48 * t78 + t72 * t75;
t40 = t48 * t75 - t78 * t72;
t39 = t53 * t74 + t77 * t85;
t38 = t53 * t77 - t74 * t85;
t31 = t41 * t77 + t46 * t74;
t30 = -t41 * t74 + t46 * t77;
t26 = t78 * t27;
t7 = pkin(5) * t23 - pkin(11) * t24 + t8;
t4 = pkin(11) * t27 + t6;
t2 = t4 * t77 + t7 * t74;
t1 = -t4 * t74 + t7 * t77;
t15 = [1, 0, 0, 0.2e1 * t71 * t109 + 0.2e1 * t51 * t73, -0.2e1 * t67 * t109 - 0.2e1 * t52 * t73, 0.2e1 * (-t51 * t67 + t52 * t71) * t69, pkin(1) ^ 2 * t62 + t51 ^ 2 + t52 ^ 2, t36 ^ 2, t35 * t111, t49 * t111, 0.2e1 * t35 * t49, t49 ^ 2, -0.2e1 * t20 * t49 + 0.2e1 * t29 * t35, 0.2e1 * t21 * t49 + 0.2e1 * t29 * t36, -0.2e1 * t10 * t28 - 0.2e1 * t11 * t27, t10 ^ 2 + t11 ^ 2 + t22 ^ 2, t24 ^ 2, t24 * t112, 0.2e1 * t24 * t27, t27 * t112, t27 ^ 2, 0.2e1 * t23 * t8 + 0.2e1 * t27 * t5, 0.2e1 * t24 * t8 - 0.2e1 * t27 * t6, t14 ^ 2, -0.2e1 * t14 * t13, 0.2e1 * t14 * t23, t13 * t112, t23 ^ 2, 0.2e1 * t1 * t23 + 0.2e1 * t13 * t3, 0.2e1 * t14 * t3 - 0.2e1 * t2 * t23; 0, 0, 0, -t97, t100, 0, -t69 * pkin(1), 0, 0, 0, 0, 0, t35 * t72 - t49 * t98, t36 * t72 + t49 * t99, -t27 * t48 + t28 * t46, -t10 * t46 + t11 * t48 + t22 * t72, 0, 0, 0, 0, 0, t23 * t46 - t27 * t40, t24 * t46 - t27 * t41, 0, 0, 0, 0, 0, t13 * t40 + t23 * t30, t14 * t40 - t23 * t31; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t46 ^ 2 + t48 ^ 2 + t72 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, -t49, t20, -t21 (-t27 * t66 - t28 * t70) * pkin(3) (t10 * t70 + t11 * t66) * pkin(3), t24 * t75, -t23 * t75 + t24 * t78, t90, t26, 0, t23 * t61 - t27 * t89 - t78 * t8, t24 * t61 - t27 * t85 + t75 * t8, t14 * t87 (-t13 * t77 - t103) * t75, t81 - t102, -t82 + t104, -t23 * t78, -t1 * t78 + t23 * t38 + (t13 * t60 + t106) * t75, t2 * t78 - t23 * t39 + (t14 * t60 + t105) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, -t99, 0 (-t46 * t70 + t48 * t66) * pkin(3), 0, 0, 0, 0, 0, -t46 * t78, t46 * t75, 0, 0, 0, 0, 0, -t30 * t78 + t40 * t93, t31 * t78 + t40 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t66 ^ 2 + t70 ^ 2) * pkin(3) ^ 2, t64, t83, 0, 0, 0, -0.2e1 * t61 * t78, t61 * t110, t65 * t64, -0.2e1 * t64 * t92, -0.2e1 * t75 * t86, t74 * t83, t78 ^ 2, 0.2e1 * t64 * t101 - 0.2e1 * t38 * t78, 0.2e1 * t60 * t64 * t77 + 0.2e1 * t39 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, t26, -t90, 0, 0, 0, 0, 0, -t82 - t104, -t81 - t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, t27, t5, -t6, t103, -t13 * t74 + t14 * t77, t94, t88, 0, -pkin(5) * t13 - pkin(11) * t94 - t105, -pkin(5) * t14 - pkin(11) * t88 + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t41, 0, 0, 0, 0, 0, -t40 * t77, t40 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t78, 0, -t89, -t85, t74 * t87 (-t63 + t65) * t75, -t91, -t86, 0, -t60 * t87 + (-pkin(5) * t75 + pkin(11) * t78) * t74, pkin(11) * t86 + (t101 - t107) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t75, 0, 0, 0, 0, 0, t86, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t63, 0.2e1 * t92, 0, 0, 0, 0.2e1 * t107, -0.2e1 * pkin(5) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t23, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t93, -t78, t38, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t77, 0, -t74 * pkin(11), -t77 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t15;
