% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t101 = cos(qJ(4));
t65 = sin(pkin(11));
t67 = cos(pkin(11));
t70 = sin(qJ(4));
t48 = t101 * t65 + t70 * t67;
t108 = -0.2e1 * t48;
t68 = cos(pkin(6));
t66 = sin(pkin(6));
t71 = sin(qJ(2));
t97 = t66 * t71;
t40 = t65 * t97 - t68 * t67;
t41 = t68 * t65 + t67 * t97;
t27 = t101 * t41 - t70 * t40;
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t73 = cos(qJ(2));
t96 = t66 * t73;
t20 = t69 * t27 + t72 * t96;
t107 = -0.2e1 * t20;
t106 = -0.2e1 * t27;
t58 = -t67 * pkin(3) - pkin(2);
t105 = 0.2e1 * t58;
t104 = pkin(1) * t71;
t103 = pkin(1) * t73;
t102 = t72 * pkin(5);
t21 = t72 * t27 - t69 * t96;
t100 = t21 * t69;
t99 = t21 * t72;
t61 = t66 ^ 2;
t98 = t61 * t73;
t95 = t68 * t71;
t26 = t101 * t40 + t70 * t41;
t94 = t69 * t26;
t47 = -t101 * t67 + t70 * t65;
t93 = t69 * t47;
t92 = t69 * t48;
t91 = t69 * t72;
t25 = t72 * t26;
t88 = pkin(9) + qJ(3);
t49 = t88 * t65;
t50 = t88 * t67;
t32 = t101 * t50 - t70 * t49;
t90 = t72 * t32;
t89 = t72 * t48;
t87 = -qJ(6) - pkin(10);
t81 = pkin(8) * t96;
t36 = t81 + (qJ(3) + t104) * t68;
t37 = (-pkin(2) * t73 - qJ(3) * t71 - pkin(1)) * t66;
t23 = t67 * t36 + t65 * t37;
t86 = t65 ^ 2 + t67 ^ 2;
t63 = t69 ^ 2;
t64 = t72 ^ 2;
t85 = t63 + t64;
t84 = qJ(6) * t48;
t83 = t47 * t108;
t82 = 0.2e1 * t96;
t80 = qJ(3) * t96;
t55 = pkin(8) * t97;
t39 = t55 + (-pkin(2) - t103) * t68;
t29 = t40 * pkin(3) + t39;
t11 = t26 * pkin(4) - t27 * pkin(10) + t29;
t22 = -t65 * t36 + t67 * t37;
t14 = -pkin(3) * t96 - t41 * pkin(9) + t22;
t18 = -t40 * pkin(9) + t23;
t9 = t101 * t18 + t70 * t14;
t7 = -pkin(10) * t96 + t9;
t3 = t72 * t11 - t69 * t7;
t30 = t47 * pkin(4) - t48 * pkin(10) + t58;
t15 = t72 * t30 - t69 * t32;
t31 = t101 * t49 + t70 * t50;
t79 = -pkin(4) * t48 - pkin(10) * t47;
t1 = t26 * pkin(5) - t21 * qJ(6) + t3;
t4 = t69 * t11 + t72 * t7;
t2 = -t20 * qJ(6) + t4;
t78 = t1 * t72 + t2 * t69;
t8 = t101 * t14 - t70 * t18;
t12 = t47 * pkin(5) - t72 * t84 + t15;
t13 = t90 + (t30 - t84) * t69;
t77 = t12 * t72 + t13 * t69;
t76 = -t22 * t65 + t23 * t67;
t51 = t87 * t69;
t52 = t87 * t72;
t75 = t72 * t51 - t69 * t52;
t6 = pkin(4) * t96 - t8;
t59 = -pkin(4) - t102;
t45 = t48 ^ 2;
t44 = pkin(1) * t95 + t81;
t43 = t68 * t103 - t55;
t42 = t72 * t47;
t24 = pkin(5) * t92 + t31;
t19 = t69 * t20;
t16 = t69 * t30 + t90;
t5 = t20 * pkin(5) + t6;
t10 = [1, 0, 0, t61 * t71 ^ 2, 0.2e1 * t71 * t98, 0.2e1 * t66 * t95, t68 * t82, t68 ^ 2, 0.2e1 * pkin(1) * t98 + 0.2e1 * t43 * t68, -0.2e1 * t61 * t104 - 0.2e1 * t44 * t68, -0.2e1 * t22 * t96 + 0.2e1 * t39 * t40, 0.2e1 * t23 * t96 + 0.2e1 * t39 * t41, -0.2e1 * t22 * t41 - 0.2e1 * t23 * t40, t22 ^ 2 + t23 ^ 2 + t39 ^ 2, t27 ^ 2, t26 * t106, t96 * t106, t26 * t82, t61 * t73 ^ 2, 0.2e1 * t29 * t26 - 0.2e1 * t8 * t96, 0.2e1 * t29 * t27 + 0.2e1 * t9 * t96, t21 ^ 2, t21 * t107, 0.2e1 * t21 * t26, t26 * t107, t26 ^ 2, 0.2e1 * t6 * t20 + 0.2e1 * t3 * t26, 0.2e1 * t6 * t21 - 0.2e1 * t4 * t26, -0.2e1 * t1 * t21 - 0.2e1 * t2 * t20, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t97, t96, t68, t43, -t44, -pkin(2) * t40 - t39 * t67 + t65 * t80, -pkin(2) * t41 + t39 * t65 + t67 * t80 (-t40 * t67 + t41 * t65) * qJ(3) + t76, -t39 * pkin(2) + t76 * qJ(3), t27 * t48, -t48 * t26 - t27 * t47, -t48 * t96, t47 * t96, 0, t58 * t26 + t29 * t47 + t31 * t96, t58 * t27 + t29 * t48 + t32 * t96, t21 * t89 (-t20 * t72 - t100) * t48, t21 * t47 + t26 * t89, -t20 * t47 - t26 * t92, t26 * t47, t15 * t26 + t31 * t20 + t3 * t47 + t6 * t92, -t16 * t26 + t31 * t21 - t4 * t47 + t6 * t89, -t12 * t21 - t13 * t20 - t78 * t48, t1 * t12 + t2 * t13 + t5 * t24; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t67, -0.2e1 * pkin(2) * t65, 0.2e1 * t86 * qJ(3), t86 * qJ(3) ^ 2 + pkin(2) ^ 2, t45, t83, 0, 0, 0, t47 * t105, t48 * t105, t64 * t45, -0.2e1 * t45 * t91, 0.2e1 * t47 * t89, t69 * t83, t47 ^ 2, 0.2e1 * t15 * t47 + 0.2e1 * t31 * t92, -0.2e1 * t16 * t47 + 0.2e1 * t31 * t89, t77 * t108, t12 ^ 2 + t13 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t41, 0, t39, 0, 0, 0, 0, 0, t26, t27, 0, 0, 0, 0, 0, t25, -t94, -t19 - t99, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t65, 0, -pkin(2), 0, 0, 0, 0, 0, t47, t48, 0, 0, 0, 0, 0, t42, -t93, -t85 * t48, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, -t96, t8, -t9, t100, -t19 + t99, t94, t25, 0, -pkin(4) * t20 - pkin(10) * t94 - t6 * t72, -pkin(4) * t21 - pkin(10) * t25 + t6 * t69, -t1 * t69 + t2 * t72 + t52 * t20 - t51 * t21, t1 * t51 - t2 * t52 + t5 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t47, 0, -t31, -t32, t69 * t89 (-t63 + t64) * t48, t93, t42, 0, -t31 * t72 + t79 * t69, t31 * t69 + t79 * t72, -t12 * t69 + t13 * t72 - t75 * t48, t12 * t51 - t13 * t52 + t24 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t63, 0.2e1 * t91, 0, 0, 0, 0.2e1 * pkin(4) * t72, -0.2e1 * pkin(4) * t69, -0.2e1 * t51 * t69 - 0.2e1 * t52 * t72, t51 ^ 2 + t52 ^ 2 + t59 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, t26, t3, -t4, -pkin(5) * t21, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -t92, t47, t15, -t16, -pkin(5) * t89, t12 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t69, 0, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t72, 0, -t69 * pkin(10), -t72 * pkin(10), -t69 * pkin(5), t51 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t10;
