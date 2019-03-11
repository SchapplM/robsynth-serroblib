% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t51 = sin(qJ(4));
t88 = t51 * pkin(3);
t40 = pkin(9) + t88;
t50 = sin(qJ(5));
t46 = t50 ^ 2;
t52 = cos(qJ(5));
t47 = t52 ^ 2;
t69 = t46 + t47;
t72 = t69 * t40;
t48 = sin(pkin(10));
t49 = cos(pkin(10));
t81 = sin(qJ(3));
t83 = cos(qJ(3));
t29 = t48 * t83 + t49 * t81;
t62 = t48 * t81 - t49 * t83;
t82 = cos(qJ(4));
t24 = t29 * t82 - t51 * t62;
t96 = -0.2e1 * t24;
t39 = -t49 * pkin(2) - pkin(1);
t25 = pkin(3) * t62 + t39;
t95 = 0.2e1 * t25;
t94 = 0.2e1 * t39;
t93 = -0.2e1 * t50;
t92 = -0.2e1 * t52;
t74 = pkin(7) + qJ(2);
t31 = t74 * t48;
t32 = t74 * t49;
t61 = -t31 * t83 - t32 * t81;
t16 = -pkin(8) * t29 + t61;
t55 = t31 * t81 - t32 * t83;
t17 = -pkin(8) * t62 - t55;
t12 = t16 * t51 + t17 * t82;
t23 = t29 * t51 + t62 * t82;
t13 = t23 * pkin(4) - t24 * pkin(9) + t25;
t5 = t12 * t52 + t13 * t50;
t91 = pkin(9) * t23;
t90 = t23 * pkin(5);
t89 = t50 * pkin(9);
t87 = t52 * pkin(9);
t11 = -t16 * t82 + t51 * t17;
t58 = pkin(5) * t50 - qJ(6) * t52;
t6 = t24 * t58 + t11;
t86 = t6 * t50;
t85 = t6 * t52;
t66 = t82 * pkin(3);
t41 = -t66 - pkin(4);
t84 = pkin(4) - t41;
t80 = t11 * t52;
t79 = t23 * t40;
t19 = t50 * t23;
t78 = t50 * t24;
t77 = t50 * t40;
t76 = t50 * t52;
t20 = t52 * t23;
t21 = t52 * t24;
t75 = t52 * t40;
t59 = pkin(5) * t52 + qJ(6) * t50;
t30 = -pkin(4) - t59;
t27 = -t66 + t30;
t73 = -t27 - t30;
t71 = t69 * pkin(9);
t70 = t48 ^ 2 + t49 ^ 2;
t68 = t23 * qJ(6);
t67 = t23 * t96;
t65 = t12 * t50 - t13 * t52;
t64 = -pkin(4) * t24 - t91;
t2 = t68 + t5;
t3 = t65 - t90;
t1 = t2 * t52 + t3 * t50;
t63 = t2 * t50 - t3 * t52;
t60 = -t24 * t30 + t91;
t57 = -t24 * t27 + t79;
t56 = t24 * t41 - t79;
t36 = 0.2e1 * t76;
t22 = t24 ^ 2;
t18 = t50 * t21;
t14 = (-t46 + t47) * t24;
t8 = t11 * t50;
t4 = [1, 0, 0, 0.2e1 * pkin(1) * t49, -0.2e1 * pkin(1) * t48, 0.2e1 * t70 * qJ(2), qJ(2) ^ 2 * t70 + pkin(1) ^ 2, t29 ^ 2, -0.2e1 * t29 * t62, 0, 0, 0, t62 * t94, t29 * t94, t22, t67, 0, 0, 0, t23 * t95, t24 * t95, t47 * t22, -0.2e1 * t22 * t76, 0.2e1 * t23 * t21, t50 * t67, t23 ^ 2, 0.2e1 * t11 * t78 - 0.2e1 * t23 * t65, 0.2e1 * t11 * t21 - 0.2e1 * t23 * t5, -0.2e1 * t23 * t3 + 0.2e1 * t6 * t78, t63 * t96, 0.2e1 * t2 * t23 - 0.2e1 * t21 * t6, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, -t49, t48, 0, -pkin(1), 0, 0, 0, 0, 0, t62, t29, 0, 0, 0, 0, 0, t23, t24, 0, 0, 0, 0, 0, t20, -t19, t20, -t69 * t24, t19, t63; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t62, 0, t61, t55, 0, 0, t24, -t23, 0, -t11, -t12, t18, t14, t19, t20, 0, t50 * t56 - t80, t52 * t56 + t8, -t50 * t57 - t85, t1, t52 * t57 - t86, t1 * t40 + t6 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t66, -0.2e1 * t88, t46, t36, 0, 0, 0, t41 * t92, 0.2e1 * t41 * t50, t27 * t92, 0.2e1 * t72, t27 * t93, t40 ^ 2 * t69 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, -t11, -t12, t18, t14, t19, t20, 0, t50 * t64 - t80, t52 * t64 + t8, -t50 * t60 - t85, t1, t52 * t60 - t86, pkin(9) * t1 + t6 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t66, -t88, t46, t36, 0, 0, 0, t84 * t52, -t84 * t50, t73 * t52, t71 + t72, t73 * t50, pkin(9) * t72 + t27 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t46, t36, 0, 0, 0, 0.2e1 * pkin(4) * t52, pkin(4) * t93, t30 * t92, 0.2e1 * t71, t30 * t93, pkin(9) ^ 2 * t69 + t30 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t78, t23, -t65, -t5, -t65 + 0.2e1 * t90, -t59 * t24, 0.2e1 * t68 + t5, -pkin(5) * t3 + qJ(6) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t50, t52, 0, t50, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t52, 0, -t77, -t75, -t77, -t58, t75, -t58 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t52, 0, -t89, -t87, -t89, -t58, t87, -t58 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t21, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
