% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t64 = sin(pkin(11));
t65 = cos(pkin(11));
t88 = sin(qJ(3));
t90 = cos(qJ(3));
t43 = t90 * t64 + t88 * t65;
t68 = sin(qJ(4));
t75 = t88 * t64 - t90 * t65;
t89 = cos(qJ(4));
t33 = t68 * t43 + t89 * t75;
t104 = -0.2e1 * t33;
t103 = 0.2e1 * t33;
t54 = -t65 * pkin(2) - pkin(1);
t38 = t75 * pkin(3) + t54;
t102 = 0.2e1 * t38;
t77 = t89 * pkin(3);
t57 = -t77 - pkin(4);
t70 = cos(qJ(5));
t93 = t70 * pkin(5);
t49 = t57 - t93;
t101 = 0.2e1 * t49;
t100 = 0.2e1 * t54;
t58 = -pkin(4) - t93;
t99 = 0.2e1 * t58;
t98 = t33 * pkin(5);
t66 = sin(qJ(6));
t97 = t66 * pkin(5);
t96 = t68 * pkin(3);
t69 = cos(qJ(6));
t95 = t69 * pkin(5);
t34 = t89 * t43 - t68 * t75;
t17 = t33 * pkin(4) - t34 * pkin(9) + t38;
t67 = sin(qJ(5));
t81 = pkin(7) + qJ(2);
t47 = t81 * t64;
t48 = t81 * t65;
t74 = -t90 * t47 - t88 * t48;
t24 = -t43 * pkin(8) + t74;
t72 = t88 * t47 - t90 * t48;
t25 = -t75 * pkin(8) - t72;
t16 = t68 * t24 + t89 * t25;
t84 = t70 * t16;
t5 = t84 + (-pkin(10) * t34 + t17) * t67;
t94 = t69 * t5;
t92 = t70 * pkin(9);
t91 = pkin(4) - t57;
t15 = -t89 * t24 + t68 * t25;
t87 = t15 * t70;
t46 = t66 * t70 + t69 * t67;
t22 = t46 * t33;
t27 = t67 * t33;
t86 = t67 * t34;
t85 = t67 * t70;
t83 = t70 * t34;
t56 = pkin(9) + t96;
t82 = t70 * t56;
t80 = t49 + t58;
t79 = t64 ^ 2 + t65 ^ 2;
t78 = t34 * t104;
t6 = -t67 * t16 + t70 * t17;
t4 = -pkin(10) * t83 + t6 + t98;
t1 = t69 * t4 - t66 * t5;
t76 = -pkin(4) * t34 - pkin(9) * t33;
t73 = -t33 * t56 + t34 * t57;
t45 = t66 * t67 - t69 * t70;
t63 = t70 ^ 2;
t62 = t67 ^ 2;
t59 = t70 * pkin(10);
t52 = 0.2e1 * t85;
t51 = t59 + t92;
t50 = (-pkin(9) - pkin(10)) * t67;
t44 = t46 ^ 2;
t41 = t59 + t82;
t40 = (-pkin(10) - t56) * t67;
t37 = t66 * t50 + t69 * t51;
t36 = t69 * t50 - t66 * t51;
t35 = -0.2e1 * t46 * t45;
t32 = t34 ^ 2;
t31 = t33 ^ 2;
t30 = t66 * t40 + t69 * t41;
t29 = t69 * t40 - t66 * t41;
t28 = t70 * t33;
t26 = t67 * t83;
t21 = t45 * t33;
t20 = (-t62 + t63) * t34;
t19 = t45 * t34;
t18 = t46 * t34;
t14 = t19 * t46;
t12 = t15 * t67;
t11 = pkin(5) * t86 + t15;
t10 = t11 * t46;
t9 = t11 * t45;
t8 = -t46 * t18 + t19 * t45;
t7 = t67 * t17 + t84;
t2 = t66 * t4 + t94;
t3 = [1, 0, 0, 0.2e1 * pkin(1) * t65, -0.2e1 * pkin(1) * t64, 0.2e1 * t79 * qJ(2), t79 * qJ(2) ^ 2 + pkin(1) ^ 2, t43 ^ 2, -0.2e1 * t43 * t75, 0, 0, 0, t75 * t100, t43 * t100, t32, t78, 0, 0, 0, t33 * t102, t34 * t102, t63 * t32, -0.2e1 * t32 * t85, t83 * t103, t67 * t78, t31, 0.2e1 * t15 * t86 + 0.2e1 * t6 * t33, 0.2e1 * t15 * t83 - 0.2e1 * t7 * t33, t19 ^ 2, 0.2e1 * t19 * t18, -t19 * t103, t18 * t104, t31, 0.2e1 * t1 * t33 + 0.2e1 * t11 * t18, -0.2e1 * t11 * t19 - 0.2e1 * t2 * t33; 0, 0, 0, -t65, t64, 0, -pkin(1), 0, 0, 0, 0, 0, t75, t43, 0, 0, 0, 0, 0, t33, t34, 0, 0, 0, 0, 0, t28, -t27, 0, 0, 0, 0, 0, -t21, -t22; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t75, 0, t74, t72, 0, 0, t34, -t33, 0, -t15, -t16, t26, t20, t27, t28, 0, t73 * t67 - t87, t70 * t73 + t12, -t14, t8, t22, -t21, 0, t49 * t18 + t29 * t33 + t9, -t49 * t19 - t30 * t33 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t77, -0.2e1 * t96, t62, t52, 0, 0, 0, -0.2e1 * t57 * t70, 0.2e1 * t57 * t67, t44, t35, 0, 0, 0, t45 * t101, t46 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t33, 0, -t15, -t16, t26, t20, t27, t28, 0, t76 * t67 - t87, t70 * t76 + t12, -t14, t8, t22, -t21, 0, t58 * t18 + t36 * t33 + t9, -t58 * t19 - t37 * t33 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t77, -t96, t62, t52, 0, 0, 0, t91 * t70, -t91 * t67, t44, t35, 0, 0, 0, t80 * t45, t80 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t62, t52, 0, 0, 0, 0.2e1 * pkin(4) * t70, -0.2e1 * pkin(4) * t67, t44, t35, 0, 0, 0, t45 * t99, t46 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t86, t33, t6, -t7, 0, 0, -t19, -t18, t33, t33 * t95 + t1, -t94 + (-t4 - t98) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t67, 0, 0, 0, 0, 0, -t45, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t70, 0, -t67 * t56, -t82, 0, 0, t46, -t45, 0, t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t70, 0, -t67 * pkin(9), -t92, 0, 0, t46, -t45, 0, t36, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t95, -0.2e1 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t18, t33, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, 0, t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, 0, t36, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t95, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
