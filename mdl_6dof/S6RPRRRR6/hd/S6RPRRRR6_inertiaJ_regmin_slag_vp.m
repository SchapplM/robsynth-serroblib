% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRR6
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t67 = sin(qJ(5));
t68 = sin(qJ(4));
t71 = cos(qJ(5));
t72 = cos(qJ(4));
t45 = t67 * t68 - t71 * t72;
t57 = -t72 * pkin(4) - pkin(3);
t34 = t45 * pkin(5) + t57;
t98 = 0.2e1 * t34;
t64 = sin(pkin(11));
t65 = cos(pkin(11));
t69 = sin(qJ(3));
t87 = cos(qJ(3));
t43 = t69 * t64 - t87 * t65;
t97 = -0.2e1 * t43;
t96 = 0.2e1 * t43;
t54 = -t65 * pkin(2) - pkin(1);
t95 = 0.2e1 * t54;
t94 = 0.2e1 * t57;
t93 = pkin(8) + pkin(9);
t92 = t43 * pkin(4);
t91 = t43 * pkin(5);
t66 = sin(qJ(6));
t90 = t66 * pkin(5);
t89 = t67 * pkin(4);
t70 = cos(qJ(6));
t58 = t70 * pkin(5);
t44 = t87 * t64 + t69 * t65;
t46 = t67 * t72 + t71 * t68;
t23 = t46 * t44;
t26 = t43 * pkin(3) - t44 * pkin(8) + t54;
t78 = pkin(7) + qJ(2);
t48 = t78 * t64;
t49 = t78 * t65;
t30 = -t69 * t48 + t87 * t49;
t80 = t72 * t30;
t10 = t80 + (-pkin(9) * t44 + t26) * t68;
t81 = t71 * t10;
t16 = t72 * t26 - t68 * t30;
t79 = t72 * t44;
t9 = -pkin(9) * t79 + t16 + t92;
t7 = t67 * t9 + t81;
t5 = -t23 * pkin(10) + t7;
t88 = t70 * t5;
t59 = t71 * pkin(4);
t28 = -t66 * t45 + t70 * t46;
t86 = t28 * t43;
t85 = t46 * t43;
t84 = t68 * t43;
t83 = t68 * t44;
t82 = t68 * t72;
t77 = t64 ^ 2 + t65 ^ 2;
t76 = t44 * t97;
t75 = t70 * t89;
t24 = t45 * t44;
t6 = -t67 * t10 + t71 * t9;
t4 = t24 * pkin(10) + t6 + t91;
t1 = t70 * t4 - t66 * t5;
t50 = t93 * t68;
t51 = t93 * t72;
t32 = -t71 * t50 - t67 * t51;
t56 = t59 + pkin(5);
t36 = t70 * t56 - t66 * t89;
t74 = -pkin(3) * t44 - pkin(8) * t43;
t2 = t66 * t4 + t88;
t29 = t87 * t48 + t69 * t49;
t33 = -t67 * t50 + t71 * t51;
t19 = pkin(4) * t83 + t29;
t63 = t72 ^ 2;
t62 = t68 ^ 2;
t41 = t44 ^ 2;
t40 = t43 ^ 2;
t38 = t72 * t43;
t37 = t66 * t56 + t75;
t31 = t45 * t43;
t27 = t70 * t45 + t66 * t46;
t22 = -t45 * pkin(10) + t33;
t21 = -t46 * pkin(10) + t32;
t18 = t27 * t43;
t17 = t68 * t26 + t80;
t15 = -t66 * t23 - t70 * t24;
t14 = t70 * t23 - t66 * t24;
t13 = t23 * pkin(5) + t19;
t12 = t66 * t21 + t70 * t22;
t11 = t70 * t21 - t66 * t22;
t3 = [1, 0, 0, 0.2e1 * pkin(1) * t65, -0.2e1 * pkin(1) * t64, 0.2e1 * t77 * qJ(2), t77 * qJ(2) ^ 2 + pkin(1) ^ 2, t41, t76, 0, 0, 0, t43 * t95, t44 * t95, t63 * t41, -0.2e1 * t41 * t82, t79 * t96, t68 * t76, t40, 0.2e1 * t16 * t43 + 0.2e1 * t29 * t83, -0.2e1 * t17 * t43 + 0.2e1 * t29 * t79, t24 ^ 2, 0.2e1 * t24 * t23, -t24 * t96, t23 * t97, t40, 0.2e1 * t19 * t23 + 0.2e1 * t6 * t43, -0.2e1 * t19 * t24 - 0.2e1 * t7 * t43, t15 ^ 2, -0.2e1 * t15 * t14, t15 * t96, t14 * t97, t40, 0.2e1 * t1 * t43 + 0.2e1 * t13 * t14, 0.2e1 * t13 * t15 - 0.2e1 * t2 * t43; 0, 0, 0, -t65, t64, 0, -pkin(1), 0, 0, 0, 0, 0, t43, t44, 0, 0, 0, 0, 0, t38, -t84, 0, 0, 0, 0, 0, -t31, -t85, 0, 0, 0, 0, 0, -t18, -t86; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43, 0, -t29, -t30, t68 * t79 (-t62 + t63) * t44, t84, t38, 0, -t29 * t72 + t74 * t68, t29 * t68 + t74 * t72, -t24 * t46, -t46 * t23 + t24 * t45, t85, -t31, 0, t19 * t45 + t57 * t23 + t32 * t43, t19 * t46 - t57 * t24 - t33 * t43, t15 * t28, -t28 * t14 - t15 * t27, t86, -t18, 0, t11 * t43 + t13 * t27 + t34 * t14, -t12 * t43 + t13 * t28 + t34 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t62, 0.2e1 * t82, 0, 0, 0, 0.2e1 * pkin(3) * t72, -0.2e1 * pkin(3) * t68, t46 ^ 2, -0.2e1 * t46 * t45, 0, 0, 0, t45 * t94, t46 * t94, t28 ^ 2, -0.2e1 * t28 * t27, 0, 0, 0, t27 * t98, t28 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t83, t43, t16, -t17, 0, 0, -t24, -t23, t43, t43 * t59 + t6, -t81 + (-t9 - t92) * t67, 0, 0, t15, -t14, t43, t36 * t43 + t1, -t37 * t43 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t68, 0, 0, 0, 0, 0, -t45, -t46, 0, 0, 0, 0, 0, -t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t72, 0, -t68 * pkin(8), -t72 * pkin(8), 0, 0, t46, -t45, 0, t32, -t33, 0, 0, t28, -t27, 0, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t59, -0.2e1 * t89, 0, 0, 0, 0, 1, 0.2e1 * t36, -0.2e1 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, t43, t6, -t7, 0, 0, t15, -t14, t43, t43 * t58 + t1, -t88 + (-t4 - t91) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t46, 0, 0, 0, 0, 0, -t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t45, 0, t32, -t33, 0, 0, t28, -t27, 0, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t59, -t89, 0, 0, 0, 0, 1, t36 + t58, -t75 + (-pkin(5) - t56) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, t43, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, 0, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t36, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t58, -t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
