% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t64 = sin(pkin(11));
t66 = cos(pkin(11));
t83 = t64 ^ 2 + t66 ^ 2;
t84 = t83 * qJ(5);
t65 = sin(pkin(10));
t67 = cos(pkin(10));
t70 = sin(qJ(3));
t72 = cos(qJ(3));
t42 = t72 * t65 + t70 * t67;
t69 = sin(qJ(4));
t79 = -t70 * t65 + t72 * t67;
t91 = cos(qJ(4));
t29 = t69 * t42 - t91 * t79;
t99 = -0.2e1 * t29;
t55 = -t67 * pkin(2) - pkin(1);
t34 = -t79 * pkin(3) + t55;
t98 = 0.2e1 * t34;
t97 = 0.2e1 * t42;
t81 = t91 * pkin(3);
t56 = -t81 - pkin(4);
t94 = t66 * pkin(5);
t43 = t56 - t94;
t96 = 0.2e1 * t43;
t54 = -pkin(4) - t94;
t95 = 0.2e1 * t54;
t93 = t69 * pkin(3);
t92 = pkin(4) - t56;
t87 = pkin(7) + qJ(2);
t45 = t87 * t65;
t47 = t87 * t67;
t80 = -t72 * t45 - t70 * t47;
t25 = -t42 * pkin(8) + t80;
t75 = t70 * t45 - t72 * t47;
t26 = t79 * pkin(8) - t75;
t18 = -t91 * t25 + t69 * t26;
t90 = t18 * t66;
t68 = sin(qJ(6));
t71 = cos(qJ(6));
t41 = t71 * t64 + t68 * t66;
t23 = t41 * t29;
t30 = t91 * t42 + t69 * t79;
t89 = t64 * t30;
t88 = t66 * t30;
t17 = t29 * pkin(4) - t30 * qJ(5) + t34;
t19 = t69 * t25 + t91 * t26;
t7 = t64 * t17 + t66 * t19;
t86 = t43 + t54;
t53 = qJ(5) + t93;
t85 = t83 * t53;
t82 = t65 ^ 2 + t67 ^ 2;
t6 = t66 * t17 - t64 * t19;
t78 = t6 * t66 + t7 * t64;
t3 = -t6 * t64 + t7 * t66;
t77 = -pkin(4) * t30 - qJ(5) * t29;
t76 = -t29 * t53 + t30 * t56;
t40 = t68 * t64 - t71 * t66;
t59 = t66 * pkin(9);
t46 = t66 * qJ(5) + t59;
t44 = (-pkin(9) - qJ(5)) * t64;
t38 = t41 ^ 2;
t37 = t66 * t53 + t59;
t36 = (-pkin(9) - t53) * t64;
t33 = t68 * t44 + t71 * t46;
t32 = t71 * t44 - t68 * t46;
t31 = -0.2e1 * t41 * t40;
t28 = t68 * t36 + t71 * t37;
t27 = t71 * t36 - t68 * t37;
t22 = t40 * t29;
t21 = t40 * t30;
t20 = t41 * t30;
t16 = t21 * t41;
t15 = t18 * t64;
t11 = pkin(5) * t89 + t18;
t10 = t11 * t41;
t9 = t11 * t40;
t8 = -t41 * t20 + t21 * t40;
t5 = -pkin(9) * t89 + t7;
t4 = t29 * pkin(5) - pkin(9) * t88 + t6;
t2 = t68 * t4 + t71 * t5;
t1 = t71 * t4 - t68 * t5;
t12 = [1, 0, 0, 0.2e1 * pkin(1) * t67, -0.2e1 * pkin(1) * t65, 0.2e1 * t82 * qJ(2), t82 * qJ(2) ^ 2 + pkin(1) ^ 2, t42 ^ 2, t79 * t97, 0, 0, 0, -0.2e1 * t55 * t79, t55 * t97, t30 ^ 2, t30 * t99, 0, 0, 0, t29 * t98, t30 * t98, 0.2e1 * t18 * t89 + 0.2e1 * t6 * t29, 0.2e1 * t18 * t88 - 0.2e1 * t7 * t29, -0.2e1 * t78 * t30, t18 ^ 2 + t6 ^ 2 + t7 ^ 2, t21 ^ 2, 0.2e1 * t21 * t20, t21 * t99, t20 * t99, t29 ^ 2, 0.2e1 * t1 * t29 + 0.2e1 * t11 * t20, -0.2e1 * t11 * t21 - 0.2e1 * t2 * t29; 0, 0, 0, -t67, t65, 0, -pkin(1), 0, 0, 0, 0, 0, -t79, t42, 0, 0, 0, 0, 0, t29, t30, t66 * t29, -t64 * t29, -t83 * t30, t78, 0, 0, 0, 0, 0, -t22, -t23; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t79, 0, t80, t75, 0, 0, t30, -t29, 0, -t18, -t19, t64 * t76 - t90, t66 * t76 + t15, t3, t18 * t56 + t3 * t53, -t16, t8, t23, -t22, 0, t43 * t20 + t27 * t29 + t9, -t43 * t21 - t28 * t29 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t81, -0.2e1 * t93, -0.2e1 * t56 * t66, 0.2e1 * t56 * t64, 0.2e1 * t85, t53 ^ 2 * t83 + t56 ^ 2, t38, t31, 0, 0, 0, t40 * t96, t41 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, -t18, -t19, t64 * t77 - t90, t66 * t77 + t15, t3, -t18 * pkin(4) + qJ(5) * t3, -t16, t8, t23, -t22, 0, t54 * t20 + t32 * t29 + t9, -t54 * t21 - t33 * t29 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t81, -t93, t92 * t66, -t92 * t64, t84 + t85, -t56 * pkin(4) + t53 * t84, t38, t31, 0, 0, 0, t86 * t40, t86 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t66, -0.2e1 * pkin(4) * t64, 0.2e1 * t84, qJ(5) ^ 2 * t83 + pkin(4) ^ 2, t38, t31, 0, 0, 0, t40 * t95, t41 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t88, 0, t18, 0, 0, 0, 0, 0, t20, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t64, 0, t56, 0, 0, 0, 0, 0, t40, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t64, 0, -pkin(4), 0, 0, 0, 0, 0, t40, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t20, t29, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, 0, t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, 0, t32, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t12;
