% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:40:41
% EndTime: 2019-05-06 22:40:46
% DurationCPUTime: 0.98s
% Computational Cost: add. (1634->135), mult. (3466->235), div. (0->0), fcn. (4119->10), ass. (0->85)
t71 = sin(pkin(11));
t72 = cos(pkin(11));
t75 = sin(qJ(4));
t79 = cos(qJ(4));
t52 = t75 * t71 - t79 * t72;
t53 = t79 * t71 + t75 * t72;
t74 = sin(qJ(5));
t78 = cos(qJ(5));
t33 = t78 * t52 + t74 * t53;
t62 = -t72 * pkin(3) - pkin(2);
t43 = t52 * pkin(4) + t62;
t24 = t33 * pkin(5) + t43;
t101 = 0.2e1 * t24;
t100 = 0.2e1 * t43;
t99 = 0.2e1 * t62;
t80 = cos(qJ(2));
t98 = -0.2e1 * t80;
t97 = 0.2e1 * t80;
t96 = pkin(7) * t71;
t73 = sin(qJ(6));
t95 = t73 * pkin(5);
t94 = t74 * pkin(4);
t76 = sin(qJ(2));
t64 = t76 * pkin(7);
t44 = t53 * t76;
t45 = t52 * t76;
t27 = t78 * t44 - t74 * t45;
t55 = -t80 * pkin(2) - t76 * qJ(3) - pkin(1);
t50 = t72 * t55;
t88 = t72 * t76;
t35 = -pkin(8) * t88 + t50 + (-pkin(3) - t96) * t80;
t90 = t80 * pkin(7);
t41 = t71 * t55 + t72 * t90;
t89 = t71 * t76;
t37 = -pkin(8) * t89 + t41;
t22 = t79 * t35 - t75 * t37;
t92 = t80 * pkin(4);
t15 = t45 * pkin(9) + t22 - t92;
t23 = t75 * t35 + t79 * t37;
t18 = -t44 * pkin(9) + t23;
t87 = t78 * t18;
t9 = t74 * t15 + t87;
t7 = -t27 * pkin(10) + t9;
t77 = cos(qJ(6));
t93 = t77 * t7;
t91 = t80 * pkin(5);
t86 = pkin(8) + qJ(3);
t54 = pkin(3) * t89 + t64;
t85 = t71 ^ 2 + t72 ^ 2;
t84 = t77 * t94;
t28 = -t74 * t44 - t78 * t45;
t8 = t78 * t15 - t74 * t18;
t6 = -t28 * pkin(10) + t8 - t91;
t1 = t77 * t6 - t73 * t7;
t56 = t86 * t71;
t57 = t86 * t72;
t38 = -t79 * t56 - t75 * t57;
t29 = -t53 * pkin(9) + t38;
t39 = -t75 * t56 + t79 * t57;
t30 = -t52 * pkin(9) + t39;
t16 = t78 * t29 - t74 * t30;
t66 = t78 * pkin(4);
t63 = t66 + pkin(5);
t46 = t77 * t63 - t73 * t94;
t36 = t44 * pkin(4) + t54;
t2 = t73 * t6 + t93;
t83 = -pkin(2) * t76 + qJ(3) * t80;
t17 = t74 * t29 + t78 * t30;
t40 = -t71 * t90 + t50;
t82 = -t40 * t71 + t41 * t72;
t70 = t80 ^ 2;
t69 = t76 ^ 2;
t65 = t77 * pkin(5);
t47 = t73 * t63 + t84;
t34 = -t74 * t52 + t78 * t53;
t21 = -t73 * t33 + t77 * t34;
t20 = t77 * t33 + t73 * t34;
t19 = t27 * pkin(5) + t36;
t14 = -t73 * t27 + t77 * t28;
t13 = t77 * t27 + t73 * t28;
t11 = -t33 * pkin(10) + t17;
t10 = -t34 * pkin(10) + t16;
t4 = t73 * t10 + t77 * t11;
t3 = t77 * t10 - t73 * t11;
t5 = [1, 0, 0, t69, t76 * t97, 0, 0, 0, pkin(1) * t97, -0.2e1 * pkin(1) * t76, -0.2e1 * t40 * t80 + 0.2e1 * t69 * t96, 0.2e1 * t69 * pkin(7) * t72 + 0.2e1 * t41 * t80, 0.2e1 * (-t40 * t72 - t41 * t71) * t76, t69 * pkin(7) ^ 2 + t40 ^ 2 + t41 ^ 2, t45 ^ 2, 0.2e1 * t45 * t44, -t45 * t98, t44 * t97, t70, -0.2e1 * t22 * t80 + 0.2e1 * t54 * t44, 0.2e1 * t23 * t80 - 0.2e1 * t54 * t45, t28 ^ 2, -0.2e1 * t28 * t27, t28 * t98, t27 * t97, t70, 0.2e1 * t36 * t27 - 0.2e1 * t8 * t80, 0.2e1 * t36 * t28 + 0.2e1 * t9 * t80, t14 ^ 2, -0.2e1 * t14 * t13, t14 * t98, t13 * t97, t70, -0.2e1 * t1 * t80 + 0.2e1 * t19 * t13, 0.2e1 * t19 * t14 + 0.2e1 * t2 * t80; 0, 0, 0, 0, 0, t76, t80, 0, -t64, -t90, -pkin(7) * t88 + t83 * t71, pkin(7) * t89 + t83 * t72, t82, -pkin(2) * t64 + t82 * qJ(3), -t45 * t53, -t53 * t44 + t45 * t52, -t53 * t80, t52 * t80, 0, -t38 * t80 + t62 * t44 + t54 * t52, t39 * t80 - t62 * t45 + t54 * t53, t28 * t34, -t34 * t27 - t28 * t33, -t34 * t80, t33 * t80, 0, -t16 * t80 + t43 * t27 + t36 * t33, t17 * t80 + t43 * t28 + t36 * t34, t14 * t21, -t21 * t13 - t14 * t20, -t21 * t80, t20 * t80, 0, t24 * t13 + t19 * t20 - t3 * t80, t24 * t14 + t19 * t21 + t4 * t80; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t72, -0.2e1 * pkin(2) * t71, 0.2e1 * t85 * qJ(3), t85 * qJ(3) ^ 2 + pkin(2) ^ 2, t53 ^ 2, -0.2e1 * t53 * t52, 0, 0, 0, t52 * t99, t53 * t99, t34 ^ 2, -0.2e1 * t34 * t33, 0, 0, 0, t33 * t100, t34 * t100, t21 ^ 2, -0.2e1 * t21 * t20, 0, 0, 0, t20 * t101, t21 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t88, 0, t64, 0, 0, 0, 0, 0, t44, -t45, 0, 0, 0, 0, 0, t27, t28, 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t71, 0, -pkin(2), 0, 0, 0, 0, 0, t52, t53, 0, 0, 0, 0, 0, t33, t34, 0, 0, 0, 0, 0, t20, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t44, -t80, t22, -t23, 0, 0, t28, -t27, -t80, -t78 * t92 + t8, -t87 + (-t15 + t92) * t74, 0, 0, t14, -t13, -t80, -t46 * t80 + t1, t47 * t80 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t52, 0, t38, -t39, 0, 0, t34, -t33, 0, t16, -t17, 0, 0, t21, -t20, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t66, -0.2e1 * t94, 0, 0, 0, 0, 1, 0.2e1 * t46, -0.2e1 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, -t80, t8, -t9, 0, 0, t14, -t13, -t80, -t77 * t91 + t1, -t93 + (-t6 + t91) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t33, 0, t16, -t17, 0, 0, t21, -t20, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t66, -t94, 0, 0, 0, 0, 1, t46 + t65, -t84 + (-pkin(5) - t63) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t65, -0.2e1 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, -t80, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t46, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t65, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
