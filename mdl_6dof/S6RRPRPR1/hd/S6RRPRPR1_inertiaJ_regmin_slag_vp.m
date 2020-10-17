% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:49:39
% EndTime: 2019-05-06 12:49:42
% DurationCPUTime: 0.83s
% Computational Cost: add. (1421->108), mult. (2706->194), div. (0->0), fcn. (3323->10), ass. (0->82)
t69 = sin(pkin(11));
t71 = cos(pkin(11));
t84 = t69 ^ 2 + t71 ^ 2;
t85 = t84 * qJ(5);
t70 = sin(pkin(10));
t72 = cos(pkin(10));
t75 = sin(qJ(2));
t77 = cos(qJ(2));
t50 = t70 * t77 + t72 * t75;
t74 = sin(qJ(4));
t80 = t70 * t75 - t72 * t77;
t92 = cos(qJ(4));
t29 = t74 * t50 + t92 * t80;
t100 = -0.2e1 * t29;
t61 = t72 * pkin(2) + pkin(3);
t95 = t70 * pkin(2);
t45 = t92 * t61 - t74 * t95;
t44 = -pkin(4) - t45;
t94 = t71 * pkin(5);
t38 = t44 - t94;
t99 = 0.2e1 * t38;
t63 = -t77 * pkin(2) - pkin(1);
t41 = t80 * pkin(3) + t63;
t98 = 0.2e1 * t41;
t62 = -pkin(4) - t94;
t97 = 0.2e1 * t62;
t96 = 0.2e1 * t77;
t93 = pkin(4) - t44;
t88 = -qJ(3) - pkin(7);
t57 = t88 * t75;
t58 = t88 * t77;
t36 = t72 * t57 + t70 * t58;
t27 = -t50 * pkin(8) + t36;
t37 = t70 * t57 - t72 * t58;
t28 = -t80 * pkin(8) + t37;
t18 = -t92 * t27 + t74 * t28;
t91 = t18 * t71;
t73 = sin(qJ(6));
t76 = cos(qJ(6));
t51 = t76 * t69 + t73 * t71;
t25 = t51 * t29;
t30 = t92 * t50 - t74 * t80;
t90 = t69 * t30;
t89 = t71 * t30;
t17 = t29 * pkin(4) - t30 * qJ(5) + t41;
t19 = t74 * t27 + t92 * t28;
t7 = t69 * t17 + t71 * t19;
t87 = t38 + t62;
t46 = -t74 * t61 - t92 * t95;
t43 = qJ(5) - t46;
t86 = t84 * t43;
t6 = t71 * t17 - t69 * t19;
t83 = t6 * t71 + t7 * t69;
t3 = -t6 * t69 + t7 * t71;
t82 = -pkin(4) * t30 - qJ(5) * t29;
t81 = -t29 * t43 + t30 * t44;
t49 = t73 * t69 - t76 * t71;
t66 = t71 * pkin(9);
t56 = t71 * qJ(5) + t66;
t55 = (-pkin(9) - qJ(5)) * t69;
t48 = t51 ^ 2;
t35 = t71 * t43 + t66;
t34 = (-pkin(9) - t43) * t69;
t33 = t73 * t55 + t76 * t56;
t32 = t76 * t55 - t73 * t56;
t31 = -0.2e1 * t51 * t49;
t24 = t49 * t29;
t23 = t73 * t34 + t76 * t35;
t22 = t76 * t34 - t73 * t35;
t21 = t49 * t30;
t20 = t51 * t30;
t16 = t18 * t69;
t14 = t21 * t51;
t11 = pkin(5) * t90 + t18;
t10 = t11 * t51;
t9 = t11 * t49;
t8 = -t51 * t20 + t21 * t49;
t5 = -pkin(9) * t90 + t7;
t4 = t29 * pkin(5) - pkin(9) * t89 + t6;
t2 = t73 * t4 + t76 * t5;
t1 = t76 * t4 - t73 * t5;
t12 = [1, 0, 0, t75 ^ 2, t75 * t96, 0, 0, 0, pkin(1) * t96, -0.2e1 * pkin(1) * t75, -0.2e1 * t36 * t50 - 0.2e1 * t37 * t80, t36 ^ 2 + t37 ^ 2 + t63 ^ 2, t30 ^ 2, t30 * t100, 0, 0, 0, t29 * t98, t30 * t98, 0.2e1 * t18 * t90 + 0.2e1 * t6 * t29, 0.2e1 * t18 * t89 - 0.2e1 * t7 * t29, -0.2e1 * t83 * t30, t18 ^ 2 + t6 ^ 2 + t7 ^ 2, t21 ^ 2, 0.2e1 * t21 * t20, t21 * t100, t20 * t100, t29 ^ 2, 0.2e1 * t1 * t29 + 0.2e1 * t11 * t20, -0.2e1 * t11 * t21 - 0.2e1 * t2 * t29; 0, 0, 0, 0, 0, t75, t77, 0, -t75 * pkin(7), -t77 * pkin(7) (-t72 * t50 - t70 * t80) * pkin(2) (t36 * t72 + t37 * t70) * pkin(2), 0, 0, t30, -t29, 0, -t18, -t19, t81 * t69 - t91, t81 * t71 + t16, t3, t18 * t44 + t3 * t43, -t14, t8, t25, -t24, 0, t38 * t20 + t22 * t29 + t9, -t38 * t21 - t23 * t29 + t10; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t70 ^ 2 + t72 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t45, 0.2e1 * t46, -0.2e1 * t44 * t71, 0.2e1 * t44 * t69, 0.2e1 * t86, t84 * t43 ^ 2 + t44 ^ 2, t48, t31, 0, 0, 0, t49 * t99, t51 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, t29, t30, t71 * t29, -t69 * t29, -t84 * t30, t83, 0, 0, 0, 0, 0, -t24, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, -t18, -t19, t82 * t69 - t91, t82 * t71 + t16, t3, -t18 * pkin(4) + t3 * qJ(5), -t14, t8, t25, -t24, 0, t62 * t20 + t32 * t29 + t9, -t62 * t21 - t33 * t29 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t45, t46, t93 * t71, -t93 * t69, t85 + t86, -t44 * pkin(4) + t43 * t85, t48, t31, 0, 0, 0, t87 * t49, t87 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t71, -0.2e1 * pkin(4) * t69, 0.2e1 * t85, t84 * qJ(5) ^ 2 + pkin(4) ^ 2, t48, t31, 0, 0, 0, t49 * t97, t51 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t89, 0, t18, 0, 0, 0, 0, 0, t20, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t69, 0, t44, 0, 0, 0, 0, 0, t49, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t69, 0, -pkin(4), 0, 0, 0, 0, 0, t49, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t20, t29, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t49, 0, t22, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t49, 0, t32, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t12;
