% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:28:07
% EndTime: 2019-05-07 18:28:10
% DurationCPUTime: 0.96s
% Computational Cost: add. (1110->157), mult. (2125->250), div. (0->0), fcn. (2244->6), ass. (0->84)
t56 = sin(qJ(4));
t59 = cos(qJ(4));
t58 = sin(qJ(2));
t60 = cos(qJ(3));
t81 = t60 * t58;
t57 = sin(qJ(3));
t84 = t57 * t58;
t26 = -t56 * t84 + t59 * t81;
t98 = -0.2e1 * t26;
t30 = t56 * t60 + t59 * t57;
t97 = -0.2e1 * t30;
t47 = -t60 * pkin(3) - pkin(2);
t96 = 0.2e1 * t47;
t95 = -0.2e1 * t58;
t61 = cos(qJ(2));
t94 = 0.2e1 * t61;
t62 = pkin(4) + pkin(5);
t93 = -pkin(9) - pkin(8);
t92 = pkin(2) * t60;
t91 = pkin(7) * t57;
t90 = t59 * pkin(3);
t89 = t61 * pkin(3);
t34 = t93 * t60;
t72 = t93 * t57;
t19 = -t56 * t34 - t59 * t72;
t88 = t19 * t61;
t20 = -t59 * t34 + t56 * t72;
t87 = t20 * t61;
t25 = t30 * t58;
t49 = t56 * pkin(3);
t42 = t49 + qJ(5);
t86 = t42 * t25;
t29 = t56 * t57 - t59 * t60;
t85 = t42 * t29;
t83 = t57 * t60;
t82 = t57 * t61;
t80 = t60 * t61;
t33 = -t61 * pkin(2) - t58 * pkin(8) - pkin(1);
t28 = t60 * t33;
t15 = -pkin(9) * t81 + t28 + (-pkin(3) - t91) * t61;
t74 = pkin(7) * t80;
t18 = t74 + (-pkin(9) * t58 + t33) * t57;
t79 = -t59 * t15 + t56 * t18;
t7 = t56 * t15 + t59 * t18;
t50 = t58 * pkin(7);
t32 = pkin(3) * t84 + t50;
t78 = qJ(5) * t25;
t77 = qJ(5) * t29;
t76 = t61 * qJ(5);
t75 = t58 * t94;
t51 = t61 * pkin(4);
t4 = t51 + t79;
t73 = -0.2e1 * t76 + t7;
t45 = pkin(4) + t90;
t66 = 0.2e1 * pkin(4);
t71 = t66 + t90;
t70 = t26 * qJ(6) - t4;
t3 = -t76 + t7;
t69 = t26 * qJ(5) - t32;
t68 = t30 * qJ(5) - t47;
t67 = (-qJ(5) - t42) * t61 + t7;
t65 = qJ(5) ^ 2;
t64 = 0.2e1 * qJ(5);
t55 = t61 ^ 2;
t54 = t60 ^ 2;
t53 = t58 ^ 2;
t52 = t57 ^ 2;
t43 = t64 + t49;
t40 = pkin(5) + t45;
t39 = t42 ^ 2;
t37 = t42 * qJ(5);
t35 = 0.2e1 * t42;
t23 = t25 * qJ(6);
t22 = t57 * t33 + t74;
t21 = -pkin(7) * t82 + t28;
t14 = t29 * pkin(4) - t68;
t11 = t29 * qJ(6) + t20;
t10 = -t30 * qJ(6) + t19;
t9 = -t62 * t29 + t68;
t8 = t25 * pkin(4) - t69;
t5 = -t62 * t25 + t69;
t2 = t23 + t3;
t1 = t61 * pkin(5) - t70;
t6 = [1, 0, 0, t53, t75, 0, 0, 0, pkin(1) * t94, pkin(1) * t95, t54 * t53, -0.2e1 * t53 * t83, t80 * t95, t57 * t75, t55, -0.2e1 * t21 * t61 + 0.2e1 * t53 * t91, 0.2e1 * t53 * pkin(7) * t60 + 0.2e1 * t22 * t61, t26 ^ 2, t25 * t98, t61 * t98, t25 * t94, t55, 0.2e1 * t32 * t25 + 0.2e1 * t61 * t79, 0.2e1 * t32 * t26 + 0.2e1 * t7 * t61, 0.2e1 * t8 * t25 + 0.2e1 * t4 * t61, -0.2e1 * t3 * t25 + 0.2e1 * t4 * t26, -0.2e1 * t8 * t26 - 0.2e1 * t3 * t61, t3 ^ 2 + t4 ^ 2 + t8 ^ 2, 0.2e1 * t1 * t61 - 0.2e1 * t5 * t25, -0.2e1 * t2 * t61 + 0.2e1 * t5 * t26, -0.2e1 * t1 * t26 + 0.2e1 * t2 * t25, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t58, t61, 0, -t50, -t61 * pkin(7), t57 * t81 (-t52 + t54) * t58, -t82, -t80, 0, -pkin(7) * t81 + (-pkin(2) * t58 + pkin(8) * t61) * t57, pkin(8) * t80 + (t91 - t92) * t58, t26 * t30, -t30 * t25 - t26 * t29, -t30 * t61, t29 * t61, 0, t47 * t25 + t32 * t29 + t88, t47 * t26 + t32 * t30 + t87, t14 * t25 + t8 * t29 + t88, t19 * t26 - t20 * t25 - t3 * t29 + t4 * t30, -t14 * t26 - t8 * t30 - t87, t8 * t14 + t4 * t19 + t3 * t20, t10 * t61 - t9 * t25 - t5 * t29, -t11 * t61 + t9 * t26 + t5 * t30, -t1 * t30 - t10 * t26 + t11 * t25 + t2 * t29, t1 * t10 + t2 * t11 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t52, 0.2e1 * t83, 0, 0, 0, 0.2e1 * t92, -0.2e1 * pkin(2) * t57, t30 ^ 2, t29 * t97, 0, 0, 0, t29 * t96, t30 * t96, 0.2e1 * t14 * t29, 0.2e1 * t19 * t30 - 0.2e1 * t20 * t29, t14 * t97, t14 ^ 2 + t19 ^ 2 + t20 ^ 2, -0.2e1 * t9 * t29, 0.2e1 * t9 * t30, -0.2e1 * t10 * t30 + 0.2e1 * t11 * t29, t10 ^ 2 + t11 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t84, -t61, t21, -t22, 0, 0, t26, -t25, -t61, -t59 * t89 - t79, t56 * t89 - t7, -t45 * t61 - t4, -t45 * t26 - t86, t67, t3 * t42 - t4 * t45 (-pkin(5) - t40) * t61 + t70, t23 + t67, t40 * t26 + t86, -t1 * t40 + t2 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t60, 0, -t57 * pkin(8), -t60 * pkin(8), 0, 0, t30, -t29, 0, -t19, -t20, -t19, -t45 * t30 - t85, t20, -t19 * t45 + t20 * t42, -t10, t11, t40 * t30 + t85, -t10 * t40 + t11 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t90, -0.2e1 * t49, 0.2e1 * t45, 0, t35, t45 ^ 2 + t39, 0.2e1 * t40, t35, 0, t40 ^ 2 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, -t61, -t79, -t7, -0.2e1 * t51 - t79, -pkin(4) * t26 - t78, t73, -t4 * pkin(4) + t3 * qJ(5) (-pkin(5) - t62) * t61 + t70, t23 + t73, t62 * t26 + t78, t2 * qJ(5) - t1 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, -t19, -t20, -t19, -pkin(4) * t30 - t77, t20, -t19 * pkin(4) + t20 * qJ(5), -t10, t11, t62 * t30 + t77, t11 * qJ(5) - t10 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t90, -t49, t71, 0, t43, t45 * pkin(4) + t37, 0.2e1 * pkin(5) + t71, t43, 0, t40 * t62 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t66, 0, t64, pkin(4) ^ 2 + t65, 0.2e1 * t62, t64, 0, t62 ^ 2 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t26, 0, t4, t61, 0, -t26, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t19, 0, 0, -t30, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t45, -1, 0, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), -1, 0, 0, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t26, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t30, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
