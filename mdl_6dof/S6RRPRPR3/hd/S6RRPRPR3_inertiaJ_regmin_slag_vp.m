% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR3
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
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:17:54
% EndTime: 2019-05-06 13:17:56
% DurationCPUTime: 0.64s
% Computational Cost: add. (1313->110), mult. (2540->217), div. (0->0), fcn. (3070->10), ass. (0->75)
t65 = sin(pkin(11));
t67 = cos(pkin(11));
t70 = sin(qJ(4));
t72 = cos(qJ(4));
t50 = -t65 * t70 + t67 * t72;
t68 = cos(pkin(10));
t61 = -t68 * pkin(2) - pkin(3);
t55 = -t72 * pkin(4) + t61;
t36 = -t50 * pkin(5) + t55;
t91 = 0.2e1 * t36;
t66 = sin(pkin(10));
t71 = sin(qJ(2));
t87 = cos(qJ(2));
t51 = t66 * t71 - t68 * t87;
t90 = -0.2e1 * t51;
t89 = 0.2e1 * t51;
t88 = pkin(4) * t65;
t86 = cos(qJ(6));
t52 = t65 * t72 + t67 * t70;
t69 = sin(qJ(6));
t32 = t69 * t50 + t86 * t52;
t85 = t32 * t51;
t84 = t70 * t51;
t53 = t66 * t87 + t68 * t71;
t83 = t70 * t53;
t82 = t70 * t72;
t56 = (-qJ(3) - pkin(7)) * t71;
t76 = t87 * pkin(7);
t57 = t87 * qJ(3) + t76;
t35 = t66 * t56 + t68 * t57;
t81 = t72 * t35;
t80 = t72 * t53;
t62 = -t87 * pkin(2) - pkin(1);
t30 = t51 * pkin(3) - t53 * pkin(8) + t62;
t18 = t72 * t30 - t70 * t35;
t79 = qJ(5) * t53;
t12 = t51 * pkin(4) - t72 * t79 + t18;
t14 = t81 + (t30 - t79) * t70;
t7 = t65 * t12 + t67 * t14;
t59 = t66 * pkin(2) + pkin(8);
t78 = qJ(5) + t59;
t46 = t78 * t70;
t47 = t78 * t72;
t29 = -t65 * t46 + t67 * t47;
t77 = 0.2e1 * t87;
t26 = -t65 * t83 + t67 * t80;
t6 = t67 * t12 - t65 * t14;
t4 = t51 * pkin(5) - t26 * pkin(9) + t6;
t25 = t52 * t53;
t5 = -t25 * pkin(9) + t7;
t1 = t86 * t4 - t69 * t5;
t28 = -t67 * t46 - t65 * t47;
t33 = -t68 * t56 + t66 * t57;
t23 = pkin(4) * t83 + t33;
t75 = -t51 * t59 + t53 * t61;
t2 = t69 * t4 + t86 * t5;
t64 = t72 ^ 2;
t63 = t70 ^ 2;
t60 = t67 * pkin(4) + pkin(5);
t49 = t53 ^ 2;
t48 = t51 ^ 2;
t45 = t72 * t51;
t43 = t69 * t60 + t86 * t88;
t42 = t86 * t60 - t69 * t88;
t31 = -t86 * t50 + t69 * t52;
t22 = t31 * t51;
t21 = t50 * pkin(9) + t29;
t20 = -t52 * pkin(9) + t28;
t19 = t70 * t30 + t81;
t17 = t25 * pkin(5) + t23;
t16 = -t69 * t25 + t86 * t26;
t15 = t86 * t25 + t69 * t26;
t9 = t69 * t20 + t86 * t21;
t8 = t86 * t20 - t69 * t21;
t3 = [1, 0, 0, t71 ^ 2, t71 * t77, 0, 0, 0, pkin(1) * t77, -0.2e1 * pkin(1) * t71, 0.2e1 * t33 * t53 - 0.2e1 * t35 * t51, t33 ^ 2 + t35 ^ 2 + t62 ^ 2, t64 * t49, -0.2e1 * t49 * t82, t80 * t89, t83 * t90, t48, 0.2e1 * t18 * t51 + 0.2e1 * t33 * t83, -0.2e1 * t19 * t51 + 0.2e1 * t33 * t80, -0.2e1 * t7 * t25 - 0.2e1 * t6 * t26, t23 ^ 2 + t6 ^ 2 + t7 ^ 2, t16 ^ 2, -0.2e1 * t16 * t15, t16 * t89, t15 * t90, t48, 0.2e1 * t1 * t51 + 0.2e1 * t17 * t15, 0.2e1 * t17 * t16 - 0.2e1 * t2 * t51; 0, 0, 0, 0, 0, t71, t87, 0, -t71 * pkin(7), -t76 (-t51 * t66 - t53 * t68) * pkin(2) (-t33 * t68 + t35 * t66) * pkin(2), t70 * t80 (-t63 + t64) * t53, t84, t45, 0, -t33 * t72 + t75 * t70, t33 * t70 + t75 * t72, -t29 * t25 - t28 * t26 + t7 * t50 - t6 * t52, t23 * t55 + t6 * t28 + t7 * t29, t16 * t32, -t32 * t15 - t16 * t31, t85, -t22, 0, t36 * t15 + t17 * t31 + t8 * t51, t36 * t16 + t17 * t32 - t9 * t51; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t66 ^ 2 + t68 ^ 2) * pkin(2) ^ 2, t63, 0.2e1 * t82, 0, 0, 0, -0.2e1 * t61 * t72, 0.2e1 * t61 * t70, -0.2e1 * t28 * t52 + 0.2e1 * t29 * t50, t28 ^ 2 + t29 ^ 2 + t55 ^ 2, t32 ^ 2, -0.2e1 * t32 * t31, 0, 0, 0, t31 * t91, t32 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, t45, -t84, -t52 * t25 - t50 * t26, t6 * t50 + t7 * t52, 0, 0, 0, 0, 0, -t22, -t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t50 + t29 * t52, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t50 ^ 2 + t52 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t83, t51, t18, -t19 (-t25 * t65 - t26 * t67) * pkin(4) (t6 * t67 + t65 * t7) * pkin(4), 0, 0, t16, -t15, t51, t42 * t51 + t1, -t43 * t51 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t72, 0, -t70 * t59, -t72 * t59 (t50 * t65 - t52 * t67) * pkin(4) (t28 * t67 + t29 * t65) * pkin(4), 0, 0, t32, -t31, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t70, 0 (t50 * t67 + t52 * t65) * pkin(4), 0, 0, 0, 0, 0, -t31, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t65 ^ 2 + t67 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t42, -0.2e1 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, t15, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, t31, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t51, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t42, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
