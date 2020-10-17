% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP13_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:13:51
% EndTime: 2019-05-06 19:13:55
% DurationCPUTime: 0.95s
% Computational Cost: add. (944->140), mult. (2127->279), div. (0->0), fcn. (2290->8), ass. (0->95)
t100 = -2 * pkin(2);
t51 = cos(pkin(6));
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t50 = sin(pkin(6));
t57 = cos(qJ(2));
t87 = t50 * t57;
t25 = t51 * t53 + t56 * t87;
t99 = -0.2e1 * t25;
t98 = 0.2e1 * t50;
t55 = cos(qJ(5));
t97 = 0.2e1 * t55;
t96 = 2 * qJ(3);
t58 = -pkin(2) - pkin(9);
t54 = sin(qJ(2));
t95 = pkin(1) * t54;
t94 = pkin(1) * t57;
t52 = sin(qJ(5));
t93 = t52 * pkin(5);
t41 = t50 * t54;
t36 = pkin(8) * t41;
t63 = -pkin(2) - t94;
t13 = pkin(3) * t41 + t36 + (-pkin(9) + t63) * t51;
t62 = -qJ(3) * t54 - pkin(1);
t19 = (t58 * t57 + t62) * t50;
t8 = t56 * t13 - t53 * t19;
t6 = -pkin(4) * t41 - t8;
t92 = t6 * t52;
t91 = t6 * t55;
t26 = t51 * t56 - t53 * t87;
t16 = t26 * t52 - t55 * t41;
t90 = t16 * t55;
t17 = t26 * t55 + t52 * t41;
t89 = t17 * t52;
t45 = t50 ^ 2;
t88 = t45 * t57;
t86 = t51 * t57;
t85 = t52 * t25;
t84 = t52 * t53;
t83 = t52 * t55;
t82 = t52 * t56;
t81 = t52 * t58;
t80 = t53 * t58;
t79 = t55 * t25;
t78 = t55 * t53;
t42 = t55 * t56;
t77 = t55 * t58;
t76 = t56 * t17;
t75 = t56 * t25;
t74 = t56 * t53;
t73 = t56 * t58;
t72 = -qJ(6) - pkin(10);
t28 = pkin(8) * t87 + t51 * t95;
t46 = t52 ^ 2;
t48 = t55 ^ 2;
t71 = t46 + t48;
t47 = t53 ^ 2;
t49 = t56 ^ 2;
t70 = -t47 - t49;
t69 = qJ(6) * t56;
t68 = 0.2e1 * t41;
t67 = -0.2e1 * t74;
t66 = t53 * t41;
t65 = t58 * t41;
t64 = t53 * t77;
t44 = t51 * qJ(3);
t22 = -t44 - t28;
t18 = pkin(3) * t87 - t22;
t11 = t25 * pkin(4) - t26 * pkin(10) + t18;
t9 = t53 * t13 + t56 * t19;
t7 = pkin(10) * t41 + t9;
t3 = t55 * t11 - t52 * t7;
t61 = -pkin(4) * t56 - pkin(10) * t53;
t1 = t25 * pkin(5) - t17 * qJ(6) + t3;
t4 = t52 * t11 + t55 * t7;
t2 = -t16 * qJ(6) + t4;
t60 = -t1 * t52 + t2 * t55;
t32 = t72 * t52;
t33 = t72 * t55;
t59 = -t32 * t52 - t33 * t55;
t43 = -t55 * pkin(5) - pkin(4);
t40 = t45 * t54 ^ 2;
t35 = t56 * t41;
t31 = t53 * pkin(4) - t56 * pkin(10) + qJ(3);
t30 = (-t58 + t93) * t56;
t29 = t55 * t31;
t27 = pkin(1) * t86 - t36;
t24 = t63 * t51 + t36;
t23 = (-pkin(2) * t57 + t62) * t50;
t21 = t52 * t31 + t64;
t20 = -t52 * t80 + t29;
t15 = t64 + (t31 - t69) * t52;
t12 = -t55 * t69 + t29 + (pkin(5) - t81) * t53;
t5 = t16 * pkin(5) + t6;
t10 = [1, 0, 0, t40, 0.2e1 * t54 * t88, t51 * t68, t86 * t98, t51 ^ 2, 0.2e1 * pkin(1) * t88 + 0.2e1 * t27 * t51, -0.2e1 * t28 * t51 - 0.2e1 * t45 * t95 (-t22 * t57 + t24 * t54) * t98, 0.2e1 * t23 * t87 + 0.2e1 * t24 * t51, -0.2e1 * t22 * t51 - 0.2e1 * t23 * t41, t22 ^ 2 + t23 ^ 2 + t24 ^ 2, t26 ^ 2, t26 * t99, t26 * t68, t41 * t99, t40, 0.2e1 * t18 * t25 + 0.2e1 * t8 * t41, 0.2e1 * t18 * t26 - 0.2e1 * t9 * t41, t17 ^ 2, -0.2e1 * t17 * t16, 0.2e1 * t17 * t25, t16 * t99, t25 ^ 2, 0.2e1 * t6 * t16 + 0.2e1 * t3 * t25, 0.2e1 * t6 * t17 - 0.2e1 * t4 * t25, -0.2e1 * t1 * t17 - 0.2e1 * t2 * t16, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t41, t87, t51, t27, -t28 (-pkin(2) * t54 + qJ(3) * t57) * t50, t36 + (t100 - t94) * t51, 0.2e1 * t44 + t28, -t24 * pkin(2) - t22 * qJ(3), t26 * t56, -t26 * t53 - t75, t35, -t66, 0, qJ(3) * t25 + t18 * t53 + t56 * t65, qJ(3) * t26 + t18 * t56 - t53 * t65, t55 * t76 (-t89 - t90) * t56, t17 * t53 + t55 * t75, -t16 * t53 - t52 * t75, t25 * t53, t20 * t25 + t3 * t53 + (-t16 * t58 + t92) * t56, -t21 * t25 - t4 * t53 + (-t17 * t58 + t91) * t56, -t12 * t17 - t15 * t16 + (-t1 * t55 - t2 * t52) * t56, t1 * t12 + t2 * t15 + t5 * t30; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t100, t96, pkin(2) ^ 2 + qJ(3) ^ 2, t49, t67, 0, 0, 0, t53 * t96, t56 * t96, t48 * t49, -0.2e1 * t49 * t83, t74 * t97, t52 * t67, t47, 0.2e1 * t20 * t53 - 0.2e1 * t49 * t81, -0.2e1 * t21 * t53 - 0.2e1 * t49 * t77, 0.2e1 * (-t12 * t55 - t15 * t52) * t56, t12 ^ 2 + t15 ^ 2 + t30 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t51, 0, t24, 0, 0, 0, 0, 0, t35, -t66, 0, 0, 0, 0, 0, -t56 * t16 - t25 * t84, -t25 * t78 - t76 (t89 - t90) * t53, -t5 * t56 + t60 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t52, t70 * t55, 0, -t30 * t56 + (-t12 * t52 + t15 * t55) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t47 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t41, t8, -t9, t89, -t52 * t16 + t17 * t55, t85, t79, 0, -pkin(4) * t16 - pkin(10) * t85 - t91, -pkin(4) * t17 - pkin(10) * t79 + t92, t33 * t16 - t32 * t17 + t60, t1 * t32 - t2 * t33 + t5 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t53, 0, t73, -t80, t52 * t42 (-t46 + t48) * t56, t84, t78, 0, t61 * t52 + t55 * t73, -t52 * t73 + t61 * t55 (-t32 * t56 + t15) * t55 + (t33 * t56 - t12) * t52, t12 * t32 - t15 * t33 + t30 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t53, 0, 0, 0, 0, 0, t42, -t82, t71 * t53, -t56 * t43 + t53 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t46, 0.2e1 * t83, 0, 0, 0, pkin(4) * t97, -0.2e1 * pkin(4) * t52, 0.2e1 * t59, t32 ^ 2 + t33 ^ 2 + t43 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t25, t3, -t4, -pkin(5) * t17, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t82, t53, t20, -t21, -pkin(5) * t42, t12 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t78, 0, -pkin(5) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t55, 0, -t52 * pkin(10), -t55 * pkin(10), -t93, t32 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t10;
