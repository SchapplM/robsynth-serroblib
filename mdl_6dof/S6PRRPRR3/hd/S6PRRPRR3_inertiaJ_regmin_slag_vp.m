% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_inertiaJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:50:24
% EndTime: 2019-05-05 04:50:27
% DurationCPUTime: 0.81s
% Computational Cost: add. (924->127), mult. (2530->264), div. (0->0), fcn. (3050->14), ass. (0->91)
t52 = sin(pkin(13));
t53 = sin(pkin(7));
t55 = cos(pkin(13));
t60 = sin(qJ(3));
t64 = cos(qJ(3));
t37 = (t52 * t64 + t55 * t60) * t53;
t56 = cos(pkin(7));
t59 = sin(qJ(5));
t63 = cos(qJ(5));
t27 = t59 * t37 - t63 * t56;
t100 = -0.2e1 * t27;
t99 = 0.2e1 * t59;
t98 = pkin(2) * t60;
t62 = cos(qJ(6));
t97 = pkin(5) * t62;
t89 = t53 * t64;
t90 = t53 * t60;
t36 = t52 * t90 - t55 * t89;
t44 = t56 * t64 * pkin(2);
t73 = pkin(9) + qJ(4);
t29 = t56 * pkin(3) - t73 * t90 + t44;
t70 = t56 * t98;
t32 = t73 * t89 + t70;
t18 = t52 * t29 + t55 * t32;
t16 = t56 * pkin(10) + t18;
t42 = (-pkin(3) * t64 - pkin(2)) * t53;
t21 = t36 * pkin(4) - t37 * pkin(10) + t42;
t8 = -t59 * t16 + t63 * t21;
t5 = -t36 * pkin(5) - t8;
t58 = sin(qJ(6));
t96 = t5 * t58;
t95 = t5 * t62;
t28 = t63 * t37 + t59 * t56;
t20 = t62 * t28 + t58 * t36;
t94 = t20 * t58;
t93 = t20 * t63;
t46 = t52 * pkin(3) + pkin(10);
t92 = t46 * t58;
t48 = t53 ^ 2;
t91 = t48 * t64;
t54 = sin(pkin(6));
t65 = cos(qJ(2));
t88 = t54 * t65;
t87 = t56 * t65;
t86 = t58 * t27;
t85 = t58 * t59;
t84 = t58 * t62;
t83 = t58 * t63;
t82 = t59 * t27;
t81 = t59 * t36;
t80 = t59 * t46;
t61 = sin(qJ(2));
t79 = t60 * t61;
t78 = t62 * t27;
t77 = t62 * t59;
t76 = t62 * t63;
t19 = t58 * t28 - t62 * t36;
t75 = t63 * t19;
t74 = t63 * t46;
t72 = 0.2e1 * t53 * t56;
t71 = t63 * t99;
t69 = t58 * t82;
t68 = t27 * t77;
t47 = -t55 * pkin(3) - pkin(4);
t17 = t55 * t29 - t52 * t32;
t9 = t63 * t16 + t59 * t21;
t15 = -t56 * pkin(4) - t17;
t57 = cos(pkin(6));
t67 = t57 * t89 + (t64 * t87 - t79) * t54;
t51 = t62 ^ 2;
t50 = t59 ^ 2;
t49 = t58 ^ 2;
t41 = -t63 * pkin(5) - t59 * pkin(11) + t47;
t40 = pkin(9) * t89 + t70;
t39 = -pkin(9) * t90 + t44;
t38 = -t53 * t88 + t57 * t56;
t34 = t63 * t36;
t26 = t58 * t41 + t62 * t74;
t25 = t62 * t41 - t58 * t74;
t23 = t57 * t90 + (t60 * t87 + t61 * t64) * t54;
t14 = t55 * t23 + t52 * t67;
t12 = t52 * t23 - t55 * t67;
t11 = t63 * t14 + t38 * t59;
t10 = t59 * t14 - t38 * t63;
t7 = t27 * pkin(5) - t28 * pkin(11) + t15;
t6 = t36 * pkin(11) + t9;
t4 = t62 * t11 + t58 * t12;
t3 = -t58 * t11 + t62 * t12;
t2 = t58 * t7 + t62 * t6;
t1 = -t58 * t6 + t62 * t7;
t13 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 ^ 2 + t14 ^ 2 + t38 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t88, -t54 * t61, 0, 0, 0, 0, 0, -t54 * t56 * t79 + ((t53 * t57 + t54 * t87) * t56 - t38 * t53) * t64, -t23 * t56 + t38 * t90, t12 * t37 - t14 * t36, -t12 * t17 + t14 * t18 + t38 * t42, 0, 0, 0, 0, 0, -t10 * t36 + t12 * t27, -t11 * t36 + t12 * t28, 0, 0, 0, 0, 0, t10 * t19 + t3 * t27, t10 * t20 - t4 * t27; 0, 1, 0, 0, t48 * t60 ^ 2, 0.2e1 * t60 * t91, t60 * t72, t64 * t72, t56 ^ 2, 0.2e1 * pkin(2) * t91 + 0.2e1 * t39 * t56, -0.2e1 * t40 * t56 - 0.2e1 * t48 * t98, -0.2e1 * t17 * t37 - 0.2e1 * t18 * t36, t17 ^ 2 + t18 ^ 2 + t42 ^ 2, t28 ^ 2, t28 * t100, 0.2e1 * t28 * t36, t36 * t100, t36 ^ 2, 0.2e1 * t15 * t27 + 0.2e1 * t8 * t36, 0.2e1 * t15 * t28 - 0.2e1 * t9 * t36, t20 ^ 2, -0.2e1 * t20 * t19, 0.2e1 * t20 * t27, t19 * t100, t27 ^ 2, 0.2e1 * t1 * t27 + 0.2e1 * t5 * t19, -0.2e1 * t2 * t27 + 0.2e1 * t5 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t23, 0 (-t12 * t55 + t14 * t52) * pkin(3), 0, 0, 0, 0, 0, -t12 * t63, t12 * t59, 0, 0, 0, 0, 0, t10 * t85 - t3 * t63, t10 * t77 + t4 * t63; 0, 0, 0, 0, 0, 0, t90, t89, t56, t39, -t40 (-t36 * t52 - t37 * t55) * pkin(3) (t17 * t55 + t18 * t52) * pkin(3), t28 * t59, t28 * t63 - t82, t81, t34, 0, -t15 * t63 + t47 * t27 - t36 * t80, t15 * t59 + t47 * t28 - t36 * t74, t20 * t77 (-t19 * t62 - t94) * t59, t68 - t93, -t69 + t75, -t27 * t63, -t1 * t63 + t25 * t27 + (t19 * t46 + t96) * t59, t2 * t63 - t26 * t27 + (t20 * t46 + t95) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t52 ^ 2 + t55 ^ 2) * pkin(3) ^ 2, t50, t71, 0, 0, 0, -0.2e1 * t47 * t63, t47 * t99, t51 * t50, -0.2e1 * t50 * t84, -0.2e1 * t59 * t76, t58 * t71, t63 ^ 2, -0.2e1 * t25 * t63 + 0.2e1 * t50 * t92, 0.2e1 * t50 * t46 * t62 + 0.2e1 * t26 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, t34, -t81, 0, 0, 0, 0, 0, -t69 - t75, -t68 - t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t10 * t62, t10 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, t36, t8, -t9, t94, -t58 * t19 + t20 * t62, t86, t78, 0, -pkin(5) * t19 - pkin(11) * t86 - t95, -pkin(5) * t20 - pkin(11) * t78 + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t63, 0, -t80, -t74, t58 * t77 (-t49 + t51) * t59, -t83, -t76, 0, -t46 * t77 + (-pkin(5) * t59 + pkin(11) * t63) * t58, pkin(11) * t76 + (t92 - t97) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t59, 0, 0, 0, 0, 0, t76, -t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t49, 0.2e1 * t84, 0, 0, 0, 0.2e1 * t97, -0.2e1 * pkin(5) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t27, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t85, -t63, t25, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t62, 0, -t58 * pkin(11), -t62 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
