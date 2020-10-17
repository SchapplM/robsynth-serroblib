% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR14_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:59:49
% EndTime: 2019-05-06 16:59:53
% DurationCPUTime: 1.10s
% Computational Cost: add. (739->133), mult. (1684->257), div. (0->0), fcn. (1806->8), ass. (0->93)
t103 = -2 * pkin(2);
t60 = cos(pkin(6));
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t59 = sin(pkin(6));
t66 = cos(qJ(2));
t91 = t59 * t66;
t29 = t60 * t65 - t62 * t91;
t102 = -0.2e1 * t29;
t34 = pkin(4) * t62 - qJ(5) * t65 + qJ(3);
t101 = -0.2e1 * t34;
t100 = 0.2e1 * t59;
t99 = 0.2e1 * qJ(3);
t98 = 0.2e1 * qJ(5);
t68 = -pkin(2) - pkin(9);
t97 = pkin(4) + pkin(10);
t63 = sin(qJ(2));
t96 = pkin(1) * t63;
t95 = pkin(1) * t66;
t28 = t60 * t62 + t65 * t91;
t46 = t59 * t63;
t61 = sin(qJ(6));
t64 = cos(qJ(6));
t18 = t28 * t61 + t46 * t64;
t94 = t18 * t64;
t93 = t29 * t61;
t53 = t59 ^ 2;
t92 = t53 * t66;
t90 = t60 * t66;
t47 = t61 * t62;
t89 = t61 * t65;
t88 = t61 * t97;
t87 = t62 * t18;
t86 = t64 * t61;
t48 = t64 * t62;
t49 = t64 * t65;
t85 = t64 * t97;
t84 = t65 * t62;
t83 = t65 * t68;
t41 = pkin(8) * t46;
t75 = -pkin(2) - t95;
t14 = pkin(3) * t46 + t41 + (-pkin(9) + t75) * t60;
t74 = -qJ(3) * t63 - pkin(1);
t20 = (t66 * t68 + t74) * t59;
t9 = t14 * t62 + t20 * t65;
t31 = pkin(8) * t91 + t60 * t96;
t55 = t62 ^ 2;
t57 = t65 ^ 2;
t40 = t55 + t57;
t82 = qJ(5) * t62;
t81 = 0.2e1 * t46;
t80 = 0.2e1 * t84;
t79 = pkin(4) * t46;
t78 = t62 * t46;
t77 = t68 * t46;
t51 = t60 * qJ(3);
t21 = -t51 - t31;
t76 = qJ(5) * t46;
t8 = t65 * t14 - t20 * t62;
t73 = t62 * t77;
t72 = t65 * t77;
t19 = pkin(3) * t91 - t21;
t6 = -t76 - t9;
t7 = -t8 - t79;
t71 = -t6 * t62 - t7 * t65;
t37 = pkin(4) * t65 + t82;
t70 = t65 * t97 + t82;
t69 = -t29 * qJ(5) + t19;
t56 = t64 ^ 2;
t54 = t61 ^ 2;
t50 = t62 * t68;
t45 = t53 * t63 ^ 2;
t38 = t65 * t46;
t36 = (pkin(5) - t68) * t65;
t35 = -pkin(5) * t62 + t50;
t33 = t40 * t68;
t32 = pkin(10) * t62 + t34;
t30 = pkin(1) * t90 - t41;
t27 = t29 ^ 2;
t26 = t29 * t65;
t25 = t29 * t64;
t23 = t60 * t75 + t41;
t22 = (-pkin(2) * t66 + t74) * t59;
t17 = -t28 * t64 + t46 * t61;
t13 = t32 * t64 + t36 * t61;
t12 = -t32 * t61 + t36 * t64;
t10 = pkin(4) * t28 + t69;
t5 = t28 * t97 + t69;
t4 = -pkin(5) * t28 - t6;
t3 = t29 * pkin(5) - t46 * t97 - t8;
t2 = t3 * t61 + t5 * t64;
t1 = t3 * t64 - t5 * t61;
t11 = [1, 0, 0, t45, 0.2e1 * t63 * t92, t60 * t81, t90 * t100, t60 ^ 2, 0.2e1 * pkin(1) * t92 + 0.2e1 * t30 * t60, -0.2e1 * t31 * t60 - 0.2e1 * t53 * t96 (-t21 * t66 + t23 * t63) * t100, 0.2e1 * t22 * t91 + 0.2e1 * t23 * t60, -0.2e1 * t21 * t60 - 0.2e1 * t22 * t46, t21 ^ 2 + t22 ^ 2 + t23 ^ 2, t27, t28 * t102, t29 * t81, -0.2e1 * t28 * t46, t45, 0.2e1 * t19 * t28 + 0.2e1 * t46 * t8, 0.2e1 * t19 * t29 - 0.2e1 * t46 * t9, 0.2e1 * t28 * t6 + 0.2e1 * t29 * t7, -0.2e1 * t10 * t28 + 0.2e1 * t46 * t7, -0.2e1 * t10 * t29 - 0.2e1 * t46 * t6, t10 ^ 2 + t6 ^ 2 + t7 ^ 2, t18 ^ 2, -0.2e1 * t18 * t17, 0.2e1 * t18 * t29, t17 * t102, t27, 0.2e1 * t1 * t29 + 0.2e1 * t17 * t4, 0.2e1 * t18 * t4 - 0.2e1 * t2 * t29; 0, 0, 0, 0, 0, t46, t91, t60, t30, -t31 (-pkin(2) * t63 + qJ(3) * t66) * t59, t41 + (t103 - t95) * t60, 0.2e1 * t51 + t31, -pkin(2) * t23 - qJ(3) * t21, t26, -t28 * t65 - t29 * t62, t38, -t78, 0, qJ(3) * t28 + t19 * t62 + t72, qJ(3) * t29 + t19 * t65 - t73 (-t29 * t68 + t7) * t65 + (-t28 * t68 + t6) * t62, -t10 * t62 - t28 * t34 - t72, -t10 * t65 - t29 * t34 + t73, t10 * t34 + t68 * t71, t61 * t87 (-t17 * t61 + t94) * t62, t18 * t65 + t29 * t47, -t17 * t65 + t29 * t48, t26, t1 * t65 + t12 * t29 + t17 * t35 - t4 * t48, -t13 * t29 + t18 * t35 - t2 * t65 + t4 * t47; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t103, t99 (pkin(2) ^ 2) + qJ(3) ^ 2, t57, -0.2e1 * t84, 0, 0, 0, t62 * t99, t65 * t99, -0.2e1 * t33, t62 * t101, t65 * t101, t40 * (t68 ^ 2) + t34 ^ 2, t54 * t55, 0.2e1 * t55 * t86, t61 * t80, t64 * t80, t57, 0.2e1 * t12 * t65 - 0.2e1 * t35 * t48, -0.2e1 * t13 * t65 + 0.2e1 * t35 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t60, 0, t23, 0, 0, 0, 0, 0, t38, -t78, -t28 * t62 - t26, -t38, t78, t71, 0, 0, 0, 0, 0, t17 * t62 - t29 * t49, t29 * t89 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, -t40, 0, 0, t33, 0, 0, 0, 0, 0, -t40 * t64, t40 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, t46, t8, -t9, -pkin(4) * t29 - qJ(5) * t28, -t8 - 0.2e1 * t79, 0.2e1 * t76 + t9, -pkin(4) * t7 - qJ(5) * t6, t94, -t17 * t64 - t18 * t61, t25, -t93, 0, qJ(5) * t17 - t29 * t85 + t4 * t61, qJ(5) * t18 + t29 * t88 + t4 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t62, 0, t83, -t50, -t37, -t83, t50, t37 * t68, t61 * t48 (-t54 + t56) * t62, t49, -t89, 0, t35 * t61 - t64 * t70, t35 * t64 + t61 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t62, 0, -t65, t62, t37, 0, 0, 0, 0, 0, t47, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t98, pkin(4) ^ 2 + qJ(5) ^ 2, t56, -0.2e1 * t86, 0, 0, 0, t61 * t98, t64 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t46, 0, t7, 0, 0, 0, 0, 0, t25, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, -t83, 0, 0, 0, 0, 0, t49, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, t29, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t48, t65, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t61, 0, -t85, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t11;
