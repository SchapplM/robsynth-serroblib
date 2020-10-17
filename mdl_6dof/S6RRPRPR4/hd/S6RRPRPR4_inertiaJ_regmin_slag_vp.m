% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:37:12
% EndTime: 2019-05-06 13:37:16
% DurationCPUTime: 0.89s
% Computational Cost: add. (1824->128), mult. (4422->274), div. (0->0), fcn. (5160->12), ass. (0->93)
t67 = sin(pkin(12));
t70 = cos(pkin(12));
t74 = sin(qJ(4));
t77 = cos(qJ(4));
t52 = t67 * t74 - t70 * t77;
t107 = t52 ^ 2;
t68 = sin(pkin(11));
t69 = sin(pkin(6));
t71 = cos(pkin(11));
t75 = sin(qJ(2));
t78 = cos(qJ(2));
t44 = (t68 * t78 + t71 * t75) * t69;
t72 = cos(pkin(6));
t35 = t44 * t74 - t72 * t77;
t36 = t44 * t77 + t72 * t74;
t22 = -t67 * t35 + t70 * t36;
t98 = t69 * t78;
t99 = t69 * t75;
t43 = t68 * t99 - t71 * t98;
t73 = sin(qJ(6));
t76 = cos(qJ(6));
t15 = t73 * t22 - t43 * t76;
t106 = -0.2e1 * t15;
t105 = -0.2e1 * t35;
t104 = 0.2e1 * t74;
t103 = pkin(1) * t75;
t58 = t72 * t78 * pkin(1);
t88 = pkin(8) + qJ(3);
t37 = t72 * pkin(2) - t88 * t99 + t58;
t85 = t72 * t103;
t40 = t88 * t98 + t85;
t26 = t68 * t37 + t71 * t40;
t24 = t72 * pkin(9) + t26;
t55 = (-pkin(2) * t78 - pkin(1)) * t69;
t28 = t43 * pkin(3) - t44 * pkin(9) + t55;
t13 = t77 * t24 + t74 * t28;
t11 = -t35 * qJ(5) + t13;
t12 = -t74 * t24 + t77 * t28;
t9 = t43 * pkin(4) - t36 * qJ(5) + t12;
t6 = t70 * t11 + t67 * t9;
t16 = t76 * t22 + t43 * t73;
t102 = t16 * t73;
t101 = t52 * t15;
t64 = t69 ^ 2;
t100 = t64 * t78;
t21 = t70 * t35 + t67 * t36;
t97 = t73 * t21;
t46 = t73 * t52;
t54 = t67 * t77 + t70 * t74;
t96 = t73 * t54;
t60 = t67 * pkin(4) + pkin(10);
t95 = t73 * t60;
t94 = t73 * t76;
t93 = t74 * t43;
t61 = t68 * pkin(2) + pkin(9);
t92 = t74 * t61;
t91 = t76 * t54;
t90 = t76 * t60;
t89 = t77 * t61;
t87 = qJ(5) + t61;
t86 = 0.2e1 * t69 * t72;
t84 = t21 * t96;
t83 = t21 * t91;
t63 = -t71 * pkin(2) - pkin(3);
t25 = t71 * t37 - t68 * t40;
t82 = t87 * t74;
t5 = -t67 * t11 + t70 * t9;
t62 = -t70 * pkin(4) - pkin(5);
t81 = -t52 * t60 + t54 * t62;
t56 = -t77 * pkin(4) + t63;
t23 = -t72 * pkin(3) - t25;
t17 = t35 * pkin(4) + t23;
t66 = t76 ^ 2;
t65 = t73 ^ 2;
t51 = t54 ^ 2;
t50 = t87 * t77;
t49 = pkin(8) * t98 + t85;
t48 = -pkin(8) * t99 + t58;
t47 = t76 * t52;
t41 = t77 * t43;
t32 = t70 * t50 - t67 * t82;
t30 = t67 * t50 + t70 * t82;
t29 = t52 * pkin(5) - t54 * pkin(10) + t56;
t20 = t76 * t21;
t19 = t73 * t29 + t76 * t32;
t18 = t76 * t29 - t73 * t32;
t14 = t16 * t52;
t7 = t21 * pkin(5) - t22 * pkin(10) + t17;
t4 = t43 * pkin(10) + t6;
t3 = -t43 * pkin(5) - t5;
t2 = t76 * t4 + t73 * t7;
t1 = -t73 * t4 + t76 * t7;
t8 = [1, 0, 0, t64 * t75 ^ 2, 0.2e1 * t75 * t100, t75 * t86, t78 * t86, t72 ^ 2, 0.2e1 * pkin(1) * t100 + 0.2e1 * t48 * t72, -0.2e1 * t64 * t103 - 0.2e1 * t49 * t72, -0.2e1 * t25 * t44 - 0.2e1 * t26 * t43, t25 ^ 2 + t26 ^ 2 + t55 ^ 2, t36 ^ 2, t36 * t105, 0.2e1 * t36 * t43, t43 * t105, t43 ^ 2, 0.2e1 * t12 * t43 + 0.2e1 * t23 * t35, -0.2e1 * t13 * t43 + 0.2e1 * t23 * t36, -0.2e1 * t6 * t21 - 0.2e1 * t5 * t22, t17 ^ 2 + t5 ^ 2 + t6 ^ 2, t16 ^ 2, t16 * t106, 0.2e1 * t16 * t21, t21 * t106, t21 ^ 2, 0.2e1 * t1 * t21 + 0.2e1 * t3 * t15, 0.2e1 * t3 * t16 - 0.2e1 * t2 * t21; 0, 0, 0, 0, 0, t99, t98, t72, t48, -t49 (-t43 * t68 - t44 * t71) * pkin(2) (t25 * t71 + t26 * t68) * pkin(2), t36 * t74, -t74 * t35 + t36 * t77, t93, t41, 0, -t23 * t77 + t63 * t35 - t43 * t92, t23 * t74 + t63 * t36 - t43 * t89, -t32 * t21 + t30 * t22 - t5 * t54 - t6 * t52, t17 * t56 - t5 * t30 + t6 * t32, t16 * t91 (-t15 * t76 - t102) * t54, t14 + t83, -t84 - t101, t21 * t52, t1 * t52 + t30 * t15 + t18 * t21 + t3 * t96, t30 * t16 - t19 * t21 - t2 * t52 + t3 * t91; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t68 ^ 2 + t71 ^ 2) * pkin(2) ^ 2, t74 ^ 2, t77 * t104, 0, 0, 0, -0.2e1 * t63 * t77, t63 * t104, 0.2e1 * t30 * t54 - 0.2e1 * t32 * t52, t30 ^ 2 + t32 ^ 2 + t56 ^ 2, t66 * t51, -0.2e1 * t51 * t94, 0.2e1 * t52 * t91, -0.2e1 * t52 * t96, t107, 0.2e1 * t18 * t52 + 0.2e1 * t30 * t96, -0.2e1 * t19 * t52 + 0.2e1 * t30 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, t41, -t93, -t54 * t21 + t52 * t22, -t5 * t52 + t6 * t54, 0, 0, 0, 0, 0, -t84 + t101, t14 - t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t52 + t32 * t54, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t51 + t107, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, t43, t12, -t13 (-t21 * t67 - t22 * t70) * pkin(4) (t5 * t70 + t6 * t67) * pkin(4), t102, -t73 * t15 + t16 * t76, t97, t20, 0, t62 * t15 - t21 * t95 - t3 * t76, t62 * t16 - t21 * t90 + t3 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t77, 0, -t92, -t89 (-t52 * t67 - t54 * t70) * pkin(4) (-t30 * t70 + t32 * t67) * pkin(4), t73 * t91 (-t65 + t66) * t54, t46, t47, 0, -t30 * t76 + t81 * t73, t30 * t73 + t81 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t74, 0 (-t52 * t70 + t54 * t67) * pkin(4), 0, 0, 0, 0, 0, -t47, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t67 ^ 2 + t70 ^ 2) * pkin(4) ^ 2, t65, 0.2e1 * t94, 0, 0, 0, -0.2e1 * t62 * t76, 0.2e1 * t62 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, t20, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, t47, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t21, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t96, t52, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t76, 0, -t95, -t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
