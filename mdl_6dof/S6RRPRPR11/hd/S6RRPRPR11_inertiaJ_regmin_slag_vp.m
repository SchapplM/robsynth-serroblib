% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:03:22
% EndTime: 2019-05-06 16:03:24
% DurationCPUTime: 0.69s
% Computational Cost: add. (815->103), mult. (1496->192), div. (0->0), fcn. (1668->8), ass. (0->80)
t66 = sin(pkin(10));
t67 = cos(pkin(10));
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t78 = t66 * t69 - t67 * t72;
t79 = t66 * t72 + t67 * t69;
t105 = (t66 * t79 - t67 * t78) * pkin(4);
t70 = sin(qJ(2));
t68 = sin(qJ(6));
t71 = cos(qJ(6));
t80 = -t68 * t78 + t71 * t79;
t104 = t80 * t70;
t86 = -t68 * t79 - t71 * t78;
t103 = t86 * t70;
t55 = t69 * pkin(4) + qJ(3);
t29 = pkin(5) * t79 + t55;
t102 = 0.2e1 * t29;
t101 = -0.2e1 * t70;
t73 = cos(qJ(2));
t100 = 0.2e1 * t73;
t99 = 0.2e1 * qJ(3);
t74 = -pkin(2) - pkin(8);
t98 = pkin(4) * t66;
t59 = t70 * pkin(7);
t48 = t70 * pkin(3) + t59;
t97 = t69 * t48;
t96 = t69 * t70;
t95 = t69 * t73;
t94 = t70 * t73;
t93 = t72 * t69;
t92 = t72 * t73;
t44 = t72 * t48;
t87 = -t70 * qJ(3) - pkin(1);
t40 = t74 * t73 + t87;
t85 = qJ(5) * t73 - t40;
t16 = t70 * pkin(4) + t85 * t69 + t44;
t18 = -t85 * t72 + t97;
t9 = t66 * t16 + t67 * t18;
t45 = (-qJ(5) + t74) * t69;
t57 = t72 * t74;
t46 = -t72 * qJ(5) + t57;
t27 = t67 * t45 + t66 * t46;
t60 = t73 * pkin(7);
t49 = t73 * pkin(3) + t60;
t63 = t70 ^ 2;
t65 = t73 ^ 2;
t91 = t63 + t65;
t90 = t73 * qJ(3);
t89 = -0.2e1 * t94;
t36 = pkin(4) * t92 + t49;
t88 = t78 ^ 2 + t79 ^ 2;
t31 = t79 * t73;
t8 = t67 * t16 - t66 * t18;
t4 = t70 * pkin(5) + t31 * pkin(9) + t8;
t30 = t78 * t73;
t5 = t30 * pkin(9) + t9;
t1 = t71 * t4 - t68 * t5;
t26 = -t66 * t45 + t67 * t46;
t2 = t68 * t4 + t71 * t5;
t84 = -t78 * t8 + t79 * t9;
t83 = -t70 * pkin(2) + t90;
t82 = -t26 * t78 + t27 * t79;
t77 = t70 * t74 + t90;
t64 = t72 ^ 2;
t62 = t69 ^ 2;
t56 = t72 * t70;
t54 = t67 * pkin(4) + pkin(5);
t47 = -t73 * pkin(2) + t87;
t33 = t68 * t54 + t71 * t98;
t32 = t71 * t54 - t68 * t98;
t25 = t72 * t40 + t97;
t24 = -t69 * t40 + t44;
t19 = -t30 * pkin(5) + t36;
t15 = t68 * t30 - t71 * t31;
t14 = -t71 * t30 - t68 * t31;
t13 = -pkin(9) * t79 + t27;
t12 = pkin(9) * t78 + t26;
t7 = t68 * t12 + t71 * t13;
t6 = t71 * t12 - t68 * t13;
t3 = [1, 0, 0, t63, 0.2e1 * t94, 0, 0, 0, pkin(1) * t100, pkin(1) * t101, 0.2e1 * t91 * pkin(7), t47 * t100, t47 * t101, t91 * pkin(7) ^ 2 + t47 ^ 2, t62 * t65, 0.2e1 * t65 * t93, t69 * t89, t72 * t89, t63, 0.2e1 * t24 * t70 + 0.2e1 * t49 * t92, -0.2e1 * t25 * t70 - 0.2e1 * t49 * t95, 0.2e1 * t9 * t30 + 0.2e1 * t8 * t31, t36 ^ 2 + t8 ^ 2 + t9 ^ 2, t15 ^ 2, -0.2e1 * t15 * t14, 0.2e1 * t15 * t70, t14 * t101, t63, 0.2e1 * t1 * t70 + 0.2e1 * t19 * t14, 0.2e1 * t19 * t15 - 0.2e1 * t2 * t70; 0, 0, 0, 0, 0, t70, t73, 0, -t59, -t60, t83, t59, t60, t83 * pkin(7), -t69 * t92 (t62 - t64) * t73, t56, -t96, 0, t49 * t69 + t77 * t72, t49 * t72 - t77 * t69, t26 * t31 + t27 * t30 - t84, t8 * t26 + t9 * t27 + t36 * t55, t15 * t86, -t14 * t86 - t15 * t80, t103, -t104, 0, t29 * t14 + t19 * t80 + t6 * t70, t29 * t15 + t19 * t86 - t7 * t70; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t99, pkin(2) ^ 2 + qJ(3) ^ 2, t64, -0.2e1 * t93, 0, 0, 0, t69 * t99, t72 * t99, -0.2e1 * t82, t26 ^ 2 + t27 ^ 2 + t55 ^ 2, t86 ^ 2, -0.2e1 * t86 * t80, 0, 0, 0, t80 * t102, t86 * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, 0, t59, 0, 0, 0, 0, 0, t56, -t96, t30 * t79 - t31 * t78, t84, 0, 0, 0, 0, 0, t103, -t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, -t88, t82, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t92, t70, t24, -t25 (t30 * t66 + t31 * t67) * pkin(4) (t66 * t9 + t67 * t8) * pkin(4), 0, 0, t15, -t14, t70, t32 * t70 + t1, -t33 * t70 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t69, 0, t57, -t69 * t74, -t105 (t26 * t67 + t27 * t66) * pkin(4), 0, 0, t86, -t80, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t69, 0, t105, 0, 0, 0, 0, 0, t86, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t66 ^ 2 + t67 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t32, -0.2e1 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, t80, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, t70, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t80, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t32, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
