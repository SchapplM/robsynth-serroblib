% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:02:06
% EndTime: 2022-01-20 12:02:09
% DurationCPUTime: 0.57s
% Computational Cost: add. (612->105), mult. (1532->146), div. (0->0), fcn. (1204->8), ass. (0->91)
t66 = sin(qJ(2));
t69 = cos(qJ(3));
t102 = t66 * t69;
t65 = sin(qJ(3));
t70 = cos(qJ(2));
t89 = qJD(3) * t69;
t111 = ((t65 * t70 + t102) * qJD(2) + t66 * t89) * pkin(1);
t110 = qJD(4) + qJD(5);
t108 = -pkin(8) - pkin(9);
t68 = cos(qJ(4));
t107 = t68 * pkin(4);
t106 = t69 * pkin(2);
t58 = t70 * pkin(1) + pkin(2);
t37 = pkin(1) * t102 + t65 * t58 + pkin(8);
t105 = -pkin(9) - t37;
t56 = t65 * pkin(2) + pkin(8);
t104 = -pkin(9) - t56;
t55 = t65 * t66 * pkin(1);
t92 = pkin(1) * qJD(2);
t84 = t70 * t92;
t21 = -t58 * t89 - t69 * t84 + (qJD(2) + qJD(3)) * t55;
t64 = sin(qJ(4));
t103 = t64 * t21;
t101 = t68 * t21;
t90 = qJD(3) * t65;
t22 = t58 * t90 + t111;
t88 = t64 * qJD(4);
t60 = pkin(4) * t88;
t15 = t60 + t22;
t63 = sin(qJ(5));
t67 = cos(qJ(5));
t42 = t63 * t68 + t67 * t64;
t24 = t110 * t42;
t36 = -t69 * t58 - pkin(3) + t55;
t33 = t36 - t107;
t41 = t63 * t64 - t67 * t68;
t100 = t15 * t41 + t33 * t24;
t23 = t110 * t41;
t99 = t15 * t42 - t33 * t23;
t61 = t68 * qJD(4);
t98 = t22 * t64 + t36 * t61;
t83 = pkin(2) * t90;
t45 = t60 + t83;
t59 = -pkin(3) - t107;
t49 = t59 - t106;
t97 = t49 * t24 + t45 * t41;
t96 = -t49 * t23 + t45 * t42;
t95 = t59 * t24 + t41 * t60;
t94 = -t59 * t23 + t42 * t60;
t57 = -pkin(3) - t106;
t93 = t57 * t61 + t64 * t83;
t91 = pkin(4) * qJD(5);
t87 = pkin(3) * t88;
t86 = pkin(3) * t61;
t85 = t66 * t92;
t82 = pkin(2) * t89;
t81 = t63 * t91;
t80 = t67 * t91;
t78 = qJD(4) * t108;
t31 = t36 * t88;
t77 = -t22 * t68 + t31;
t76 = qJD(4) * t105;
t75 = qJD(4) * t104;
t74 = t64 * t82;
t73 = t68 * t82;
t46 = t57 * t88;
t71 = -t68 * t83 + t46;
t62 = t68 * pkin(9);
t54 = 0.2e1 * t64 * t61;
t51 = t68 * pkin(8) + t62;
t50 = t108 * t64;
t44 = t68 * t78;
t43 = t64 * t78;
t40 = 0.2e1 * (-t64 ^ 2 + t68 ^ 2) * qJD(4);
t39 = t68 * t56 + t62;
t38 = t104 * t64;
t30 = t68 * t75 - t74;
t29 = t64 * t75 + t73;
t26 = t68 * t37 + t62;
t25 = t105 * t64;
t14 = -0.2e1 * t42 * t23;
t9 = -t63 * t43 + t67 * t44 + (-t50 * t63 - t51 * t67) * qJD(5);
t8 = -t67 * t43 - t63 * t44 + (-t50 * t67 + t51 * t63) * qJD(5);
t7 = t68 * t76 + t103;
t6 = t64 * t76 - t101;
t5 = 0.2e1 * t23 * t41 - 0.2e1 * t42 * t24;
t4 = -t63 * t29 + t67 * t30 + (-t38 * t63 - t39 * t67) * qJD(5);
t3 = -t67 * t29 - t63 * t30 + (-t38 * t67 + t39 * t63) * qJD(5);
t2 = -t63 * t6 + t67 * t7 + (-t25 * t63 - t26 * t67) * qJD(5);
t1 = -t67 * t6 - t63 * t7 + (-t25 * t67 + t26 * t63) * qJD(5);
t10 = [0, 0, 0, 0, -0.2e1 * t85, -0.2e1 * t84, 0, -0.2e1 * t22, 0.2e1 * t21, t54, t40, 0, 0, 0, 0.2e1 * t77, 0.2e1 * t98, t14, t5, 0, 0, 0, 0.2e1 * t100, 0.2e1 * t99; 0, 0, 0, 0, -t85, -t84, 0, (-pkin(2) - t58) * t90 - t111, t21 - t82, t54, t40, 0, 0, 0, t31 + t46 + (-t22 - t83) * t68, t93 + t98, t14, t5, 0, 0, 0, t97 + t100, t96 + t99; 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t83, -0.2e1 * t82, t54, t40, 0, 0, 0, 0.2e1 * t71, 0.2e1 * t93, t14, t5, 0, 0, 0, 0.2e1 * t97, 0.2e1 * t96; 0, 0, 0, 0, 0, 0, 0, -t22, t21, t54, t40, 0, 0, 0, t77 - t87, -t86 + t98, t14, t5, 0, 0, 0, t95 + t100, t94 + t99; 0, 0, 0, 0, 0, 0, 0, -t83, -t82, t54, t40, 0, 0, 0, t71 - t87, -t86 + t93, t14, t5, 0, 0, 0, t95 + t97, t94 + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t40, 0, 0, 0, -0.2e1 * t87, -0.2e1 * t86, t14, t5, 0, 0, 0, 0.2e1 * t95, 0.2e1 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t88, 0, -t37 * t61 + t103, t37 * t88 + t101, 0, 0, -t23, -t24, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t88, 0, -t56 * t61 - t74, t56 * t88 - t73, 0, 0, -t23, -t24, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t88, 0, -pkin(8) * t61, pkin(8) * t88, 0, 0, -t23, -t24, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t81, -0.2e1 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
