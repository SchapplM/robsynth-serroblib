% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR15_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:08
% EndTime: 2019-12-31 20:43:12
% DurationCPUTime: 0.91s
% Computational Cost: add. (736->140), mult. (1771->264), div. (0->0), fcn. (1430->6), ass. (0->108)
t121 = pkin(3) + pkin(6);
t62 = sin(qJ(2));
t108 = t62 * qJ(3);
t65 = cos(qJ(2));
t66 = -pkin(2) - pkin(7);
t72 = -t65 * t66 + t108;
t29 = -pkin(1) - t72;
t43 = t121 * t62;
t61 = sin(qJ(4));
t35 = t61 * t43;
t64 = cos(qJ(4));
t112 = t64 * t29 + t35;
t59 = t65 ^ 2;
t82 = qJD(2) * (t62 ^ 2 - t59);
t56 = t61 ^ 2;
t111 = -t64 ^ 2 + t56;
t81 = t111 * qJD(4);
t125 = qJD(4) + qJD(5);
t44 = t121 * t65;
t101 = t44 * qJD(4);
t99 = t65 * qJD(3);
t124 = t72 * qJD(2) - t101 - t99;
t123 = t125 * t65;
t122 = 0.2e1 * qJD(3);
t120 = t62 * pkin(4);
t105 = qJD(4) * t64;
t106 = qJD(4) * t61;
t109 = qJ(3) * t65;
t100 = t62 * qJD(2);
t80 = pkin(2) * t100 - t62 * qJD(3);
t19 = (pkin(7) * t62 - t109) * qJD(2) + t80;
t53 = t65 * qJD(2);
t51 = pkin(6) * t53;
t38 = pkin(3) * t53 + t51;
t6 = -t43 * t105 + t29 * t106 - t64 * t19 - t61 * t38;
t88 = t64 * t100;
t104 = qJD(4) * t65;
t92 = t61 * t104;
t69 = t88 + t92;
t5 = t69 * pkin(8) - t6;
t63 = cos(qJ(5));
t119 = t63 * t5;
t118 = pkin(8) - t66;
t113 = t64 * t65;
t14 = -pkin(8) * t113 + t112;
t60 = sin(qJ(5));
t117 = t60 * t14;
t116 = t60 * t61;
t115 = t63 * t14;
t114 = t63 * t64;
t107 = qJD(2) * t44;
t103 = qJD(4) * t66;
t102 = qJD(5) * t60;
t98 = qJ(3) * qJD(4);
t97 = -0.2e1 * pkin(1) * qJD(2);
t96 = pkin(4) * t53;
t95 = pkin(4) * t102;
t94 = qJD(5) * t63 * pkin(4);
t93 = pkin(6) * t100;
t91 = t64 * t104;
t90 = t61 * t100;
t89 = t62 * t53;
t87 = t61 * t105;
t83 = -t61 * t19 + t64 * t38;
t84 = pkin(8) * t65 - t29;
t4 = (-pkin(8) * t61 * t62 + pkin(4) * t65) * qJD(2) + (t84 * t64 - t35) * qJD(4) + t83;
t86 = t63 * t4 - t60 * t5;
t40 = t118 * t64;
t36 = t64 * t43;
t13 = t84 * t61 + t120 + t36;
t85 = -t13 - t120;
t79 = t61 * t88;
t78 = -t65 * pkin(2) - t108;
t77 = t63 * t13 - t117;
t76 = t60 * t13 + t115;
t39 = t118 * t61;
t75 = -t63 * t39 - t60 * t40;
t74 = -t60 * t39 + t63 * t40;
t33 = t60 * t64 + t63 * t61;
t73 = -t114 + t116;
t16 = -t61 * t102 - t60 * t106 + t125 * t114;
t70 = -t16 * t62 - t33 * t53;
t37 = t121 * t100;
t68 = -t37 + (-t62 * t66 - t109) * qJD(4);
t67 = t78 * qJD(2) + t99;
t50 = t61 * pkin(4) + qJ(3);
t47 = pkin(4) * t105 + qJD(3);
t46 = 0.2e1 * t89;
t41 = -pkin(1) + t78;
t32 = qJD(4) * t40;
t31 = t118 * t106;
t26 = -t62 * t105 - t61 * t53;
t25 = -t62 * t106 + t64 * t53;
t24 = pkin(4) * t113 + t44;
t23 = -qJ(3) * t53 + t80;
t21 = t33 * t65;
t20 = -t63 * t113 + t65 * t116;
t17 = -pkin(4) * t92 + (-pkin(4) * t64 - t121) * t100;
t15 = t125 * t33;
t12 = -t15 * t62 - t53 * t73;
t11 = t33 * t123 - t60 * t90 + t63 * t88;
t10 = t33 * t100 + t73 * t123;
t9 = -t75 * qJD(5) + t63 * t31 + t60 * t32;
t8 = t74 * qJD(5) - t60 * t31 + t63 * t32;
t7 = -t112 * qJD(4) + t83;
t2 = -t76 * qJD(5) + t86;
t1 = -t77 * qJD(5) - t60 * t4 - t119;
t3 = [0, 0, 0, t46, -0.2e1 * t82, 0, 0, 0, t62 * t97, t65 * t97, 0, -0.2e1 * t41 * t100 + 0.2e1 * t23 * t65, -0.2e1 * t23 * t62 - 0.2e1 * t41 * t53, 0.2e1 * t41 * t23, -0.2e1 * t56 * t89 + 0.2e1 * t59 * t87, -0.2e1 * t59 * t81 - 0.4e1 * t65 * t79, 0.2e1 * t61 * t82 - 0.2e1 * t62 * t91, 0.2e1 * t62 * t92 + 0.2e1 * t64 * t82, t46, 0.2e1 * (-t64 * t107 + t7) * t62 + 0.2e1 * ((-t61 * t29 + t36) * qJD(2) - t37 * t64 - t61 * t101) * t65, 0.2e1 * (t61 * t107 + t6) * t62 + 0.2e1 * (-t112 * qJD(2) - t64 * t101 + t37 * t61) * t65, -0.2e1 * t21 * t10, 0.2e1 * t10 * t20 - 0.2e1 * t21 * t11, 0.2e1 * t10 * t62 - 0.2e1 * t21 * t53, 0.2e1 * t11 * t62 + 0.2e1 * t20 * t53, t46, -0.2e1 * t24 * t11 - 0.2e1 * t17 * t20 + 0.2e1 * t2 * t62 + 0.2e1 * t77 * t53, 0.2e1 * t1 * t62 + 0.2e1 * t24 * t10 - 0.2e1 * t17 * t21 - 0.2e1 * t76 * t53; 0, 0, 0, 0, 0, t53, -t100, 0, -t51, t93, t67, t51, -t93, t67 * pkin(6), t65 * t81 + t79, -t111 * t100 + 0.4e1 * t65 * t87, t25, t26, 0, -t124 * t64 + t68 * t61, t124 * t61 + t68 * t64, -t10 * t73 + t21 * t15, -t10 * t33 - t11 * t73 - t15 * t20 + t21 * t16, t12, t70, 0, -t50 * t11 + t24 * t16 + t17 * t33 - t47 * t20 - t74 * t53 + t9 * t62, t50 * t10 - t24 * t15 - t17 * t73 - t47 * t21 - t75 * t53 + t8 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, qJ(3) * t122, -0.2e1 * t87, 0.2e1 * t81, 0, 0, 0, 0.2e1 * qJD(3) * t61 + 0.2e1 * t64 * t98, 0.2e1 * qJD(3) * t64 - 0.2e1 * t61 * t98, 0.2e1 * t73 * t15, 0.2e1 * t15 * t33 + 0.2e1 * t16 * t73, 0, 0, 0, 0.2e1 * t50 * t16 + 0.2e1 * t47 * t33, -0.2e1 * t50 * t15 - 0.2e1 * t47 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, t51, 0, 0, 0, 0, 0, t25, t26, 0, 0, 0, 0, 0, t12, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90 - t91, t69, t53, t7, t6, 0, 0, t10, t11, t53, t63 * t96 + (t85 * t60 - t115) * qJD(5) + t86, -t119 + (-t4 - t96) * t60 + (t85 * t63 + t117) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t105, 0, -t61 * t103, -t64 * t103, 0, 0, -t15, -t16, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t105, 0, 0, 0, 0, 0, -t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t95, -0.2e1 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, t53, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
