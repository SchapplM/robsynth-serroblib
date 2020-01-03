% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR7
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
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:40
% EndTime: 2019-12-31 20:15:43
% DurationCPUTime: 0.58s
% Computational Cost: add. (789->147), mult. (1256->192), div. (0->0), fcn. (753->6), ass. (0->106)
t73 = sin(qJ(5));
t102 = qJD(5) * t73;
t74 = sin(qJ(4));
t104 = qJD(4) * t74;
t126 = -t74 * t102 - t73 * t104;
t67 = t74 * pkin(4);
t64 = qJ(3) + t67;
t69 = qJD(4) + qJD(5);
t125 = qJD(5) - t69;
t70 = qJD(1) + qJD(2);
t68 = t70 ^ 2;
t79 = -pkin(2) - pkin(7);
t78 = cos(qJ(2));
t93 = -t78 * pkin(1) - pkin(2);
t62 = -pkin(7) + t93;
t124 = -pkin(8) + t62;
t123 = -pkin(8) + t79;
t76 = cos(qJ(5));
t77 = cos(qJ(4));
t114 = t76 * t77;
t87 = t69 * t114;
t20 = t87 + t126;
t103 = qJD(4) * t77;
t66 = pkin(4) * t103;
t57 = qJD(3) + t66;
t106 = pkin(1) * qJD(2);
t92 = qJD(1) * t106;
t61 = t78 * t92;
t25 = t57 * t70 + t61;
t107 = pkin(1) * qJD(1);
t75 = sin(qJ(2));
t97 = t75 * t107;
t31 = t64 * t70 + t97;
t45 = t73 * t77 + t76 * t74;
t122 = t31 * t20 + t25 * t45;
t19 = t69 * t45;
t116 = t73 * t74;
t46 = t114 - t116;
t121 = -t31 * t19 + t25 * t46;
t17 = t20 * t69;
t100 = t70 * t114;
t28 = t70 * t116 - t100;
t29 = t45 * t70;
t120 = t28 * t29;
t119 = t69 * t75;
t118 = t70 * t74;
t117 = t70 * t77;
t96 = t78 * t107;
t88 = qJD(3) - t96;
t32 = t79 * t70 + t88;
t21 = -pkin(8) * t118 + t74 * t32;
t115 = t76 * t21;
t80 = qJD(4) ^ 2;
t113 = t80 * t74;
t112 = t80 * t77;
t43 = t70 * qJD(3) + t61;
t105 = t70 * qJ(3);
t47 = t97 + t105;
t111 = t47 * t103 + t43 * t74;
t110 = t126 * t70;
t109 = -t68 - t80;
t108 = t74 ^ 2 - t77 ^ 2;
t101 = pkin(4) * t117;
t99 = t75 * t106;
t98 = t78 * t106;
t39 = t124 * t77;
t52 = t123 * t77;
t22 = -pkin(8) * t117 + t77 * t32;
t18 = qJD(4) * pkin(4) + t22;
t91 = -pkin(4) * t69 - t18;
t90 = pkin(8) * t70 - t32;
t63 = t75 * pkin(1) + qJ(3);
t89 = t75 * t92;
t56 = qJD(3) + t98;
t85 = t56 * t70 - t62 * t80;
t54 = t77 * t89;
t11 = t90 * t104 + t54;
t12 = -t90 * t103 + t74 * t89;
t84 = t76 * t11 - t73 * t12 + t31 * t28;
t83 = t63 * t70 + t99;
t82 = -t47 * t70 + t89;
t81 = t21 * t102 + (-t21 * t69 - t11) * t73 + t31 * t29;
t9 = t19 * t70;
t65 = pkin(8) * t104;
t53 = t63 + t67;
t51 = t123 * t74;
t48 = -0.2e1 * t103 * t118;
t44 = -t70 * pkin(2) + t88;
t42 = t56 + t66;
t41 = qJD(4) * t52;
t40 = -t79 * t104 + t65;
t38 = t124 * t74;
t36 = t43 * t77;
t34 = (qJD(1) + t70) * t99;
t33 = (qJD(2) - t70) * t97;
t27 = 0.2e1 * t108 * t70 * qJD(4);
t24 = qJD(4) * t39 + t74 * t99;
t23 = -t62 * t104 + t77 * t99 + t65;
t16 = t19 * t69;
t10 = t70 * t87 + t110;
t5 = t28 ^ 2 - t29 ^ 2;
t4 = -t110 + (-t100 - t28) * t69;
t3 = t29 * t69 - t9;
t2 = t28 * t19 - t9 * t46;
t1 = -t46 * t10 + t19 * t29 + t28 * t20 + t9 * t45;
t6 = [0, 0, 0, 0, -t34, -t70 * t98 - t61, t34, t61 + (qJD(3) + t56) * t70, t43 * t63 + t47 * t56 + (t93 * qJD(1) + t44) * t99, t48, t27, -t113, -t112, 0, t83 * t103 + t85 * t74 + t111, t36 + t85 * t77 + (-t47 - t83) * t104, t2, t1, -t16, -t17, 0, t42 * t29 + t53 * t10 + (t76 * t23 - t73 * t24 + (-t38 * t76 - t39 * t73) * qJD(5)) * t69 + t122, -t42 * t28 - t53 * t9 - (t73 * t23 + t76 * t24 + (-t38 * t73 + t39 * t76) * qJD(5)) * t69 + t121; 0, 0, 0, 0, -t33, t70 * t96 - t61, t33, t61 + (0.2e1 * qJD(3) - t96) * t70, t43 * qJ(3) + t47 * qJD(3) + (-t47 * t78 + (-pkin(2) * qJD(2) - t44) * t75) * t107, t48, t27, -t113, -t112, 0, -t97 * t103 - t79 * t113 + (qJ(3) * t103 + t88 * t74) * t70 + t111, t36 + (t88 * t70 - t79 * t80) * t77 + (-t47 + t97 - t105) * t104, t2, t1, -t16, -t17, 0, t57 * t29 + t64 * t10 + (t76 * t40 - t73 * t41 + (-t51 * t76 - t52 * t73) * qJD(5)) * t69 + (-t46 * t119 - t78 * t29) * t107 + t122, -t57 * t28 - t64 * t9 - (t73 * t40 + t76 * t41 + (-t51 * t73 + t52 * t76) * qJD(5)) * t69 + (t45 * t119 + t78 * t28) * t107 + t121; 0, 0, 0, 0, 0, 0, 0, -t68, t82, 0, 0, 0, 0, 0, t109 * t74, t109 * t77, 0, 0, 0, 0, 0, -t70 * t29 - t16, t70 * t28 - t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t68 * t74, -t108 * t68, 0, 0, 0, -t47 * t117 + t54, -t82 * t74, -t120, t5, t3, t4, 0, -t29 * t101 - (-t73 * t22 - t115) * t69 + (t91 * t73 - t115) * qJD(5) + t84, t28 * t101 + (t91 * qJD(5) + t22 * t69 - t12) * t76 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, t5, t3, t4, 0, t84 + t125 * (-t73 * t18 - t115), (-t125 * t18 - t12) * t76 + t81;];
tauc_reg = t6;
