% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:51
% EndTime: 2021-01-15 15:04:56
% DurationCPUTime: 1.14s
% Computational Cost: add. (899->209), mult. (1862->278), div. (0->0), fcn. (1234->6), ass. (0->128)
t66 = sin(pkin(8));
t69 = sin(qJ(4));
t105 = qJD(2) * qJD(4);
t70 = cos(qJ(4));
t94 = t70 * t105;
t75 = qJDD(2) * t69 + t94;
t155 = t75 * t66;
t63 = pkin(7) + qJ(2);
t59 = sin(t63);
t60 = cos(t63);
t96 = g(1) * t59 - g(2) * t60;
t115 = qJDD(2) * pkin(2);
t154 = t115 + t96;
t67 = cos(pkin(8));
t114 = t67 * qJD(2);
t48 = -qJD(4) + t114;
t153 = qJD(4) + t48;
t106 = qJD(2) * qJD(3);
t108 = qJ(3) * qJDD(2);
t77 = t106 + t108;
t152 = -t66 * (-qJ(5) - pkin(6)) + (t70 * pkin(4) + pkin(3)) * t67;
t134 = t67 * t69;
t21 = t59 * t134 + t60 * t70;
t23 = -t60 * t134 + t59 * t70;
t151 = -g(1) * t23 + g(2) * t21;
t143 = g(3) * t69;
t150 = t66 * t143 + t151;
t149 = pkin(4) * t155 + qJDD(5);
t147 = pkin(4) * t69;
t109 = t67 * qJDD(2);
t46 = -qJDD(4) + t109;
t142 = t46 * pkin(4);
t39 = -t67 * pkin(3) - t66 * pkin(6) - pkin(2);
t28 = t39 * qJD(2) + qJD(3);
t112 = qJ(3) * qJD(2);
t37 = t66 * qJD(1) + t67 * t112;
t93 = t70 * t28 - t69 * t37;
t123 = qJD(2) * t66;
t97 = qJ(5) * t123;
t10 = -t70 * t97 + t93;
t3 = -t48 * pkin(4) + t10;
t141 = -t10 + t3;
t57 = t67 * qJDD(1);
t25 = t77 * t66 - t57;
t140 = t25 * t66;
t139 = t46 * t67;
t137 = t60 * t69;
t61 = t66 ^ 2;
t71 = qJD(2) ^ 2;
t136 = t61 * t71;
t133 = t67 * t70;
t117 = qJD(4) * t70;
t120 = qJD(3) * t67;
t132 = t39 * t117 + t70 * t120;
t47 = qJ(3) * t133;
t131 = t69 * t39 + t47;
t130 = t60 * pkin(2) + t59 * qJ(3);
t129 = t67 ^ 2 + t61;
t64 = t69 ^ 2;
t65 = t70 ^ 2;
t128 = -t64 - t65;
t127 = t64 - t65;
t126 = qJ(3) * t69;
t125 = qJ(5) * t66;
t35 = (qJ(3) + t147) * t66;
t124 = qJD(2) * t35;
t122 = qJD(2) * t69;
t121 = qJD(2) * t70;
t119 = qJD(4) * t37;
t118 = qJD(4) * t69;
t116 = qJD(5) * t66;
t58 = t67 * qJD(1);
t16 = qJD(5) - t58 + t124;
t113 = qJD(5) + t16;
t110 = qJDD(2) * t70;
t107 = qJ(5) * qJDD(2);
t104 = qJD(2) * qJD(5);
t103 = t69 * t136;
t102 = t67 * t126;
t101 = t70 * t125;
t100 = t66 * t122;
t99 = t48 * t118;
t36 = t66 * t112 - t58;
t98 = t36 * t123;
t95 = t69 * t105;
t26 = t66 * qJDD(1) + t77 * t67;
t27 = t39 * qJDD(2) + qJDD(3);
t92 = t28 * t117 - t37 * t118 + t70 * t26 + t69 * t27;
t91 = t16 + t124;
t90 = -qJD(4) * t28 - t26;
t89 = t46 - t109;
t88 = t46 + t109;
t87 = t66 * t95;
t86 = -g(1) * t21 - g(2) * t23;
t22 = -t59 * t133 + t137;
t24 = t60 * t133 + t59 * t69;
t85 = -g(1) * t22 - g(2) * t24;
t84 = g(1) * t60 + g(2) * t59;
t80 = -t69 * t28 - t70 * t37;
t11 = -t69 * t97 - t80;
t83 = t11 * t70 - t3 * t69;
t82 = -t11 * t69 - t3 * t70;
t81 = t26 * t67 + t140;
t79 = t36 * t66 + t37 * t67;
t12 = t25 + t149;
t32 = (pkin(4) * t117 + qJD(3)) * t66;
t76 = qJD(2) * t32 + qJDD(2) * t35 + t12;
t19 = t70 * t27;
t74 = qJ(5) * t87 + t90 * t69 + t19;
t73 = -t48 ^ 2 - t136;
t72 = g(3) * t66 * t70 + g(1) * t24 - g(2) * t22 - t92;
t56 = qJDD(3) - t115;
t52 = t60 * qJ(3);
t44 = t66 * t110;
t38 = t67 * t87;
t34 = t70 * t39;
t31 = t66 * t48 * t121;
t14 = -t69 * t125 + t131;
t13 = -t101 + t34 + (-pkin(4) - t126) * t67;
t9 = t69 * t46 + t73 * t70;
t8 = -t70 * t46 + t73 * t69;
t7 = -t69 * t120 - t70 * t116 + (-t47 + (-t39 + t125) * t69) * qJD(4);
t6 = -t69 * t116 + (-t101 - t102) * qJD(4) + t132;
t5 = (t89 * t69 + (t48 - t114) * t117) * t66;
t4 = t38 + (t89 * t70 - t99) * t66;
t2 = (-t75 * qJ(5) - t69 * t104) * t66 + t92;
t1 = -t142 + (-t119 + (-t104 - t107) * t66) * t70 + t74;
t15 = [qJDD(1) - g(3), 0, 0, 0, 0, 0, -t25 * t67 + t26 * t66 - g(3), 0, 0, 0, 0, 0, t5, t4, t5, t4, 0, -t12 * t67 - g(3) + (t82 * qJD(4) - t1 * t69 + t2 * t70) * t66; 0, qJDD(2), t96, t84, (-t56 + t154) * t67, t77 * t129 + t81 - t84, -t56 * pkin(2) - g(1) * (-t59 * pkin(2) + t52) - g(2) * t130 + t79 * qJD(3) + t81 * qJ(3), (qJDD(2) * t65 - 0.2e1 * t69 * t94) * t61, 0.2e1 * (t127 * t105 - t69 * t110) * t61, t38 + (-t88 * t70 + t99) * t66, (t88 * t69 + (t48 + t114) * t117) * t66, t139, -t19 * t67 - t34 * t46 + ((qJD(2) * t61 + t48 * t67) * qJ(3) + t79) * t117 + (-(-qJD(4) * t39 - t120) * t48 - t90 * t67 + t61 * t106 + t140 + (t61 * qJDD(2) + t139) * qJ(3)) * t69 + t85, (-qJD(4) * t102 + t132) * t48 + t131 * t46 + t92 * t67 + (-t36 * t118 + t25 * t70) * t66 + (t70 * t106 + (-t95 + t110) * qJ(3)) * t61 + t86, -t1 * t67 - t13 * t46 - t7 * t48 + (t91 * t117 + t76 * t69) * t66 + t85, t14 * t46 + t2 * t67 + t6 * t48 + (-t91 * t118 + t76 * t70) * t66 + t86, ((-qJD(4) * t11 - qJDD(2) * t13 - t1 + (-qJD(4) * t14 - t7) * qJD(2)) * t70 + (qJD(4) * t3 - qJDD(2) * t14 - t2 + (qJD(4) * t13 - t6) * qJD(2)) * t69 + t96) * t66, t2 * t14 + t11 * t6 + t1 * t13 + t3 * t7 + t12 * t35 + t16 * t32 - g(1) * (pkin(4) * t137 + t52) - g(2) * (t152 * t60 + t130) + (-g(1) * (-pkin(2) - t152) - g(2) * t147) * t59; 0, 0, 0, 0, -t109, -t129 * t71, -t79 * qJD(2) + qJDD(3) - t154, 0, 0, 0, 0, 0, t8, t9, t8, t9, t128 * t66 * qJDD(2), t1 * t70 + t2 * t69 + t83 * qJD(4) + (-t16 * t66 - t83 * t67) * qJD(2) - t96; 0, 0, 0, 0, 0, 0, 0, t70 * t103, -t127 * t136, -t153 * t100 + t44, -t31 - t155, -t46, t153 * t80 - t69 * t26 - t70 * t98 + t150 + t19, -t93 * t48 + t69 * t98 + t72, -0.2e1 * t142 - t11 * t48 + (-pkin(4) * t103 - t119 + (-t113 * qJD(2) - t107) * t66) * t70 + t74 + t150, -t65 * pkin(4) * t136 - t10 * t48 + (t69 * t107 + (qJ(5) * t117 + t113 * t69) * qJD(2)) * t66 + t72, (-pkin(4) * t110 + (pkin(4) * qJD(4) - t141) * t122) * t66, t141 * t11 + (t1 + (-t16 * t121 + t143) * t66 + t151) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 + t155, t44 + (-qJD(4) + t48) * t100, t128 * t136, g(3) * t67 - t57 + (t108 + (qJD(3) - t82) * qJD(2) - t84) * t66 + t149;];
tau_reg = t15;
