% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:06
% EndTime: 2019-12-31 17:10:10
% DurationCPUTime: 1.39s
% Computational Cost: add. (1245->260), mult. (2935->376), div. (0->0), fcn. (2043->10), ass. (0->137)
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t102 = sin(pkin(7));
t105 = sin(qJ(2));
t154 = qJD(1) * t105;
t139 = t102 * t154;
t103 = cos(pkin(7));
t147 = t103 * qJD(2);
t68 = t139 - t147;
t138 = t103 * t154;
t152 = qJD(2) * t102;
t70 = t138 + t152;
t127 = t104 * t68 - t107 * t70;
t108 = cos(qJ(2));
t153 = qJD(1) * t108;
t88 = -qJD(4) + t153;
t176 = t127 * t88;
t146 = qJD(1) * qJD(2);
t135 = t105 * t146;
t98 = t108 * qJDD(1);
t117 = t135 - t98;
t106 = sin(qJ(1));
t109 = cos(qJ(1));
t133 = g(1) * t109 + g(2) * t106;
t168 = g(3) * t108;
t136 = t108 * t146;
t144 = t105 * qJDD(1);
t91 = pkin(5) * t144;
t55 = -qJDD(2) * pkin(2) + pkin(5) * t136 + qJDD(3) + t91;
t175 = -t133 * t105 + t168 + t55;
t174 = -qJD(4) - t88;
t132 = g(1) * t106 - g(2) * t109;
t118 = t136 + t144;
t94 = t103 * qJDD(2);
t38 = t118 * t102 - t94;
t145 = qJDD(2) * t102;
t39 = t118 * t103 + t145;
t6 = -t127 * qJD(4) + t104 * t39 + t107 * t38;
t173 = pkin(5) * t68;
t172 = pkin(5) * t70;
t169 = g(3) * t105;
t163 = t104 * t70;
t21 = t107 * t68 + t163;
t167 = t21 * t88;
t166 = pkin(6) + qJ(3);
t125 = pkin(2) * t108 + qJ(3) * t105 + pkin(1);
t130 = pkin(2) * t105 - qJ(3) * t108;
t57 = t130 * qJD(2) - qJD(3) * t105;
t20 = t57 * qJD(1) - t125 * qJDD(1);
t46 = -t117 * pkin(5) + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t9 = t102 * t20 + t103 * t46;
t162 = t102 * t104;
t73 = -t107 * t103 + t162;
t119 = t73 * t108;
t165 = qJD(1) * t119 - t73 * qJD(4);
t74 = t102 * t107 + t103 * t104;
t120 = t74 * t108;
t164 = -qJD(1) * t120 + t74 * qJD(4);
t63 = t125 * qJD(1);
t93 = pkin(5) * t153;
t83 = qJD(2) * qJ(3) + t93;
t27 = -t102 * t63 + t103 * t83;
t76 = t130 * qJD(1);
t35 = pkin(5) * t139 + t103 * t76;
t151 = qJD(2) * t105;
t142 = pkin(5) * t151;
t28 = t102 * t142 + t103 * t57;
t159 = t103 * t108;
t41 = pkin(5) * t159 - t102 * t125;
t161 = t102 * t108;
t160 = t103 * t105;
t158 = t106 * t108;
t157 = t108 * t109;
t78 = -qJD(2) * pkin(2) + pkin(5) * t154 + qJD(3);
t156 = qJD(3) - t78;
t100 = t105 ^ 2;
t155 = -t108 ^ 2 + t100;
t150 = qJD(2) * t108;
t149 = qJD(4) * t105;
t148 = qJD(4) * t107;
t143 = pkin(3) * t153;
t141 = pkin(3) * t102 + pkin(5);
t8 = -t102 * t46 + t103 * t20;
t4 = t117 * pkin(3) - pkin(6) * t39 + t8;
t7 = -pkin(6) * t38 + t9;
t140 = -t104 * t7 + t107 * t4;
t137 = qJ(3) * t98;
t26 = -t102 * t83 - t103 * t63;
t131 = t104 * t4 + t107 * t7;
t10 = -pkin(6) * t70 - t143 + t26;
t11 = -pkin(6) * t68 + t27;
t1 = t10 * t107 - t104 * t11;
t2 = t10 * t104 + t107 * t11;
t67 = t103 * t125;
t25 = -pkin(6) * t160 - t67 + (-pkin(5) * t102 - pkin(3)) * t108;
t30 = -pkin(6) * t102 * t105 + t41;
t129 = -t104 * t30 + t107 * t25;
t128 = t104 * t25 + t107 * t30;
t126 = pkin(3) * t105 - pkin(6) * t159;
t124 = -0.2e1 * pkin(1) * t146 - pkin(5) * qJDD(2);
t82 = t166 * t103;
t123 = t126 * qJD(1) + t102 * qJD(3) + qJD(4) * t82 + t35;
t121 = -pkin(5) * t160 - pkin(6) * t161;
t61 = t102 * t76;
t81 = t166 * t102;
t122 = -t121 * qJD(1) + t103 * qJD(3) - qJD(4) * t81 - t61;
t5 = -qJD(4) * t163 - t104 * t38 + t107 * t39 - t68 * t148;
t111 = qJD(1) ^ 2;
t116 = pkin(1) * t111 + t133;
t114 = -t133 * t108 - t169;
t110 = qJD(2) ^ 2;
t113 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t110 + t132;
t99 = pkin(7) + qJ(4);
t97 = cos(t99);
t96 = sin(t99);
t90 = -pkin(3) * t103 - pkin(2);
t77 = t141 * t105;
t72 = qJDD(4) + t117;
t65 = t141 * t150;
t64 = t102 * t143 + t93;
t54 = t73 * t105;
t53 = t74 * t105;
t50 = t106 * t96 + t97 * t157;
t49 = t106 * t97 - t96 * t157;
t48 = t109 * t96 - t97 * t158;
t47 = t109 * t97 + t96 * t158;
t44 = t102 * t57;
t40 = -pkin(5) * t161 - t67;
t36 = -pkin(5) * t138 + t61;
t34 = pkin(3) * t68 + t78;
t29 = -t103 * t142 + t44;
t19 = t121 * qJD(2) + t44;
t15 = qJD(2) * t120 + t148 * t160 - t149 * t162;
t14 = -qJD(2) * t119 - t74 * t149;
t13 = pkin(3) * t38 + t55;
t12 = t126 * qJD(2) + t28;
t3 = [qJDD(1), t132, t133, qJDD(1) * t100 + 0.2e1 * t108 * t135, 0.2e1 * t105 * t98 - 0.2e1 * t155 * t146, qJDD(2) * t105 + t108 * t110, qJDD(2) * t108 - t105 * t110, 0, t124 * t105 + t113 * t108, -t113 * t105 + t124 * t108, -t133 * t102 + (pkin(5) * t38 + t55 * t102 + (qJD(1) * t40 + t26) * qJD(2)) * t105 + (-t28 * qJD(1) - t40 * qJDD(1) - t8 + t132 * t103 + (t102 * t78 + t173) * qJD(2)) * t108, -t133 * t103 + (pkin(5) * t39 + t55 * t103 + (-qJD(1) * t41 - t27) * qJD(2)) * t105 + (t29 * qJD(1) + t41 * qJDD(1) + t9 - t132 * t102 + (t103 * t78 + t172) * qJD(2)) * t108, -t28 * t70 - t29 * t68 - t38 * t41 - t39 * t40 + (-t102 * t27 - t103 * t26) * t150 + (-t102 * t9 - t103 * t8 + t132) * t105, t26 * t28 + t27 * t29 + t8 * t40 + t9 * t41 + (t55 * t105 + t78 * t150 - t133) * pkin(5) + t132 * t125, -t127 * t14 - t5 * t54, t127 * t15 - t14 * t21 - t5 * t53 + t54 * t6, -t108 * t5 - t127 * t151 - t14 * t88 - t54 * t72, t108 * t6 + t15 * t88 - t21 * t151 - t53 * t72, -t108 * t72 - t88 * t151, -(-t104 * t19 + t107 * t12) * t88 + t129 * t72 - t140 * t108 + t1 * t151 + t65 * t21 + t77 * t6 + t13 * t53 + t34 * t15 - g(1) * t48 - g(2) * t50 + (t2 * t108 + t128 * t88) * qJD(4), (t104 * t12 + t107 * t19) * t88 - t128 * t72 + t131 * t108 - t2 * t151 - t65 * t127 + t77 * t5 - t13 * t54 + t34 * t14 - g(1) * t47 - g(2) * t49 + (t1 * t108 + t129 * t88) * qJD(4); 0, 0, 0, -t105 * t111 * t108, t155 * t111, t144, t98, qJDD(2), t116 * t105 - t168 - t91, t169 + (-pkin(5) * qJDD(1) + t116) * t108, t102 * t137 - pkin(2) * t38 - t175 * t103 + ((-qJ(3) * t152 - t26) * t105 + (t156 * t102 - t173 + t35) * t108) * qJD(1), t103 * t137 - pkin(2) * t39 + t175 * t102 + ((-qJ(3) * t147 + t27) * t105 + (t156 * t103 - t172 - t36) * t108) * qJD(1), t35 * t70 + t36 * t68 + (-qJ(3) * t38 - qJD(3) * t68 + t26 * t153 + t9) * t103 + (qJ(3) * t39 + qJD(3) * t70 + t27 * t153 - t8) * t102 + t114, -t78 * t93 - t26 * t35 - t27 * t36 + (-t26 * t102 + t27 * t103) * qJD(3) - t175 * pkin(2) + (-t8 * t102 + t9 * t103 + t114) * qJ(3), -t127 * t165 + t5 * t74, t127 * t164 - t165 * t21 - t5 * t73 - t6 * t74, t127 * t154 - t165 * t88 + t72 * t74, t21 * t154 + t164 * t88 - t72 * t73, t88 * t154, (-t104 * t82 - t107 * t81) * t72 + t90 * t6 + t13 * t73 - t64 * t21 - t97 * t168 + (t122 * t104 + t123 * t107) * t88 + t164 * t34 + (-t1 * qJD(1) + t133 * t97) * t105, -(-t104 * t81 + t107 * t82) * t72 + t90 * t5 + t13 * t74 + t64 * t127 + t96 * t168 + (-t123 * t104 + t122 * t107) * t88 + t165 * t34 + (t2 * qJD(1) - t133 * t96) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * t144 - t94 + (-t70 + t152) * t153, t103 * t144 + t145 + (t68 + t147) * t153, -t68 ^ 2 - t70 ^ 2, t26 * t70 + t27 * t68 + t175, 0, 0, 0, 0, 0, t6 + t176, t5 + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127 * t21, t127 ^ 2 - t21 ^ 2, t5 - t167, -t6 + t176, t72, -g(1) * t49 + g(2) * t47 + t127 * t34 + t96 * t169 + t174 * t2 + t140, g(1) * t50 - g(2) * t48 + t174 * t1 + t97 * t169 + t21 * t34 - t131;];
tau_reg = t3;
