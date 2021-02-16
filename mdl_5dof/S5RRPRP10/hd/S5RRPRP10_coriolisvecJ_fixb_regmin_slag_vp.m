% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:02:45
% EndTime: 2021-01-15 21:02:54
% DurationCPUTime: 1.68s
% Computational Cost: add. (1598->256), mult. (3661->353), div. (0->0), fcn. (2022->4), ass. (0->148)
t93 = sin(qJ(4));
t140 = t93 * qJD(2);
t96 = cos(qJ(2));
t146 = qJD(1) * t96;
t95 = cos(qJ(4));
t53 = t95 * t146 + t140;
t94 = sin(qJ(2));
t139 = t94 * qJD(1);
t80 = qJD(4) + t139;
t164 = t53 * t80;
t134 = qJD(1) * qJD(2);
t126 = t94 * t134;
t23 = t53 * qJD(4) - t93 * t126;
t181 = t23 - t164;
t180 = t23 + t164;
t171 = pkin(3) + pkin(6);
t137 = t95 * qJD(2);
t55 = -t93 * t146 + t137;
t124 = -t94 * qJ(3) - pkin(1);
t97 = -pkin(2) - pkin(7);
t50 = t97 * t96 + t124;
t30 = t50 * qJD(1);
t82 = pkin(6) * t139;
t174 = qJD(3) + t82;
t135 = pkin(3) * t139 + t174;
t33 = t97 * qJD(2) + t135;
t12 = t95 * t30 + t93 * t33;
t115 = pkin(7) * t94 - qJ(3) * t96;
t138 = t94 * qJD(3);
t103 = t115 * qJD(2) - t138;
t79 = pkin(2) * t126;
t20 = t103 * qJD(1) + t79;
t125 = t96 * t134;
t78 = pkin(6) * t125;
t49 = pkin(3) * t125 + t78;
t123 = -t93 * t20 + t95 * t49;
t102 = -t12 * qJD(4) + t123;
t101 = t23 * qJ(5) + t102;
t116 = pkin(4) * t125;
t1 = -t55 * qJD(5) + t101 + t116;
t7 = -t53 * qJ(5) + t12;
t168 = t7 * t80;
t142 = qJD(4) * t95;
t143 = qJD(4) * t93;
t119 = -t33 * t142 + t30 * t143 - t95 * t20 - t93 * t49;
t24 = t55 * qJD(4) - t95 * t126;
t107 = t24 * qJ(5) + t119;
t2 = -t53 * qJD(5) - t107;
t11 = -t93 * t30 + t95 * t33;
t6 = -t55 * qJ(5) + t11;
t5 = t80 * pkin(4) + t6;
t179 = -(t5 * t80 - t2) * t93 + (t1 + t168) * t95;
t178 = -0.2e1 * t134;
t70 = t171 * t94;
t151 = t95 * t50 + t93 * t70;
t163 = t55 * t80;
t176 = -t24 + t163;
t175 = t24 + t163;
t172 = t55 ^ 2;
t170 = t5 - t6;
t169 = pkin(4) * t96;
t86 = pkin(2) * t139;
t39 = t115 * qJD(1) + t86;
t83 = pkin(6) * t146;
t62 = pkin(3) * t146 + t83;
t122 = -t93 * t39 + t95 * t62;
t136 = t95 * qJD(5);
t147 = qJ(5) - t97;
t160 = t93 * t94;
t167 = t147 * t143 - t136 - (-qJ(5) * t160 + t169) * qJD(1) - t122;
t145 = qJD(2) * t94;
t61 = t171 * t145;
t89 = qJD(2) * qJD(3);
t36 = -qJD(1) * t61 + t89;
t166 = t36 * t93;
t165 = t36 * t95;
t162 = t55 * t96;
t161 = t80 * t97;
t159 = t94 * t95;
t158 = t95 * t23;
t99 = qJD(1) ^ 2;
t157 = t96 * t99;
t98 = qJD(2) ^ 2;
t156 = t98 * t94;
t155 = t98 * t96;
t153 = t95 * t39 + t93 * t62;
t66 = t147 * t95;
t154 = t95 * qJ(5) * t139 + qJD(4) * t66 + t93 * qJD(5) + t153;
t127 = -pkin(4) * t95 - pkin(3);
t152 = pkin(4) * t142 - t127 * t139 + t174;
t71 = t171 * t96;
t91 = t94 ^ 2;
t92 = t96 ^ 2;
t150 = t91 - t92;
t149 = qJ(5) * t96;
t148 = qJD(2) * pkin(2);
t68 = -t96 * pkin(2) + t124;
t44 = qJD(1) * t68;
t144 = qJD(2) * t96;
t141 = qJD(4) * t96;
t90 = qJD(2) * qJ(3);
t133 = t80 * t159;
t132 = t94 * t157;
t43 = t90 + t62;
t131 = t93 * t141;
t130 = t80 * t142;
t129 = t95 * t141;
t121 = -t53 * pkin(4) - qJD(5);
t120 = -t50 + t149;
t118 = pkin(1) * t178;
t117 = qJD(3) - t148;
t114 = -qJD(1) * t92 + t80 * t94;
t113 = -0.2e1 * qJD(2) * t44;
t112 = t80 * t93;
t105 = -t96 * t90 - t138;
t28 = t105 * qJD(1) + t79;
t85 = pkin(2) * t145;
t41 = t105 + t85;
t111 = pkin(6) * t98 + qJD(1) * t41 + t28;
t110 = t97 * t144 + t43 * t94;
t14 = t24 * pkin(4) + t36;
t18 = -t121 + t43;
t109 = -t14 * t93 - t18 * t142;
t108 = t14 * t95 - t18 * t143;
t26 = t85 + t103;
t63 = t171 * t144;
t106 = t70 * t142 - t50 * t143 + t95 * t26 + t93 * t63;
t64 = pkin(6) * t126 - t89;
t67 = t117 + t82;
t69 = -t83 - t90;
t100 = -t64 * t96 + (t67 * t96 + (t69 + t83) * t94) * qJD(2);
t81 = t93 * pkin(4) + qJ(3);
t73 = t95 * t125;
t65 = t147 * t93;
t59 = -qJ(3) * t146 + t86;
t58 = t95 * t70;
t52 = t53 ^ 2;
t48 = t95 * t63;
t42 = t95 * t169 + t71;
t32 = t44 * t139;
t19 = -pkin(4) * t131 + (-pkin(6) + t127) * t145;
t16 = -t95 * t149 + t151;
t15 = t94 * pkin(4) + t120 * t93 + t58;
t10 = -t130 - qJD(2) * t55 + (-t96 * t140 - t133) * qJD(1);
t9 = -qJD(2) * t53 - t80 * t112 + t73;
t4 = -t96 * t136 + (t94 * t137 + t131) * qJ(5) + t106;
t3 = pkin(4) * t144 + t48 + t120 * t142 + (-qJ(5) * t145 - qJD(4) * t70 + qJD(5) * t96 - t26) * t93;
t8 = [0, 0, 0, 0.2e1 * t94 * t125, t150 * t178, t155, -t156, 0, -pkin(6) * t155 + t94 * t118, pkin(6) * t156 + t96 * t118, t100, t111 * t96 + t94 * t113, -t111 * t94 + t96 * t113, t100 * pkin(6) + t28 * t68 + t44 * t41, t23 * t93 * t96 + (t94 * t140 - t129) * t55, (-t53 * t93 + t55 * t95) * t145 + (t158 + t93 * t24 + (t53 * t95 + t55 * t93) * qJD(4)) * t96, -t80 * t129 - t23 * t94 + (t114 * t93 + t162) * qJD(2), t80 * t131 - t24 * t94 + (t114 * t95 - t53 * t96) * qJD(2), (t80 + t139) * t144, (-t93 * t26 + t48) * t80 - t61 * t53 + t71 * t24 + (-t43 * t137 + t123) * t94 + (-t12 * t94 - t151 * t80) * qJD(4) + (-t43 * t143 + t165 + ((-t93 * t50 + t58) * qJD(1) + t11) * qJD(2)) * t96, -t106 * t80 - t61 * t55 - t71 * t23 + (t43 * t140 + t119) * t94 + (-t43 * t142 - t166 + (-t151 * qJD(1) - t12) * qJD(2)) * t96, t19 * t53 + t42 * t24 + t3 * t80 + (-t137 * t18 + t1) * t94 + ((qJD(1) * t15 + t5) * qJD(2) + t108) * t96, t19 * t55 - t42 * t23 - t4 * t80 + (t140 * t18 - t2) * t94 + ((-qJD(1) * t16 - t7) * qJD(2) + t109) * t96, t15 * t23 - t16 * t24 - t3 * t55 - t4 * t53 + (-t5 * t93 + t7 * t95) * t145 + (t1 * t93 - t2 * t95 + (t5 * t95 + t7 * t93) * qJD(4)) * t96, t1 * t15 + t14 * t42 + t2 * t16 + t18 * t19 + t5 * t3 + t7 * t4; 0, 0, 0, -t132, t150 * t99, 0, 0, 0, t99 * pkin(1) * t94, pkin(1) * t157, ((-t69 - t90) * t94 + (t117 - t67) * t96) * qJD(1), -t59 * t146 + t32, 0.2e1 * t89 + (t44 * t96 + t59 * t94) * qJD(1), -t64 * qJ(3) - t69 * qJD(3) - t44 * t59 + (-t69 * t94 + (-t67 - t148) * t96) * qJD(1) * pkin(6), -t55 * t112 - t158, -t175 * t95 + t180 * t93, -t80 * t143 + t73 + (-t80 * t160 - t162) * qJD(1), -t130 + (-t133 + (t53 - t140) * t96) * qJD(1), -t80 * t146, qJ(3) * t24 + t166 - t122 * t80 + t135 * t53 + (-t93 * t161 + t43 * t95) * qJD(4) + (-t11 * t96 + t110 * t95) * qJD(1), -qJ(3) * t23 + t165 + t153 * t80 + t135 * t55 + (-t161 * t95 - t43 * t93) * qJD(4) + (-t110 * t93 + t12 * t96) * qJD(1), t81 * t24 + t167 * t80 + t152 * t53 + (t18 * t159 + (-qJD(2) * t66 - t5) * t96) * qJD(1) - t109, -t81 * t23 + t154 * t80 + t152 * t55 + (-t18 * t160 + (qJD(2) * t65 + t7) * t96) * qJD(1) + t108, t154 * t53 - t167 * t55 - t66 * t23 + t65 * t24 - t179, -t1 * t66 + t14 * t81 + t152 * t18 - t154 * t7 + t167 * t5 - t2 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, -t91 * t99 - t98, t69 * qJD(2) + t32 + t78, 0, 0, 0, 0, 0, t9, t10, t9, t10, t176 * t93 + t181 * t95, -t18 * qJD(2) + t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t53, -t52 + t172, -t181, t176, t125, t12 * t80 - t43 * t55 + t102, t11 * t80 + t43 * t53 + t119, 0.2e1 * t116 + t168 + (t121 - t18) * t55 + t101, -t172 * pkin(4) + t6 * t80 + (qJD(5) + t18) * t53 + t107, t23 * pkin(4) - t170 * t53, t170 * t7 + (-t18 * t55 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, -t180, -t52 - t172, t5 * t55 + t7 * t53 + t14;];
tauc_reg = t8;
