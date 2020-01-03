% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:31
% EndTime: 2019-12-31 20:01:37
% DurationCPUTime: 1.64s
% Computational Cost: add. (2727->277), mult. (6956->374), div. (0->0), fcn. (4811->6), ass. (0->134)
t101 = sin(pkin(8));
t104 = sin(qJ(2));
t141 = qJD(1) * t104;
t102 = cos(pkin(8));
t106 = cos(qJ(2));
t146 = t102 * t106;
t73 = qJD(1) * t146 - t101 * t141;
t69 = qJD(4) - t73;
t159 = -qJ(3) - pkin(6);
t90 = t159 * t106;
t87 = qJD(1) * t90;
t78 = t101 * t87;
t153 = qJD(2) * pkin(2);
t89 = t159 * t104;
t86 = qJD(1) * t89;
t81 = t86 + t153;
t48 = t102 * t81 + t78;
t43 = -qJD(2) * pkin(3) - t48;
t103 = sin(qJ(4));
t105 = cos(qJ(4));
t138 = t105 * qJD(2);
t84 = t101 * t106 + t102 * t104;
t75 = t84 * qJD(1);
t56 = t103 * t75 - t138;
t58 = qJD(2) * t103 + t105 * t75;
t14 = pkin(4) * t56 - qJ(5) * t58 + t43;
t74 = t84 * qJD(2);
t67 = qJD(1) * t74;
t95 = pkin(2) * t101 + pkin(7);
t160 = t95 * t67;
t172 = t69 * t14 - t160;
t137 = qJD(1) * qJD(2);
t130 = t106 * t137;
t131 = t104 * t137;
t115 = -t101 * t131 + t102 * t130;
t30 = qJD(4) * t58 + t103 * t115;
t171 = t58 ^ 2;
t170 = t69 ^ 2;
t169 = pkin(4) * t67;
t132 = -pkin(2) * t106 - pkin(1);
t119 = t132 * qJD(1);
t88 = qJD(3) + t119;
t32 = -t73 * pkin(3) - t75 * pkin(7) + t88;
t152 = t102 * t87;
t49 = t101 * t81 - t152;
t44 = qJD(2) * pkin(7) + t49;
t13 = t103 * t32 + t105 * t44;
t8 = qJ(5) * t69 + t13;
t168 = t69 * t8;
t167 = t13 * t69;
t129 = qJD(2) * t159;
t71 = -t104 * qJD(3) + t106 * t129;
t109 = qJD(1) * t71;
t70 = t106 * qJD(3) + t104 * t129;
t65 = t70 * qJD(1);
t27 = t101 * t65 - t102 * t109;
t166 = t27 * t84;
t83 = t101 * t104 - t146;
t77 = t83 * qJD(2);
t165 = t43 * t77;
t164 = t56 * t73;
t163 = t58 * t56;
t128 = t58 * t69;
t162 = t58 * t75;
t161 = t75 * t56;
t118 = pkin(4) * t103 - qJ(5) * t105;
t50 = t101 * t86 - t152;
t158 = t103 * qJD(5) - t69 * t118 + t50;
t139 = qJD(4) * t105;
t157 = -t103 * t30 - t56 * t139;
t39 = pkin(2) * t141 + pkin(3) * t75 - pkin(7) * t73;
t51 = t102 * t86 + t78;
t156 = t103 * t39 + t105 * t51;
t47 = pkin(3) * t83 - pkin(7) * t84 + t132;
t55 = t101 * t89 - t102 * t90;
t155 = t103 * t47 + t105 * t55;
t154 = qJ(5) * t67;
t62 = t103 * t67;
t151 = t103 * t69;
t63 = t105 * t67;
t140 = qJD(4) * t103;
t29 = -qJD(4) * t138 - t105 * t115 + t75 * t140;
t150 = t29 * t103;
t149 = -t104 ^ 2 + t106 ^ 2;
t148 = qJD(4) * t84;
t147 = qJD(4) * t95;
t108 = qJD(1) ^ 2;
t145 = t106 * t108;
t107 = qJD(2) ^ 2;
t144 = t107 * t104;
t143 = t107 * t106;
t12 = -t103 * t44 + t105 * t32;
t142 = qJD(5) - t12;
t136 = t104 * t153;
t135 = t69 * t147;
t134 = t84 * t140;
t133 = t84 * t139;
t96 = -pkin(2) * t102 - pkin(3);
t37 = t101 * t70 - t102 * t71;
t54 = -t101 * t90 - t102 * t89;
t28 = t101 * t109 + t102 * t65;
t93 = pkin(2) * t131;
t31 = t67 * pkin(3) - t115 * pkin(7) + t93;
t127 = t103 * t28 - t105 * t31 + t44 * t139 + t32 * t140;
t125 = -0.2e1 * pkin(1) * t137;
t5 = pkin(4) * t30 + qJ(5) * t29 - qJD(5) * t58 + t27;
t124 = -t5 - t135;
t123 = -t14 * t77 + t5 * t84;
t7 = -pkin(4) * t69 + t142;
t122 = -t103 * t8 + t105 * t7;
t121 = -t55 * t67 + t166;
t120 = t67 * t84 - t69 * t77;
t117 = t62 + (-t105 * t73 + t139) * t69;
t116 = -t69 * t140 + t73 * t151 + t63;
t114 = t14 * t58 + t127;
t113 = t103 * t31 + t105 * t28 + t32 * t139 - t44 * t140;
t38 = t101 * t71 + t102 * t70;
t40 = pkin(3) * t74 + pkin(7) * t77 + t136;
t112 = t103 * t40 + t105 * t38 + t47 * t139 - t55 * t140;
t111 = t69 * t43 - t160;
t82 = -pkin(4) * t105 - qJ(5) * t103 + t96;
t19 = pkin(4) * t58 + qJ(5) * t56;
t18 = t118 * t84 + t54;
t16 = -pkin(4) * t83 + t103 * t55 - t105 * t47;
t15 = qJ(5) * t83 + t155;
t11 = t56 * t69 - t29;
t10 = -pkin(4) * t75 + t103 * t51 - t105 * t39;
t9 = qJ(5) * t75 + t156;
t6 = (-pkin(4) * t77 + qJ(5) * t148) * t103 + (qJ(5) * t77 + (pkin(4) * qJD(4) - qJD(5)) * t84) * t105 + t37;
t4 = -t74 * pkin(4) + t155 * qJD(4) + t103 * t38 - t105 * t40;
t3 = qJ(5) * t74 + qJD(5) * t83 + t112;
t2 = t127 - t169;
t1 = qJD(5) * t69 + t113 + t154;
t17 = [0, 0, 0, 0.2e1 * t104 * t130, 0.2e1 * t149 * t137, t143, -t144, 0, -pkin(6) * t143 + t104 * t125, pkin(6) * t144 + t106 * t125, t54 * t115 - t28 * t83 + t37 * t75 + t38 * t73 + t48 * t77 - t49 * t74 + t121, t27 * t54 + t28 * t55 - t48 * t37 + t49 * t38 + (t88 + t119) * t136, -t58 * t134 + (-t29 * t84 - t58 * t77) * t105, -(-t103 * t58 - t105 * t56) * t77 + (t150 - t105 * t30 + (t103 * t56 - t105 * t58) * qJD(4)) * t84, t120 * t105 - t69 * t134 - t29 * t83 + t58 * t74, -t120 * t103 - t69 * t133 - t30 * t83 - t56 * t74, t67 * t83 + t69 * t74, -t127 * t83 + t12 * t74 + t37 * t56 + t54 * t30 + ((-qJD(4) * t55 + t40) * t69 + t47 * t67 + t43 * t148) * t105 + ((-qJD(4) * t47 - t38) * t69 - t165 + t121) * t103, -t112 * t69 - t155 * t67 - t113 * t83 - t13 * t74 + t37 * t58 - t54 * t29 - t43 * t134 + (-t165 + t166) * t105, t123 * t103 + t14 * t133 - t16 * t67 + t18 * t30 - t2 * t83 - t4 * t69 + t56 * t6 - t7 * t74, -t15 * t30 - t16 * t29 - t3 * t56 + t4 * t58 - t122 * t77 + (-t1 * t103 + t105 * t2 + (-t103 * t7 - t105 * t8) * qJD(4)) * t84, t1 * t83 - t123 * t105 + t14 * t134 + t15 * t67 + t18 * t29 + t3 * t69 - t58 * t6 + t74 * t8, t1 * t15 + t14 * t6 + t16 * t2 + t18 * t5 + t3 * t8 + t4 * t7; 0, 0, 0, -t104 * t145, -t149 * t108, 0, 0, 0, t108 * pkin(1) * t104, pkin(1) * t145, (t49 - t50) * t75 + (-t51 + t48) * t73 + (-t101 * t67 - t102 * t115) * pkin(2), t48 * t50 - t49 * t51 + (t101 * t28 - t102 * t27 - t88 * t141) * pkin(2), t105 * t128 - t150, (-t29 + t164) * t105 - t58 * t151 + t157, t117 - t162, t116 + t161, -t69 * t75, -t12 * t75 + t96 * t30 - t50 * t56 + (-t27 + (-t39 - t147) * t69) * t105 + (t51 * t69 + t111) * t103, -t96 * t29 + t156 * t69 + t13 * t75 - t50 * t58 + (t27 + t135) * t103 + t111 * t105, t10 * t69 + t172 * t103 + t124 * t105 - t158 * t56 + t30 * t82 + t7 * t75, -t10 * t58 + t56 * t9 + (-t30 * t95 - t7 * t73 + t1 + (t58 * t95 + t7) * qJD(4)) * t105 + (-t29 * t95 + t73 * t8 + t2 + (t56 * t95 - t8) * qJD(4)) * t103, t124 * t103 - t172 * t105 + t158 * t58 + t29 * t82 - t69 * t9 - t75 * t8, -t10 * t7 + t5 * t82 - t8 * t9 - t158 * t14 + (t122 * qJD(4) + t1 * t105 + t103 * t2) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 ^ 2 - t75 ^ 2, t48 * t75 - t49 * t73 + t93, 0, 0, 0, 0, 0, t116 - t161, -t170 * t105 - t162 - t62, -t151 * t69 - t161 + t63, (t29 + t164) * t105 + t103 * t128 + t157, t117 + t162, -t14 * t75 + (-t2 + t168) * t105 + (t69 * t7 + t1) * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, -t56 ^ 2 + t171, t11, t128 - t30, t67, -t43 * t58 - t127 + t167, t12 * t69 + t43 * t56 - t113, -t19 * t56 - t114 + t167 + 0.2e1 * t169, pkin(4) * t29 - qJ(5) * t30 + (-t13 + t8) * t58 + (t7 - t142) * t56, 0.2e1 * t154 - t14 * t56 + t19 * t58 + (0.2e1 * qJD(5) - t12) * t69 + t113, -t2 * pkin(4) + t1 * qJ(5) - t7 * t13 - t14 * t19 + t142 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t75 + t163, t11, -t170 - t171, t114 - t168 - t169;];
tauc_reg = t17;
