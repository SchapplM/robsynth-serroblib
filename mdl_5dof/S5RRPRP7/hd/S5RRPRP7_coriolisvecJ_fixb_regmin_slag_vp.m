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
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:42
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 20:41:23
% EndTime: 2021-01-15 20:41:31
% DurationCPUTime: 1.66s
% Computational Cost: add. (2795->292), mult. (7150->393), div. (0->0), fcn. (4931->6), ass. (0->139)
t103 = sin(pkin(8));
t105 = sin(qJ(2));
t144 = qJD(1) * t105;
t107 = cos(qJ(2));
t154 = cos(pkin(8));
t130 = t154 * t107;
t93 = qJD(1) * t130;
t73 = t103 * t144 - t93;
t70 = qJD(4) + t73;
t140 = qJD(1) * qJD(2);
t178 = -0.2e1 * t140;
t163 = -qJ(3) - pkin(6);
t91 = t163 * t107;
t88 = qJD(1) * t91;
t79 = t103 * t88;
t157 = qJD(2) * pkin(2);
t90 = t163 * t105;
t87 = qJD(1) * t90;
t82 = t87 + t157;
t48 = t154 * t82 + t79;
t43 = -qJD(2) * pkin(3) - t48;
t104 = sin(qJ(4));
t106 = cos(qJ(4));
t141 = t106 * qJD(2);
t85 = t103 * t107 + t154 * t105;
t153 = qJD(1) * t85;
t56 = t104 * t153 - t141;
t58 = qJD(2) * t104 + t106 * t153;
t14 = t56 * pkin(4) - t58 * qJ(5) + t43;
t75 = t85 * qJD(2);
t68 = qJD(1) * t75;
t96 = pkin(2) * t103 + pkin(7);
t164 = t96 * t68;
t177 = t14 * t70 - t164;
t133 = t105 * t140;
t92 = t103 * t133;
t112 = qJD(2) * t93 - t92;
t30 = qJD(4) * t58 + t104 * t112;
t176 = t58 ^ 2;
t175 = t70 ^ 2;
t174 = pkin(4) * t68;
t99 = -pkin(2) * t107 - pkin(1);
t152 = qJD(1) * t99;
t89 = qJD(3) + t152;
t32 = t73 * pkin(3) - pkin(7) * t153 + t89;
t134 = t154 * t88;
t49 = t103 * t82 - t134;
t44 = qJD(2) * pkin(7) + t49;
t13 = t104 * t32 + t106 * t44;
t8 = qJ(5) * t70 + t13;
t173 = t70 * t8;
t172 = pkin(2) * t105;
t171 = t13 * t70;
t132 = qJD(2) * t163;
t71 = t107 * qJD(3) + t105 * t132;
t65 = t71 * qJD(1);
t111 = -t105 * qJD(3) + t107 * t132;
t66 = t111 * qJD(1);
t27 = t103 * t65 - t154 * t66;
t170 = t27 * t85;
t113 = -t103 * t105 + t130;
t78 = t113 * qJD(2);
t169 = t43 * t78;
t168 = t56 * t73;
t167 = t58 * t56;
t131 = t58 * t70;
t166 = t58 * t153;
t165 = t153 * t56;
t120 = pkin(4) * t104 - qJ(5) * t106;
t50 = t103 * t87 - t134;
t162 = t104 * qJD(5) - t70 * t120 + t50;
t142 = qJD(4) * t106;
t161 = -t104 * t30 - t56 * t142;
t139 = pkin(2) * t144;
t39 = pkin(3) * t153 + pkin(7) * t73 + t139;
t51 = t154 * t87 + t79;
t160 = t104 * t39 + t106 * t51;
t47 = -pkin(3) * t113 - pkin(7) * t85 + t99;
t55 = t103 * t90 - t154 * t91;
t159 = t104 * t47 + t106 * t55;
t158 = qJ(5) * t68;
t62 = t104 * t68;
t156 = t104 * t70;
t63 = t106 * t68;
t143 = qJD(4) * t104;
t29 = -qJD(4) * t141 - t106 * t112 + t143 * t153;
t155 = t29 * t104;
t151 = qJD(4) * t85;
t150 = qJD(4) * t96;
t109 = qJD(1) ^ 2;
t149 = t107 * t109;
t108 = qJD(2) ^ 2;
t148 = t108 * t105;
t147 = t108 * t107;
t12 = -t104 * t44 + t106 * t32;
t146 = qJD(5) - t12;
t145 = t105 ^ 2 - t107 ^ 2;
t138 = t105 * t157;
t137 = t70 * t150;
t136 = t85 * t143;
t135 = t85 * t142;
t37 = t103 * t71 - t154 * t111;
t54 = -t103 * t91 - t154 * t90;
t28 = t103 * t66 + t154 * t65;
t94 = pkin(2) * t133;
t31 = t68 * pkin(3) - t112 * pkin(7) + t94;
t129 = t104 * t28 - t106 * t31 + t44 * t142 + t32 * t143;
t128 = 0.2e1 * t153;
t126 = pkin(1) * t178;
t5 = pkin(4) * t30 + qJ(5) * t29 - qJD(5) * t58 + t27;
t125 = -t5 - t137;
t97 = -t154 * pkin(2) - pkin(3);
t124 = t14 * t78 + t5 * t85;
t7 = -pkin(4) * t70 + t146;
t123 = -t104 * t8 + t106 * t7;
t122 = -t55 * t68 + t170;
t121 = t68 * t85 + t70 * t78;
t119 = t62 + (t106 * t73 + t142) * t70;
t118 = -t70 * t143 - t73 * t156 + t63;
t117 = t14 * t58 + t129;
t116 = t104 * t31 + t106 * t28 + t32 * t142 - t44 * t143;
t38 = t103 * t111 + t154 * t71;
t40 = pkin(3) * t75 - pkin(7) * t78 + t138;
t115 = t104 * t40 + t106 * t38 + t47 * t142 - t55 * t143;
t114 = t70 * t43 - t164;
t83 = -t106 * pkin(4) - t104 * qJ(5) + t97;
t19 = pkin(4) * t58 + qJ(5) * t56;
t18 = t120 * t85 + t54;
t16 = pkin(4) * t113 + t104 * t55 - t106 * t47;
t15 = -qJ(5) * t113 + t159;
t11 = t56 * t70 - t29;
t10 = -pkin(4) * t153 + t104 * t51 - t106 * t39;
t9 = qJ(5) * t153 + t160;
t6 = (pkin(4) * t78 + qJ(5) * t151) * t104 + (-qJ(5) * t78 + (pkin(4) * qJD(4) - qJD(5)) * t85) * t106 + t37;
t4 = -t75 * pkin(4) + t159 * qJD(4) + t104 * t38 - t106 * t40;
t3 = qJ(5) * t75 - qJD(5) * t113 + t115;
t2 = t129 - t174;
t1 = qJD(5) * t70 + t116 + t158;
t17 = [0, 0, 0, 0.2e1 * t107 * t133, t145 * t178, t147, -t148, 0, -pkin(6) * t147 + t105 * t126, pkin(6) * t148 + t107 * t126, t68 * t99 + t75 * t89 + (-t37 + (-qJD(1) * t113 + t73) * t172) * qJD(2), t89 * t78 - t99 * t92 + (t128 * t172 + t99 * t93 - t38) * qJD(2), t112 * t54 + t113 * t28 + t153 * t37 - t38 * t73 - t48 * t78 - t49 * t75 + t122, t27 * t54 + t28 * t55 - t37 * t48 + t38 * t49 + (t89 + t152) * t138, -t58 * t136 + (-t29 * t85 + t58 * t78) * t106, (-t104 * t58 - t106 * t56) * t78 + (t155 - t106 * t30 + (t104 * t56 - t106 * t58) * qJD(4)) * t85, t106 * t121 + t113 * t29 - t70 * t136 + t58 * t75, -t104 * t121 + t113 * t30 - t135 * t70 - t56 * t75, -t113 * t68 + t70 * t75, t129 * t113 + t12 * t75 + t37 * t56 + t54 * t30 + ((-qJD(4) * t55 + t40) * t70 + t47 * t68 + t43 * t151) * t106 + ((-qJD(4) * t47 - t38) * t70 + t169 + t122) * t104, -t115 * t70 - t159 * t68 + t116 * t113 - t13 * t75 + t37 * t58 - t54 * t29 - t43 * t136 + (t169 + t170) * t106, t104 * t124 + t113 * t2 + t135 * t14 - t16 * t68 + t18 * t30 - t4 * t70 + t56 * t6 - t7 * t75, -t15 * t30 - t16 * t29 - t3 * t56 + t4 * t58 + t123 * t78 + (-t1 * t104 + t106 * t2 + (-t104 * t7 - t106 * t8) * qJD(4)) * t85, -t1 * t113 - t106 * t124 + t136 * t14 + t15 * t68 + t18 * t29 + t3 * t70 - t58 * t6 + t75 * t8, t1 * t15 + t14 * t6 + t16 * t2 + t18 * t5 + t3 * t8 + t4 * t7; 0, 0, 0, -t105 * t149, t145 * t109, 0, 0, 0, t109 * pkin(1) * t105, pkin(1) * t149, qJD(2) * t50 - t73 * t139 - t153 * t89 - t27, t51 * qJD(2) - t139 * t153 + t89 * t73 - t28, (t49 - t50) * t153 + (t51 - t48) * t73 + (-t103 * t68 - t154 * t112) * pkin(2), t48 * t50 - t49 * t51 + (t103 * t28 - t89 * t144 - t154 * t27) * pkin(2), t106 * t131 - t155, (-t29 - t168) * t106 - t58 * t156 + t161, t119 - t166, t118 + t165, -t70 * t153, -t12 * t153 + t97 * t30 - t50 * t56 + (-t27 + (-t39 - t150) * t70) * t106 + (t51 * t70 + t114) * t104, -t97 * t29 + t160 * t70 + t13 * t153 - t50 * t58 + (t27 + t137) * t104 + t114 * t106, t10 * t70 + t177 * t104 + t125 * t106 + t153 * t7 - t162 * t56 + t30 * t83, -t10 * t58 + t56 * t9 + (-t30 * t96 + t7 * t73 + t1 + (t58 * t96 + t7) * qJD(4)) * t106 + (-t29 * t96 - t73 * t8 + t2 + (t56 * t96 - t8) * qJD(4)) * t104, t125 * t104 - t177 * t106 - t153 * t8 + t162 * t58 + t29 * t83 - t70 * t9, -t10 * t7 + t5 * t83 - t8 * t9 - t162 * t14 + (qJD(4) * t123 + t1 * t106 + t104 * t2) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128 * qJD(2), -t92 + (t93 - t73) * qJD(2), -t153 ^ 2 - t73 ^ 2, t153 * t48 + t49 * t73 + t94, 0, 0, 0, 0, 0, t118 - t165, -t106 * t175 - t166 - t62, -t156 * t70 - t165 + t63, (t29 - t168) * t106 + t104 * t131 + t161, t119 + t166, -t14 * t153 + (-t2 + t173) * t106 + (t7 * t70 + t1) * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, -t56 ^ 2 + t176, t11, t131 - t30, t68, -t43 * t58 - t129 + t171, t12 * t70 + t43 * t56 - t116, -t19 * t56 - t117 + t171 + 0.2e1 * t174, pkin(4) * t29 - qJ(5) * t30 + (-t13 + t8) * t58 + (t7 - t146) * t56, 0.2e1 * t158 - t14 * t56 + t19 * t58 + (0.2e1 * qJD(5) - t12) * t70 + t116, -pkin(4) * t2 + qJ(5) * t1 - t13 * t7 - t14 * t19 + t146 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t153 + t167, t11, -t175 - t176, t117 - t173 - t174;];
tauc_reg = t17;
