% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:22
% EndTime: 2019-03-09 02:11:26
% DurationCPUTime: 1.42s
% Computational Cost: add. (1780->270), mult. (3424->358), div. (0->0), fcn. (1833->4), ass. (0->139)
t74 = cos(qJ(4));
t136 = qJD(1) * t74;
t73 = cos(qJ(5));
t110 = t73 * t136;
t71 = sin(qJ(5));
t128 = t71 * qJD(4);
t48 = t110 + t128;
t72 = sin(qJ(4));
t127 = t72 * qJD(1);
t61 = qJD(5) + t127;
t151 = t48 * t61;
t132 = qJD(5) * t48;
t122 = qJD(1) * qJD(4);
t108 = t72 * t122;
t56 = t71 * t108;
t21 = -t56 + t132;
t172 = t21 - t151;
t130 = qJD(5) * t72;
t102 = qJD(1) + t130;
t126 = t73 * qJD(4);
t131 = qJD(5) * t71;
t112 = t74 * t131;
t81 = t72 * t126 + t112;
t20 = t81 * qJD(1) - qJD(5) * t126;
t143 = t74 * t20;
t149 = t61 * t71;
t89 = t61 * t73;
t171 = qJD(4) * ((-t48 + t110) * t72 + t74 * t89) - t102 * t149 - t143;
t46 = t71 * t136 - t126;
t170 = t46 * t89;
t107 = t74 * t122;
t101 = pkin(5) * t107;
t129 = qJD(5) * t73;
t123 = qJD(1) * qJD(2);
t133 = qJD(4) * t74;
t63 = qJD(1) * qJ(2) + qJD(3);
t58 = -pkin(7) * qJD(1) + t63;
t31 = t72 * t123 + t58 * t133;
t70 = pkin(1) + qJ(3);
t52 = t72 * pkin(4) - t74 * pkin(8) + t70;
t32 = t52 * qJD(1) - qJD(2);
t99 = pkin(4) * t74 + pkin(8) * t72;
t44 = t99 * qJD(4) + qJD(3);
t33 = t44 * qJD(1);
t51 = t72 * t58;
t39 = qJD(4) * pkin(8) + t51;
t106 = -t39 * t129 - t32 * t131 - t71 * t31 + t73 * t33;
t2 = -t101 - t106;
t12 = t71 * t32 + t73 * t39;
t8 = t61 * qJ(6) + t12;
t169 = -t61 * t8 + t2;
t168 = qJD(1) * t70;
t146 = t72 * t73;
t119 = t61 * t146;
t150 = t48 * t74;
t111 = t61 * t129;
t97 = t71 * t107 + t111;
t166 = (t119 + t150) * qJD(1) + t97;
t164 = t48 ^ 2;
t134 = qJD(4) * t72;
t30 = -t74 * t123 + t58 * t134;
t3 = t21 * pkin(5) + t20 * qJ(6) - t48 * qJD(6) + t30;
t163 = t3 * t71;
t162 = t3 * t73;
t142 = t74 * t58;
t40 = -qJD(4) * pkin(4) - t142;
t10 = t46 * pkin(5) - t48 * qJ(6) + t40;
t160 = t10 * t48;
t159 = t20 * t71;
t158 = t20 * t72;
t157 = t21 * t72;
t156 = t30 * t71;
t155 = t30 * t73;
t154 = t40 * t71;
t153 = t40 * t73;
t152 = t48 * t46;
t69 = -pkin(7) + qJ(2);
t148 = t69 * t71;
t147 = t71 * t72;
t50 = t99 * qJD(1);
t145 = t73 * t50;
t144 = t73 * t52;
t93 = pkin(5) * t71 - qJ(6) * t73;
t141 = t71 * qJD(6) - t61 * t93 + t51;
t140 = t73 * t142 + t71 * t50;
t139 = t69 * t146 + t71 * t52;
t68 = t74 ^ 2;
t138 = t72 ^ 2 - t68;
t75 = qJD(4) ^ 2;
t76 = qJD(1) ^ 2;
t137 = -t75 - t76;
t135 = qJD(2) * t72;
t59 = -qJD(2) + t168;
t125 = qJD(2) - t59;
t11 = t73 * t32 - t71 * t39;
t124 = qJD(6) - t11;
t121 = pkin(8) * t149;
t120 = pkin(8) * t89;
t117 = t61 * t147;
t116 = -t32 * t129 - t73 * t31 - t71 * t33;
t115 = pkin(8) * t133;
t114 = t74 * t126;
t113 = t61 * t131;
t65 = 0.2e1 * t123;
t109 = 0.2e1 * qJD(3) * qJD(1);
t105 = t69 * t114 + t52 * t129 + t73 * t135 + t71 * t44;
t104 = t125 * qJD(1);
t103 = -0.2e1 * t107;
t100 = qJ(6) * t107;
t7 = -t61 * pkin(5) + t124;
t96 = t7 * t73 - t71 * t8;
t95 = t7 * t71 + t73 * t8;
t94 = t73 * pkin(5) + t71 * qJ(6);
t92 = qJD(2) + t59 + t168;
t91 = qJD(1) * t68 - t61 * t72;
t88 = -t69 * t75 + t109;
t86 = -t69 + t93;
t85 = -t10 * t72 + t115;
t84 = t40 * t72 - t115;
t83 = t12 * t61 + t106;
t82 = t39 * t131 + t116;
t1 = t61 * qJD(6) + t100 - t82;
t80 = -t95 * qJD(5) - t1 * t71 + t2 * t73;
t79 = t96 * qJD(5) + t1 * t73 + t2 * t71;
t78 = t113 + t46 * t136 + (-t114 + t117) * qJD(1);
t77 = t46 * t134 + (-t21 - t56) * t74 + (-t102 * t73 - t74 * t128) * t61;
t53 = -pkin(4) - t94;
t28 = t86 * t74;
t19 = t20 * t73;
t18 = -t144 + (-pkin(5) + t148) * t72;
t17 = t72 * qJ(6) + t139;
t15 = t48 * pkin(5) + t46 * qJ(6);
t14 = -t145 + (-pkin(5) * qJD(1) + t58 * t71) * t74;
t13 = qJ(6) * t136 + t140;
t9 = t46 * t61 - t20;
t6 = -t86 * t134 + (t94 * qJD(5) - qJD(6) * t73 - qJD(2)) * t74;
t5 = -pkin(5) * t133 + (t69 * t130 - t44) * t73 + (qJD(5) * t52 + t69 * t133 + t135) * t71;
t4 = qJ(6) * t133 + (-t69 * t131 + qJD(6)) * t72 + t105;
t16 = [0, 0, 0, 0, t65, qJ(2) * t65, t65, t109, t63 * qJD(2) + t59 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t70) * qJD(1), t72 * t103, 0.2e1 * t138 * t122, -t75 * t72, -t75 * t74, 0, t92 * t133 + t88 * t72, -t92 * t134 + t88 * t74, -t73 * t143 - t81 * t48 (t46 * t73 + t48 * t71) * t134 + (t159 - t21 * t73 + (t46 * t71 - t48 * t73) * qJD(5)) * t74, -t61 * t112 - t158 + (t91 * t73 + t150) * qJD(4), -t74 * t111 - t157 + (-t74 * t46 - t91 * t71) * qJD(4) (t61 + t127) * t133 (-t52 * t131 + t73 * t44) * t61 + ((-qJD(2) * t71 - t69 * t129) * t61 + (t69 * t46 - t154) * qJD(4) + t106) * t72 + (t40 * t129 - qJD(2) * t46 - t69 * t21 + t156 + (-t61 * t148 + (-t69 * t147 + t144) * qJD(1) + t11) * qJD(4)) * t74, -t105 * t61 + ((t61 * t69 + t39) * t131 + (t48 * t69 - t153) * qJD(4) + t116) * t72 + (-t40 * t131 - qJD(2) * t48 + t69 * t20 + t155 + (-t139 * qJD(1) - t12) * qJD(4)) * t74, t28 * t21 + t6 * t46 - t5 * t61 + (-t10 * t128 - t2) * t72 + (t10 * t129 + t163 + (-qJD(1) * t18 - t7) * qJD(4)) * t74, -t96 * t134 - t17 * t21 - t18 * t20 - t4 * t46 + t5 * t48 + t80 * t74, t28 * t20 + t4 * t61 - t6 * t48 + (t10 * t126 + t1) * t72 + (t10 * t131 - t162 + (qJD(1) * t17 + t8) * qJD(4)) * t74, t1 * t17 + t10 * t6 + t18 * t2 + t28 * t3 + t4 * t8 + t5 * t7; 0, 0, 0, 0, -t76, -t76 * qJ(2), -t76, 0 (-qJD(3) - t63) * qJD(1), 0, 0, 0, 0, 0, t103, 0.2e1 * t108, 0, 0, 0, 0, 0, t78, t166, t78, t172 * t71 + t170 - t19, -t166 (t10 * t74 - t72 * t95) * qJD(1) + t80; 0, 0, 0, 0, 0, 0, 0, -t76, t104, 0, 0, 0, 0, 0, t137 * t72, t137 * t74, 0, 0, 0, 0, 0, t77, -t171, t77 (t102 * t48 - t46 * t133 - t157) * t73 + (t102 * t46 + t48 * t133 - t158) * t71, t171, t96 * qJD(1) + (qJD(4) * t95 - t3) * t74 + (qJD(4) * t10 + t79) * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t76 * t72, -t138 * t76, 0, 0, 0, t74 * t104, -t125 * t127, t48 * t89 - t159, -t19 - t170 + (-t21 - t151) * t71 (t119 - t150) * qJD(1) + t97, -t113 + (-t117 + (t46 + t126) * t74) * qJD(1), -t61 * t136, -pkin(4) * t21 - t155 - (-t71 * t142 + t145) * t61 - t46 * t51 + (-t120 + t154) * qJD(5) + (-t11 * t74 + t84 * t71) * qJD(1), pkin(4) * t20 + t156 + t140 * t61 - t48 * t51 + (t121 + t153) * qJD(5) + (t12 * t74 + t84 * t73) * qJD(1), t14 * t61 + t53 * t21 - t162 - t141 * t46 + (t10 * t71 - t120) * qJD(5) + (t7 * t74 - t85 * t71) * qJD(1), t13 * t46 - t14 * t48 + (t1 + t61 * t7 + (-t21 + t132) * pkin(8)) * t73 + ((qJD(5) * t46 - t20) * pkin(8) + t169) * t71, -t13 * t61 + t53 * t20 - t163 + t141 * t48 + (-t10 * t73 - t121) * qJD(5) + (t85 * t73 - t74 * t8) * qJD(1), pkin(8) * t79 - t10 * t141 - t8 * t13 - t7 * t14 + t3 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, -t46 ^ 2 + t164, t9, -t172, t107, -t40 * t48 + t83, t11 * t61 + t40 * t46 + t82, -t15 * t46 + 0.2e1 * t101 - t160 + t83, pkin(5) * t20 - t21 * qJ(6) + (-t12 + t8) * t48 + (t7 - t124) * t46, 0.2e1 * t100 - t10 * t46 + t15 * t48 + (0.2e1 * qJD(6) - t11) * t61 - t82, -t2 * pkin(5) + t1 * qJ(6) - t10 * t15 - t7 * t12 + t124 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107 + t152, t9, -t61 ^ 2 - t164, t160 + t169;];
tauc_reg  = t16;
