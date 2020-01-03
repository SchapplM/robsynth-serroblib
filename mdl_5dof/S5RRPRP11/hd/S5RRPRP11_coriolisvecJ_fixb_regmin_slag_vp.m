% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP11
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:04
% EndTime: 2019-12-31 20:14:09
% DurationCPUTime: 1.48s
% Computational Cost: add. (1594->267), mult. (3596->361), div. (0->0), fcn. (1939->4), ass. (0->143)
t167 = pkin(3) + pkin(6);
t88 = cos(qJ(4));
t135 = t88 * qJD(2);
t89 = cos(qJ(2));
t143 = qJD(1) * t89;
t86 = sin(qJ(4));
t51 = -t86 * t143 + t135;
t133 = qJD(1) * qJD(2);
t176 = -0.2e1 * t133;
t121 = t89 * t133;
t116 = pkin(4) * t121;
t139 = qJD(4) * t88;
t140 = qJD(4) * t86;
t87 = sin(qJ(2));
t122 = t87 * t133;
t73 = pkin(2) * t122;
t111 = pkin(7) * t87 - qJ(3) * t89;
t136 = t87 * qJD(3);
t95 = t111 * qJD(2) - t136;
t21 = t95 * qJD(1) + t73;
t120 = -t87 * qJ(3) - pkin(1);
t90 = -pkin(2) - pkin(7);
t47 = t90 * t89 + t120;
t32 = t47 * qJD(1);
t137 = t87 * qJD(1);
t75 = pkin(6) * t137;
t171 = qJD(3) + t75;
t134 = pkin(3) * t137 + t171;
t34 = t90 * qJD(2) + t134;
t72 = pkin(6) * t121;
t46 = pkin(3) * t121 + t72;
t119 = -t32 * t139 - t34 * t140 - t86 * t21 + t88 * t46;
t2 = -t116 - t119;
t10 = t88 * t32 + t86 * t34;
t74 = qJD(4) + t137;
t7 = t74 * qJ(5) + t10;
t175 = t7 * t74 - t2;
t65 = t167 * t87;
t174 = t88 * t47 + t86 * t65;
t138 = t86 * qJD(2);
t49 = t88 * t143 + t138;
t159 = t49 * t74;
t24 = t49 * qJD(4) - t86 * t122;
t173 = t24 - t159;
t157 = t51 * t74;
t25 = t51 * qJD(4) - t88 * t122;
t172 = -t25 + t157;
t125 = t74 * t139;
t153 = t87 * t88;
t132 = t74 * t153;
t170 = qJD(1) * (t89 * t138 + t132) + qJD(2) * t51 + t125;
t142 = qJD(2) * t87;
t78 = pkin(2) * t142;
t29 = t78 + t95;
t141 = qJD(2) * t89;
t58 = t167 * t141;
t169 = -qJD(4) * t174 - t86 * t29 + t88 * t58;
t168 = t51 ^ 2;
t56 = t167 * t142;
t82 = qJD(2) * qJD(3);
t36 = -qJD(1) * t56 + t82;
t4 = t25 * pkin(4) + t24 * qJ(5) - t51 * qJD(5) + t36;
t166 = t4 * t86;
t165 = t4 * t88;
t76 = pkin(6) * t143;
t57 = pkin(3) * t143 + t76;
t83 = qJD(2) * qJ(3);
t42 = t83 + t57;
t12 = t49 * pkin(4) - t51 * qJ(5) + t42;
t163 = t12 * t51;
t162 = t24 * t88;
t161 = t36 * t86;
t160 = t36 * t88;
t158 = t51 * t49;
t156 = t51 * t89;
t155 = t74 * t87;
t154 = t74 * t90;
t92 = qJD(1) ^ 2;
t152 = t89 * t92;
t91 = qJD(2) ^ 2;
t151 = t91 * t87;
t150 = t91 * t89;
t113 = pkin(4) * t88 + qJ(5) * t86;
t103 = -pkin(3) - t113;
t149 = -t113 * qJD(4) + t88 * qJD(5) + t103 * t137 - t171;
t79 = pkin(2) * t137;
t39 = t111 * qJD(1) + t79;
t148 = t88 * t39 + t86 * t57;
t66 = t167 * t89;
t84 = t87 ^ 2;
t85 = t89 ^ 2;
t146 = t84 - t85;
t145 = qJD(2) * pkin(2);
t9 = -t86 * t32 + t88 * t34;
t144 = qJD(5) - t9;
t63 = -t89 * pkin(2) + t120;
t43 = qJD(1) * t63;
t131 = t86 * t154;
t130 = t87 * t152;
t129 = t88 * t154;
t128 = -t34 * t139 - t88 * t21 - t86 * t46;
t127 = t90 * t141;
t126 = t74 * t140;
t124 = t89 * t139;
t118 = pkin(1) * t176;
t117 = qJD(3) - t145;
t68 = t88 * t121;
t115 = qJ(5) * t121;
t6 = -t74 * pkin(4) + t144;
t114 = t6 * t86 + t7 * t88;
t112 = -t86 * pkin(4) + t88 * qJ(5);
t109 = -t86 * t39 + t88 * t57;
t107 = -t86 * t47 + t88 * t65;
t106 = -qJD(1) * t85 + t155;
t105 = -0.2e1 * qJD(2) * t43;
t104 = t74 * t86;
t97 = -t89 * t83 - t136;
t30 = t97 * qJD(1) + t73;
t41 = t78 + t97;
t102 = pkin(6) * t91 + qJD(1) * t41 + t30;
t101 = t10 * t74 + t119;
t99 = t32 * t140 + t128;
t98 = t65 * t139 - t47 * t140 + t88 * t29 + t86 * t58;
t94 = -qJD(2) * t49 - t74 * t104 + t68;
t59 = pkin(6) * t122 - t82;
t61 = t117 + t75;
t64 = -t76 - t83;
t93 = -t59 * t89 + (t61 * t89 + (t64 + t76) * t87) * qJD(2);
t62 = qJ(3) - t112;
t60 = t90 * t68;
t54 = -qJ(3) * t143 + t79;
t33 = t43 * t137;
t28 = t113 * t89 + t66;
t19 = t51 * pkin(4) + t49 * qJ(5);
t16 = -t87 * pkin(4) - t107;
t15 = t87 * qJ(5) + t174;
t14 = -pkin(4) * t143 - t109;
t13 = qJ(5) * t143 + t148;
t8 = (t112 * qJD(4) + qJD(5) * t86) * t89 + (-pkin(6) + t103) * t142;
t5 = -pkin(4) * t141 - t169;
t3 = qJ(5) * t141 + t87 * qJD(5) + t98;
t1 = t74 * qJD(5) + t115 - t99;
t11 = [0, 0, 0, 0.2e1 * t87 * t121, t146 * t176, t150, -t151, 0, -pkin(6) * t150 + t87 * t118, pkin(6) * t151 + t89 * t118, t93, t102 * t89 + t87 * t105, -t102 * t87 + t105 * t89, pkin(6) * t93 + t30 * t63 + t43 * t41, t24 * t86 * t89 + (t138 * t87 - t124) * t51, (-t49 * t86 + t51 * t88) * t142 + (t162 + t86 * t25 + (t49 * t88 + t51 * t86) * qJD(4)) * t89, -t74 * t124 - t24 * t87 + (t106 * t86 + t156) * qJD(2), t89 * t126 - t25 * t87 + (t106 * t88 - t49 * t89) * qJD(2), (t74 + t137) * t141, t169 * t74 - t56 * t49 + t66 * t25 + (-t135 * t42 + t119) * t87 + (-t42 * t140 + t160 + (qJD(1) * t107 + t9) * qJD(2)) * t89, -t98 * t74 - t56 * t51 - t66 * t24 + ((qJD(2) * t42 + qJD(4) * t32) * t86 + t128) * t87 + (-t42 * t139 - t161 + (-qJD(1) * t174 - t10) * qJD(2)) * t89, t28 * t25 + t8 * t49 - t5 * t74 + (-t12 * t135 - t2) * t87 + (-t12 * t140 + t165 + (-qJD(1) * t16 - t6) * qJD(2)) * t89, -t15 * t25 - t16 * t24 - t3 * t49 + t5 * t51 + t114 * t142 + (-t1 * t88 - t2 * t86 + (-t6 * t88 + t7 * t86) * qJD(4)) * t89, t28 * t24 + t3 * t74 - t8 * t51 + (-t12 * t138 + t1) * t87 + (t12 * t139 + t166 + (qJD(1) * t15 + t7) * qJD(2)) * t89, t1 * t15 + t12 * t8 + t2 * t16 + t4 * t28 + t7 * t3 + t6 * t5; 0, 0, 0, -t130, t146 * t92, 0, 0, 0, t92 * pkin(1) * t87, pkin(1) * t152, ((-t64 - t83) * t87 + (t117 - t61) * t89) * qJD(1), -t54 * t143 + t33, 0.2e1 * t82 + (t43 * t89 + t54 * t87) * qJD(1), -t59 * qJ(3) - t64 * qJD(3) - t43 * t54 + (-t64 * t87 + (-t61 - t145) * t89) * qJD(1) * pkin(6), -t104 * t51 - t162, (-t25 - t157) * t88 + (t24 + t159) * t86, -t126 + t68 + (-t155 * t86 - t156) * qJD(1), -t125 + (-t132 + (t49 - t138) * t89) * qJD(1), -t74 * t143, t60 + qJ(3) * t25 + t161 - t109 * t74 + t134 * t49 + (t42 * t88 - t131) * qJD(4) + (t153 * t42 - t9 * t89) * qJD(1), -qJ(3) * t24 + t160 + t148 * t74 + t134 * t51 + (-t42 * t86 - t129) * qJD(4) + (t10 * t89 + (-t42 * t87 - t127) * t86) * qJD(1), t14 * t74 + t62 * t25 + t166 + t60 - t149 * t49 + (t12 * t88 - t131) * qJD(4) + (t12 * t153 + t6 * t89) * qJD(1), t13 * t49 - t14 * t51 + (-t7 * t137 + t24 * t90 + t2 + (-t49 * t90 - t7) * qJD(4)) * t88 + (-t6 * t137 - t25 * t90 - t1 + (t51 * t90 - t6) * qJD(4)) * t86, -t13 * t74 + t62 * t24 - t165 + t149 * t51 + (t12 * t86 + t129) * qJD(4) + (-t7 * t89 + (t12 * t87 + t127) * t86) * qJD(1), -t7 * t13 - t6 * t14 + t4 * t62 - t149 * t12 + (qJD(4) * t114 + t1 * t86 - t2 * t88) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t84 * t92 - t91, t64 * qJD(2) + t33 + t72, 0, 0, 0, 0, 0, t94, -t170, t94, t172 * t86 + t173 * t88, t170, -t12 * qJD(2) + t175 * t88 + (t6 * t74 + t1) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, -t49 ^ 2 + t168, -t173, t172, t121, -t42 * t51 + t101, t42 * t49 + t9 * t74 + t99, -t19 * t49 + t101 + 0.2e1 * t116 - t163, pkin(4) * t24 - t25 * qJ(5) + (-t10 + t7) * t51 + (t6 - t144) * t49, 0.2e1 * t115 - t12 * t49 + t19 * t51 + (0.2e1 * qJD(5) - t9) * t74 - t99, -t2 * pkin(4) + t1 * qJ(5) - t6 * t10 - t12 * t19 + t144 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121 + t158, -t173, -t74 ^ 2 - t168, t163 - t175;];
tauc_reg = t11;
