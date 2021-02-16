% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:51
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:50:11
% EndTime: 2021-01-16 00:50:23
% DurationCPUTime: 2.55s
% Computational Cost: add. (3155->309), mult. (8404->450), div. (0->0), fcn. (7103->12), ass. (0->172)
t114 = cos(pkin(6));
t101 = qJD(1) * t114 + qJD(2);
t109 = sin(pkin(12));
t111 = sin(pkin(6));
t117 = sin(qJ(3));
t120 = cos(qJ(3));
t112 = cos(pkin(12));
t113 = cos(pkin(7));
t181 = t112 * t113;
t126 = (t109 * t120 + t117 * t181) * t111;
t110 = sin(pkin(7));
t184 = t110 * t117;
t55 = qJD(1) * t126 + t101 * t184;
t116 = sin(qJ(4));
t119 = cos(qJ(4));
t139 = pkin(4) * t116 - pkin(10) * t119;
t93 = t139 * qJD(4);
t225 = t55 - t93;
t183 = t110 * t120;
t216 = (-t109 * t117 + t120 * t181) * t111;
t224 = t114 * t183 + t216;
t223 = qJD(1) * t216 + t101 * t183;
t118 = cos(qJ(5));
t170 = qJD(4) * t116;
t115 = sin(qJ(5));
t180 = t115 * t119;
t212 = pkin(9) * t115;
t176 = qJD(1) * t111;
t158 = t112 * t176;
t142 = t113 * t158;
t159 = t109 * t176;
t54 = -t117 * t159 + t120 * (t101 * t110 + t142);
t222 = t225 * t118 - t170 * t212 - t54 * t180;
t167 = qJD(5) * t118;
t178 = t118 * t119;
t96 = -pkin(4) * t119 - pkin(10) * t116 - pkin(3);
t221 = -t225 * t115 + t96 * t167 - t54 * t178;
t53 = qJD(3) * pkin(9) + t55;
t71 = t101 * t113 - t110 * t158;
t220 = -t116 * t53 + t119 * t71;
t174 = qJD(3) * t117;
t157 = t110 * t174;
t172 = qJD(3) * t120;
t51 = t101 * t157 + t142 * t174 + t159 * t172;
t219 = qJD(3) * t55 - t51;
t173 = qJD(3) * t119;
t102 = -qJD(5) + t173;
t165 = t118 * qJD(4);
t175 = qJD(3) * t116;
t87 = t115 * t175 - t165;
t200 = t102 * t87;
t151 = t119 * t165;
t168 = qJD(5) * t115;
t154 = t116 * t168;
t163 = qJD(4) * qJD(5);
t68 = -t118 * t163 + (-t151 + t154) * qJD(3);
t218 = -t68 + t200;
t171 = qJD(4) * t115;
t89 = t118 * t175 + t171;
t199 = t102 * t89;
t153 = t116 * t167;
t169 = qJD(4) * t119;
t129 = t115 * t169 + t153;
t69 = qJD(3) * t129 + t115 * t163;
t217 = t69 - t199;
t215 = t89 ^ 2;
t214 = pkin(5) * t87;
t37 = t116 * t71 + t119 * t53;
t34 = qJD(4) * pkin(10) + t37;
t45 = qJD(3) * t96 - t54;
t13 = -t115 * t34 + t118 * t45;
t10 = -qJ(6) * t89 + t13;
t9 = -pkin(5) * t102 + t10;
t213 = t10 - t9;
t211 = t54 * t87;
t210 = t54 * t89;
t209 = -qJ(6) - pkin(10);
t134 = pkin(5) * t116 - qJ(6) * t178;
t145 = qJD(5) * t209;
t92 = t139 * qJD(3);
t146 = -t115 * t220 + t118 * t92;
t208 = qJD(3) * t134 + qJD(6) * t115 - t118 * t145 + t146;
t152 = t115 * t173;
t166 = qJD(6) * t118;
t204 = t115 * t92 + t118 * t220;
t207 = -qJ(6) * t152 - t115 * t145 - t166 + t204;
t103 = pkin(9) * t178;
t187 = qJ(6) * t116;
t206 = t116 * t166 - t134 * qJD(4) - (-t103 + (-t96 + t187) * t115) * qJD(5) + t222;
t179 = t116 * t118;
t205 = -(-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t179 - (-qJD(6) * t116 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t119) * t115 - t221;
t202 = qJD(3) * pkin(3);
t197 = t115 * t45;
t14 = t118 * t34 + t197;
t11 = -qJ(6) * t87 + t14;
t201 = t102 * t11;
t50 = t223 * qJD(3);
t17 = t116 * t50 + t53 * t169 + t71 * t170;
t12 = pkin(5) * t69 + t17;
t198 = t115 * t12;
t196 = t115 * t68;
t195 = t115 * t89;
t194 = t118 * t12;
t192 = t17 * t115;
t191 = t17 * t118;
t33 = -qJD(4) * pkin(4) - t220;
t190 = t33 * t115;
t188 = t115 * t96 + t103;
t185 = t102 * t118;
t122 = qJD(3) ^ 2;
t182 = t110 * t122;
t107 = t116 ^ 2;
t177 = -t119 ^ 2 + t107;
t164 = qJD(3) * qJD(4);
t160 = t117 * t182;
t156 = t110 * t172;
t155 = t102 * t168;
t149 = t116 * t164;
t148 = -qJD(6) - t214;
t16 = qJD(4) * t220 + t119 * t50;
t41 = qJD(3) * t93 + t51;
t147 = t115 * t16 - t118 * t41;
t144 = -t115 * t41 - t118 * t16 - t45 * t167 + t34 * t168;
t143 = pkin(5) * t149;
t141 = t116 * t156;
t140 = t119 * t156;
t60 = t114 * t184 + t126;
t74 = -t110 * t111 * t112 + t113 * t114;
t44 = t116 * t74 + t119 * t60;
t23 = -t115 * t224 + t118 * t44;
t22 = -t115 * t44 - t118 * t224;
t43 = t116 * t60 - t119 * t74;
t137 = qJD(3) * t107 - t102 * t119;
t121 = qJD(4) ^ 2;
t136 = pkin(9) * t121 - t219;
t52 = -t54 - t202;
t135 = qJD(4) * (t52 + t54 - t202);
t133 = qJ(6) * t69 + t144;
t76 = t113 * t116 + t119 * t184;
t63 = -t115 * t76 - t118 * t183;
t132 = t115 * t183 - t118 * t76;
t75 = -t113 * t119 + t116 * t184;
t124 = -qJD(5) * t14 - t147;
t123 = qJ(6) * t68 + t124;
t105 = -pkin(5) * t118 - pkin(4);
t98 = t209 * t118;
t97 = t209 * t115;
t94 = (pkin(5) * t115 + pkin(9)) * t116;
t86 = t118 * t96;
t84 = t87 ^ 2;
t70 = pkin(5) * t129 + pkin(9) * t169;
t65 = -t115 * t187 + t188;
t62 = qJD(4) * t76 + t141;
t61 = -qJD(4) * t75 + t140;
t58 = -qJ(6) * t179 + t86 + (-pkin(5) - t212) * t119;
t57 = t60 * qJD(3);
t56 = t224 * qJD(3);
t32 = qJD(5) * t132 - t115 * t61 + t118 * t157;
t31 = qJD(5) * t63 + t115 * t157 + t118 * t61;
t27 = pkin(5) * t152 + t37;
t24 = -t148 + t33;
t21 = qJD(4) * t44 + t116 * t56;
t20 = -qJD(4) * t43 + t119 * t56;
t8 = -t102 * t32 + t63 * t149 + t62 * t87 + t69 * t75;
t7 = t102 * t31 + t132 * t149 + t62 * t89 - t68 * t75;
t6 = -qJD(5) * t23 - t115 * t20 + t118 * t57;
t5 = qJD(5) * t22 + t115 * t57 + t118 * t20;
t4 = -qJD(6) * t87 - t133;
t3 = -qJD(6) * t89 + t123 + t143;
t2 = -t102 * t6 + t22 * t149 + t21 * t87 + t43 * t69;
t1 = t102 * t5 - t23 * t149 + t21 * t89 - t43 * t68;
t15 = [0, 0, 0, -t57 * qJD(3), -t56 * qJD(3), 0, 0, 0, 0, 0, -qJD(4) * t21 + (-t119 * t57 - t170 * t224) * qJD(3), -qJD(4) * t20 + (t116 * t57 - t169 * t224) * qJD(3), 0, 0, 0, 0, 0, t2, t1, t2, t1, t22 * t68 - t23 * t69 - t5 * t87 - t6 * t89, t11 * t5 + t12 * t43 + t21 * t24 + t22 * t3 + t23 * t4 + t6 * t9; 0, 0, 0, -t160, -t120 * t182, 0, 0, 0, 0, 0, -t119 * t160 + (-t62 - t141) * qJD(4), t116 * t160 + (-t61 - t140) * qJD(4), 0, 0, 0, 0, 0, t8, t7, t8, t7, t132 * t69 - t31 * t87 - t32 * t89 + t63 * t68, t11 * t31 + t12 * t75 - t132 * t4 + t24 * t62 + t3 * t63 + t32 * t9; 0, 0, 0, t219, (-t223 + t54) * qJD(3), 0.2e1 * t119 * t149, -0.2e1 * t177 * t164, t121 * t119, -t121 * t116, 0, t116 * t135 - t119 * t136, t116 * t136 + t119 * t135, t89 * t151 + (-t118 * t68 - t89 * t168) * t116, (-t118 * t87 - t195) * t169 + (t196 - t118 * t69 + (t115 * t87 - t118 * t89) * qJD(5)) * t116, t102 * t154 + t119 * t68 + (t116 * t89 + t118 * t137) * qJD(4), t102 * t153 + t119 * t69 + (-t115 * t137 - t116 * t87) * qJD(4), (-t102 - t173) * t170, (t96 * t168 + t222) * t102 + ((pkin(9) * t87 + t190) * qJD(4) + (t197 + (pkin(9) * t102 + t34) * t118) * qJD(5) + t147) * t119 + (t33 * t167 + pkin(9) * t69 + t192 - t211 + ((-pkin(9) * t180 + t86) * qJD(3) + t13) * qJD(4)) * t116, t221 * t102 + (t33 * t165 + (qJD(4) * t89 - t155) * pkin(9) - t144) * t119 + (-t33 * t168 - pkin(9) * t68 + t191 - t210 + (-pkin(9) * t185 - t188 * qJD(3) - t14) * qJD(4)) * t116, t69 * t94 + t70 * t87 + (t24 * t171 - t3) * t119 + t206 * t102 + (t24 * t167 + t198 - t211 + (qJD(3) * t58 + t9) * qJD(4)) * t116, -t68 * t94 + t70 * t89 + (t165 * t24 + t4) * t119 - t205 * t102 + (-t24 * t168 + t194 - t210 + (-qJD(3) * t65 - t11) * qJD(4)) * t116, t58 * t68 - t65 * t69 + t206 * t89 + t205 * t87 + (-t11 * t115 - t118 * t9) * t169 + (-t115 * t4 - t118 * t3 + (-t11 * t118 + t115 * t9) * qJD(5)) * t116, t12 * t94 + t3 * t58 + t4 * t65 - t206 * t9 + (-t116 * t54 + t70) * t24 - t205 * t11; 0, 0, 0, 0, 0, -t116 * t122 * t119, t177 * t122, 0, 0, 0, qJD(4) * t37 - t52 * t175 - t17, (-qJD(3) * t52 - t50) * t119, -t89 * t185 - t196, -t217 * t115 + t218 * t118, -t102 * t167 + (t102 * t178 + (-t89 + t171) * t116) * qJD(3), t155 + (-t102 * t180 + (t87 + t165) * t116) * qJD(3), t102 * t175, -pkin(4) * t69 - t191 + t146 * t102 - t37 * t87 + (pkin(10) * t185 + t190) * qJD(5) + (-t13 * t116 + (-pkin(10) * t170 - t119 * t33) * t115) * qJD(3), pkin(4) * t68 + t192 - t204 * t102 - t37 * t89 + (-pkin(10) * t102 * t115 + t118 * t33) * qJD(5) + (-t33 * t178 + (-pkin(10) * t165 + t14) * t116) * qJD(3), t105 * t69 - t194 - t27 * t87 + t208 * t102 + (t24 + t214) * t168 + (-t24 * t180 + (qJD(4) * t97 - t9) * t116) * qJD(3), -t105 * t68 + t198 - t27 * t89 - t207 * t102 + (pkin(5) * t195 + t118 * t24) * qJD(5) + (-t24 * t178 + (qJD(4) * t98 + t11) * t116) * qJD(3), t68 * t97 + t69 * t98 + t208 * t89 + t207 * t87 + (t102 * t9 + t4) * t118 + (-t3 + t201) * t115, t105 * t12 + t3 * t97 - t4 * t98 - t208 * t9 + (pkin(5) * t168 - t27) * t24 - t207 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t87, -t84 + t215, -t68 - t200, -t199 - t69, t149, -t102 * t14 - t33 * t89 + t124, -t102 * t13 + t33 * t87 + t144, 0.2e1 * t143 - t201 + (t148 - t24) * t89 + t123, -pkin(5) * t215 - t10 * t102 + (qJD(6) + t24) * t87 + t133, pkin(5) * t68 + t213 * t87, -t213 * t11 + (-t24 * t89 + t3) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217, t218, -t84 - t215, t11 * t87 + t89 * t9 + t12;];
tauc_reg = t15;
