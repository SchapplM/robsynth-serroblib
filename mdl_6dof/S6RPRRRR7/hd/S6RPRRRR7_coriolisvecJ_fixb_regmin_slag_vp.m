% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x34]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:17:55
% EndTime: 2019-03-09 07:18:04
% DurationCPUTime: 3.15s
% Computational Cost: add. (5757->315), mult. (12680->426), div. (0->0), fcn. (9342->8), ass. (0->190)
t159 = cos(qJ(6));
t213 = qJD(6) * t159;
t157 = sin(qJ(4));
t161 = cos(qJ(4));
t162 = cos(qJ(3));
t220 = qJD(1) * t162;
t158 = sin(qJ(3));
t221 = qJD(1) * t158;
t118 = t157 * t221 - t161 * t220;
t119 = -t157 * t220 - t161 * t221;
t156 = sin(qJ(5));
t160 = cos(qJ(5));
t80 = t118 * t156 + t160 * t119;
t267 = t159 * t80;
t274 = t213 - t267;
t129 = pkin(3) * t221 + qJD(1) * qJ(2);
t94 = -pkin(4) * t119 + t129;
t249 = t80 * t94;
t216 = qJD(5) * t156;
t163 = -pkin(1) - pkin(7);
t132 = t163 * qJD(1) + qJD(2);
t109 = -pkin(8) * t221 + t132 * t158;
t100 = t161 * t109;
t110 = -pkin(8) * t220 + t162 * t132;
t102 = qJD(3) * pkin(3) + t110;
t182 = -t102 * t157 - t100;
t198 = pkin(8) * qJD(1) - t132;
t219 = qJD(3) * t158;
t103 = t198 * t219;
t218 = qJD(3) * t162;
t104 = t198 * t218;
t195 = t161 * t103 + t157 * t104;
t170 = qJD(4) * t182 + t195;
t123 = t157 * t162 + t158 * t161;
t150 = qJD(3) + qJD(4);
t92 = t150 * t123;
t84 = t92 * qJD(1);
t31 = t84 * pkin(9) + t170;
t256 = pkin(9) * t119;
t57 = -t182 + t256;
t201 = t156 * t31 - t57 * t216;
t228 = t161 * t162;
t178 = t157 * t158 - t228;
t261 = qJD(1) * t178;
t166 = t150 * t261;
t217 = qJD(4) * t157;
t196 = t157 * t103 - t109 * t217;
t262 = (qJD(4) * t102 - t104) * t161;
t30 = pkin(9) * t166 + t196 + t262;
t113 = t118 * pkin(9);
t99 = t157 * t109;
t194 = t161 * t102 - t99;
t56 = t113 + t194;
t54 = pkin(4) * t150 + t56;
t3 = (qJD(5) * t54 + t30) * t160 + t201;
t273 = -t249 - t3;
t181 = t160 * t118 - t119 * t156;
t250 = t181 * t94;
t202 = t156 * t30 - t160 * t31;
t237 = t160 * t57;
t28 = t156 * t54 + t237;
t4 = qJD(5) * t28 + t202;
t272 = t250 - t4;
t148 = qJD(5) + t150;
t180 = t160 * t123 - t156 * t178;
t93 = t150 * t228 - t157 * t219 - t158 * t217;
t45 = qJD(5) * t180 + t156 * t93 + t160 * t92;
t271 = t148 * t45;
t224 = -qJD(6) + t80;
t270 = qJD(6) + t224;
t155 = sin(qJ(6));
t214 = qJD(6) * t155;
t215 = qJD(5) * t160;
t39 = t118 * t216 + t119 * t215 + t156 * t166 - t160 * t84;
t18 = t148 * t213 + t159 * t39 + t181 * t214;
t66 = t148 * t155 - t159 * t181;
t88 = t123 * t156 + t160 * t178;
t269 = -t88 * t18 - t45 * t66;
t16 = t18 * t155;
t7 = t274 * t66 + t16;
t171 = qJD(5) * t181 + t156 * t84 + t160 * t166;
t36 = t155 * t171;
t69 = t224 * t213;
t244 = -t36 - t69;
t6 = t181 * t66 + t224 * t267 + t244;
t266 = t155 * t224;
t38 = t159 * t171;
t64 = -t159 * t148 - t155 * t181;
t5 = -t181 * t64 - t224 * t266 - t38;
t19 = qJD(6) * t66 + t155 * t39;
t1 = -t155 * t19 + t18 * t159 + t266 * t66 - t274 * t64;
t240 = t156 * t57;
t27 = t160 * t54 - t240;
t25 = -pkin(5) * t148 - t27;
t255 = t25 * t80;
t248 = t181 * t80;
t251 = t224 * t181;
t22 = t181 ^ 2 - t80 ^ 2;
t20 = -t148 * t80 + t39;
t26 = pkin(10) * t148 + t28;
t41 = -pkin(5) * t80 + pkin(10) * t181 + t94;
t13 = t155 * t41 + t159 * t26;
t186 = -t13 * t181 + t4 * t155 + t25 * t213;
t183 = t155 * t26 - t159 * t41;
t207 = -t183 * t181 + t25 * t214;
t21 = -t148 * t181 + t171;
t53 = -pkin(5) * t181 - pkin(10) * t80;
t167 = -qJD(5) * t88 - t156 * t92 + t160 * t93;
t263 = t167 * t148;
t247 = pkin(8) - t163;
t127 = t247 * t158;
t128 = t247 * t162;
t232 = t128 * t161;
t70 = pkin(9) * t178 + t127 * t157 - t232;
t179 = t127 * t161 + t128 * t157;
t71 = -pkin(9) * t123 - t179;
t51 = t156 * t70 + t160 * t71;
t140 = t158 * pkin(3) + qJ(2);
t105 = pkin(4) * t123 + t140;
t52 = pkin(5) * t180 + pkin(10) * t88 + t105;
t121 = t247 * t219;
t122 = qJD(3) * t128;
t175 = -qJD(4) * t232 + t157 * t121 - t161 * t122 + t127 * t217;
t48 = -pkin(9) * t93 + t175;
t168 = qJD(4) * t179 + t161 * t121 + t157 * t122;
t49 = t92 * pkin(9) + t168;
t50 = t156 * t71 - t160 * t70;
t8 = -qJD(5) * t50 + t156 * t49 + t160 * t48;
t260 = -t180 * (qJD(6) * t41 + t3) + (qJD(6) * t52 + t8) * t224 - t25 * t45 - t4 * t88 + t51 * t171;
t151 = qJD(1) * qJD(2);
t259 = 0.2e1 * t151;
t257 = pkin(4) * t118;
t254 = t25 * t88;
t253 = t171 * t88;
t252 = t52 * t171;
t143 = pkin(3) * t161 + pkin(4);
t230 = t156 * t157;
t189 = -t110 * t157 - t100;
t58 = t189 - t256;
t235 = t161 * t110 - t99;
t59 = t113 + t235;
t245 = t156 * t58 + t160 * t59 - t143 * t215 - (-t157 * t216 + (t160 * t161 - t230) * qJD(4)) * pkin(3);
t229 = t157 * t160;
t243 = -t156 * t59 + t160 * t58 + t143 * t216 + (t157 * t215 + (t156 * t161 + t229) * qJD(4)) * pkin(3);
t236 = t93 * t150;
t233 = t118 * t119;
t231 = t129 * t118;
t164 = qJD(3) ^ 2;
t227 = t164 * t158;
t226 = t164 * t162;
t165 = qJD(1) ^ 2;
t225 = t165 * qJ(2);
t212 = qJD(1) * qJD(3);
t206 = t162 * t212;
t126 = pkin(3) * t206 + t151;
t223 = t158 ^ 2 - t162 ^ 2;
t222 = -t164 - t165;
t133 = pkin(3) * t218 + qJD(2);
t210 = 0.2e1 * qJD(1);
t146 = pkin(3) * t220;
t209 = t88 * t214;
t208 = -pkin(4) * t148 - t54;
t205 = -pkin(3) * t150 - t102;
t112 = pkin(3) * t229 + t143 * t156 + pkin(10);
t47 = -t257 + t53;
t191 = qJD(6) * t112 + t146 + t47;
t141 = pkin(4) * t156 + pkin(10);
t190 = qJD(6) * t141 + t47;
t188 = qJD(6) * t180 + qJD(1);
t32 = t156 * t56 + t237;
t185 = pkin(4) * t216 - t32;
t73 = pkin(4) * t93 + t133;
t184 = t224 * t45 + t253;
t177 = -t129 * t119 - t196;
t174 = t112 * t171 - t224 * t245 - t255;
t33 = t160 * t56 - t240;
t169 = t141 * t171 - t255 - (-pkin(4) * t215 + t33) * t224;
t67 = -pkin(4) * t166 + t126;
t142 = -pkin(4) * t160 - pkin(5);
t111 = pkin(3) * t230 - t143 * t160 - pkin(5);
t97 = t146 - t257;
t85 = t92 * t150;
t68 = t118 ^ 2 - t119 ^ 2;
t61 = -t118 * t150 + t166;
t60 = -t119 * t150 - t84;
t14 = pkin(5) * t167 + pkin(10) * t45 + t73;
t11 = -pkin(5) * t171 - t39 * pkin(10) + t67;
t10 = t159 * t11;
t9 = qJD(5) * t51 + t156 * t48 - t160 * t49;
t2 = [0, 0, 0, 0, t259, qJ(2) * t259, -0.2e1 * t158 * t206, 0.2e1 * t223 * t212, -t227, -t226, 0, -t163 * t227 + (qJ(2) * t218 + qJD(2) * t158) * t210, -t163 * t226 + (-qJ(2) * t219 + qJD(2) * t162) * t210, t118 * t92 + t178 * t84, t118 * t93 - t92 * t119 + t84 * t123 - t166 * t178, -t85, -t236, 0, -t133 * t119 + t126 * t123 + t129 * t93 + (-t140 * t261 + t168) * t150, -t133 * t118 - t126 * t178 - t129 * t92 - t140 * t84 - t150 * t175, t181 * t45 - t39 * t88, t167 * t181 - t180 * t39 - t45 * t80 - t253, -t271, -t263, 0, -t105 * t171 - t148 * t9 + t167 * t94 + t180 * t67 - t73 * t80, t105 * t39 - t148 * t8 - t181 * t73 - t45 * t94 - t67 * t88, t269 * t159 + t66 * t209 -(-t155 * t66 - t159 * t64) * t45 - (-t16 - t159 * t19 + (t155 * t64 - t159 * t66) * qJD(6)) * t88, t159 * t184 + t167 * t66 + t18 * t180 - t209 * t224, -t155 * t184 - t167 * t64 - t180 * t19 - t69 * t88, -t167 * t224 - t171 * t180, t10 * t180 - t183 * t167 + t50 * t19 + t9 * t64 + (-t14 * t224 - t252 + (-t180 * t26 + t224 * t51 - t254) * qJD(6)) * t159 + t260 * t155, -t13 * t167 + t50 * t18 + t9 * t66 + ((-qJD(6) * t51 + t14) * t224 + t252 - (-qJD(6) * t26 + t11) * t180 + qJD(6) * t254) * t155 + t260 * t159; 0, 0, 0, 0, -t165, -t225, 0, 0, 0, 0, 0, t222 * t158, t222 * t162, 0, 0, 0, 0, 0, qJD(1) * t119 - t85, qJD(1) * t118 - t236, 0, 0, 0, 0, 0, qJD(1) * t80 - t271, qJD(1) * t181 - t263, 0, 0, 0, 0, 0, t180 * t36 + t88 * t19 + t45 * t64 - (-t155 * t167 - t159 * t188) * t224, t180 * t38 - (t155 * t188 - t159 * t167) * t224 - t269; 0, 0, 0, 0, 0, 0, t162 * t165 * t158, -t223 * t165, 0, 0, 0, -t162 * t225, t158 * t225, t233, t68, t60, t61, 0, t119 * t146 + t231 - t189 * t150 + (t157 * t205 - t100) * qJD(4) + t195, t118 * t146 + t235 * t150 + (qJD(4) * t205 + t104) * t161 + t177, t248, t22, t20, t21, 0, -t243 * t148 + t80 * t97 + t272, t245 * t148 + t181 * t97 + t273, t7, t1, t6, t5, -t251, t111 * t19 + t243 * t64 + (t191 * t224 - t4) * t159 + t174 * t155 + t207, t111 * t18 + t174 * t159 - t191 * t266 + t243 * t66 + t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, t68, t60, t61, 0, -t150 * t182 + t170 + t231, t150 * t194 + t177 - t262, t248, t22, t20, t21, 0, -t80 * t257 + t148 * t32 + t250 + (t156 * t208 - t237) * qJD(5) - t202, -t181 * t257 + t148 * t33 - t249 + (qJD(5) * t208 - t30) * t160 - t201, t7, t1, t6, t5, -t251, t142 * t19 + t185 * t64 + (t190 * t224 - t4) * t159 + t169 * t155 + t207, t142 * t18 + t159 * t169 + t185 * t66 - t190 * t266 + t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, t22, t20, t21, 0, t28 * t148 + t272, t27 * t148 + t273, t7, t1, t6, t5, -t251, -pkin(5) * t19 - t4 * t159 + (-t155 * t27 + t159 * t53) * t224 - t28 * t64 - t155 * t255 - t244 * pkin(10) + t207, -pkin(5) * t18 - (t155 * t53 + t159 * t27) * t224 - t28 * t66 - t25 * t267 + (-t214 * t224 + t38) * pkin(10) + t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t64, -t64 ^ 2 + t66 ^ 2, -t224 * t64 + t18, -t224 * t66 - t19, -t171, -t270 * t13 - t155 * t3 - t25 * t66 + t10, -t155 * t11 - t159 * t3 + t270 * t183 + t25 * t64;];
tauc_reg  = t2;
