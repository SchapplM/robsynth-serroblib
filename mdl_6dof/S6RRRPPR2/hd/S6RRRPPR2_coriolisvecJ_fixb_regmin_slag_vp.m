% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:43
% EndTime: 2019-03-09 15:26:50
% DurationCPUTime: 3.02s
% Computational Cost: add. (5882->344), mult. (15382->446), div. (0->0), fcn. (11421->8), ass. (0->208)
t168 = cos(qJ(3));
t169 = cos(qJ(2));
t222 = qJD(1) * t169;
t212 = t168 * t222;
t165 = sin(qJ(3));
t166 = sin(qJ(2));
t223 = qJD(1) * t166;
t213 = t165 * t223;
t119 = -t212 + t213;
t121 = -t165 * t222 - t168 * t223;
t162 = sin(pkin(10));
t163 = cos(pkin(10));
t188 = -t119 * t162 - t163 * t121;
t260 = qJD(6) + t188;
t270 = qJD(6) - t260;
t164 = sin(qJ(6));
t266 = t164 * t260;
t167 = cos(qJ(6));
t158 = qJD(2) + qJD(3);
t217 = qJD(1) * qJD(2);
t209 = t169 * t217;
t97 = qJD(3) * t212 - t158 * t213 + t168 * t209;
t133 = t165 * t169 + t166 * t168;
t103 = t158 * t133;
t98 = t103 * qJD(1);
t65 = -t162 * t98 + t163 * t97;
t63 = t167 * t65;
t269 = -t260 * t266 + t63;
t93 = t119 * t163 - t121 * t162;
t256 = pkin(5) * t93;
t268 = t158 * t93;
t78 = t158 * t164 - t167 * t93;
t267 = t260 * t78;
t147 = t162 * t165 * pkin(2);
t220 = qJD(3) * t168;
t258 = pkin(7) + pkin(8);
t140 = t258 * t169;
t136 = qJD(1) * t140;
t126 = t168 * t136;
t139 = t258 * t166;
t134 = qJD(1) * t139;
t198 = t134 * t165 - t126;
t235 = qJ(4) * t119;
t81 = t198 + t235;
t116 = t121 * qJ(4);
t122 = t165 * t136;
t225 = -t168 * t134 - t122;
t82 = t116 + t225;
t236 = -qJD(3) * t147 - t162 * t81 + (pkin(2) * t220 - t82) * t163;
t265 = t188 ^ 2;
t264 = -0.2e1 * t217;
t257 = pkin(5) * t188;
t239 = qJD(5) + t236;
t232 = t163 * t165;
t237 = -t162 * t82 + t163 * t81 + (t162 * t168 + t232) * qJD(3) * pkin(2);
t248 = qJD(2) * pkin(2);
t128 = -t134 + t248;
t214 = qJD(2) * t258;
t196 = qJD(1) * t214;
t129 = t166 * t196;
t263 = (qJD(3) * t128 - t129) * t168;
t187 = -t128 * t165 - t126;
t77 = -t187 - t235;
t72 = t162 * t77;
t201 = t168 * t128 - t122;
t76 = t116 + t201;
t49 = t163 * t76 - t72;
t227 = qJD(5) - t49;
t261 = qJD(1) * t133;
t130 = t169 * t196;
t221 = qJD(3) * t165;
t199 = -t165 * t130 - t136 * t221;
t39 = -qJ(4) * t98 - qJD(4) * t119 + t199 + t263;
t200 = t165 * t129 - t168 * t130;
t175 = qJD(3) * t187 + t200;
t40 = -qJ(4) * t97 + qJD(4) * t121 + t175;
t17 = t162 * t39 - t163 * t40;
t154 = -pkin(2) * t169 - pkin(1);
t138 = t154 * qJD(1);
t104 = pkin(3) * t119 + qJD(4) + t138;
t173 = -qJ(5) * t188 + t104;
t51 = pkin(4) * t93 + t173;
t184 = t188 * t51 + t17;
t259 = pkin(4) + pkin(9);
t255 = pkin(3) * t121;
t89 = -qJ(4) * t133 - t139 * t168 - t140 * t165;
t132 = t165 * t166 - t168 * t169;
t186 = t139 * t165 - t140 * t168;
t90 = -qJ(4) * t132 - t186;
t59 = t162 * t90 - t163 * t89;
t254 = t17 * t59;
t100 = t163 * t132 + t133 * t162;
t101 = -t132 * t162 + t133 * t163;
t185 = pkin(3) * t132 + t154;
t178 = -qJ(5) * t101 + t185;
t41 = t259 * t100 + t178;
t253 = t41 * t65;
t246 = t163 * t77;
t48 = t162 * t76 + t246;
t252 = t48 * t188;
t251 = t78 * t93;
t80 = t158 * t167 + t164 * t93;
t250 = t80 * t93;
t249 = t260 * t93;
t18 = t162 * t40 + t163 * t39;
t69 = pkin(3) * t158 + t76;
t47 = t162 * t69 + t246;
t45 = -qJ(5) * t158 - t47;
t26 = -t45 - t256;
t247 = t100 * t26;
t245 = t164 * t65;
t244 = t164 * t80;
t218 = qJD(6) * t167;
t219 = qJD(6) * t164;
t64 = t162 * t97 + t163 * t98;
t34 = -t158 * t219 + t164 * t64 + t93 * t218;
t242 = t167 * t34;
t241 = t167 * t260;
t240 = t257 + t239;
t238 = t256 + t237;
t234 = t121 * t119;
t233 = t138 * t121;
t171 = qJD(1) ^ 2;
t231 = t169 * t171;
t170 = qJD(2) ^ 2;
t230 = t170 * t166;
t229 = t170 * t169;
t228 = t257 + t227;
t224 = t166 ^ 2 - t169 ^ 2;
t216 = t158 * qJD(5) + t18;
t156 = t166 * t248;
t155 = pkin(2) * t223;
t4 = -pkin(5) * t64 + t216;
t46 = t163 * t69 - t72;
t195 = qJD(5) - t46;
t24 = -t259 * t158 + t195 + t257;
t29 = t259 * t93 + t173;
t9 = t164 * t24 + t167 * t29;
t215 = t4 * t167 - t9 * t93;
t151 = -pkin(3) * t163 - pkin(4);
t211 = pkin(3) * t98 + qJD(2) * t155;
t210 = t237 * t188;
t208 = -pkin(2) * t158 - t128;
t207 = pkin(3) * t103 + t156;
t135 = t166 * t214;
t137 = t169 * t214;
t181 = -t168 * t135 - t165 * t137 - t139 * t220 - t140 * t221;
t52 = -qJ(4) * t103 - qJD(4) * t132 + t181;
t102 = t158 * t132;
t174 = qJD(3) * t186 + t135 * t165 - t168 * t137;
t53 = qJ(4) * t102 - qJD(4) * t133 + t174;
t20 = t162 * t52 - t163 * t53;
t206 = t260 * t26;
t205 = t260 ^ 2;
t203 = pkin(1) * t264;
t153 = pkin(2) * t168 + pkin(3);
t114 = t153 * t163 - t147;
t109 = -pkin(4) - t114;
t106 = -pkin(9) + t109;
t57 = pkin(4) * t188 + qJ(5) * t93 - t255;
t56 = t155 + t57;
t87 = t188 * pkin(9);
t202 = -qJD(6) * t106 + t56 + t87;
t44 = -pkin(4) * t158 + t195;
t194 = -t188 * t45 + t44 * t93;
t193 = t188 * t47 - t46 * t93;
t192 = -t93 ^ 2 - t265;
t190 = t164 * t29 - t167 * t24;
t191 = t4 * t164 - t190 * t93 + (t167 * t188 + t218) * t26;
t21 = t162 * t53 + t163 * t52;
t60 = t162 * t89 + t163 * t90;
t189 = -t51 * t93 + t216;
t115 = pkin(2) * t232 + t153 * t162;
t183 = -t219 * t260 + t63;
t182 = t138 * t119 - t199;
t42 = pkin(5) * t101 + t59;
t66 = -t102 * t162 + t163 * t103;
t180 = t4 * t100 + t26 * t66 - t42 * t65;
t179 = -t241 * t260 - t245;
t177 = -qJ(5) * t65 - qJD(5) * t188 + t211;
t67 = -t102 * t163 - t103 * t162;
t176 = -qJ(5) * t67 - qJD(5) * t101 + t207;
t19 = pkin(4) * t64 + t177;
t172 = t101 * t17 + t188 * t20 - t21 * t93 + t59 * t65 - t60 * t64;
t150 = pkin(3) * t162 + qJ(5);
t148 = -pkin(9) + t151;
t108 = qJ(5) + t115;
t83 = -t119 ^ 2 + t121 ^ 2;
t71 = (-t121 - t261) * t158;
t70 = t119 * t158 + t97;
t62 = t167 * t64;
t58 = pkin(4) * t100 + t178;
t43 = -pkin(5) * t100 + t60;
t35 = qJD(6) * t80 - t62;
t33 = t57 + t87;
t27 = t48 - t256;
t22 = pkin(4) * t66 + t176;
t15 = t259 * t66 + t176;
t14 = -t266 * t80 + t242;
t13 = -pkin(5) * t66 + t21;
t12 = pkin(5) * t67 + t20;
t11 = t179 - t251;
t10 = t250 + t269;
t7 = t259 * t64 + t177;
t6 = pkin(5) * t65 + t17;
t5 = t167 * t6;
t1 = (-t260 * t80 - t35) * t167 + (-t34 + t267) * t164;
t2 = [0, 0, 0, 0.2e1 * t166 * t209, t224 * t264, t229, -t230, 0, -pkin(7) * t229 + t166 * t203, pkin(7) * t230 + t169 * t203, t102 * t121 + t133 * t97, t102 * t119 + t103 * t121 - t132 * t97 - t133 * t98, -t102 * t158, -t103 * t158, 0, t154 * t98 + t138 * t103 + t174 * t158 + (qJD(1) * t132 + t119) * t156, t154 * t97 - t138 * t102 - t181 * t158 + (-t121 + t261) * t156, -t100 * t18 - t46 * t67 - t47 * t66 + t172, t104 * t207 + t18 * t60 + t185 * t211 - t46 * t20 + t47 * t21 + t254, -t100 * t216 + t44 * t67 + t45 * t66 + t172, -t100 * t19 + t158 * t20 - t22 * t93 - t51 * t66 - t58 * t64, -t101 * t19 + t158 * t21 - t188 * t22 - t51 * t67 - t58 * t65, t19 * t58 + t20 * t44 - t21 * t45 + t216 * t60 + t22 * t51 + t254, t66 * t244 + (t164 * t34 + t218 * t80) * t100 (-t164 * t78 + t167 * t80) * t66 + (-t164 * t35 + t242 + (-t167 * t78 - t244) * qJD(6)) * t100, t66 * t266 + t101 * t34 + t67 * t80 + (t218 * t260 + t245) * t100, t100 * t183 - t101 * t35 + t241 * t66 - t67 * t78, t101 * t65 + t260 * t67, t5 * t101 + t13 * t78 + t43 * t35 - t190 * t67 + (-t7 * t101 - t15 * t260 - t253) * t164 + (t12 * t260 - t180) * t167 + ((-t164 * t42 - t167 * t41) * t260 - t9 * t101 + t164 * t247) * qJD(6), t13 * t80 + t43 * t34 - t9 * t67 + (-(qJD(6) * t42 + t15) * t260 - t253 - (qJD(6) * t24 + t7) * t101 + qJD(6) * t247) * t167 + (-(-qJD(6) * t41 + t12) * t260 - (-qJD(6) * t29 + t6) * t101 + t180) * t164; 0, 0, 0, -t166 * t231, t224 * t171, 0, 0, 0, t171 * pkin(1) * t166, pkin(1) * t231, -t234, t83, t70, t71, 0, -t119 * t155 + t233 - t198 * t158 + (t165 * t208 - t126) * qJD(3) + t200, t121 * t155 + t225 * t158 + (qJD(3) * t208 + t129) * t168 + t182, -t114 * t65 - t115 * t64 - t236 * t93 + t193 + t210, t18 * t115 - t17 * t114 - t104 * (t155 - t255) + t236 * t47 - t237 * t46, -t108 * t64 + t109 * t65 - t239 * t93 + t194 + t210, t158 * t237 + t56 * t93 + t184, t158 * t239 + t188 * t56 + t189, t108 * t216 + t109 * t17 + t237 * t44 - t239 * t45 - t51 * t56, t14, t1, t10, t11, t249, t106 * t63 + t108 * t35 + t240 * t78 + (t164 * t202 + t167 * t238) * t260 + t191, t108 * t34 + t240 * t80 + t202 * t241 + (-t106 * t65 - t238 * t260 - t206) * t164 + t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, t83, t70, t71, 0, -t158 * t187 + t175 + t233, t158 * t201 + t182 - t263, -t252 + t49 * t93 + (-t162 * t64 - t163 * t65) * pkin(3) + t193, t46 * t48 - t47 * t49 + (t104 * t121 + t162 * t18 - t163 * t17) * pkin(3), -t150 * t64 + t151 * t65 - t227 * t93 + t194 - t252, -t158 * t48 + t57 * t93 + t184, t158 * t227 + t188 * t57 + t189, t150 * t216 + t151 * t17 - t227 * t45 - t44 * t48 - t51 * t57, t14, t1, t10, t11, t249, t150 * t35 - (-t164 * t33 + t167 * t27) * t260 + t228 * t78 + t183 * t148 + t191, t150 * t34 + t228 * t80 + (-qJD(6) * t148 + t33) * t241 + (-t148 * t65 + t260 * t27 - t206) * t164 + t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, t188 * t46 + t47 * t93 + t211, t192, -t158 * t188 - t64, -t65 + t268, -t188 * t44 - t45 * t93 + t19, 0, 0, 0, 0, 0, t179 + t251, t250 - t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 + t268, -t188 * t93, -t158 ^ 2 - t265, t158 * t45 + t184, 0, 0, 0, 0, 0, -t158 * t78 - t164 * t205 + t63, -t158 * t80 - t167 * t205 - t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t78, -t78 ^ 2 + t80 ^ 2, t34 + t267, -t270 * t80 + t62, t65, -t164 * t7 - t26 * t80 - t270 * t9 + t5, -t164 * t6 - t167 * t7 + t190 * t270 + t26 * t78;];
tauc_reg  = t2;
