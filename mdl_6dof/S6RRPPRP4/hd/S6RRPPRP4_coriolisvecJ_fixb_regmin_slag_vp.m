% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:41
% EndTime: 2019-03-09 08:39:52
% DurationCPUTime: 3.51s
% Computational Cost: add. (3901->415), mult. (9530->540), div. (0->0), fcn. (6286->6), ass. (0->200)
t165 = sin(pkin(9));
t171 = cos(qJ(2));
t216 = qJD(1) * qJD(2);
t205 = t171 * t216;
t140 = t165 * t205;
t168 = sin(qJ(5));
t170 = cos(qJ(5));
t166 = cos(pkin(9));
t198 = t166 * t205;
t169 = sin(qJ(2));
t227 = qJD(1) * t169;
t208 = t165 * t227;
t219 = t166 * qJD(2);
t117 = -t208 + t219;
t207 = t166 * t227;
t221 = t165 * qJD(2);
t118 = t207 + t221;
t64 = t168 * t117 - t170 * t118;
t176 = qJD(5) * t64 + t170 * t140 - t168 * t198;
t218 = t171 * qJD(1);
t152 = qJD(5) + t218;
t249 = t64 * t152;
t273 = t176 - t249;
t262 = t64 ^ 2;
t120 = t168 * t165 + t170 * t166;
t104 = t120 * qJD(5);
t181 = t120 * t171;
t248 = qJD(1) * t181 + t104;
t210 = t165 * t218;
t222 = qJD(5) * t170;
t223 = qJD(5) * t168;
t238 = t168 * t166;
t247 = t165 * t222 - t166 * t223 + t170 * t210 - t218 * t238;
t121 = t170 * t165 - t238;
t97 = t121 * t169;
t188 = t170 * t117 + t168 * t118;
t272 = t188 ^ 2;
t271 = -0.2e1 * t216;
t204 = t169 * t216;
t150 = pkin(5) * t204;
t239 = t166 * t171;
t215 = pkin(8) * t239;
t261 = -pkin(3) - pkin(4);
t155 = pkin(7) * t227;
t128 = (qJD(3) - t155) * qJD(2);
t193 = pkin(2) * t169 - qJ(3) * t171;
t102 = qJD(2) * t193 - t169 * qJD(3);
t93 = t102 * qJD(1);
t55 = -t165 * t128 + t166 * t93;
t26 = (t261 * t169 - t215) * t216 - t55;
t56 = t166 * t128 + t165 * t93;
t213 = qJ(4) * t204 + t56;
t27 = (pkin(8) * t221 - qJD(4)) * t218 + t213;
t237 = t169 * qJ(3);
t131 = -t171 * pkin(2) - pkin(1) - t237;
t110 = t131 * qJD(1);
t156 = pkin(7) * t218;
t138 = qJD(2) * qJ(3) + t156;
t71 = t166 * t110 - t165 * t138;
t58 = pkin(3) * t218 + qJD(4) - t71;
t28 = pkin(4) * t218 - t118 * pkin(8) + t58;
t206 = qJ(4) * t218;
t72 = t165 * t110 + t166 * t138;
t60 = -t206 + t72;
t35 = -t117 * pkin(8) + t60;
t203 = t168 * t27 - t170 * t26 + t35 * t222 + t28 * t223;
t2 = t150 + t203;
t11 = t168 * t28 + t170 * t35;
t8 = t152 * qJ(6) + t11;
t270 = t8 * t152 - t2;
t212 = -pkin(7) * t165 - pkin(3);
t178 = -t215 + (-pkin(4) + t212) * t169;
t124 = t193 * qJD(1);
t241 = t166 * t124;
t46 = qJD(1) * t178 - t241;
t106 = t165 * t124;
t153 = qJ(4) * t227;
t240 = t166 * t169;
t243 = t165 * t171;
t182 = -pkin(7) * t240 + pkin(8) * t243;
t59 = qJD(1) * t182 + t106 + t153;
t257 = -pkin(8) + qJ(3);
t135 = t257 * t165;
t136 = t257 * t166;
t81 = t168 * t135 + t170 * t136;
t269 = qJD(3) * t121 - qJD(5) * t81 + t168 * t59 - t170 * t46;
t250 = t188 * t152;
t29 = t117 * t222 + t118 * t223 - t168 * t140 - t170 * t198;
t268 = t29 - t250;
t187 = t170 * t135 - t168 * t136;
t267 = -qJD(3) * t120 - qJD(5) * t187 + t168 * t46 + t170 * t59;
t220 = t165 * qJD(4);
t230 = t166 * t206 - t156;
t200 = -t261 * t210 + t220 - t230;
t115 = t118 ^ 2;
t266 = -t117 ^ 2 - t115;
t217 = t171 * qJD(4);
t226 = qJD(2) * t169;
t265 = qJ(4) * t226 - t217;
t235 = t170 * t171;
t264 = qJD(1) * (-t152 * t235 + t168 * t226) + t118 * t64 - t152 * t222;
t148 = pkin(7) * t243;
t162 = t171 * pkin(3);
t61 = t171 * pkin(4) + t148 + t162 + (-pkin(8) * t169 - t131) * t166;
t244 = t165 * t169;
t149 = pkin(7) * t239;
t90 = t165 * t131 + t149;
t83 = -t171 * qJ(4) + t90;
t70 = pkin(8) * t244 + t83;
t189 = t168 * t61 + t170 * t70;
t242 = t166 * t102;
t39 = qJD(2) * t178 - t242;
t94 = t165 * t102;
t40 = qJD(2) * t182 + t265 + t94;
t263 = -qJD(5) * t189 - t168 * t40 + t170 * t39;
t167 = qJD(2) * pkin(2);
t130 = qJD(3) + t155 - t167;
t52 = -t117 * pkin(3) - t118 * qJ(4) + t130;
t33 = t117 * pkin(4) - t52;
t9 = pkin(5) * t188 + qJ(6) * t64 + t33;
t260 = t9 * t64;
t259 = t64 * t188;
t256 = -qJ(6) * t227 + t267;
t255 = pkin(5) * t227 + t269;
t254 = -pkin(5) * t247 - qJ(6) * t248 + t121 * qJD(6) - t200;
t252 = t169 * t58;
t251 = t169 * t60;
t246 = qJD(2) * t187;
t245 = qJD(2) * t81;
t173 = qJD(1) ^ 2;
t234 = t171 * t173;
t172 = qJD(2) ^ 2;
t233 = t172 * t169;
t232 = t172 * t171;
t10 = -t168 * t35 + t170 * t28;
t231 = qJD(6) - t10;
t209 = t171 * t219;
t229 = -qJ(4) * t209 - qJD(4) * t240;
t164 = t171 ^ 2;
t228 = t169 ^ 2 - t164;
t225 = qJD(2) * t171;
t224 = qJD(3) * t118;
t129 = -t166 * pkin(3) - t165 * qJ(4) - pkin(2);
t214 = pkin(7) * t226;
t211 = pkin(3) * t165 + pkin(7);
t202 = -t130 - t167;
t201 = pkin(1) * t271;
t88 = pkin(3) * t210 - t230;
t199 = t88 + t220;
t89 = t166 * t131 - t148;
t107 = t166 * pkin(4) - t129;
t151 = pkin(7) * t205;
t47 = pkin(3) * t140 - qJ(4) * t198 - t118 * qJD(4) + t151;
t197 = qJ(6) * t204;
t196 = t261 * t165 - pkin(7);
t195 = t212 * t169;
t194 = -t118 * t218 + t140;
t190 = -t168 * t70 + t170 * t61;
t78 = -t166 * t214 + t94;
t86 = -pkin(7) * t207 + t106;
t186 = t11 * t152 - t203;
t184 = -t168 * t26 - t170 * t27 - t28 * t222 + t35 * t223;
t183 = t168 * t39 + t170 * t40 + t61 * t222 - t70 * t223;
t146 = qJ(4) * t240;
t82 = t169 * t196 + t146;
t177 = t29 + t250;
t57 = t196 * t225 - t229;
t36 = -pkin(4) * t140 - t47;
t175 = -t152 * t223 - t118 * t188 + (-t152 * t168 * t171 - t170 * t226) * qJD(1);
t174 = t176 + t249;
t3 = -pkin(5) * t176 + t29 * qJ(6) + qJD(6) * t64 + t36;
t139 = qJD(3) * t210;
t100 = t117 * t218;
t99 = t166 * qJD(3) * t117;
t98 = t120 * t169;
t95 = t211 * t169 - t146;
t85 = pkin(7) * t208 + t241;
t84 = t162 - t89;
t79 = -t100 + t198;
t77 = t165 * t214 + t242;
t76 = qJD(1) * t195 - t241;
t75 = t153 + t86;
t74 = t211 * t225 + t229;
t63 = qJD(2) * t195 - t242;
t54 = t78 + t265;
t51 = qJD(2) * t181 + qJD(5) * t97;
t50 = t104 * t169 + t168 * t209 - t221 * t235;
t45 = t120 * pkin(5) - t121 * qJ(6) + t107;
t42 = -pkin(3) * t204 - t55;
t34 = -qJD(1) * t217 + t213;
t21 = -t97 * pkin(5) - t98 * qJ(6) + t82;
t18 = -pkin(5) * t64 + qJ(6) * t188;
t16 = -t171 * pkin(5) - t190;
t15 = t171 * qJ(6) + t189;
t7 = -t152 * pkin(5) + t231;
t6 = t50 * pkin(5) - t51 * qJ(6) - t98 * qJD(6) + t57;
t5 = pkin(5) * t226 - t263;
t4 = -qJ(6) * t226 + t171 * qJD(6) + t183;
t1 = t152 * qJD(6) - t184 - t197;
t12 = [0, 0, 0, 0.2e1 * t171 * t204, t228 * t271, t232, -t233, 0, -pkin(7) * t232 + t169 * t201, pkin(7) * t233 + t171 * t201 (-qJD(1) * t77 - t55) * t171 + ((-pkin(7) * t117 + t130 * t165) * t171 + (t71 + (t89 + 0.2e1 * t148) * qJD(1)) * t169) * qJD(2) (qJD(1) * t78 + t56) * t171 + ((pkin(7) * t118 + t130 * t166) * t171 + (-t72 + (-t90 + 0.2e1 * t149) * qJD(1)) * t169) * qJD(2), t78 * t117 - t77 * t118 + (-t165 * t56 - t166 * t55) * t169 + (-t165 * t72 - t166 * t71 + (-t165 * t90 - t166 * t89) * qJD(1)) * t225, t55 * t89 + t56 * t90 + t71 * t77 + t72 * t78 + (t130 + t155) * pkin(7) * t225, t47 * t244 - t74 * t117 + (qJD(1) * t63 + t42) * t171 + (t52 * t243 - t252 + (-t169 * t84 + t95 * t243) * qJD(1)) * qJD(2), t54 * t117 + t63 * t118 + (-t165 * t34 + t166 * t42) * t169 + (-t165 * t60 + t166 * t58 + (-t165 * t83 + t166 * t84) * qJD(1)) * t225, -t47 * t240 - t74 * t118 + (-qJD(1) * t54 - t34) * t171 + (-t52 * t239 + t251 + (t169 * t83 - t95 * t239) * qJD(1)) * qJD(2), t34 * t83 + t42 * t84 + t47 * t95 + t52 * t74 + t60 * t54 + t58 * t63, -t29 * t98 - t51 * t64, t176 * t98 - t188 * t51 - t29 * t97 + t50 * t64, t51 * t152 - t29 * t171 + (-qJD(1) * t98 + t64) * t226, -t50 * t152 + t176 * t171 + (-qJD(1) * t97 + t188) * t226 (-t152 - t218) * t226, t263 * t152 - t203 * t171 + t57 * t188 - t82 * t176 - t36 * t97 + t33 * t50 + (-qJD(1) * t190 - t10) * t226, -t183 * t152 + t184 * t171 - t57 * t64 - t82 * t29 + t36 * t98 + t33 * t51 + (qJD(1) * t189 + t11) * t226, -t5 * t152 - t2 * t171 - t21 * t176 - t3 * t97 + t9 * t50 + t6 * t188 + (qJD(1) * t16 + t7) * t226, t1 * t97 + t15 * t176 - t16 * t29 - t188 * t4 + t2 * t98 - t5 * t64 - t8 * t50 + t7 * t51, t1 * t171 + t4 * t152 + t21 * t29 - t3 * t98 - t9 * t51 + t6 * t64 + (-qJD(1) * t15 - t8) * t226, t1 * t15 + t2 * t16 + t3 * t21 + t8 * t4 + t7 * t5 + t9 * t6; 0, 0, 0, -t169 * t234, t228 * t173, 0, 0, 0, t173 * pkin(1) * t169, pkin(1) * t234, t139 + ((-qJ(3) * t221 - t71) * t169 + (t85 + t202 * t165 + (t117 - t219) * pkin(7)) * t171) * qJD(1) ((-qJ(3) * t219 + t72) * t169 + (-t86 + (-t118 + t221) * pkin(7) + (qJD(3) + t202) * t166) * t171) * qJD(1), -t86 * t117 + t85 * t118 + t99 + (t71 * t218 + t56) * t166 + (t72 * t218 + t224 - t55) * t165, -t71 * t85 - t72 * t86 + (-t165 * t71 + t166 * t72) * qJD(3) + (-t55 * t165 + t56 * t166) * qJ(3) + t202 * t156, -t47 * t166 + t139 + t199 * t117 + (t252 - t171 * t76 + (-t171 * t52 + (t129 * t171 - t237) * qJD(2)) * t165) * qJD(1), -t75 * t117 - t76 * t118 + t99 + (-t58 * t218 + t34) * t166 + (t60 * t218 + t224 + t42) * t165, -t47 * t165 + t199 * t118 + (-t251 + t171 * t75 + (qJ(3) * t226 + (-qJD(2) * t129 - qJD(3) + t52) * t171) * t166) * qJD(1), t47 * t129 - t52 * t88 - t58 * t76 - t60 * t75 + (qJ(3) * t34 + qJD(3) * t60) * t166 + (qJ(3) * t42 + qJD(3) * t58 - qJD(4) * t52) * t165, -t29 * t121 + t248 * t64, t29 * t120 + t121 * t176 + t188 * t248 + t247 * t64, -t248 * t152 + (-qJD(2) * t121 - t64) * t227, -t247 * t152 + (qJD(2) * t120 - t188) * t227, t152 * t227, -t107 * t176 + t36 * t120 + t200 * t188 + t247 * t33 + t269 * t152 + (t10 - t246) * t227, -t107 * t29 + t36 * t121 - t200 * t64 - t248 * t33 + t267 * t152 + (-t11 + t245) * t227, t3 * t120 - t45 * t176 + t247 * t9 - t254 * t188 + t255 * t152 + (-t7 - t246) * t227, -t1 * t120 + t2 * t121 + t176 * t81 + t187 * t29 + t188 * t256 - t247 * t8 - t248 * t7 + t255 * t64, -t3 * t121 + t45 * t29 + t248 * t9 - t254 * t64 - t256 * t152 + (t8 - t245) * t227, t1 * t81 - t187 * t2 - t254 * t9 - t255 * t7 - t256 * t8 + t3 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, t79, t266, -t72 * t117 + t71 * t118 + t151, t194, t266, -t79, -t60 * t117 - t58 * t118 + t47, 0, 0, 0, 0, 0, t174, t177, t174, t272 + t262, -t177, -t188 * t8 - t7 * t64 - t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118 * t117 - t204, t100 + t198, -t164 * t173 - t115, t52 * t118 + (-pkin(3) * t226 + t171 * t60) * qJD(1) - t55, 0, 0, 0, 0, 0, t175, t264, t175, t273 * t168 + t268 * t170, -t264, -t9 * t118 + t270 * t170 + (t152 * t7 + t1) * t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t259, t262 - t272, -t268, t273, -t204, t33 * t64 + t186, t10 * t152 + t188 * t33 + t184, -t18 * t188 - 0.2e1 * t150 + t186 + t260, pkin(5) * t29 + t176 * qJ(6) - (-t11 + t8) * t64 + (t7 - t231) * t188, -0.2e1 * t197 - t18 * t64 - t9 * t188 + (0.2e1 * qJD(6) - t10) * t152 - t184, -t2 * pkin(5) + t1 * qJ(6) - t7 * t11 - t9 * t18 + t231 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204 - t259, -t268, -t152 ^ 2 - t262, -t260 - t270;];
tauc_reg  = t12;
