% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:25:40
% EndTime: 2019-12-05 18:25:50
% DurationCPUTime: 3.77s
% Computational Cost: add. (3832->395), mult. (8607->514), div. (0->0), fcn. (5590->10), ass. (0->230)
t157 = cos(qJ(1));
t154 = sin(qJ(1));
t280 = g(2) * t154;
t104 = g(1) * t157 + t280;
t148 = qJ(2) + qJ(4);
t142 = sin(t148);
t134 = g(3) * t142;
t143 = cos(t148);
t214 = t104 * t143 + t134;
t158 = pkin(2) + pkin(1);
t296 = t158 * qJD(2);
t152 = sin(qJ(4));
t156 = cos(qJ(2));
t267 = pkin(3) + qJ(3);
t102 = t267 * t156;
t153 = sin(qJ(2));
t205 = qJD(2) * t267;
t75 = t156 * qJD(3) - t153 * t205;
t51 = t75 * qJD(1) + qJDD(1) * t102;
t259 = t152 * t51;
t283 = cos(qJ(4));
t101 = t267 * t153;
t149 = qJDD(2) * pkin(1);
t172 = -t153 * qJD(3) - t156 * t205;
t46 = qJDD(2) * pkin(2) + t172 * qJD(1) - qJDD(1) * t101 + t149;
t43 = t283 * t46;
t91 = qJD(1) * t101;
t170 = -t91 + t296;
t92 = qJD(1) * t102;
t219 = t283 * t92;
t48 = t152 * t170 + t219;
t14 = t48 * qJD(4) + t259 - t43;
t279 = g(3) * t143;
t295 = t14 + t279;
t145 = g(2) * t157;
t282 = g(1) * t154;
t287 = -t145 + t282;
t141 = t156 * qJDD(1);
t215 = pkin(1) * t141;
t228 = qJD(1) * qJD(2);
t208 = t153 * t228;
t223 = pkin(1) * t208 + qJDD(3);
t84 = -t215 + t223;
t294 = t287 - t84;
t251 = t142 * t157;
t292 = g(1) * t251 + t142 * t280;
t224 = qJD(2) + qJD(4);
t151 = sin(qJ(5));
t155 = cos(qJ(5));
t90 = t152 * t156 + t283 * t153;
t82 = t90 * qJD(1);
t175 = -t151 * t224 - t155 * t82;
t212 = t283 * t156;
t120 = qJD(1) * t212;
t235 = qJD(1) * t153;
t211 = t152 * t235;
t80 = -t120 + t211;
t77 = qJD(5) + t80;
t201 = t151 * t77;
t291 = t175 * t201;
t290 = t80 * t224;
t289 = t156 * t158;
t288 = t287 * t142;
t41 = t224 * pkin(4) + t48;
t87 = -qJD(1) * t289 + qJD(3);
t55 = -t82 * pkin(4) + t87;
t21 = -t151 * t41 + t155 * t55;
t22 = t151 * t55 + t155 * t41;
t286 = -t151 * t21 + t155 * t22;
t222 = qJDD(2) + qJDD(4);
t250 = t152 * t153;
t187 = t224 * t250;
t204 = qJDD(1) * t283;
t196 = -t224 * t120 - t152 * t141 - t153 * t204;
t37 = qJD(1) * t187 + t196;
t19 = -t175 * qJD(5) - t151 * t37 - t155 * t222;
t285 = t82 ^ 2;
t147 = t156 ^ 2;
t284 = 0.2e1 * t147;
t278 = g(3) * t156;
t277 = t21 * t77;
t276 = t22 * t77;
t164 = t283 * t170;
t232 = qJD(4) * t152;
t13 = qJD(4) * t164 + t152 * t46 - t92 * t232 + t283 * t51;
t11 = t222 * pkin(4) + t13;
t59 = pkin(2) * t208 - qJDD(1) * t289 + t223;
t23 = t37 * pkin(4) + t59;
t3 = qJD(5) * t21 + t155 * t11 + t151 * t23;
t2 = t3 * t155;
t258 = t152 * t92;
t47 = -t164 + t258;
t52 = -t152 * t91 + t219;
t275 = t47 * t52;
t274 = t47 * t80;
t197 = t155 * t224;
t60 = t151 * t82 - t197;
t273 = t60 * t80;
t272 = t60 * t82;
t271 = t175 * t60;
t270 = t175 * t82;
t269 = t77 * t82;
t268 = t82 * t80;
t230 = qJD(5) * t155;
t266 = -t151 * t19 - t60 * t230;
t265 = qJD(2) * pkin(1);
t140 = t153 * qJDD(1);
t185 = t152 * t140 - t156 * t204;
t57 = t224 * t90;
t38 = t57 * qJD(1) + t185;
t36 = qJDD(5) + t38;
t263 = t151 * t36;
t262 = t151 * t60;
t261 = t151 * t80;
t260 = t152 * t47;
t256 = t155 * t175;
t255 = t155 * t80;
t231 = qJD(5) * t151;
t18 = -qJD(5) * t197 - t151 * t222 + t155 * t37 + t82 * t231;
t254 = t18 * t151;
t253 = t19 * t155;
t252 = qJD(5) * t77;
t249 = t153 * t154;
t248 = t153 * t157;
t160 = qJD(1) ^ 2;
t247 = t153 * t160;
t246 = t154 * t151;
t245 = t154 * t155;
t244 = t154 * t156;
t243 = t156 * t157;
t241 = t156 * t160;
t240 = t157 * t267;
t239 = t157 * t151;
t238 = t157 * t155;
t93 = t158 * t235;
t94 = t153 * t296;
t146 = t153 ^ 2;
t237 = t146 - t147;
t236 = t146 + t284;
t234 = qJD(1) * t156;
t220 = pkin(1) * t234;
t114 = qJD(3) - t220;
t229 = -qJD(3) + t114;
t227 = qJD(1) * qJD(3);
t226 = qJDD(2) * t156;
t225 = t147 * qJDD(1);
t221 = pkin(4) * t252;
t118 = t152 * t158 + pkin(4);
t218 = t118 * t252;
t217 = t90 * t231;
t216 = t90 * t230;
t44 = t47 * t231;
t45 = t47 * t230;
t213 = g(1) * t243 + g(2) * t244 + g(3) * t153;
t209 = qJD(4) * t283;
t207 = t156 * t228;
t159 = qJD(2) ^ 2;
t203 = -qJ(3) * t159 - t84;
t200 = t155 * t77;
t199 = t154 * t267 + t157 * t289;
t195 = t158 * t209;
t194 = t153 * t207;
t193 = g(1) * t248 + g(2) * t249 - t278;
t192 = -pkin(4) * t36 + t274;
t176 = -t283 * t101 - t152 * t102;
t64 = -t152 * t101 + t283 * t102;
t27 = t64 * qJD(4) + t152 * t75 - t283 * t172;
t190 = -t14 * t176 + t47 * t27;
t56 = -qJD(2) * t212 - t156 * t209 + t187;
t189 = t14 * t90 - t47 * t56;
t188 = t36 * t90 - t56 * t77;
t186 = -0.2e1 * t194;
t184 = t151 * t22 + t155 * t21;
t66 = -t90 * pkin(4) - t289;
t33 = t151 * t66 + t155 * t64;
t32 = -t151 * t64 + t155 * t66;
t183 = t292 * t155 - t21 * t82 + t44;
t182 = t295 * t151 + t22 * t82 + t45;
t181 = t155 * t36 + (-t231 - t261) * t77;
t180 = -t21 * t255 - t22 * t261 + t2 - t214;
t179 = -t47 * t82 - t287;
t178 = -qJD(5) * t55 - t11 + t134;
t177 = t104 * t142;
t174 = -qJ(3) * qJDD(1) + (-qJD(3) - t114) * qJD(1);
t173 = -t87 * t82 - t279 + t292 + t43;
t100 = -qJ(3) * t235 + t265;
t171 = -t100 * qJD(2) * t156 - t104;
t169 = t87 * t80 - t13 + t214;
t20 = t155 * t23;
t4 = -qJD(5) * t22 - t151 * t11 + t20;
t167 = -t184 * qJD(5) - t4 * t151;
t166 = -t118 * t36 - t77 * t195 + t274;
t165 = -t224 * t211 - t196;
t163 = t167 + t2;
t162 = t104 * t153 - t283 * t14 - t278;
t161 = qJ(3) ^ 2;
t130 = g(1) * t244;
t129 = g(2) * t248;
t119 = t153 * t241;
t99 = t237 * t160;
t98 = -t159 * t153 + t226;
t97 = qJDD(2) * t153 + t159 * t156;
t89 = -t212 + t250;
t86 = t186 + t225;
t85 = t146 * qJDD(1) + 0.2e1 * t194;
t78 = t80 ^ 2;
t74 = t143 * t238 + t246;
t73 = -t143 * t239 + t245;
t72 = -t143 * t245 + t239;
t71 = t143 * t246 + t238;
t67 = 0.2e1 * t153 * t141 - 0.2e1 * t237 * t228;
t65 = -t153 * t227 + t149 + (-t207 - t140) * qJ(3);
t58 = t80 * pkin(4) + t93;
t53 = -t283 * t91 - t258;
t40 = t56 * pkin(4) + t94;
t39 = -t78 + t285;
t31 = pkin(4) * t261 - t155 * t47;
t30 = pkin(4) * t255 + t151 * t47;
t29 = t151 * t58 + t155 * t53;
t28 = -t151 * t53 + t155 * t58;
t26 = t176 * qJD(4) + t152 * t172 + t283 * t75;
t24 = t165 + t290;
t10 = t77 * t200 + t263 + t270;
t9 = t181 + t272;
t8 = t60 * t201 - t253;
t7 = -t175 * t200 - t254;
t6 = -t33 * qJD(5) - t151 * t26 + t155 * t40;
t5 = t32 * qJD(5) + t151 * t40 + t155 * t26;
t1 = (-t18 - t273) * t155 + t291 + t266;
t12 = [0, 0, 0, 0, 0, qJDD(1), t287, t104, 0, 0, t85, t67, t97, t86, t98, 0, -g(2) * t243 + t130, -g(1) * t249 + t129, -t104, 0, t85, t67, t97, t86, t98, 0, pkin(1) * t225 + t130 + (t203 - t145) * t156 + (-qJ(3) * qJDD(2) + (-0.2e1 * t220 + t229) * qJD(2)) * t153, -qJ(3) * t226 + t129 + (-t203 - t215 - t282) * t153 + (t237 * qJD(1) * pkin(1) + t229 * t156) * qJD(2), -t65 * t153 + t236 * t227 + (t236 * qJDD(1) + t186) * qJ(3) + t171, t161 * t225 + t294 * t156 * pkin(1) + (t227 * t284 + t171) * qJ(3) + (-t65 * qJ(3) - t100 * qJD(3) + (pkin(1) * t114 - 0.2e1 * t161 * t234) * qJD(2)) * t153, -t37 * t90 - t56 * t82, t37 * t89 - t38 * t90 + t56 * t80 - t57 * t82, t222 * t90 - t224 * t56, t38 * t89 + t57 * t80, -t222 * t89 - t224 * t57, 0, t143 * t287 + t176 * t222 - t224 * t27 - t289 * t38 + t87 * t57 + t59 * t89 + t94 * t80, -t222 * t64 - t224 * t26 + t289 * t37 - t87 * t56 + t59 * t90 + t94 * t82 - t288, -t13 * t89 + t176 * t37 - t26 * t80 + t27 * t82 - t38 * t64 - t48 * t57 - t104 + t189, t13 * t64 + t48 * t26 - t59 * t289 + t87 * t94 - g(1) * (-t154 * t289 + t240) - g(2) * t199 + t190, t175 * t217 + (t175 * t56 - t18 * t90) * t155, (-t151 * t175 + t155 * t60) * t56 + (t254 - t253 + (t256 + t262) * qJD(5)) * t90, t155 * t188 - t175 * t57 - t18 * t89 - t217 * t77, t60 * t216 + (t19 * t90 - t56 * t60) * t151, -t151 * t188 - t19 * t89 - t216 * t77 - t60 * t57, t36 * t89 + t57 * t77, -g(1) * t72 - g(2) * t74 + t151 * t189 - t176 * t19 + t21 * t57 + t27 * t60 + t32 * t36 + t4 * t89 + t45 * t90 + t6 * t77, -g(1) * t71 - g(2) * t73 + t155 * t189 - t175 * t27 + t176 * t18 - t22 * t57 - t3 * t89 - t33 * t36 - t44 * t90 - t5 * t77, t32 * t18 - t33 * t19 - t5 * t60 + t6 * t175 + t184 * t56 + t288 + (-qJD(5) * t286 - t3 * t151 - t4 * t155) * t90, t3 * t33 + t22 * t5 + t4 * t32 + t21 * t6 - g(1) * (t240 + (-pkin(4) * t142 - t289) * t154) - g(2) * (pkin(4) * t251 + t199) + t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t99, t140, t119, t141, qJDD(2), t193, t213, 0, 0, -t119, t99, t140, t119, t141, qJDD(2), 0.2e1 * t149 + (pkin(1) * t241 + t174) * t153 + t193, -t146 * t160 * pkin(1) + t174 * t156 + t213, -pkin(1) * t140 + (qJ(3) * t247 + (t100 - t265) * qJD(1)) * t156, (qJ(3) * qJD(1) * t100 + t161 * t247) * t156 + (-t278 + t65 + (-qJD(1) * t114 + t104) * t153) * pkin(1), t268, t39, t24, -t268, -t185, t222, t283 * t158 * t222 - t92 * t209 - t93 * t80 + t52 * t224 + (-t51 + (-t158 * t224 - t170) * qJD(4)) * t152 + t173, -t93 * t82 + t53 * t224 + (-t152 * t222 - t209 * t224) * t158 + t169, (t48 - t52) * t82 + (t47 + t53) * t80 + (t283 * t37 - t152 * t38 + (t152 * t82 - t283 * t80) * qJD(4)) * t158, -t275 - t48 * t53 - t87 * t93 + (t13 * t152 + (t283 * t48 + t260) * qJD(4) + t162) * t158, t7, t1, t10, t8, t9, -t269, -t28 * t77 - t52 * t60 + (-t283 * t19 + t60 * t232) * t158 + (-t295 - t218) * t155 + t166 * t151 + t183, t29 * t77 + t52 * t175 + (-t175 * t232 + t283 * t18) * t158 + t166 * t155 + (-t177 + t218) * t151 + t182, -t28 * t175 + t29 * t60 + (-t60 * t195 - t118 * t19 + (-t118 * t175 - t21) * qJD(5)) * t155 + (-t175 * t195 - t118 * t18 - t4 + (t118 * t60 - t22) * qJD(5)) * t151 + t180, -t21 * t28 - t22 * t29 - t275 - t214 * pkin(4) + ((t283 * t286 + t260) * qJD(4) + t162) * t158 + t163 * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141 + 0.2e1 * t208, t140 + 0.2e1 * t207, (-t146 - t147) * t160, -qJ(3) * t147 * t160 + t100 * t235 - t294, 0, 0, 0, 0, 0, 0, t82 * t224 + t38, t165 - t290, -t78 - t285, t48 * t80 + t179 + t59, 0, 0, 0, 0, 0, 0, t181 - t272, -t155 * t77 ^ 2 - t263 + t270, (t18 - t273) * t155 - t291 + t266, (t4 + t276) * t155 + (t3 - t277) * t151 + t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, t39, t24, -t268, -t185, t222, t48 * qJD(2) + t173 - t259, -t224 * t47 + t169, 0, 0, t7, t1, t10, t8, t9, -t269, -t30 * t77 - t48 * t60 + t192 * t151 + (-t295 - t221) * t155 + t183, t31 * t77 + t48 * t175 + t192 * t155 + (-t177 + t221) * t151 + t182, -t30 * t175 + t31 * t60 + (-t254 - t253 + (-t256 + t262) * qJD(5)) * pkin(4) + t167 + t180, -t21 * t30 - t22 * t31 - t47 * t48 + (t163 - t214) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t271, t175 ^ 2 - t60 ^ 2, t60 * t77 - t18, t271, -t175 * t77 - t19, t36, -g(1) * t73 + g(2) * t71 + t151 * t178 + t175 * t47 - t230 * t41 + t20 + t276, g(1) * t74 - g(2) * t72 + t277 + t47 * t60 + (qJD(5) * t41 - t23) * t151 + t178 * t155, 0, 0;];
tau_reg = t12;
