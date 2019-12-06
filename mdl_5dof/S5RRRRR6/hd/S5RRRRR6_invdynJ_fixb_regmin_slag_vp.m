% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:38
% EndTime: 2019-12-05 19:00:45
% DurationCPUTime: 2.44s
% Computational Cost: add. (3550->296), mult. (5394->388), div. (0->0), fcn. (3946->16), ass. (0->208)
t178 = qJ(3) + qJ(4);
t169 = qJ(5) + t178;
t154 = sin(t169);
t155 = cos(t169);
t180 = sin(qJ(5));
t185 = cos(qJ(5));
t175 = qJD(1) + qJD(2);
t186 = cos(qJ(4));
t187 = cos(qJ(3));
t262 = t186 * t187;
t241 = t175 * t262;
t181 = sin(qJ(4));
t182 = sin(qJ(3));
t267 = t181 * t182;
t242 = t175 * t267;
t92 = -t241 + t242;
t263 = t182 * t186;
t114 = t181 * t187 + t263;
t94 = t114 * t175;
t212 = t180 * t92 - t185 * t94;
t179 = qJ(1) + qJ(2);
t166 = sin(t179);
t168 = cos(t179);
t234 = -g(2) * t166 + g(3) * t168;
t172 = qJDD(3) + qJDD(4);
t183 = sin(qJ(2));
t288 = pkin(1) * t183;
t245 = qJD(1) * t288;
t125 = pkin(7) * t175 + t245;
t235 = pkin(8) * t175 + t125;
t84 = t235 * t187;
t77 = t186 * t84;
t83 = t235 * t182;
t78 = qJD(3) * pkin(3) - t83;
t211 = -t181 * t78 - t77;
t252 = qJD(3) * t187;
t236 = t175 * t252;
t173 = qJDD(1) + qJDD(2);
t264 = t182 * t173;
t247 = qJDD(1) * t183;
t188 = cos(qJ(2));
t254 = qJD(2) * t188;
t97 = t173 * pkin(7) + (qJD(1) * t254 + t247) * pkin(1);
t41 = -t125 * t252 + qJDD(3) * pkin(3) - t182 * t97 + (-t236 - t264) * pkin(8);
t253 = qJD(3) * t182;
t237 = t175 * t253;
t261 = t187 * t173;
t42 = -t125 * t253 + t187 * t97 + (-t237 + t261) * pkin(8);
t200 = t211 * qJD(4) - t181 * t42 + t186 * t41;
t174 = qJD(3) + qJD(4);
t37 = qJD(4) * t241 + t173 * t263 - t174 * t242 + t181 * t261 + t186 * t236;
t5 = pkin(4) * t172 - pkin(9) * t37 + t200;
t159 = -pkin(3) * t187 - pkin(2);
t256 = qJD(1) * t188;
t243 = pkin(1) * t256;
t95 = t159 * t175 - t243;
t57 = t92 * pkin(4) + t95;
t251 = qJD(4) * t181;
t292 = (qJD(4) * t78 + t42) * t186 + t181 * t41 - t84 * t251;
t216 = t181 * t264 - t186 * t261;
t67 = t174 * t114;
t38 = t175 * t67 + t216;
t6 = -t38 * pkin(9) + t292;
t299 = -g(1) * t155 + t234 * t154 - t180 * t6 + t185 * t5 + t57 * t212;
t249 = qJD(5) * t180;
t273 = t155 * t168;
t274 = t155 * t166;
t290 = pkin(9) * t92;
t31 = -t211 - t290;
t54 = -t180 * t94 - t185 * t92;
t298 = g(1) * t154 - g(2) * t274 + g(3) * t273 + t31 * t249 - t57 * t54;
t276 = t185 * t31;
t75 = t181 * t84;
t230 = t186 * t78 - t75;
t90 = t94 * pkin(9);
t30 = t230 - t90;
t28 = pkin(4) * t174 + t30;
t215 = -t180 * t28 - t276;
t297 = t215 * qJD(5) + t299;
t246 = -qJD(4) - qJD(5);
t164 = qJD(3) - t246;
t296 = (-t164 * t31 - t5) * t180 + t298;
t282 = t212 * t54;
t295 = g(2) * t168 + g(3) * t166;
t9 = t212 ^ 2 - t54 ^ 2;
t248 = qJD(5) * t185;
t10 = -t180 * t38 + t185 * t37 - t92 * t248 - t94 * t249;
t7 = -t164 * t54 + t10;
t199 = t212 * qJD(5) - t180 * t37 - t185 * t38;
t8 = -t164 * t212 + t199;
t113 = -t262 + t267;
t291 = -pkin(7) - pkin(8);
t239 = qJD(3) * t291;
t118 = t182 * t239;
t119 = t187 * t239;
t142 = t291 * t182;
t170 = t187 * pkin(8);
t143 = pkin(7) * t187 + t170;
t250 = qJD(4) * t186;
t294 = -t113 * t243 - t186 * t118 - t181 * t119 - t142 * t250 + t143 * t251;
t259 = t181 * t142 + t186 * t143;
t293 = -t259 * qJD(4) + t114 * t243 - t181 * t118 + t186 * t119;
t156 = pkin(7) + t288;
t278 = -pkin(8) - t156;
t110 = t278 * t182;
t111 = t156 * t187 + t170;
t260 = t181 * t110 + t186 * t111;
t66 = t174 * t113;
t289 = t66 * pkin(9);
t287 = pkin(1) * t188;
t286 = pkin(2) * t173;
t285 = pkin(2) * t175;
t284 = pkin(9) * t114;
t279 = t94 * t92;
t277 = -t186 * t83 - t75;
t126 = -t243 - t285;
t255 = qJD(2) * t183;
t162 = pkin(1) * t255;
t258 = -qJD(1) * t162 + qJDD(1) * t287;
t96 = -t258 - t286;
t275 = t126 * t252 + t96 * t182;
t163 = qJDD(5) + t172;
t272 = t163 * t181;
t271 = t163 * t185;
t167 = cos(t178);
t270 = t166 * t167;
t269 = t167 * t168;
t268 = t175 * t182;
t266 = t181 * t185;
t176 = t182 ^ 2;
t257 = -t187 ^ 2 + t176;
t244 = pkin(1) * t254;
t161 = pkin(3) * t253;
t240 = t126 * t253 + t295 * t187;
t238 = t175 * t255;
t60 = pkin(4) * t67 + t161;
t233 = -qJD(5) * t28 - t6;
t229 = t181 * t83 - t77;
t228 = qJD(3) * t278;
t226 = t186 * t110 - t111 * t181;
t225 = t186 * t142 - t143 * t181;
t64 = -t113 * t180 + t114 * t185;
t22 = t64 * qJD(5) - t180 * t66 + t185 * t67;
t61 = pkin(3) * t237 + t159 * t173 - t258;
t23 = pkin(4) * t38 + t61;
t63 = t185 * t113 + t114 * t180;
t224 = g(2) * t273 + g(3) * t274 + t57 * t22 + t23 * t63;
t223 = g(2) * t269 + g(3) * t270 + t61 * t113 + t95 * t67;
t222 = t175 * t245;
t221 = t258 + t295;
t220 = t60 - t245;
t58 = t225 - t284;
t65 = t67 * pkin(9);
t219 = -qJD(5) * t58 + t294 + t65;
t109 = t113 * pkin(9);
t59 = -t109 + t259;
t218 = qJD(5) * t59 - t289 - t293;
t47 = t226 - t284;
t48 = -t109 + t260;
t214 = -t180 * t48 + t185 * t47;
t213 = t180 * t47 + t185 * t48;
t91 = pkin(4) * t113 + t159;
t79 = t182 * t228 + t187 * t244;
t80 = -t182 * t244 + t187 * t228;
t209 = t110 * t250 - t111 * t251 + t181 * t80 + t186 * t79;
t207 = -t161 + t245;
t206 = -t126 * t175 + t234 - t97;
t21 = -t63 * qJD(5) - t180 * t67 - t185 * t66;
t205 = -t154 * t295 + t57 * t21 + t23 * t64;
t165 = sin(t178);
t204 = t61 * t114 - t165 * t295 - t95 * t66;
t190 = qJD(3) ^ 2;
t203 = -pkin(7) * t190 + t222 + t286;
t158 = -pkin(2) - t287;
t202 = -pkin(1) * t238 - t156 * t190 - t158 * t173;
t198 = -pkin(7) * qJDD(3) + (t243 - t285) * qJD(3);
t197 = -t260 * qJD(4) - t181 * t79 + t186 * t80;
t195 = -qJDD(3) * t156 + (t158 * t175 - t244) * qJD(3);
t193 = g(1) * t165 - g(2) * t270 + g(3) * t269 + t95 * t92 - t292;
t191 = -g(1) * t167 + t234 * t165 - t95 * t94 + t200;
t189 = cos(qJ(1));
t184 = sin(qJ(1));
t171 = t175 ^ 2;
t157 = pkin(3) * t186 + pkin(4);
t135 = t159 - t287;
t134 = qJDD(3) * t187 - t182 * t190;
t133 = qJDD(3) * t182 + t187 * t190;
t120 = t162 + t161;
t98 = t173 * t176 + 0.2e1 * t182 * t236;
t82 = t91 - t287;
t69 = -0.2e1 * t257 * t175 * qJD(3) + 0.2e1 * t182 * t261;
t68 = pkin(3) * t268 + pkin(4) * t94;
t56 = t162 + t60;
t44 = -t113 * t172 - t174 * t67;
t43 = t114 * t172 - t174 * t66;
t40 = -t92 ^ 2 + t94 ^ 2;
t33 = -t90 + t277;
t32 = t229 + t290;
t26 = t174 * t92 + t37;
t20 = t197 + t289;
t19 = t209 - t65;
t18 = t114 * t37 - t66 * t94;
t13 = -t163 * t63 - t164 * t22;
t12 = t163 * t64 + t164 * t21;
t3 = -t113 * t37 - t114 * t38 + t66 * t92 - t67 * t94;
t2 = t10 * t64 - t21 * t212;
t1 = -t10 * t63 + t199 * t64 + t21 * t54 + t212 * t22;
t4 = [qJDD(1), g(2) * t189 + g(3) * t184, -g(2) * t184 + g(3) * t189, t173, (t173 * t188 - t238) * pkin(1) + t221, ((-qJDD(1) - t173) * t183 + (-qJD(1) - t175) * t254) * pkin(1) + t234, t98, t69, t133, t134, 0, t195 * t182 + (t202 - t96) * t187 + t240, t195 * t187 + (-t202 - t295) * t182 + t275, t18, t3, t43, t44, 0, t120 * t92 + t135 * t38 + t226 * t172 + t197 * t174 + t223, t120 * t94 + t135 * t37 - t260 * t172 - t209 * t174 + t204, t2, t1, t12, t13, 0, -t56 * t54 - t82 * t199 + (-t213 * qJD(5) - t180 * t19 + t185 * t20) * t164 + t214 * t163 + t224, -t56 * t212 + t82 * t10 - (qJD(5) * t214 + t180 * t20 + t185 * t19) * t164 - t213 * t163 + t205; 0, 0, 0, t173, t221 + t222, (-t247 + (-qJD(2) + t175) * t256) * pkin(1) + t234, t98, t69, t133, t134, 0, t198 * t182 + (t203 - t96) * t187 + t240, t198 * t187 + (-t203 - t295) * t182 + t275, t18, t3, t43, t44, 0, t159 * t38 + t225 * t172 + t293 * t174 - t207 * t92 + t223, t159 * t37 - t259 * t172 + t294 * t174 - t207 * t94 + t204, t2, t1, t12, t13, 0, -t91 * t199 + (-t180 * t59 + t185 * t58) * t163 - t220 * t54 + (t219 * t180 - t218 * t185) * t164 + t224, t91 * t10 - (t180 * t58 + t185 * t59) * t163 - t220 * t212 + (t180 * t218 + t185 * t219) * t164 + t205; 0, 0, 0, 0, 0, 0, -t182 * t171 * t187, t257 * t171, t264, t261, qJDD(3), -g(1) * t187 + t206 * t182, g(1) * t182 + t206 * t187, t279, t40, t26, -t216, t172, -t229 * t174 + (t172 * t186 - t174 * t251 - t92 * t268) * pkin(3) + t191, t277 * t174 + (-t172 * t181 - t174 * t250 - t94 * t268) * pkin(3) + t193, t282, t9, t7, t8, t163, t157 * t271 + t68 * t54 - (-t180 * t33 + t185 * t32) * t164 + (-t180 * t272 + (-t180 * t186 - t266) * t164 * qJD(4)) * pkin(3) + ((-pkin(3) * t266 - t157 * t180) * t164 + t215) * qJD(5) + t299, t68 * t212 + (-t157 * t163 - t5 + (-pkin(3) * t181 * t246 + t32) * t164) * t180 + (-pkin(3) * t272 + (-pkin(3) * t250 - qJD(5) * t157 + t33) * t164 + t233) * t185 + t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t279, t40, t26, -t216, t172, -t211 * t174 + t191, t230 * t174 + t193, t282, t9, t7, t8, t163, -(-t180 * t30 - t276) * t164 + (-t164 * t249 + t54 * t94 + t271) * pkin(4) + t297, (t164 * t30 + t233) * t185 + (-t163 * t180 - t164 * t248 + t212 * t94) * pkin(4) + t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t282, t9, t7, t8, t163, -t215 * t164 + t297, (-t6 + (-qJD(5) + t164) * t28) * t185 + t296;];
tau_reg = t4;
