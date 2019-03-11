% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:39:01
% EndTime: 2019-03-09 03:39:09
% DurationCPUTime: 3.72s
% Computational Cost: add. (4780->394), mult. (10541->531), div. (0->0), fcn. (7924->18), ass. (0->226)
t179 = sin(pkin(11));
t186 = sin(qJ(3));
t252 = qJD(1) * t186;
t181 = cos(pkin(11));
t190 = cos(qJ(3));
t258 = t181 * t190;
t129 = qJD(1) * t258 - t179 * t252;
t297 = qJD(5) + qJD(6);
t307 = t129 - t297;
t125 = qJD(5) - t129;
t156 = pkin(3) * t179 + pkin(8);
t169 = t190 * qJDD(2);
t180 = sin(pkin(10));
t157 = pkin(1) * t180 + pkin(7);
t147 = t157 * qJDD(1);
t200 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t147;
t255 = qJ(4) + t157;
t227 = t255 * qJD(1);
t212 = t227 * qJD(3);
t49 = qJDD(3) * pkin(3) - t200 * t186 - t190 * t212 + t169;
t53 = (qJDD(2) - t212) * t186 + t200 * t190;
t21 = -t179 * t53 + t181 * t49;
t19 = -qJDD(3) * pkin(4) - t21;
t174 = qJ(3) + pkin(11);
t163 = sin(t174);
t165 = cos(t174);
t175 = qJ(1) + pkin(10);
t164 = sin(t175);
t166 = cos(t175);
t223 = g(1) * t166 + g(2) * t164;
t202 = -g(3) * t165 + t223 * t163;
t306 = -qJD(5) * t156 * t125 - t19 + t202;
t185 = sin(qJ(5));
t188 = cos(qJ(6));
t184 = sin(qJ(6));
t189 = cos(qJ(5));
t257 = t184 * t189;
t141 = t185 * t188 + t257;
t283 = t307 * t141;
t139 = t179 * t190 + t181 * t186;
t131 = t139 * qJD(1);
t247 = t189 * qJD(3);
t110 = t131 * t185 - t247;
t112 = qJD(3) * t185 + t131 * t189;
t214 = t110 * t184 - t188 * t112;
t55 = t188 * t110 + t112 * t184;
t305 = t214 * t55;
t251 = qJD(5) * t185;
t271 = t129 * t185;
t304 = t251 - t271;
t303 = t214 ^ 2 - t55 ^ 2;
t123 = qJD(6) + t125;
t248 = qJD(6) * t188;
t249 = qJD(6) * t184;
t246 = qJD(1) * qJD(3);
t235 = t190 * t246;
t236 = t186 * t246;
t94 = qJDD(1) * t139 - t179 * t236 + t181 * t235;
t42 = qJD(5) * t247 + t185 * qJDD(3) - t131 * t251 + t189 * t94;
t43 = qJD(5) * t112 - t189 * qJDD(3) + t185 * t94;
t9 = -t110 * t248 - t112 * t249 - t184 * t43 + t188 * t42;
t302 = t123 * t55 + t9;
t178 = qJ(5) + qJ(6);
t172 = sin(t178);
t263 = t166 * t172;
t173 = cos(t178);
t266 = t164 * t173;
t106 = -t165 * t266 + t263;
t262 = t166 * t173;
t267 = t164 * t172;
t108 = t165 * t262 + t267;
t113 = t190 * qJD(2) - t227 * t186;
t282 = qJD(3) * pkin(3);
t104 = t113 + t282;
t114 = t186 * qJD(2) + t227 * t190;
t259 = t181 * t114;
t52 = t179 * t104 + t259;
t46 = qJD(3) * pkin(8) + t52;
t162 = pkin(3) * t190 + pkin(2);
t182 = cos(pkin(10));
t293 = pkin(1) * t182;
t213 = -t162 - t293;
t128 = qJD(1) * t213 + qJD(4);
t67 = -t129 * pkin(4) - t131 * pkin(8) + t128;
t27 = t185 * t67 + t189 * t46;
t15 = -pkin(9) * t110 + t27;
t13 = t15 * t249;
t292 = g(3) * t163;
t101 = t179 * t114;
t51 = t104 * t181 - t101;
t45 = -qJD(3) * pkin(4) - t51;
t34 = pkin(5) * t110 + t45;
t301 = g(1) * t108 - g(2) * t106 + t173 * t292 + t34 * t55 + t13;
t105 = t165 * t267 + t262;
t107 = -t165 * t263 + t266;
t22 = t179 * t49 + t181 * t53;
t20 = qJDD(3) * pkin(8) + t22;
t199 = pkin(3) * t236 + qJDD(1) * t213 + qJDD(4);
t244 = t190 * qJDD(1);
t245 = t186 * qJDD(1);
t217 = -t179 * t245 + t181 * t244;
t93 = -qJD(3) * t131 + t217;
t32 = -t93 * pkin(4) - t94 * pkin(8) + t199;
t31 = t189 * t32;
t130 = t139 * qJD(3);
t92 = qJD(1) * t130 + qJDD(5) - t217;
t2 = t92 * pkin(5) - t42 * pkin(9) - qJD(5) * t27 - t185 * t20 + t31;
t250 = qJD(5) * t189;
t210 = t185 * t32 + t189 * t20 + t67 * t250 - t46 * t251;
t3 = -pkin(9) * t43 + t210;
t240 = -t184 * t3 + t188 * t2;
t26 = -t185 * t46 + t189 * t67;
t14 = -pkin(9) * t112 + t26;
t12 = pkin(5) * t125 + t14;
t279 = t15 * t188;
t5 = t12 * t184 + t279;
t300 = -g(1) * t107 + g(2) * t105 - t5 * qJD(6) + t172 * t292 + t34 * t214 + t240;
t198 = qJD(6) * t214 - t184 * t42 - t188 * t43;
t299 = -t123 * t214 + t198;
t85 = t141 * t139;
t140 = t184 * t185 - t188 * t189;
t284 = t307 * t140;
t88 = qJDD(6) + t92;
t296 = -t284 * t123 - t141 * t88;
t239 = t139 * t251;
t138 = t179 * t186 - t258;
t133 = t138 * qJD(3);
t270 = t133 * t189;
t207 = t239 + t270;
t268 = t139 * t189;
t295 = -t125 * t207 + t92 * t268;
t294 = -t130 * t214 + t9 * t138;
t290 = g(3) * t190;
t289 = pkin(9) + t156;
t256 = t185 * t133;
t269 = t139 * t185;
t24 = -t133 * t257 - t184 * t239 - t249 * t269 + (t268 * t297 - t256) * t188;
t288 = -t24 * t123 - t85 * t88;
t287 = t112 * t130 + t42 * t138;
t60 = t113 * t181 - t101;
t79 = pkin(3) * t252 + pkin(4) * t131 - pkin(8) * t129;
t286 = t185 * t79 + t189 * t60;
t134 = t255 * t186;
t135 = t255 * t190;
t90 = -t134 * t179 + t135 * t181;
t82 = t189 * t90;
t83 = pkin(4) * t138 - pkin(8) * t139 + t213;
t285 = t185 * t83 + t82;
t281 = t131 * t55;
t278 = t185 * t92;
t277 = t42 * t185;
t276 = t214 * t131;
t275 = t110 * t125;
t274 = t110 * t131;
t273 = t112 * t125;
t272 = t112 * t131;
t265 = t164 * t185;
t264 = t164 * t189;
t261 = t166 * t185;
t260 = t166 * t189;
t254 = qJDD(2) - g(3);
t176 = t186 ^ 2;
t253 = -t190 ^ 2 + t176;
t159 = -pkin(2) - t293;
t150 = qJD(1) * t159;
t242 = t186 * t282;
t158 = -pkin(3) * t181 - pkin(4);
t238 = t139 * t250;
t234 = qJD(6) * t12 + t3;
t232 = qJD(5) * t289;
t231 = -qJD(5) * t67 - t20;
t59 = t113 * t179 + t259;
t229 = qJD(3) * t255;
t115 = t190 * qJD(4) - t186 * t229;
t116 = -t186 * qJD(4) - t190 * t229;
t65 = t179 * t115 - t181 * t116;
t89 = t181 * t134 + t135 * t179;
t228 = t125 * t189;
t226 = t283 * t123 - t140 * t88;
t225 = pkin(5) * t304 - t59;
t224 = -t46 * t250 + t31;
t222 = -g(1) * t164 + g(2) * t166;
t187 = sin(qJ(1));
t191 = cos(qJ(1));
t221 = g(1) * t187 - g(2) * t191;
t136 = t289 * t185;
t220 = -pkin(9) * t271 + qJD(6) * t136 + t185 * t232 + t286;
t137 = t289 * t189;
t72 = t189 * t79;
t219 = pkin(5) * t131 + qJD(6) * t137 - t185 * t60 + t72 + (-pkin(9) * t129 + t232) * t189;
t23 = t140 * t133 - t297 * t85;
t86 = t140 * t139;
t218 = -t123 * t23 + t86 * t88;
t216 = -t130 * t55 + t138 * t198;
t215 = -t130 * t110 - t138 * t43;
t211 = -t125 * t304 + t189 * t92;
t66 = t115 * t181 + t116 * t179;
t80 = pkin(4) * t130 + pkin(8) * t133 + t242;
t209 = t185 * t80 + t189 * t66 + t83 * t250 - t90 * t251;
t208 = t238 - t256;
t205 = t125 * t45 - t156 * t92;
t204 = -qJD(1) * t150 - t147 + t223;
t203 = 0.2e1 * t150 * qJD(3) - qJDD(3) * t157;
t192 = qJD(3) ^ 2;
t197 = -0.2e1 * qJDD(1) * t159 - t157 * t192 - t222;
t196 = -t125 * t208 - t92 * t269;
t193 = qJD(1) ^ 2;
t183 = -qJ(4) - pkin(7);
t146 = qJDD(3) * t190 - t186 * t192;
t145 = qJDD(3) * t186 + t190 * t192;
t144 = -pkin(5) * t189 + t158;
t120 = t165 * t260 + t265;
t119 = -t165 * t261 + t264;
t118 = -t165 * t264 + t261;
t117 = t165 * t265 + t260;
t78 = t189 * t83;
t73 = t189 * t80;
t64 = pkin(5) * t269 + t89;
t33 = pkin(5) * t208 + t65;
t29 = -pkin(9) * t269 + t285;
t28 = pkin(5) * t138 - pkin(9) * t268 - t185 * t90 + t78;
t11 = pkin(5) * t43 + t19;
t7 = -pkin(9) * t208 + t209;
t6 = pkin(9) * t270 + t130 * pkin(5) - t185 * t66 + t73 + (-t82 + (pkin(9) * t139 - t83) * t185) * qJD(5);
t4 = t12 * t188 - t15 * t184;
t1 = [qJDD(1), t221, g(1) * t191 + g(2) * t187 (t221 + (t180 ^ 2 + t182 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t176 + 0.2e1 * t186 * t235, 0.2e1 * t186 * t244 - 0.2e1 * t253 * t246, t145, t146, 0, t186 * t203 + t190 * t197, -t186 * t197 + t190 * t203, t129 * t66 - t130 * t52 + t131 * t65 + t133 * t51 - t138 * t22 - t139 * t21 + t89 * t94 + t90 * t93 - t223, t22 * t90 + t52 * t66 - t21 * t89 - t51 * t65 + t128 * t242 - g(1) * (-pkin(1) * t187 - t162 * t164 - t166 * t183) - g(2) * (pkin(1) * t191 + t162 * t166 - t164 * t183) + t199 * t213, -t112 * t207 + t268 * t42 -(-t110 * t189 - t112 * t185) * t133 + (-t277 - t189 * t43 + (t110 * t185 - t112 * t189) * qJD(5)) * t139, t287 + t295, t196 + t215, t125 * t130 + t138 * t92 (-t250 * t90 + t73) * t125 + t78 * t92 + t224 * t138 + t26 * t130 + t65 * t110 + t89 * t43 + t45 * t238 - g(1) * t118 - g(2) * t120 + ((-qJD(5) * t83 - t66) * t125 - t90 * t92 + t231 * t138 + t19 * t139 - t45 * t133) * t185, -t209 * t125 - t285 * t92 - t210 * t138 - t27 * t130 + t65 * t112 + t89 * t42 - t45 * t270 - g(1) * t117 - g(2) * t119 + (t19 * t189 - t251 * t45) * t139, -t214 * t23 - t86 * t9, -t198 * t86 + t214 * t24 - t23 * t55 - t85 * t9, -t218 + t294, t216 + t288, t123 * t130 + t138 * t88 (-t184 * t7 + t188 * t6) * t123 + (-t184 * t29 + t188 * t28) * t88 + t240 * t138 + t4 * t130 + t33 * t55 - t64 * t198 + t11 * t85 + t34 * t24 - g(1) * t106 - g(2) * t108 + ((-t184 * t28 - t188 * t29) * t123 - t5 * t138) * qJD(6), -g(1) * t105 - g(2) * t107 - t11 * t86 + t13 * t138 - t5 * t130 + t34 * t23 - t33 * t214 + t64 * t9 + (-(-qJD(6) * t29 + t6) * t123 - t28 * t88 - t2 * t138) * t184 + (-(qJD(6) * t28 + t7) * t123 - t29 * t88 - t234 * t138) * t188; 0, 0, 0, t254, 0, 0, 0, 0, 0, t146, -t145, -t129 * t133 + t130 * t131 + t138 * t94 + t139 * t93, -t130 * t51 - t133 * t52 - t138 * t21 + t139 * t22 - g(3), 0, 0, 0, 0, 0, t196 - t215, t287 - t295, 0, 0, 0, 0, 0, -t216 + t288, t218 + t294; 0, 0, 0, 0, -t186 * t193 * t190, t253 * t193, t245, t244, qJDD(3), t186 * t204 + t169 - t290, -t254 * t186 + t204 * t190 (t52 - t59) * t131 + (t51 - t60) * t129 + (t179 * t93 - t181 * t94) * pkin(3), t51 * t59 - t52 * t60 + (-t290 + t179 * t22 + t181 * t21 + (-qJD(1) * t128 + t223) * t186) * pkin(3), t112 * t228 + t277 (t42 - t275) * t189 + (-t43 - t273) * t185, t125 * t228 - t272 + t278, t211 + t274, -t125 * t131, -t59 * t110 - t72 * t125 - t26 * t131 + t158 * t43 + (t60 * t125 + t205) * t185 + t306 * t189, -t59 * t112 + t286 * t125 + t27 * t131 + t158 * t42 - t185 * t306 + t205 * t189, t9 * t141 - t214 * t284, -t9 * t140 + t141 * t198 - t214 * t283 - t284 * t55, t276 - t296, t226 + t281, -t123 * t131 (-t136 * t188 - t137 * t184) * t88 - t144 * t198 + t11 * t140 - t4 * t131 + t225 * t55 - t283 * t34 + (t184 * t220 - t188 * t219) * t123 + t202 * t173 -(-t136 * t184 + t137 * t188) * t88 + t144 * t9 + t11 * t141 + t5 * t131 - t225 * t214 + t284 * t34 + (t184 * t219 + t188 * t220) * t123 - t202 * t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129 ^ 2 - t131 ^ 2, -t52 * t129 + t51 * t131 + t199 + t222, 0, 0, 0, 0, 0, t211 - t274, -t125 ^ 2 * t189 - t272 - t278, 0, 0, 0, 0, 0, t226 - t281, t276 + t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112 * t110, -t110 ^ 2 + t112 ^ 2, t42 + t275, t273 - t43, t92, -g(1) * t119 + g(2) * t117 - t45 * t112 + t27 * t125 + (t231 + t292) * t185 + t224, g(1) * t120 - g(2) * t118 + t110 * t45 + t125 * t26 + t189 * t292 - t210, -t305, t303, t302, t299, t88 -(-t14 * t184 - t279) * t123 + (-t112 * t55 - t123 * t249 + t188 * t88) * pkin(5) + t300 (-t123 * t15 - t2) * t184 + (t123 * t14 - t234) * t188 + (t112 * t214 - t123 * t248 - t184 * t88) * pkin(5) + t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t305, t303, t302, t299, t88, t5 * t123 + t300, t4 * t123 - t184 * t2 - t188 * t234 + t301;];
tau_reg  = t1;
