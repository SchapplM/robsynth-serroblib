% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR13
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR13_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:33:55
% EndTime: 2019-12-31 20:34:07
% DurationCPUTime: 4.32s
% Computational Cost: add. (4194->418), mult. (9698->578), div. (0->0), fcn. (7360->14), ass. (0->207)
t200 = cos(qJ(2));
t253 = t200 * qJD(1);
t173 = -qJD(4) + t253;
t168 = -qJD(5) + t173;
t194 = sin(qJ(5));
t198 = cos(qJ(5));
t192 = sin(pkin(9));
t196 = sin(qJ(2));
t264 = qJD(1) * t196;
t245 = t192 * t264;
t193 = cos(pkin(9));
t255 = t193 * qJD(2);
t143 = t245 - t255;
t244 = t193 * t264;
t256 = t192 * qJD(2);
t145 = t244 + t256;
t195 = sin(qJ(4));
t199 = cos(qJ(4));
t79 = t195 * t143 - t199 * t145;
t80 = t199 * t143 + t195 * t145;
t293 = t194 * t79 - t198 * t80;
t304 = t168 * t293;
t228 = pkin(2) * t196 - qJ(3) * t200;
t154 = t228 * qJD(1);
t134 = t192 * t154;
t276 = t193 * t196;
t277 = t192 * t200;
t218 = -pkin(6) * t276 - pkin(7) * t277;
t86 = qJD(1) * t218 + t134;
t303 = t193 * qJD(3) - t86;
t249 = t196 * qJDD(1);
t178 = pkin(6) * t249;
t251 = qJD(1) * qJD(2);
t241 = t200 * t251;
t127 = -qJDD(2) * pkin(2) + pkin(6) * t241 + qJDD(3) + t178;
t197 = sin(qJ(1));
t201 = cos(qJ(1));
t230 = g(1) * t201 + g(2) * t197;
t285 = g(3) * t200;
t212 = t196 * t230 - t285;
t292 = t127 - t212;
t227 = t194 * t80 + t198 * t79;
t302 = t227 * t293;
t301 = t79 * t173;
t300 = t80 * t173;
t299 = t168 * t227;
t274 = t195 * t192;
t151 = -t199 * t193 + t274;
t216 = t200 * t151;
t268 = qJD(1) * t216 - t151 * qJD(4);
t152 = t199 * t192 + t195 * t193;
t217 = t200 * t152;
t267 = -qJD(1) * t217 + t152 * qJD(4);
t298 = t227 ^ 2 - t293 ^ 2;
t257 = qJD(5) * t198;
t258 = qJD(5) * t194;
t182 = t193 * qJDD(2);
t215 = t241 + t249;
t103 = t192 * t215 - t182;
t250 = t192 * qJDD(2);
t104 = t193 * t215 + t250;
t259 = qJD(4) * t199;
t261 = qJD(4) * t195;
t30 = -t195 * t103 + t199 * t104 - t143 * t259 - t145 * t261;
t31 = -qJD(4) * t79 + t199 * t103 + t195 * t104;
t6 = -t194 * t31 + t198 * t30 - t80 * t257 + t258 * t79;
t297 = t6 + t304;
t189 = pkin(9) + qJ(4);
t187 = qJ(5) + t189;
t176 = cos(t187);
t175 = sin(t187);
t272 = t201 * t175;
t273 = t197 * t200;
t106 = -t176 * t273 + t272;
t271 = t201 * t176;
t108 = t197 * t175 + t200 * t271;
t247 = pkin(3) * t253;
t223 = t200 * pkin(2) + t196 * qJ(3) + pkin(1);
t136 = t223 * qJD(1);
t180 = pkin(6) * t253;
t163 = qJD(2) * qJ(3) + t180;
t88 = -t193 * t136 - t192 * t163;
t48 = -t145 * pkin(7) - t247 + t88;
t89 = -t192 * t136 + t193 * t163;
t51 = -t143 * pkin(7) + t89;
t19 = t195 * t48 + t199 * t51;
t13 = -t80 * pkin(8) + t19;
t11 = t13 * t258;
t286 = g(3) * t196;
t156 = -qJD(2) * pkin(2) + pkin(6) * t264 + qJD(3);
t98 = t143 * pkin(3) + t156;
t44 = t80 * pkin(4) + t98;
t296 = g(1) * t108 - g(2) * t106 + t176 * t286 - t293 * t44 + t11;
t105 = t175 * t273 + t271;
t107 = t197 * t176 - t200 * t272;
t186 = t200 * qJDD(1);
t240 = t196 * t251;
t214 = t240 - t186;
t150 = qJDD(4) + t214;
t117 = -pkin(6) * t214 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t129 = qJD(2) * t228 - t196 * qJD(3);
t78 = qJD(1) * t129 - qJDD(1) * t223;
t40 = -t192 * t117 + t193 * t78;
t23 = pkin(3) * t214 - t104 * pkin(7) + t40;
t41 = t193 * t117 + t192 * t78;
t32 = -t103 * pkin(7) + t41;
t237 = -t195 * t32 + t199 * t23;
t206 = -qJD(4) * t19 + t237;
t2 = t150 * pkin(4) - t30 * pkin(8) + t206;
t221 = t195 * t23 + t199 * t32 + t48 * t259 - t261 * t51;
t3 = -t31 * pkin(8) + t221;
t246 = -t194 * t3 + t198 * t2;
t18 = -t195 * t51 + t199 * t48;
t12 = pkin(8) * t79 + t18;
t10 = -t173 * pkin(4) + t12;
t279 = t198 * t13;
t5 = t194 * t10 + t279;
t295 = -g(1) * t107 + g(2) * t105 - qJD(5) * t5 + t175 * t286 + t44 * t227 + t246;
t205 = qJD(5) * t227 - t194 * t30 - t198 * t31;
t294 = t205 + t299;
t284 = pkin(7) + qJ(3);
t161 = t284 * t192;
t162 = t284 * t193;
t266 = -t195 * t161 + t199 * t162;
t225 = t192 * qJD(3) + qJD(4) * t162;
t275 = t193 * t200;
t224 = pkin(3) * t196 - pkin(7) * t275;
t99 = pkin(6) * t245 + t193 * t154;
t68 = qJD(1) * t224 + t99;
t291 = -t161 * t259 + t303 * t199 + (-t225 - t68) * t195;
t229 = g(1) * t197 - g(2) * t201;
t290 = pkin(6) * t143;
t289 = pkin(6) * t145;
t84 = t198 * t151 + t194 * t152;
t283 = -qJD(5) * t84 - t194 * t267 + t198 * t268;
t85 = -t194 * t151 + t198 * t152;
t282 = qJD(5) * t85 + t194 * t268 + t198 * t267;
t142 = t193 * t223;
t87 = -pkin(7) * t276 - t142 + (-pkin(6) * t192 - pkin(3)) * t200;
t110 = pkin(6) * t275 - t192 * t223;
t278 = t192 * t196;
t94 = -pkin(7) * t278 + t110;
t280 = t195 * t87 + t199 * t94;
t184 = sin(t189);
t270 = t201 * t184;
t185 = cos(t189);
t269 = t201 * t185;
t263 = qJD(2) * t196;
t248 = pkin(6) * t263;
t92 = t193 * t129 + t192 * t248;
t138 = t192 * t247 + t180;
t262 = qJD(2) * t200;
t139 = t200 * pkin(3) * t256 + pkin(6) * t262;
t155 = pkin(3) * t278 + t196 * pkin(6);
t190 = t196 ^ 2;
t265 = -t200 ^ 2 + t190;
t260 = qJD(4) * t196;
t252 = qJD(3) - t156;
t177 = -t193 * pkin(3) - pkin(2);
t243 = pkin(4) * t267 - t138;
t242 = qJ(3) * t186;
t239 = qJD(5) * t10 + t3;
t57 = qJD(2) * t224 + t92;
t115 = t192 * t129;
t70 = qJD(2) * t218 + t115;
t236 = -t195 * t70 + t199 * t57;
t235 = -t195 * t94 + t199 * t87;
t233 = -t199 * t161 - t195 * t162;
t60 = t199 * t68;
t62 = -t151 * pkin(8) + t266;
t232 = pkin(4) * t264 + pkin(8) * t268 + t152 * qJD(3) + t266 * qJD(4) + qJD(5) * t62 - t195 * t86 + t60;
t61 = -t152 * pkin(8) + t233;
t231 = -pkin(8) * t267 + qJD(5) * t61 + t291;
t125 = t152 * t196;
t126 = t151 * t196;
t55 = t198 * t125 - t194 * t126;
t56 = -t194 * t125 - t198 * t126;
t222 = -0.2e1 * pkin(1) * t251 - pkin(6) * qJDD(2);
t220 = t195 * t57 + t199 * t70 + t87 * t259 - t261 * t94;
t203 = qJD(1) ^ 2;
t213 = pkin(1) * t203 + t230;
t58 = t103 * pkin(3) + t127;
t210 = -t200 * t230 - t286;
t202 = qJD(2) ^ 2;
t209 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t202 + t229;
t137 = qJDD(5) + t150;
t121 = t197 * t184 + t200 * t269;
t120 = t197 * t185 - t200 * t270;
t119 = -t185 * t273 + t270;
t118 = t184 * t273 + t269;
t114 = t151 * pkin(4) + t177;
t109 = -pkin(6) * t277 - t142;
t100 = -pkin(6) * t244 + t134;
t93 = -t193 * t248 + t115;
t91 = t125 * pkin(4) + t155;
t64 = qJD(2) * t217 + t259 * t276 - t260 * t274;
t63 = -qJD(2) * t216 - t152 * t260;
t45 = t64 * pkin(4) + t139;
t29 = -t125 * pkin(8) + t280;
t24 = -t200 * pkin(4) + t126 * pkin(8) + t235;
t16 = t31 * pkin(4) + t58;
t15 = qJD(5) * t56 + t194 * t63 + t198 * t64;
t14 = -qJD(5) * t55 - t194 * t64 + t198 * t63;
t9 = -t64 * pkin(8) + t220;
t8 = pkin(4) * t263 - t63 * pkin(8) - qJD(4) * t280 + t236;
t4 = t198 * t10 - t194 * t13;
t1 = [qJDD(1), t229, t230, t190 * qJDD(1) + 0.2e1 * t200 * t240, 0.2e1 * t186 * t196 - 0.2e1 * t251 * t265, qJDD(2) * t196 + t202 * t200, qJDD(2) * t200 - t202 * t196, 0, t196 * t222 + t200 * t209, -t196 * t209 + t200 * t222, -t230 * t192 + (pkin(6) * t103 + t127 * t192 + (qJD(1) * t109 + t88) * qJD(2)) * t196 + (-t92 * qJD(1) - t109 * qJDD(1) - t40 + t229 * t193 + (t156 * t192 + t290) * qJD(2)) * t200, -t230 * t193 + (pkin(6) * t104 + t127 * t193 + (-qJD(1) * t110 - t89) * qJD(2)) * t196 + (t93 * qJD(1) + t110 * qJDD(1) + t41 - t229 * t192 + (t156 * t193 + t289) * qJD(2)) * t200, -t110 * t103 - t109 * t104 - t93 * t143 - t92 * t145 + (-t192 * t89 - t193 * t88) * t262 + (-t192 * t41 - t193 * t40 + t229) * t196, t40 * t109 + t41 * t110 + t88 * t92 + t89 * t93 + (t127 * t196 + t156 * t262 - t230) * pkin(6) + t229 * t223, -t30 * t126 - t63 * t79, -t30 * t125 + t126 * t31 - t63 * t80 + t64 * t79, -t126 * t150 - t63 * t173 - t30 * t200 - t263 * t79, -t125 * t150 + t64 * t173 + t31 * t200 - t263 * t80, -t150 * t200 - t173 * t263, -t236 * t173 + t235 * t150 - t237 * t200 + t18 * t263 + t139 * t80 + t155 * t31 + t58 * t125 + t98 * t64 - g(1) * t119 - g(2) * t121 + (t173 * t280 + t19 * t200) * qJD(4), -g(1) * t118 - g(2) * t120 - t58 * t126 - t139 * t79 - t150 * t280 + t155 * t30 + t173 * t220 - t19 * t263 + t200 * t221 + t98 * t63, -t14 * t227 + t6 * t56, t14 * t293 + t15 * t227 + t205 * t56 - t6 * t55, t56 * t137 - t14 * t168 - t6 * t200 - t227 * t263, -t55 * t137 + t15 * t168 - t200 * t205 + t263 * t293, -t137 * t200 - t168 * t263, -(-t194 * t9 + t198 * t8) * t168 + (-t194 * t29 + t198 * t24) * t137 - t246 * t200 + t4 * t263 - t45 * t293 - t91 * t205 + t16 * t55 + t44 * t15 - g(1) * t106 - g(2) * t108 + (-(-t194 * t24 - t198 * t29) * t168 + t5 * t200) * qJD(5), -t5 * t263 - g(1) * t105 - g(2) * t107 - t11 * t200 + t44 * t14 + t16 * t56 - t45 * t227 + t91 * t6 + ((-qJD(5) * t29 + t8) * t168 - t24 * t137 + t2 * t200) * t194 + ((qJD(5) * t24 + t9) * t168 - t29 * t137 + t239 * t200) * t198; 0, 0, 0, -t196 * t203 * t200, t265 * t203, t249, t186, qJDD(2), t196 * t213 - t178 - t285, t286 + (-pkin(6) * qJDD(1) + t213) * t200, t192 * t242 - pkin(2) * t103 - t292 * t193 + ((-qJ(3) * t256 - t88) * t196 + (t192 * t252 - t290 + t99) * t200) * qJD(1), t193 * t242 - pkin(2) * t104 + t292 * t192 + ((-qJ(3) * t255 + t89) * t196 + (t193 * t252 - t100 - t289) * t200) * qJD(1), t100 * t143 + t99 * t145 + (-qJ(3) * t103 - qJD(3) * t143 + t253 * t88 + t41) * t193 + (qJ(3) * t104 + qJD(3) * t145 + t253 * t89 - t40) * t192 + t210, -t156 * t180 - t89 * t100 - t88 * t99 + (-t88 * t192 + t89 * t193) * qJD(3) - t292 * pkin(2) + (-t40 * t192 + t41 * t193 + t210) * qJ(3), t30 * t152 - t268 * t79, -t30 * t151 - t152 * t31 + t267 * t79 - t268 * t80, t152 * t150 - t173 * t268 + t264 * t79, -t151 * t150 + t173 * t267 + t264 * t80, t173 * t264, t233 * t150 + t177 * t31 + t58 * t151 - t18 * t264 - t138 * t80 + t267 * t98 + (t60 + t225 * t199 + (-qJD(4) * t161 + t303) * t195) * t173 + t212 * t185, t138 * t79 - t266 * t150 + t58 * t152 + t291 * t173 + t177 * t30 - t184 * t212 + t19 * t264 + t268 * t98, -t227 * t283 + t6 * t85, t205 * t85 + t227 * t282 + t283 * t293 - t6 * t84, t85 * t137 - t168 * t283 + t227 * t264, -t84 * t137 + t168 * t282 - t264 * t293, t168 * t264, (-t194 * t62 + t198 * t61) * t137 - t114 * t205 + t16 * t84 - t4 * t264 + t282 * t44 - t243 * t293 + (t194 * t231 + t198 * t232) * t168 + t212 * t176, -(t194 * t61 + t198 * t62) * t137 + t114 * t6 + t16 * t85 + t5 * t264 + t283 * t44 - t243 * t227 + (-t194 * t232 + t198 * t231) * t168 - t212 * t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192 * t249 - t182 + (-t145 + t256) * t253, t193 * t249 + t250 + (t143 + t255) * t253, -t143 ^ 2 - t145 ^ 2, t89 * t143 + t88 * t145 + t292, 0, 0, 0, 0, 0, t31 + t301, t30 + t300, 0, 0, 0, 0, 0, -t205 + t299, t6 - t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 * t80, t79 ^ 2 - t80 ^ 2, t30 - t300, -t31 + t301, t150, -g(1) * t120 + g(2) * t118 - t19 * t173 + t184 * t286 + t79 * t98 + t206, g(1) * t121 - g(2) * t119 - t18 * t173 + t185 * t286 + t98 * t80 - t221, t302, t298, t297, t294, t137, (-t194 * t12 - t279) * t168 + (t198 * t137 + t168 * t258 - t293 * t79) * pkin(4) + t295, (t13 * t168 - t2) * t194 + (-t12 * t168 - t239) * t198 + (-t194 * t137 + t168 * t257 - t227 * t79) * pkin(4) + t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t298, t297, t294, t137, -t5 * t168 + t295, -t4 * t168 - t194 * t2 - t198 * t239 + t296;];
tau_reg = t1;
