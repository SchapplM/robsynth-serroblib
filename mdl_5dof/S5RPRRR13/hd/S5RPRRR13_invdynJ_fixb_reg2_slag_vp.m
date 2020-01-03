% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR13
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR13_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:33
% EndTime: 2019-12-31 19:15:43
% DurationCPUTime: 4.84s
% Computational Cost: add. (4918->465), mult. (9713->612), div. (0->0), fcn. (6112->10), ass. (0->232)
t160 = sin(qJ(3));
t248 = qJD(1) * t160;
t134 = qJD(4) + t248;
t164 = cos(qJ(3));
t295 = g(3) * t160;
t165 = cos(qJ(1));
t152 = g(2) * t165;
t161 = sin(qJ(1));
t153 = g(1) * t161;
t306 = t153 - t152;
t180 = t306 * t164 - t295;
t167 = -pkin(1) - pkin(6);
t126 = t167 * qJDD(1) + qJDD(2);
t129 = t167 * qJD(1) + qJD(2);
t244 = qJD(3) * t160;
t277 = qJDD(3) * pkin(3);
t59 = -t164 * t126 + t129 * t244 - t277;
t177 = -t59 - t180;
t316 = qJD(4) * pkin(7) * t134 - t177;
t159 = sin(qJ(4));
t247 = qJD(1) * t164;
t220 = t159 * t247;
t163 = cos(qJ(4));
t237 = t163 * qJD(3);
t103 = t220 - t237;
t245 = qJD(3) * t159;
t105 = t163 * t247 + t245;
t158 = sin(qJ(5));
t162 = cos(qJ(5));
t192 = t103 * t158 - t162 * t105;
t48 = t162 * t103 + t105 * t158;
t293 = t48 * t192;
t221 = t159 * t248;
t301 = pkin(8) + pkin(7);
t222 = qJD(4) * t301;
t199 = pkin(3) * t164 + pkin(7) * t160;
t109 = t199 * qJD(1);
t260 = t163 * t164;
t56 = t159 * t109 + t129 * t260;
t315 = pkin(8) * t221 + t159 * t222 + t56;
t265 = t160 * t163;
t230 = pkin(8) * t265;
t268 = t159 * t164;
t55 = t163 * t109 - t129 * t268;
t314 = -t163 * t222 - (pkin(4) * t164 + t230) * qJD(1) - t55;
t232 = t164 * qJDD(1);
t313 = qJD(3) * qJD(4) + t232;
t312 = t192 ^ 2 - t48 ^ 2;
t240 = qJD(4) * t164;
t216 = t163 * t240;
t225 = qJD(1) * t216 + t313 * t159;
t187 = t163 * qJDD(3) - t225;
t236 = qJD(1) * qJD(3);
t212 = t160 * t236;
t174 = t159 * t212 + t187;
t238 = qJD(5) * t162;
t239 = qJD(5) * t158;
t218 = t160 * t237;
t182 = t159 * t240 + t218;
t42 = qJD(1) * t182 - t159 * qJDD(3) - t313 * t163;
t10 = t103 * t238 + t105 * t239 - t158 * t174 + t162 * t42;
t128 = qJD(5) + t134;
t311 = t128 * t48 - t10;
t198 = pkin(3) * t160 - pkin(7) * t164;
t113 = qJ(2) + t198;
t86 = t113 * qJD(1);
t112 = t160 * t129;
t93 = qJD(3) * pkin(7) + t112;
t40 = -t159 * t93 + t163 * t86;
t33 = -pkin(8) * t105 + t40;
t27 = pkin(4) * t134 + t33;
t41 = t159 * t86 + t163 * t93;
t34 = -pkin(8) * t103 + t41;
t211 = t164 * t236;
t233 = t160 * qJDD(1);
t100 = qJDD(4) + t211 + t233;
t101 = qJD(3) * t199 + qJD(2);
t46 = qJD(1) * t101 + qJDD(1) * t113;
t44 = t163 * t46;
t243 = qJD(3) * t164;
t60 = qJDD(3) * pkin(7) + t126 * t160 + t129 * t243;
t15 = -qJD(4) * t41 - t159 * t60 + t44;
t6 = pkin(4) * t100 + pkin(8) * t42 + t15;
t241 = qJD(4) * t163;
t242 = qJD(4) * t159;
t14 = t159 * t46 + t163 * t60 + t86 * t241 - t242 * t93;
t7 = pkin(8) * t174 + t14;
t1 = (qJD(5) * t27 + t7) * t162 + t158 * t6 - t34 * t239;
t157 = qJ(4) + qJ(5);
t147 = cos(t157);
t294 = g(3) * t164;
t286 = qJD(3) * pkin(3);
t94 = -t129 * t164 - t286;
t58 = pkin(4) * t103 + t94;
t146 = sin(t157);
t266 = t160 * t161;
t71 = t146 * t165 + t147 * t266;
t264 = t160 * t165;
t73 = -t146 * t161 + t147 * t264;
t310 = g(1) * t71 - g(2) * t73 + t147 * t294 + t48 * t58 - t1;
t282 = t162 * t34;
t9 = t158 * t27 + t282;
t2 = -qJD(5) * t9 - t158 * t7 + t162 * t6;
t70 = -t146 * t266 + t147 * t165;
t72 = t146 * t264 + t147 * t161;
t309 = -g(1) * t70 - g(2) * t72 + t146 * t294 + t192 * t58 + t2;
t175 = qJD(5) * t192 + t158 * t42 + t162 * t174;
t308 = -t128 * t192 + t175;
t307 = -t134 * t40 + t14;
t261 = t162 * t163;
t270 = t158 * t159;
t107 = -t261 + t270;
t77 = t107 * t160;
t108 = t158 * t163 + t159 * t162;
t75 = t108 * t160;
t155 = t160 ^ 2;
t156 = t164 ^ 2;
t250 = t155 + t156;
t206 = t250 * t126;
t305 = -pkin(4) * t159 + t167;
t259 = t163 * t165;
t87 = -t159 * t266 + t259;
t262 = t161 * t163;
t89 = t159 * t264 + t262;
t304 = -g(1) * t87 - g(2) * t89;
t231 = qJD(4) + qJD(5);
t203 = qJD(4) * t160 + qJD(1);
t215 = t159 * t243;
t303 = t163 * t203 + t215;
t302 = 0.2e1 * qJ(2);
t297 = pkin(7) * t100;
t296 = g(2) * t161;
t78 = t107 * t164;
t85 = t108 * qJD(1);
t292 = -qJD(3) * t78 - t231 * t75 - t85;
t291 = t107 * qJD(1) - t108 * t243 + t231 * t77;
t121 = t301 * t159;
t122 = t301 * t163;
t61 = -t121 * t162 - t122 * t158;
t290 = qJD(5) * t61 + t314 * t158 - t315 * t162;
t62 = -t121 * t158 + t122 * t162;
t289 = -qJD(5) * t62 + t315 * t158 + t314 * t162;
t288 = t158 * t221 - t162 * t241 - t163 * t238 + t231 * t270 - t248 * t261;
t54 = t231 * t108;
t287 = t160 * t85 + t54;
t284 = t134 * t41;
t283 = t158 * t34;
t281 = t163 * t40;
t280 = t42 * t159;
t263 = t160 * t167;
t127 = t163 * t263;
t65 = t159 * t113 + t127;
t279 = pkin(1) * qJDD(1);
t169 = qJD(1) ^ 2;
t278 = qJ(2) * t169;
t276 = t100 * t163;
t275 = t103 * t134;
t274 = t103 * t159;
t273 = t105 * t103;
t272 = t105 * t163;
t271 = t134 * t159;
t269 = t159 * t100;
t267 = t159 * t165;
t258 = t164 * t165;
t257 = t164 * t301;
t168 = qJD(3) ^ 2;
t256 = t167 * t168;
t255 = g(1) * t258 + t164 * t296;
t228 = 0.2e1 * qJD(1) * qJD(2);
t254 = (qJDD(1) * qJ(2) + t228) * qJ(2);
t253 = t165 * pkin(1) + t161 * qJ(2);
t251 = t155 - t156;
t249 = -t168 - t169;
t246 = qJD(3) * t105;
t234 = qJDD(3) * t160;
t227 = t159 * t263;
t226 = t164 * t169 * t160;
t223 = t165 * pkin(6) + t253;
t219 = t159 * t244;
t214 = t164 * t237;
t210 = -t159 * t167 + pkin(4);
t209 = -g(2) * t264 + t294;
t205 = t250 * qJDD(1);
t204 = qJDD(2) - t279;
t202 = g(2) * t223;
t201 = t160 * t211;
t200 = -t112 + (t221 + t242) * pkin(4);
t197 = g(1) * t165 + t296;
t195 = -t306 - t278;
t99 = t163 * t113;
t45 = -pkin(8) * t260 + t160 * t210 + t99;
t57 = -pkin(8) * t268 + t65;
t19 = -t158 * t57 + t162 * t45;
t20 = t158 * t45 + t162 * t57;
t194 = -t159 * t41 - t281;
t193 = t159 * t40 - t163 * t41;
t141 = pkin(4) * t163 + pkin(3);
t190 = t141 * t160 - t257;
t189 = qJDD(1) * t302 + t228;
t188 = -t126 + t278 + t153;
t186 = t134 * t241 + t269;
t185 = -t134 * t242 + t276;
t184 = -g(1) * t266 - t209;
t183 = qJDD(3) * t167 + t236 * t302;
t181 = -t216 + t219;
t31 = -qJD(4) * t227 + t159 * t101 + t113 * t241 + t167 * t214;
t178 = t159 * t203 - t214;
t176 = t189 - t197;
t173 = qJD(4) * t194 + t14 * t163 - t15 * t159;
t172 = t174 * t163;
t149 = t165 * qJ(2);
t144 = qJDD(3) * t164;
t102 = t305 * t164;
t96 = qJDD(5) + t100;
t90 = -t159 * t161 + t160 * t259;
t88 = t160 * t262 + t267;
t81 = t163 * t101;
t76 = t108 * t164;
t64 = t99 - t227;
t63 = -pkin(4) * t181 + t167 * t244;
t32 = -t65 * qJD(4) - t167 * t215 + t81;
t26 = -t239 * t268 + (t231 * t260 - t219) * t162 - t182 * t158;
t24 = -t158 * t219 + t162 * t218 + t164 * t54;
t22 = -pkin(4) * t174 + t59;
t21 = pkin(8) * t181 + t31;
t18 = t81 + (-t127 + (pkin(8) * t164 - t113) * t159) * qJD(4) + (t164 * t210 + t230) * qJD(3);
t13 = t162 * t33 - t283;
t12 = -t158 * t33 - t282;
t8 = t162 * t27 - t283;
t4 = -qJD(5) * t20 - t158 * t21 + t162 * t18;
t3 = qJD(5) * t19 + t158 * t18 + t162 * t21;
t5 = [0, 0, 0, 0, 0, qJDD(1), t306, t197, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t306 - 0.2e1 * t279, t176, -t204 * pkin(1) - g(1) * (-pkin(1) * t161 + t149) - g(2) * t253 + t254, qJDD(1) * t156 - 0.2e1 * t201, -0.2e1 * t160 * t232 + 0.2e1 * t236 * t251, -t160 * t168 + t144, qJDD(1) * t155 + 0.2e1 * t201, -t164 * t168 - t234, 0, t183 * t164 + (t176 - t256) * t160, -t183 * t160 + (t189 - t256) * t164 - t255, -t167 * t205 - t206 + t306, -g(1) * (t167 * t161 + t149) - t202 + t167 * t206 + t254, -t105 * t182 - t42 * t260, (t103 * t163 + t105 * t159) * t244 + (t172 + t280 + (-t272 + t274) * qJD(4)) * t164, (-t134 * t237 - t42) * t160 + (t185 + t246) * t164, -t103 * t181 - t174 * t268, ((t134 + t248) * t245 + t187) * t160 + (-t103 * qJD(3) - t186) * t164, t100 * t160 + t134 * t243, -g(1) * t90 - g(2) * t88 + t64 * t100 + t32 * t134 + (t15 + (t103 * t167 - t159 * t94) * qJD(3)) * t160 + (t40 * qJD(3) + t59 * t159 + t167 * t174 + t241 * t94) * t164, g(1) * t89 - g(2) * t87 - t100 * t65 - t134 * t31 + (-t14 + (t105 * t167 - t163 * t94) * qJD(3)) * t160 + (-qJD(3) * t41 + t59 * t163 + t167 * t42 - t242 * t94) * t164, -t31 * t103 + t65 * t187 - t32 * t105 + t64 * t42 + (t281 + (qJD(1) * t65 + t41) * t159) * t244 + (qJD(4) * t193 - t14 * t159 - t15 * t163) * t164 + t255, t14 * t65 + t41 * t31 + t15 * t64 + t40 * t32 - g(1) * (pkin(3) * t264 - pkin(7) * t258 + t149) - t202 + (-t164 * t59 + t244 * t94) * t167 + (-g(1) * t167 - g(2) * t198) * t161, t10 * t78 + t192 * t24, t10 * t76 - t175 * t78 + t192 * t26 + t24 * t48, -t10 * t160 - t128 * t24 - t192 * t243 - t78 * t96, -t175 * t76 + t26 * t48, -t128 * t26 + t160 * t175 - t243 * t48 - t76 * t96, t128 * t243 + t160 * t96, -g(1) * t73 - g(2) * t71 + t102 * t175 + t128 * t4 + t160 * t2 + t19 * t96 + t22 * t76 + t243 * t8 + t26 * t58 + t48 * t63, g(1) * t72 - g(2) * t70 - t1 * t160 + t10 * t102 - t128 * t3 - t192 * t63 - t20 * t96 - t22 * t78 - t24 * t58 - t243 * t9, -t1 * t76 + t10 * t19 + t175 * t20 + t192 * t4 + t2 * t78 + t24 * t8 - t26 * t9 - t3 * t48 + t255, t1 * t20 + t9 * t3 + t2 * t19 + t8 * t4 - t22 * t102 + t58 * t63 - g(1) * (t141 * t264 - t165 * t257 + t149) - g(2) * (pkin(4) * t267 + t223) + (-g(1) * t305 - g(2) * t190) * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t169, t195 + t204, 0, 0, 0, 0, 0, 0, t160 * t249 + t144, t164 * t249 - t234, -t205, t206 + t195, 0, 0, 0, 0, 0, 0, t164 * t187 + (-t269 + (t103 + t220) * qJD(3)) * t160 - t303 * t134, t164 * t42 + (t246 - t276) * t160 + t178 * t134, (t187 * t163 + (t163 * t212 - t42) * t159) * t160 + t303 * t105 + t178 * t103, t194 * qJD(1) + (-qJD(3) * t193 - t59) * t164 + (qJD(3) * t94 + t173) * t160 - t306, 0, 0, 0, 0, 0, 0, t128 * t291 + t164 * t175 + t244 * t48 - t75 * t96, t10 * t164 - t128 * t292 - t192 * t244 + t77 * t96, -t10 * t75 - t175 * t77 + t192 * t291 - t292 * t48, -t1 * t77 - t164 * t22 - t2 * t75 + t58 * t244 + t291 * t8 + t292 * t9 - t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, -t251 * t169, t232, -t226, -t233, qJDD(3), t295 + (-t188 + t152) * t164, t160 * t188 + t209, 0, 0, t134 * t272 - t280, (-t42 - t275) * t163 + (-t105 * qJD(4) + (-t105 + t245) * t248 + t187) * t159, (-t105 * t164 + t134 * t265) * qJD(1) + t186, t103 * t271 + t172, (t103 * t164 - t160 * t271) * qJD(1) + t185, -t134 * t247, -pkin(3) * t225 - t55 * t134 - t40 * t247 - t103 * t112 + (-t297 + t94 * qJD(4) + (t94 + t286) * t248) * t159 + (t277 - t316) * t163, t41 * t247 - t105 * t112 + pkin(3) * t42 + t134 * t56 + (t134 * t94 - t297) * t163 + t316 * t159, t56 * t103 + t55 * t105 + t307 * t163 + (-t15 - t284) * t159 + (-t280 + t172 + (t272 + t274) * qJD(4)) * pkin(7) + t184, -t94 * t112 - t40 * t55 - t41 * t56 + t177 * pkin(3) + (-t160 * t306 + t173 - t294) * pkin(7), -t10 * t108 + t192 * t288, t10 * t107 + t108 * t175 + t192 * t287 + t288 * t48, t108 * t96 - t288 * t128 + t192 * t247, -t107 * t175 + t287 * t48, -t107 * t96 - t287 * t128 + t48 * t247, -t128 * t247, t107 * t22 + t128 * t289 + t141 * t175 - t147 * t180 + t200 * t48 - t247 * t8 + t287 * t58 + t61 * t96, t10 * t141 + t108 * t22 - t128 * t290 + t146 * t180 - t192 * t200 + t247 * t9 - t288 * t58 - t62 * t96, -t1 * t107 + t10 * t61 - t108 * t2 + t175 * t62 + t192 * t289 - t287 * t9 + t288 * t8 - t290 * t48 + t184, g(3) * t190 + t1 * t62 - t22 * t141 + t2 * t61 + t200 * t58 + t289 * t8 + t290 * t9 - t306 * (t141 * t164 + t160 * t301); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, -t103 ^ 2 + t105 ^ 2, t275 - t42, -t273, t105 * t134 + t174, t100, -t93 * t241 - t105 * t94 + t284 + t44 + (-qJD(4) * t86 + t294 - t60) * t159 + t304, g(1) * t88 - g(2) * t90 + g(3) * t260 + t103 * t94 - t307, 0, 0, -t293, t312, t311, t293, t308, t96, -t12 * t128 + (-t105 * t48 - t128 * t239 + t162 * t96) * pkin(4) + t309, t128 * t13 + (t105 * t192 - t128 * t238 - t158 * t96) * pkin(4) + t310, -t12 * t192 + t13 * t48 - t192 * t9 - t48 * t8 + (t10 * t162 + t175 * t158 + (-t158 * t192 - t162 * t48) * qJD(5)) * pkin(4), -t8 * t12 - t9 * t13 + (t1 * t158 + t2 * t162 - t58 * t105 + g(3) * t268 + (-t158 * t8 + t162 * t9) * qJD(5) + t304) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t293, t312, t311, t293, t308, t96, t128 * t9 + t309, t128 * t8 + t310, 0, 0;];
tau_reg = t5;
