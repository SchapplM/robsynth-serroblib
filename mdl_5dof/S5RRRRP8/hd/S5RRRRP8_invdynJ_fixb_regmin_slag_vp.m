% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:11
% EndTime: 2019-12-31 22:02:21
% DurationCPUTime: 3.63s
% Computational Cost: add. (3829->397), mult. (8656->539), div. (0->0), fcn. (5964->10), ass. (0->207)
t176 = sin(qJ(3));
t290 = pkin(7) + pkin(8);
t232 = qJD(3) * t290;
t177 = sin(qJ(2));
t179 = cos(qJ(3));
t259 = t177 * t179;
t180 = cos(qJ(2));
t260 = t176 * t180;
t210 = pkin(2) * t177 - pkin(7) * t180;
t120 = t210 * qJD(1);
t98 = t176 * t120;
t312 = -t176 * t232 - t98 - (-pkin(6) * t259 - pkin(8) * t260) * qJD(1);
t257 = t179 * t180;
t203 = pkin(3) * t177 - pkin(8) * t257;
t250 = qJD(1) * t177;
t225 = t176 * t250;
t253 = pkin(6) * t225 + t179 * t120;
t311 = qJD(1) * t203 + t179 * t232 + t253;
t244 = qJD(2) * t179;
t112 = -t225 + t244;
t246 = qJD(2) * t176;
t113 = t179 * t250 + t246;
t175 = sin(qJ(4));
t286 = cos(qJ(4));
t199 = t175 * t112 + t286 * t113;
t291 = t199 ^ 2;
t59 = -t286 * t112 + t113 * t175;
t57 = t59 ^ 2;
t310 = -t57 + t291;
t262 = t175 * t176;
t198 = t286 * t179 - t262;
t249 = qJD(1) * t180;
t296 = qJD(3) + qJD(4);
t223 = t286 * qJD(4);
t297 = t286 * qJD(3) + t223;
t272 = -t297 * t179 + t198 * t249 + t296 * t262;
t115 = t175 * t179 + t286 * t176;
t72 = t296 * t115;
t271 = -t115 * t249 + t72;
t309 = t199 * t59;
t308 = t59 * qJ(5);
t243 = qJD(2) * t180;
t226 = t176 * t243;
t240 = qJD(3) * t179;
t307 = t177 * t240 + t226;
t237 = qJD(1) * qJD(2);
t221 = t180 * t237;
t235 = t177 * qJDD(1);
t306 = qJD(2) * qJD(3) + t221 + t235;
t152 = -qJD(3) + t249;
t140 = -qJD(4) + t152;
t241 = qJD(3) * t177;
t220 = qJD(1) * t241;
t212 = t306 * t176 + t179 * t220;
t192 = t179 * qJDD(2) - t212;
t239 = qJD(4) * t175;
t52 = (qJDD(2) - t220) * t176 + t306 * t179;
t16 = -t112 * t223 + t113 * t239 - t175 * t192 - t286 * t52;
t305 = -t140 * t59 - t16;
t174 = qJ(3) + qJ(4);
t168 = cos(t174);
t166 = t180 * qJDD(1);
t298 = -t177 * t237 + t166;
t109 = qJDD(3) - t298;
t123 = t210 * qJD(2);
t127 = -pkin(2) * t180 - pkin(7) * t177 - pkin(1);
t73 = qJD(1) * t123 + qJDD(1) * t127;
t65 = t179 * t73;
t105 = t127 * qJD(1);
t163 = pkin(6) * t249;
t132 = qJD(2) * pkin(7) + t163;
t69 = t105 * t176 + t132 * t179;
t89 = t298 * pkin(6) + qJDD(2) * pkin(7);
t11 = pkin(3) * t109 - pkin(8) * t52 - qJD(3) * t69 - t176 * t89 + t65;
t242 = qJD(3) * t176;
t195 = t105 * t240 - t132 * t242 + t176 * t73 + t179 * t89;
t15 = pkin(8) * t192 + t195;
t68 = t179 * t105 - t132 * t176;
t42 = -pkin(8) * t113 + t68;
t37 = -pkin(3) * t152 + t42;
t43 = pkin(8) * t112 + t69;
t219 = -t175 * t11 - t286 * t15 - t37 * t223 + t43 * t239;
t281 = g(3) * t177;
t131 = -qJD(2) * pkin(2) + pkin(6) * t250;
t75 = -pkin(3) * t112 + t131;
t167 = sin(t174);
t181 = cos(qJ(1));
t178 = sin(qJ(1));
t258 = t178 * t180;
t82 = t167 * t181 - t168 * t258;
t256 = t180 * t181;
t84 = t167 * t178 + t168 * t256;
t304 = g(1) * t84 - g(2) * t82 + t168 * t281 + t59 * t75 + t219;
t33 = pkin(4) * t59 + qJD(5) + t75;
t303 = t199 * t33;
t302 = qJ(5) * t199;
t285 = pkin(3) * t176;
t301 = pkin(3) * t242 - t249 * t285 - t163;
t133 = t290 * t176;
t134 = t290 * t179;
t254 = -t175 * t133 + t286 * t134;
t300 = t254 * qJD(4) + t312 * t175 + t311 * t286;
t299 = -t133 * t223 - t134 * t239 - t311 * t175 + t312 * t286;
t209 = g(1) * t181 + g(2) * t178;
t280 = g(3) * t180;
t190 = -t177 * t209 + t280;
t161 = pkin(6) * t235;
t268 = qJDD(2) * pkin(2);
t90 = pkin(6) * t221 + t161 - t268;
t295 = -qJD(3) * pkin(7) * t152 + t190 + t90;
t263 = t168 * t181;
t81 = t167 * t258 + t263;
t264 = t168 * t178;
t83 = -t167 * t256 + t264;
t294 = -g(1) * t83 + g(2) * t81 + t167 * t281;
t41 = t286 * t43;
t19 = t175 * t37 + t41;
t188 = -qJD(4) * t19 + t286 * t11 - t175 * t15;
t293 = -t75 * t199 + t188 + t294;
t17 = qJD(4) * t199 + t175 * t52 - t286 * t192;
t292 = -t140 * t199 - t17;
t39 = t175 * t43;
t18 = t286 * t37 - t39;
t9 = t18 - t302;
t6 = -pkin(4) * t140 + t9;
t289 = -t9 + t6;
t284 = pkin(6) * t176;
t170 = t179 * pkin(3);
t279 = pkin(2) + t170;
t125 = pkin(4) * t167 + t285;
t278 = pkin(6) + t125;
t277 = -t271 * qJ(5) + qJD(5) * t198 + t299;
t276 = -pkin(4) * t250 + t272 * qJ(5) - t115 * qJD(5) - t300;
t275 = t286 * t42 - t39;
t111 = t179 * t127;
t67 = -pkin(8) * t259 + t111 + (-pkin(3) - t284) * t180;
t154 = pkin(6) * t257;
t252 = t176 * t127 + t154;
t261 = t176 * t177;
t74 = -pkin(8) * t261 + t252;
t273 = t175 * t67 + t286 * t74;
t270 = t52 * t176;
t269 = t176 * t123 + t127 * t240;
t267 = t112 * t152;
t266 = t113 * t152;
t265 = t113 * t179;
t245 = qJD(2) * t177;
t255 = t179 * t123 + t245 * t284;
t126 = pkin(4) * t168 + t170;
t124 = pkin(3) * t261 + t177 * pkin(6);
t172 = t177 ^ 2;
t251 = -t180 ^ 2 + t172;
t248 = qJD(2) * t112;
t247 = qJD(2) * t113;
t238 = t131 * qJD(3);
t76 = t307 * pkin(3) + pkin(6) * t243;
t231 = t152 * t244;
t230 = t152 * t242;
t229 = t152 * t240;
t228 = t176 * t241;
t218 = -t175 * t42 - t41;
t215 = -t175 * t74 + t286 * t67;
t214 = -qJD(3) * t105 - t89;
t213 = -t286 * t133 - t134 * t175;
t211 = t286 * t243;
t208 = g(1) * t178 - g(2) * t181;
t207 = t132 * t240 - t65;
t206 = -pkin(7) * t109 + t238;
t119 = pkin(2) + t126;
t171 = -qJ(5) - t290;
t205 = t119 * t180 - t171 * t177;
t202 = pkin(1) + t205;
t201 = -0.2e1 * pkin(1) * t237 - pkin(6) * qJDD(2);
t28 = t203 * qJD(2) + (-t154 + (pkin(8) * t177 - t127) * t176) * qJD(3) + t255;
t30 = -t307 * pkin(8) + (-t177 * t244 - t180 * t242) * pkin(6) + t269;
t200 = t175 * t28 + t67 * t223 - t239 * t74 + t286 * t30;
t197 = t109 * t176 - t229;
t196 = t179 * t109 + t230;
t183 = qJD(1) ^ 2;
t193 = pkin(1) * t183 + t209;
t182 = qJD(2) ^ 2;
t189 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t182 + t208;
t187 = -qJD(4) * t273 - t175 * t30 + t286 * t28;
t36 = -pkin(3) * t192 + t90;
t5 = t17 * pkin(4) + qJDD(5) + t36;
t159 = t286 * pkin(3) + pkin(4);
t104 = qJDD(4) + t109;
t96 = t176 * t178 + t179 * t256;
t95 = -t176 * t256 + t178 * t179;
t94 = t176 * t181 - t178 * t257;
t93 = t176 * t258 + t179 * t181;
t87 = t198 * t177;
t86 = t115 * t177;
t47 = qJ(5) * t198 + t254;
t46 = -qJ(5) * t115 + t213;
t32 = t176 * t211 - t175 * t228 - t239 * t261 + (t175 * t243 + t297 * t177) * t179;
t31 = t175 * t226 + t177 * t72 - t179 * t211;
t25 = -qJ(5) * t86 + t273;
t24 = -pkin(4) * t180 - qJ(5) * t87 + t215;
t13 = t275 - t302;
t12 = t218 + t308;
t10 = t19 - t308;
t4 = -qJ(5) * t32 - qJD(5) * t86 + t200;
t3 = pkin(4) * t245 + t31 * qJ(5) - t87 * qJD(5) + t187;
t2 = -qJ(5) * t17 - qJD(5) * t59 - t219;
t1 = t104 * pkin(4) + t16 * qJ(5) - qJD(5) * t199 + t188;
t7 = [qJDD(1), t208, t209, qJDD(1) * t172 + 0.2e1 * t177 * t221, 0.2e1 * t166 * t177 - 0.2e1 * t237 * t251, qJDD(2) * t177 + t180 * t182, qJDD(2) * t180 - t177 * t182, 0, t177 * t201 + t180 * t189, -t177 * t189 + t180 * t201, t52 * t259 + (t179 * t243 - t228) * t113, (t112 * t179 - t113 * t176) * t243 + (t179 * t192 - t270 + (-t112 * t176 - t265) * qJD(3)) * t177, (-t52 - t231) * t180 + (t196 + t247) * t177, (t152 * t246 - t192) * t180 + (-t197 + t248) * t177, -t109 * t180 - t152 * t245, -(-t127 * t242 + t255) * t152 + t111 * t109 - g(1) * t94 - g(2) * t96 + ((t229 - t248) * pkin(6) + (-pkin(6) * t109 + qJD(2) * t131 - t214) * t176 + t207) * t180 + (-pkin(6) * t192 + t68 * qJD(2) + t90 * t176 + t179 * t238) * t177, t269 * t152 - t252 * t109 - g(1) * t93 - g(2) * t95 + (t131 * t244 + (-t230 + t247) * pkin(6) + t195) * t180 + (-t176 * t238 - t69 * qJD(2) + t90 * t179 + (t52 - t231) * pkin(6)) * t177, -t16 * t87 - t199 * t31, t16 * t86 - t17 * t87 - t199 * t32 + t31 * t59, t104 * t87 + t140 * t31 + t16 * t180 + t199 * t245, -t104 * t86 + t140 * t32 + t17 * t180 - t245 * t59, -t104 * t180 - t140 * t245, -g(1) * t82 - g(2) * t84 + t104 * t215 + t124 * t17 - t140 * t187 + t18 * t245 - t180 * t188 + t75 * t32 + t36 * t86 + t76 * t59, -g(1) * t81 - g(2) * t83 - t104 * t273 - t124 * t16 + t140 * t200 - t180 * t219 - t19 * t245 + t199 * t76 - t75 * t31 + t36 * t87, -t1 * t87 - t10 * t32 + t16 * t24 - t17 * t25 + t177 * t208 - t199 * t3 - t2 * t86 + t31 * t6 - t4 * t59, t2 * t25 + t10 * t4 + t1 * t24 + t6 * t3 + t5 * (pkin(4) * t86 + t124) + t33 * (pkin(4) * t32 + t76) + (-g(1) * t278 - g(2) * t202) * t181 + (g(1) * t202 - g(2) * t278) * t178; 0, 0, 0, -t177 * t183 * t180, t251 * t183, t235, t166, qJDD(2), t177 * t193 - t161 - t280, t281 + (-pkin(6) * qJDD(1) + t193) * t180, -t152 * t265 + t270, (t52 - t267) * t179 + (t192 + t266) * t176, (-t113 * t177 + t152 * t257) * qJD(1) + t197, (-t112 * t177 - t152 * t260) * qJD(1) + t196, t152 * t250, -pkin(2) * t212 + t253 * t152 + t206 * t176 + (-t68 * t177 + (pkin(6) * t112 - t131 * t176) * t180) * qJD(1) + (t268 - t295) * t179, -pkin(2) * t52 - t98 * t152 + t206 * t179 + (-t131 * t257 + t69 * t177 + (-t113 * t180 + t152 * t259) * pkin(6)) * qJD(1) + t295 * t176, -t115 * t16 - t199 * t272, -t115 * t17 - t16 * t198 - t199 * t271 + t272 * t59, t104 * t115 + t140 * t272 - t199 * t250, t104 * t198 + t140 * t271 + t250 * t59, t140 * t250, t104 * t213 - t36 * t198 - t279 * t17 - t168 * t280 - t18 * t250 + t271 * t75 + t301 * t59 + (g(1) * t263 + g(2) * t264) * t177 + t300 * t140, -t254 * t104 + t36 * t115 + t299 * t140 + t16 * t279 + t190 * t167 + t19 * t250 + t301 * t199 - t272 * t75, -t1 * t115 - t10 * t271 + t16 * t46 - t17 * t47 - t180 * t209 + t198 * t2 - t199 * t276 + t272 * t6 - t277 * t59 - t281, t2 * t47 + t1 * t46 + t5 * (-pkin(4) * t198 - t279) - g(3) * t205 + t276 * t6 + (pkin(4) * t271 + t301) * t33 + t277 * t10 + t209 * (t119 * t177 + t171 * t180); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113 * t112, -t112 ^ 2 + t113 ^ 2, t52 + t267, t192 - t266, t109, -g(1) * t95 + g(2) * t93 - t113 * t131 - t152 * t69 + (t214 + t281) * t176 - t207, g(1) * t96 - g(2) * t94 + g(3) * t259 - t112 * t131 - t152 * t68 - t195, t309, t310, t305, t292, t104, t218 * t140 + (t286 * t104 - t113 * t59 + t140 * t239) * pkin(3) + t293, -t275 * t140 + (-t175 * t104 - t113 * t199 + t140 * t223) * pkin(3) + t304, t10 * t199 + t12 * t199 + t13 * t59 + t159 * t16 - t6 * t59 + (-t17 * t175 + (t175 * t199 - t286 * t59) * qJD(4)) * pkin(3), t1 * t159 - t10 * t13 - t6 * t12 - pkin(4) * t303 - g(1) * (-t125 * t256 + t126 * t178) - g(2) * (-t125 * t258 - t126 * t181) + t125 * t281 + (-t33 * t113 + t2 * t175 + (t286 * t10 - t175 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t309, t310, t305, t292, t104, -t19 * t140 + t293, -t140 * t18 + t304, pkin(4) * t16 - t289 * t59, t289 * t10 + (t1 + t294 - t303) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 - t291, t10 * t59 + t199 * t6 + t190 + t5;];
tau_reg = t7;
