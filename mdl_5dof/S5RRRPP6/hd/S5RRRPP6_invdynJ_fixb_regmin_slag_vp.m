% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:02
% EndTime: 2019-12-31 21:02:12
% DurationCPUTime: 3.87s
% Computational Cost: add. (4042->422), mult. (9087->551), div. (0->0), fcn. (6016->10), ass. (0->210)
t176 = sin(qJ(3));
t236 = qJD(3) * t176;
t180 = cos(qJ(2));
t242 = qJD(1) * t180;
t303 = -t176 * t242 + t236;
t142 = -qJD(3) + t242;
t169 = g(3) * t180;
t177 = sin(qJ(2));
t178 = sin(qJ(1));
t181 = cos(qJ(1));
t207 = g(1) * t181 + g(2) * t178;
t297 = -t207 * t177 + t169;
t230 = t177 * qJDD(1);
t155 = pkin(6) * t230;
t231 = qJD(1) * qJD(2);
t214 = t180 * t231;
t94 = -qJDD(2) * pkin(2) + pkin(6) * t214 + t155;
t302 = -qJD(3) * pkin(7) * t142 + t297 + t94;
t173 = sin(pkin(8));
t174 = cos(pkin(8));
t179 = cos(qJ(3));
t112 = t173 * t179 + t174 * t176;
t97 = t112 * qJD(3);
t274 = t112 * t242 - t97;
t234 = qJD(3) * t179;
t266 = t174 * t179;
t273 = t303 * t173 - t174 * t234 + t242 * t266;
t243 = qJD(1) * t177;
t223 = t176 * t243;
t233 = t179 * qJD(2);
t116 = t223 - t233;
t239 = qJD(2) * t176;
t118 = t179 * t243 + t239;
t64 = t174 * t116 + t118 * t173;
t301 = t142 * t64;
t235 = qJD(3) * t177;
t300 = -qJD(1) * t235 + qJDD(2);
t59 = ((qJD(3) + t242) * qJD(2) + t230) * t176 - t300 * t179;
t202 = -t116 * t173 + t174 * t118;
t299 = t202 ^ 2;
t279 = qJ(4) + pkin(7);
t212 = qJD(3) * t279;
t189 = -t176 * qJD(4) - t179 * t212;
t256 = t179 * t180;
t288 = pkin(3) * t177;
t200 = -qJ(4) * t256 + t288;
t208 = pkin(2) * t177 - pkin(7) * t180;
t119 = t208 * qJD(1);
t248 = pkin(6) * t223 + t179 * t119;
t57 = t200 * qJD(1) + t248;
t103 = t176 * t119;
t260 = t177 * t179;
t263 = t176 * t180;
t68 = t103 + (-pkin(6) * t260 - qJ(4) * t263) * qJD(1);
t232 = t179 * qJD(4);
t92 = -t176 * t212 + t232;
t277 = (t189 - t57) * t174 + (t68 - t92) * t173;
t163 = t180 * qJDD(1);
t296 = -t177 * t231 + t163;
t157 = pkin(6) * t242;
t295 = t303 * pkin(3) - t157;
t58 = qJD(3) * t233 + (t214 + t230) * t179 + t300 * t176;
t28 = t173 * t58 + t174 * t59;
t29 = -t173 * t59 + t174 * t58;
t294 = t28 * pkin(4) - t29 * qJ(5) - t202 * qJD(5);
t170 = qJ(3) + pkin(8);
t159 = sin(t170);
t110 = qJDD(3) - t296;
t120 = t208 * qJD(2);
t124 = -pkin(2) * t180 - pkin(7) * t177 - pkin(1);
t75 = qJD(1) * t120 + t124 * qJDD(1);
t70 = t179 * t75;
t108 = t124 * qJD(1);
t127 = qJD(2) * pkin(7) + t157;
t74 = t108 * t176 + t127 * t179;
t93 = t296 * pkin(6) + qJDD(2) * pkin(7);
t12 = t110 * pkin(3) - t58 * qJ(4) - t74 * qJD(3) - t118 * qJD(4) - t176 * t93 + t70;
t193 = t108 * t234 - t127 * t236 + t176 * t75 + t179 * t93;
t15 = -qJ(4) * t59 - qJD(4) * t116 + t193;
t3 = t174 * t12 - t173 * t15;
t226 = -qJDD(5) + t3;
t126 = -qJD(2) * pkin(2) + pkin(6) * t243;
t79 = pkin(3) * t116 + qJD(4) + t126;
t24 = pkin(4) * t64 - qJ(5) * t202 + t79;
t281 = g(3) * t177;
t160 = cos(t170);
t257 = t178 * t180;
t84 = t159 * t257 + t160 * t181;
t253 = t181 * t159;
t86 = -t178 * t160 + t180 * t253;
t293 = g(1) * t86 + g(2) * t84 + t159 * t281 - t24 * t202 + t226;
t292 = t207 * t180 + t281;
t290 = -0.2e1 * pkin(1);
t4 = t173 * t12 + t174 * t15;
t287 = pkin(4) * t110;
t286 = pkin(6) * t176;
t285 = g(1) * t178;
t282 = g(2) * t181;
t50 = -qJ(4) * t116 + t74;
t47 = t174 * t50;
t73 = t179 * t108 - t127 * t176;
t49 = -qJ(4) * t118 + t73;
t22 = t173 * t49 + t47;
t280 = t22 * t202;
t147 = pkin(6) * t256;
t238 = qJD(2) * t177;
t249 = t179 * t120 + t238 * t286;
t33 = -t177 * t232 + t200 * qJD(2) + (-t147 + (qJ(4) * t177 - t124) * t176) * qJD(3) + t249;
t250 = t176 * t120 + t124 * t234;
t38 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t260 + (-qJD(4) * t177 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t180) * t176 + t250;
t11 = t173 * t33 + t174 * t38;
t278 = -t274 * pkin(4) + t273 * qJ(5) - qJD(5) * t112 + t295;
t45 = -pkin(3) * t142 + t49;
t21 = t173 * t45 + t47;
t276 = pkin(4) * t243 - t277;
t32 = t173 * t57 + t174 * t68;
t26 = qJ(5) * t243 + t32;
t54 = t173 * t189 + t174 * t92;
t275 = t54 - t26;
t114 = t179 * t124;
t71 = -qJ(4) * t260 + t114 + (-pkin(3) - t286) * t180;
t246 = t176 * t124 + t147;
t265 = t176 * t177;
t76 = -qJ(4) * t265 + t246;
t40 = t173 * t71 + t174 * t76;
t272 = t173 * t50;
t270 = t58 * t176;
t269 = t116 * t142;
t268 = t118 * t142;
t267 = t118 * t179;
t264 = t176 * t178;
t262 = t176 * t181;
t261 = t177 * t178;
t259 = t177 * t181;
t258 = t178 * t179;
t255 = t179 * t181;
t154 = pkin(3) * t179 + pkin(2);
t128 = t180 * t154;
t254 = t180 * t181;
t23 = t174 * t49 - t272;
t251 = qJD(5) - t23;
t144 = pkin(3) * t265;
t245 = t177 * pkin(6) + t144;
t171 = t177 ^ 2;
t244 = -t180 ^ 2 + t171;
t241 = qJD(2) * t116;
t240 = qJD(2) * t118;
t237 = qJD(2) * t180;
t229 = t110 * qJ(5) + t4;
t225 = t176 * t254;
t220 = t176 * t237;
t224 = pkin(3) * t220 + pkin(6) * t237 + t234 * t288;
t221 = t142 * t233;
t219 = t180 * t233;
t218 = t142 * t236;
t217 = t142 * t234;
t216 = t279 * t176;
t211 = -qJD(3) * t108 - t93;
t149 = g(1) * t261;
t209 = -g(2) * t259 + t149;
t206 = t127 * t234 - t70;
t205 = -t64 ^ 2 - t299;
t204 = pkin(4) * t160 + qJ(5) * t159;
t10 = -t173 * t38 + t174 * t33;
t20 = t174 * t45 - t272;
t39 = -t173 * t76 + t174 * t71;
t203 = -pkin(7) * t110 + qJD(3) * t126;
t111 = t173 * t176 - t266;
t196 = -pkin(6) * qJDD(2) + t231 * t290;
t99 = t176 * t257 + t255;
t195 = t176 * t110 - t217;
t194 = t179 * t110 + t218;
t183 = qJD(1) ^ 2;
t192 = pkin(1) * t183 + t207;
t182 = qJD(2) ^ 2;
t191 = pkin(6) * t182 + qJDD(1) * t290 + t282;
t190 = t181 * pkin(1) + pkin(3) * t264 + t178 * pkin(6) + t154 * t254 + t259 * t279;
t44 = t59 * pkin(3) + qJDD(4) + t94;
t188 = -t279 * t261 + pkin(3) * t262 + t181 * pkin(6) + (-pkin(1) - t128) * t178;
t125 = t279 * t179;
t77 = t125 * t173 + t174 * t216;
t78 = t174 * t125 - t173 * t216;
t185 = -t78 * t28 + t29 * t77 - t54 * t64 - t292;
t184 = t44 + t297;
t152 = -pkin(3) * t174 - pkin(4);
t150 = pkin(3) * t173 + qJ(5);
t146 = pkin(3) * t258;
t102 = t179 * t254 + t264;
t101 = -t225 + t258;
t100 = -t178 * t256 + t262;
t90 = -t173 * t265 + t174 * t260;
t89 = t112 * t177;
t87 = t159 * t178 + t160 * t254;
t85 = t160 * t257 - t253;
t61 = pkin(4) * t111 - qJ(5) * t112 - t154;
t52 = t173 * t220 - t174 * t219 + t177 * t97;
t51 = t111 * t235 - t112 * t237;
t46 = pkin(4) * t89 - qJ(5) * t90 + t245;
t37 = pkin(4) * t180 - t39;
t36 = -qJ(5) * t180 + t40;
t25 = pkin(3) * t118 + pkin(4) * t202 + qJ(5) * t64;
t18 = -qJ(5) * t142 + t21;
t17 = pkin(4) * t142 + qJD(5) - t20;
t16 = -pkin(4) * t51 + qJ(5) * t52 - qJD(5) * t90 + t224;
t7 = -pkin(4) * t238 - t10;
t6 = qJ(5) * t238 - qJD(5) * t180 + t11;
t5 = t44 + t294;
t2 = -t226 - t287;
t1 = -qJD(5) * t142 + t229;
t8 = [qJDD(1), -t282 + t285, t207, qJDD(1) * t171 + 0.2e1 * t177 * t214, 0.2e1 * t177 * t163 - 0.2e1 * t244 * t231, qJDD(2) * t177 + t180 * t182, qJDD(2) * t180 - t177 * t182, 0, t196 * t177 + (-t191 + t285) * t180, t177 * t191 + t180 * t196 - t149, t58 * t260 + (-t176 * t235 + t219) * t118, (-t116 * t179 - t118 * t176) * t237 + (-t270 - t179 * t59 + (t116 * t176 - t267) * qJD(3)) * t177, (-t58 - t221) * t180 + (t194 + t240) * t177, (t142 * t239 + t59) * t180 + (-t195 - t241) * t177, -t110 * t180 - t142 * t238, -(-t124 * t236 + t249) * t142 + t114 * t110 - g(1) * t100 - g(2) * t102 + ((t217 + t241) * pkin(6) + (-pkin(6) * t110 + qJD(2) * t126 - t211) * t176 + t206) * t180 + (pkin(6) * t59 + qJD(2) * t73 + t126 * t234 + t94 * t176) * t177, t250 * t142 - t246 * t110 - g(1) * t99 - g(2) * t101 + (t126 * t233 + (-t218 + t240) * pkin(6) + t193) * t180 + (-t126 * t236 - t74 * qJD(2) + t94 * t179 + (t58 - t221) * pkin(6)) * t177, -t10 * t202 - t11 * t64 + t20 * t52 + t21 * t51 - t28 * t40 - t29 * t39 - t3 * t90 - t4 * t89 + t209, -g(1) * t188 - g(2) * t190 + t20 * t10 + t21 * t11 + t224 * t79 + t245 * t44 + t3 * t39 + t4 * t40, g(1) * t85 - g(2) * t87 - t110 * t37 + t142 * t7 + t16 * t64 - t17 * t238 + t180 * t2 - t24 * t51 + t28 * t46 + t5 * t89, -t1 * t89 - t17 * t52 + t18 * t51 + t2 * t90 + t202 * t7 - t28 * t36 + t29 * t37 - t6 * t64 + t209, g(1) * t84 - g(2) * t86 - t1 * t180 + t110 * t36 - t142 * t6 - t16 * t202 + t18 * t238 + t24 * t52 - t29 * t46 - t5 * t90, t1 * t36 + t18 * t6 + t5 * t46 + t24 * t16 + t2 * t37 + t17 * t7 - g(1) * (-t85 * pkin(4) - t84 * qJ(5) + t188) - g(2) * (pkin(4) * t87 + qJ(5) * t86 + t190); 0, 0, 0, -t177 * t183 * t180, t244 * t183, t230, t163, qJDD(2), t177 * t192 - t155 - t169, t281 + (-pkin(6) * qJDD(1) + t192) * t180, -t142 * t267 + t270, (t58 + t269) * t179 + (-t59 + t268) * t176, (-t118 * t177 + t142 * t256) * qJD(1) + t195, (t116 * t177 - t142 * t263) * qJD(1) + t194, t142 * t243, -pkin(2) * t59 + t248 * t142 + t203 * t176 + (-t73 * t177 + (-pkin(6) * t116 - t126 * t176) * t180) * qJD(1) - t302 * t179, -pkin(2) * t58 - t103 * t142 + t203 * t179 + (-t126 * t256 + t74 * t177 + (-t118 * t180 + t142 * t260) * pkin(6)) * qJD(1) + t302 * t176, -t111 * t4 - t112 * t3 + t273 * t20 - t202 * t277 + t274 * t21 + t32 * t64 + t185, t4 * t78 - t3 * t77 - t44 * t154 - g(3) * (t177 * t279 + t128) + t295 * t79 + (t54 - t32) * t21 + t277 * t20 + t207 * (t154 * t177 - t180 * t279), -t77 * t110 + t5 * t111 + t276 * t142 - t160 * t297 + t17 * t243 - t274 * t24 + t278 * t64 + t61 * t28, -t1 * t111 + t112 * t2 - t273 * t17 + t274 * t18 + t202 * t276 + t26 * t64 + t185, t78 * t110 - t5 * t112 - t275 * t142 - t159 * t297 - t18 * t243 - t202 * t278 + t273 * t24 - t61 * t29, -g(3) * t128 + t1 * t78 + t2 * t77 + t5 * t61 + t278 * t24 + t275 * t18 + t276 * t17 + (-g(3) * t204 - t207 * t279) * t180 + (-g(3) * t279 + t207 * (t154 + t204)) * t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118 * t116, -t116 ^ 2 + t118 ^ 2, t58 - t269, -t268 - t59, t110, -g(1) * t101 + g(2) * t99 - t126 * t118 - t74 * t142 + (t211 + t281) * t176 - t206, g(1) * t102 - g(2) * t100 + g(3) * t260 + t116 * t126 - t142 * t73 - t193, t21 * t202 - t280 + (-t173 * t28 - t174 * t29) * pkin(3) + (-t20 + t23) * t64, -g(1) * t146 + t20 * t22 - t21 * t23 + (g(2) * t255 - t79 * t118 + t4 * t173 + t3 * t174 + t292 * t176) * pkin(3), -t22 * t142 - t25 * t64 + (pkin(4) - t152) * t110 + t293, -t150 * t28 + t152 * t29 + t18 * t202 - t280 + (t17 - t251) * t64, -t160 * t281 - g(1) * t87 - g(2) * t85 + t150 * t110 - t24 * t64 + t25 * t202 + (-0.2e1 * qJD(5) + t23) * t142 + t229, t1 * t150 + t2 * t152 - t24 * t25 - t17 * t22 - g(1) * (-pkin(3) * t225 - pkin(4) * t86 + qJ(5) * t87 + t146) - g(2) * (-pkin(3) * t99 - t84 * pkin(4) + t85 * qJ(5)) - g(3) * (-t144 + (-pkin(4) * t159 + qJ(5) * t160) * t177) + t251 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, t20 * t202 + t21 * t64 + t184, -t142 * t202 + t28, t205, -t29 - t301, -t17 * t202 + t18 * t64 + t184 + t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202 * t64 - t110, t29 - t301, -t142 ^ 2 - t299, t142 * t18 - t287 - t293;];
tau_reg = t8;
