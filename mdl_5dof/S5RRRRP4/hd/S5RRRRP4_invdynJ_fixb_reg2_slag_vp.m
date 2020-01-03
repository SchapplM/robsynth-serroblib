% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRP4
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:18
% EndTime: 2019-12-31 21:51:23
% DurationCPUTime: 2.66s
% Computational Cost: add. (4629->370), mult. (7052->427), div. (0->0), fcn. (4468->12), ass. (0->216)
t176 = sin(qJ(3));
t182 = -pkin(8) - pkin(7);
t240 = qJD(3) * t182;
t105 = t176 * t240;
t179 = cos(qJ(3));
t295 = cos(qJ(4));
t238 = t295 * t179;
t175 = sin(qJ(4));
t267 = t175 * t176;
t203 = t238 - t267;
t130 = t182 * t176;
t164 = t179 * pkin(8);
t131 = t179 * pkin(7) + t164;
t204 = t295 * t130 - t175 * t131;
t180 = cos(qJ(2));
t255 = qJD(1) * t180;
t245 = pkin(1) * t255;
t266 = t175 * t179;
t282 = t204 * qJD(4) + t295 * t105 - t203 * t245 + t240 * t266;
t171 = t176 ^ 2;
t172 = t179 ^ 2;
t256 = t171 + t172;
t168 = qJDD(1) + qJDD(2);
t177 = sin(qJ(2));
t249 = qJDD(1) * t177;
t253 = qJD(2) * t180;
t93 = t168 * pkin(7) + (qJD(1) * t253 + t249) * pkin(1);
t230 = t256 * t93;
t174 = qJ(1) + qJ(2);
t163 = cos(t174);
t293 = g(2) * t163;
t280 = pkin(1) * qJD(1);
t246 = t177 * t280;
t287 = t180 * pkin(1);
t258 = -qJD(2) * t246 + qJDD(1) * t287;
t291 = t168 * pkin(2);
t92 = -t258 - t291;
t311 = t92 + t293;
t170 = qJD(1) + qJD(2);
t303 = t256 * t180;
t310 = t170 * t303;
t219 = qJD(3) * t238;
t233 = qJD(4) * t295;
t309 = -t179 * t233 - t219;
t107 = t170 * pkin(7) + t246;
t290 = t170 * pkin(2);
t108 = -t245 - t290;
t308 = t107 * t303 + t108 * t177;
t169 = qJD(3) + qJD(4);
t208 = t169 * t267;
t239 = t295 * t176;
t264 = t179 * t168;
t222 = -t168 * t239 + t309 * t170 - t175 * t264;
t34 = t170 * t208 + t222;
t265 = t176 * t168;
t210 = -t168 * t238 + t175 * t265;
t103 = t239 + t266;
t70 = t169 * t103;
t35 = t70 * t170 + t210;
t288 = t179 * pkin(3);
t152 = pkin(2) + t288;
t252 = qJD(3) * t176;
t236 = t170 * t252;
t63 = pkin(3) * t236 - t152 * t168 - t258;
t89 = t103 * t170;
t10 = t35 * pkin(4) + t34 * qJ(5) - t89 * qJD(5) + t63;
t243 = t170 * t267;
t87 = -t170 * t238 + t243;
t91 = -t152 * t170 - t245;
t40 = t87 * pkin(4) - t89 * qJ(5) + t91;
t69 = t208 + t309;
t307 = -t10 * t103 + t40 * t69;
t306 = -t10 * t203 + t40 * t70;
t305 = -t203 * t63 + t91 * t70;
t304 = t63 * t103 - t91 * t69;
t173 = qJ(3) + qJ(4);
t160 = sin(t173);
t162 = cos(t173);
t261 = t162 * pkin(4) + t160 * qJ(5);
t161 = sin(t174);
t302 = g(1) * t163 + g(2) * t161;
t167 = qJDD(3) + qJDD(4);
t155 = t167 * qJ(5);
t158 = t169 * qJD(5);
t301 = t155 + t158;
t159 = t167 * pkin(4);
t300 = qJDD(5) - t159;
t250 = qJD(4) * t175;
t251 = qJD(3) * t179;
t235 = t170 * t251;
t45 = -t107 * t251 + qJDD(3) * pkin(3) - t176 * t93 + (-t235 - t265) * pkin(8);
t50 = -t107 * t252 + t179 * t93 + (-t236 + t264) * pkin(8);
t232 = pkin(8) * t170 + t107;
t80 = t232 * t176;
t76 = qJD(3) * pkin(3) - t80;
t81 = t232 * t179;
t12 = t175 * t45 + t76 * t233 - t81 * t250 + t295 * t50;
t229 = t175 * t50 + t81 * t233 + t76 * t250 - t295 * t45;
t279 = t175 * t81;
t48 = t295 * t76 - t279;
t244 = t295 * t81;
t49 = t175 * t76 + t244;
t299 = t103 * t229 + t12 * t203 + t48 * t69 - t49 * t70;
t11 = t229 + t300;
t263 = qJD(5) - t48;
t36 = -t169 * pkin(4) + t263;
t39 = t169 * qJ(5) + t49;
t9 = t12 + t301;
t298 = t11 * t103 + t203 * t9 - t36 * t69 - t39 * t70;
t272 = t160 * t163;
t273 = t160 * t161;
t262 = g(1) * t273 - g(2) * t272;
t78 = t175 * t130 + t295 * t131;
t297 = -t78 * t167 - t282 * t169 - t262;
t296 = t89 ^ 2;
t147 = g(1) * t161;
t178 = sin(qJ(1));
t294 = g(1) * t178;
t292 = g(3) * t179;
t289 = t178 * pkin(1);
t286 = t40 * t87;
t285 = t89 * t87;
t284 = t91 * t87;
t150 = t177 * pkin(1) + pkin(7);
t283 = -pkin(8) - t150;
t281 = t78 * qJD(4) - t103 * t245 + t175 * t105 - t182 * t219;
t278 = t49 * t169;
t276 = t108 * t252 + t179 * t147;
t52 = -t295 * t80 - t279;
t275 = pkin(3) * t233 + qJD(5) - t52;
t271 = t161 * t162;
t270 = t162 * t163;
t269 = t163 * t182;
t268 = t170 * t176;
t260 = t163 * pkin(2) + t161 * pkin(7);
t257 = t171 - t172;
t254 = qJD(2) * t177;
t248 = t108 * t251 + t311 * t176;
t247 = pkin(1) * t253;
t156 = pkin(3) * t252;
t41 = t87 ^ 2 - t296;
t166 = t170 ^ 2;
t242 = t176 * t166 * t179;
t112 = t163 * t152;
t241 = pkin(4) * t270 + qJ(5) * t272 + t112;
t237 = t170 * t254;
t231 = qJD(3) * t283;
t227 = t256 * t168;
t225 = -t161 * t182 + t112;
t224 = -t302 + t230;
t223 = t170 * t246;
t221 = t176 * t235;
t32 = t70 * pkin(4) + t69 * qJ(5) - t103 * qJD(5) + t156;
t218 = -t32 + t246;
t51 = -t175 * t80 + t244;
t217 = pkin(3) * t250 - t51;
t216 = g(1) * t271 - g(2) * t270;
t214 = g(1) * (-t161 * pkin(2) + t163 * pkin(7));
t213 = -pkin(3) * t176 - pkin(4) * t160;
t61 = t89 * pkin(4) + t87 * qJ(5);
t181 = cos(qJ(1));
t211 = -g(2) * t181 + t294;
t209 = -t203 * t35 + t87 * t70;
t55 = -t167 * t203 + t169 * t70;
t207 = -t161 * t152 - t269;
t206 = t147 + t258 - t293;
t100 = t283 * t176;
t101 = t179 * t150 + t164;
t205 = t295 * t100 - t175 * t101;
t65 = t175 * t100 + t295 * t101;
t202 = -g(1) * t270 - g(2) * t271 - g(3) * t160 + t12;
t201 = -t246 + t156;
t193 = -t176 * t247 + t179 * t231;
t79 = t176 * t231 + t179 * t247;
t19 = t205 * qJD(4) + t175 * t193 + t295 * t79;
t200 = -t65 * t167 - t19 * t169 - t262;
t199 = g(1) * t272 + g(2) * t273 - g(3) * t162 - t229;
t66 = -pkin(4) * t203 - t103 * qJ(5) - t152;
t198 = -t108 * t170 + t302 - t93;
t20 = t65 * qJD(4) + t175 * t79 - t295 * t193;
t197 = -t19 * t87 + t20 * t89 + t205 * t34 - t65 * t35 - t302;
t196 = t48 * t169 - t202;
t183 = qJD(3) ^ 2;
t195 = pkin(7) * t183 - t223 - t291;
t194 = t167 * t205 - t20 * t169 + t216;
t2 = -t103 * t35 - t203 * t34 + t69 * t87 - t89 * t70;
t153 = -pkin(2) - t287;
t192 = pkin(1) * t237 + t150 * t183 + t153 * t168;
t191 = -t91 * t89 + t199;
t190 = -pkin(7) * qJDD(3) + (t245 - t290) * qJD(3);
t189 = t167 * t204 - t281 * t169 + t216;
t188 = -qJDD(3) * t150 + (t153 * t170 - t247) * qJD(3);
t187 = (-g(1) * (-t152 - t261) + g(2) * t182) * t161;
t186 = t204 * t34 + t281 * t89 - t282 * t87 - t78 * t35 - t302;
t185 = t40 * t89 - t199 + t300;
t165 = t181 * pkin(1);
t157 = pkin(1) * t254;
t151 = -t295 * pkin(3) - pkin(4);
t140 = t175 * pkin(3) + qJ(5);
t125 = -t152 - t287;
t124 = t176 * qJDD(3) + t183 * t179;
t123 = qJDD(3) * t179 - t183 * t176;
t111 = qJ(5) * t270;
t109 = qJ(5) * t271;
t106 = t157 + t156;
t95 = t172 * t168 - 0.2e1 * t221;
t94 = t171 * t168 + 0.2e1 * t221;
t71 = -0.2e1 * t257 * t170 * qJD(3) + 0.2e1 * t176 * t264;
t62 = t66 - t287;
t56 = pkin(3) * t268 + t61;
t53 = t103 * t167 - t69 * t169;
t28 = -t89 * t169 + t35;
t27 = -t222 + (-t243 + t87) * t169;
t25 = t157 + t32;
t14 = -t34 * t103 - t89 * t69;
t1 = [0, 0, 0, 0, 0, qJDD(1), t211, g(1) * t181 + g(2) * t178, 0, 0, 0, 0, 0, 0, 0, t168, (t168 * t180 - t237) * pkin(1) + t206, ((-qJDD(1) - t168) * t177 + (-qJD(1) - t170) * t253) * pkin(1) + t302, 0, (t211 + (t177 ^ 2 + t180 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t94, t71, t124, t95, t123, 0, t188 * t176 + (-t192 - t311) * t179 + t276, t188 * t179 + (t192 - t147) * t176 + t248, pkin(1) * qJD(2) * t310 + t150 * t227 + t224, t92 * t153 - t214 - g(2) * (t165 + t260) + t150 * t230 + (t308 * qJD(2) + t294) * pkin(1), t14, t2, t53, t209, -t55, 0, t106 * t87 + t125 * t35 + t194 + t305, t106 * t89 - t125 * t34 + t200 + t304, t197 + t299, t12 * t65 + t49 * t19 - t229 * t205 - t48 * t20 + t63 * t125 + t91 * t106 - g(1) * (t207 - t289) - g(2) * (t165 + t225), t14, t53, -t2, 0, t55, t209, t25 * t87 + t62 * t35 + t194 + t306, t197 + t298, -t25 * t89 + t62 * t34 - t200 + t307, t9 * t65 + t39 * t19 + t10 * t62 + t40 * t25 - t11 * t205 + t36 * t20 - g(1) * (-t269 - t289) - g(2) * (t165 + t241) + t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t206 + t223, (-t249 + (-qJD(2) + t170) * t255) * pkin(1) + t302, 0, 0, t94, t71, t124, t95, t123, 0, t190 * t176 + (-t195 - t311) * t179 + t276, t190 * t179 + (t195 - t147) * t176 + t248, pkin(7) * t227 - t280 * t310 + t224, -t92 * pkin(2) + pkin(7) * t230 - g(2) * t260 - t308 * t280 - t214, t14, t2, t53, t209, -t55, 0, -t152 * t35 + t201 * t87 + t189 + t305, t152 * t34 + t201 * t89 + t297 + t304, t186 + t299, -g(1) * t207 - g(2) * t225 + t12 * t78 - t63 * t152 + t201 * t91 - t204 * t229 - t281 * t48 + t282 * t49, t14, t53, -t2, 0, t55, t209, -t218 * t87 + t66 * t35 + t189 + t306, t186 + t298, t218 * t89 + t66 * t34 - t297 + t307, g(1) * t269 - g(2) * t241 + t10 * t66 - t11 * t204 - t218 * t40 + t281 * t36 + t282 * t39 + t9 * t78 + t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t242, t257 * t166, t265, t242, t264, qJDD(3), t176 * t198 - t292, g(3) * t176 + t179 * t198, 0, 0, t285, -t41, t27, -t285, -t210, t167, t51 * t169 + (t295 * t167 - t169 * t250 - t87 * t268) * pkin(3) + t191, t52 * t169 + t284 + (-t167 * t175 - t169 * t233 - t268 * t89) * pkin(3) - t202, (t49 - t51) * t89 + (-t48 + t52) * t87 + (t295 * t34 - t175 * t35 + (t175 * t89 - t295 * t87) * qJD(4)) * pkin(3), t48 * t51 - t49 * t52 + (-t295 * t229 - t292 + t12 * t175 + (-t175 * t48 + t295 * t49) * qJD(4) + (-t170 * t91 + t302) * t176) * pkin(3), t285, t27, t41, t167, t28, -t285, -t151 * t167 - t169 * t217 - t56 * t87 - t185, -t140 * t35 - t151 * t34 + (t217 + t39) * t89 + (t36 - t275) * t87, t140 * t167 + t275 * t169 + t56 * t89 + t202 - t286 + t301, t9 * t140 + t11 * t151 - t40 * t56 - g(1) * (t163 * t213 + t111) - g(2) * (t161 * t213 + t109) - g(3) * (t261 + t288) + t275 * t39 + t217 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, -t41, t27, -t285, -t210, t167, t191 + t278, t196 + t284, 0, 0, t285, t27, t41, t167, t28, -t285, -t61 * t87 + t159 - t185 + t278, pkin(4) * t34 - t35 * qJ(5) + (t39 - t49) * t89 + (t36 - t263) * t87, t61 * t89 + 0.2e1 * t155 + 0.2e1 * t158 - t196 - t286, t9 * qJ(5) - t11 * pkin(4) - t40 * t61 - t36 * t49 - g(1) * (-pkin(4) * t272 + t111) - g(2) * (-pkin(4) * t273 + t109) - g(3) * t261 + t263 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167 + t285, t27, -t169 ^ 2 - t296, -t39 * t169 + t185;];
tau_reg = t1;
