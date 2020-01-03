% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:45
% EndTime: 2019-12-31 21:11:50
% DurationCPUTime: 2.66s
% Computational Cost: add. (3273->377), mult. (4876->438), div. (0->0), fcn. (2818->10), ass. (0->221)
t173 = sin(qJ(3));
t169 = t173 ^ 2;
t177 = cos(qJ(3));
t170 = t177 ^ 2;
t255 = t169 + t170;
t163 = qJDD(1) + qJDD(2);
t174 = sin(qJ(2));
t246 = qJDD(1) * t174;
t178 = cos(qJ(2));
t252 = qJD(2) * t178;
t77 = t163 * pkin(7) + (qJD(1) * t252 + t246) * pkin(1);
t312 = t77 * t255;
t172 = sin(qJ(5));
t176 = cos(qJ(5));
t247 = qJD(5) * t176;
t248 = qJD(5) * t172;
t249 = qJD(3) * t177;
t250 = qJD(3) * t173;
t43 = t172 * t249 + t173 * t247 - t176 * t250 - t177 * t248;
t296 = pkin(3) + pkin(4);
t165 = qJD(1) + qJD(2);
t268 = t165 * t173;
t282 = pkin(1) * qJD(1);
t242 = t174 * t282;
t105 = t165 * pkin(7) + t242;
t89 = t173 * t105;
t311 = -pkin(8) * t268 + qJD(4) + t89;
t42 = -t296 * qJD(3) + t311;
t168 = qJD(3) * qJ(4);
t90 = t177 * t105;
t58 = -t177 * t165 * pkin(8) + t90;
t50 = t168 + t58;
t21 = t172 * t42 + t176 * t50;
t137 = t173 * t163;
t232 = t165 * t249;
t64 = t173 * t77;
t238 = t105 * t249 + qJDD(4) + t64;
t22 = -t296 * qJDD(3) + (-t232 - t137) * pkin(8) + t238;
t233 = t165 * t250;
t264 = t177 * t163;
t166 = qJDD(3) * qJ(4);
t167 = qJD(3) * qJD(4);
t66 = t177 * t77;
t35 = -t105 * t250 + t166 + t167 + t66;
t25 = (t233 - t264) * pkin(8) + t35;
t5 = -t21 * qJD(5) - t172 * t25 + t176 * t22;
t164 = qJD(3) - qJD(5);
t98 = t173 * t172 + t177 * t176;
t194 = t98 * qJD(5);
t205 = t177 * t172 - t173 * t176;
t23 = t205 * t163 + t165 * t194 - t172 * t233 - t176 * t232;
t73 = t98 * t165;
t304 = t73 * t164 + t23;
t156 = t173 * qJ(4);
t228 = pkin(2) + t156;
t193 = t296 * t177 + t228;
t254 = qJD(1) * t178;
t241 = pkin(1) * t254;
t40 = t165 * t193 + t241;
t171 = qJ(1) + qJ(2);
t154 = sin(t171);
t68 = t205 * t154;
t155 = cos(t171);
t70 = t205 * t155;
t75 = t205 * t165;
t310 = g(1) * t70 + g(2) * t68 + g(3) * t98 + t40 * t75 + t5;
t266 = t170 * t163;
t267 = t169 * t163;
t308 = t266 + t267;
t290 = t165 * pkin(2);
t106 = -t241 - t290;
t224 = t178 * t255;
t307 = t105 * t224 + t106 * t174;
t275 = qJDD(3) * pkin(3);
t36 = t238 - t275;
t306 = t36 * t173 + t35 * t177;
t24 = t98 * t163 + t43 * t165;
t303 = -t75 * t164 + t24;
t300 = t296 * t173;
t299 = g(1) * t155 + g(2) * t154;
t223 = -qJD(3) * pkin(3) + qJD(4);
t72 = t223 + t89;
t81 = t90 + t168;
t298 = t72 * t249 - t81 * t250 + t306;
t251 = qJD(3) * t165;
t256 = -t169 + t170;
t297 = 0.2e1 * t173 * t264 + 0.2e1 * t256 * t251;
t295 = pkin(7) - pkin(8);
t181 = qJD(3) ^ 2;
t294 = pkin(7) * t181;
t143 = g(1) * t154;
t175 = sin(qJ(1));
t293 = g(1) * t175;
t292 = g(2) * t155;
t291 = t163 * pkin(2);
t158 = t177 * pkin(3);
t157 = t177 * pkin(4);
t289 = t178 * pkin(1);
t288 = t75 * t73;
t147 = t174 * pkin(1) + pkin(7);
t287 = -pkin(8) + t147;
t150 = pkin(8) * t250;
t100 = -pkin(7) * t250 + t150;
t122 = t295 * t177;
t101 = qJD(3) * t122;
t121 = t295 * t173;
t51 = t176 * t121 - t172 * t122;
t286 = t51 * qJD(5) + t176 * t100 + t172 * t101 - t98 * t241;
t52 = t172 * t121 + t176 * t122;
t285 = -t52 * qJD(5) - t172 * t100 + t176 * t101 + t205 * t241;
t107 = -t172 * qJ(4) - t176 * t296;
t284 = t107 * qJD(5) - t172 * t58 + t311 * t176;
t108 = t176 * qJ(4) - t172 * t296;
t283 = -t108 * qJD(5) - t311 * t172 - t176 * t58;
t20 = -t172 * t50 + t176 * t42;
t281 = t20 * t164;
t280 = t21 * t164;
t277 = pkin(7) * qJDD(3);
t276 = qJ(4) * t177;
t273 = t147 * t181;
t272 = t154 * t173;
t271 = t154 * t177;
t270 = t155 * t173;
t269 = t155 * t177;
t261 = g(1) * t272 - g(2) * t270;
t260 = qJD(2) * t242 - qJDD(1) * t289;
t259 = t155 * pkin(2) + t154 * pkin(7);
t257 = t158 + t156;
t253 = qJD(2) * t174;
t152 = t173 * qJD(4);
t245 = qJDD(3) * t147;
t244 = pkin(1) * t252;
t243 = pkin(1) * t253;
t240 = t73 ^ 2 - t75 ^ 2;
t131 = g(1) * t271;
t219 = t165 * t242;
t237 = t177 * t219 + t241 * t250 + t131;
t84 = pkin(3) * t250 - qJ(4) * t249 - t152;
t236 = t157 + t257;
t235 = -g(1) * t270 - g(2) * t272 + g(3) * t177;
t148 = -pkin(2) - t289;
t234 = t165 * t253;
t76 = t260 - t291;
t229 = -t76 - t292;
t95 = t287 * t177;
t140 = t155 * pkin(7);
t227 = -t155 * pkin(8) + t140;
t222 = t106 * t249 + t76 * t173 - t261;
t221 = t164 ^ 2;
t218 = -t64 - t235;
t217 = pkin(3) * t269 + t155 * t156 + t259;
t216 = t173 * t232;
t59 = -pkin(4) * t250 - t84;
t215 = t59 + t242;
t214 = g(1) * (-t154 * pkin(2) + t140);
t213 = -t291 + t294;
t179 = cos(qJ(1));
t211 = -g(2) * t179 + t293;
t4 = t172 * t22 + t176 * t25 + t42 * t247 - t50 * t248;
t44 = t98 * qJD(3) - t194;
t210 = -t20 * t44 + t205 * t5 - t21 * t43 - t4 * t98 + t299;
t209 = pkin(3) * t173 - t276;
t94 = t287 * t173;
t37 = -t172 * t95 + t176 * t94;
t38 = t172 * t94 + t176 * t95;
t208 = t81 * t173 - t72 * t177;
t207 = t173 * t72 + t177 * t81;
t159 = t179 * pkin(1);
t206 = t159 + t217;
t91 = t148 - t257;
t204 = t143 - t260 - t292;
t203 = -t228 - t158;
t202 = t255 * t165 * t244 + t308 * t147 - t299;
t114 = -pkin(2) - t257;
t26 = (t209 * qJD(3) - t152) * t165 + t203 * t163 + t260;
t201 = -t114 * t163 - t26 - t294;
t200 = t165 * t91 - t244;
t197 = t276 - t300;
t14 = (t197 * qJD(3) + t152) * t165 + t193 * t163 - t260;
t69 = t98 * t154;
t71 = t98 * t155;
t199 = g(1) * t69 - g(2) * t71 + t14 * t98 + t40 * t43;
t198 = -g(1) * t68 + g(2) * t70 - t14 * t205 + t40 * t44;
t60 = t84 + t243;
t191 = -t163 * t91 - t165 * t60 - t26 - t273;
t190 = pkin(1) * t234 + t148 * t163 + t273;
t189 = -g(1) * t140 - t203 * t143;
t188 = -t165 * t224 * t282 + t308 * pkin(7) - t299;
t187 = -t245 + (t148 * t165 - t244) * qJD(3);
t186 = -t208 * qJD(3) + t306;
t184 = (-g(1) * (t203 - t157) + g(2) * pkin(8)) * t154;
t183 = -g(1) * t71 - g(2) * t69 + g(3) * t205 - t40 * t73 + t4;
t162 = qJDD(3) - qJDD(5);
t161 = t165 ^ 2;
t126 = pkin(4) * t269;
t125 = qJ(4) * t269;
t123 = qJ(4) * t271;
t120 = t173 * t161 * t177;
t116 = t173 * qJDD(3) + t181 * t177;
t115 = qJDD(3) * t177 - t181 * t173;
t93 = t256 * t161;
t92 = pkin(2) + t236;
t87 = t106 * t250;
t85 = t209 * t165;
t80 = -0.2e1 * t216 + t266;
t79 = 0.2e1 * t216 + t267;
t78 = t157 - t91;
t56 = t197 * t165;
t54 = qJD(3) * t95 + t173 * t244;
t53 = -t147 * t250 + t177 * t244 + t150;
t48 = t203 * t165 - t241;
t47 = t59 - t243;
t41 = t48 * t250;
t32 = t162 * t205 - t44 * t164;
t31 = t98 * t162 + t43 * t164;
t11 = -t38 * qJD(5) - t172 * t53 + t176 * t54;
t10 = t37 * qJD(5) + t172 * t54 + t176 * t53;
t7 = t24 * t98 + t73 * t43;
t6 = t205 * t23 - t75 * t44;
t1 = t205 * t24 + t23 * t98 + t75 * t43 - t44 * t73;
t2 = [0, 0, 0, 0, 0, qJDD(1), t211, g(1) * t179 + g(2) * t175, 0, 0, 0, 0, 0, 0, 0, t163, (t163 * t178 - t234) * pkin(1) + t204, ((-qJDD(1) - t163) * t174 + (-qJD(1) - t165) * t252) * pkin(1) + t299, 0, (t211 + (t174 ^ 2 + t178 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t79, t297, t116, t80, t115, 0, t131 + t87 + t187 * t173 + (-t190 + t229) * t177, t173 * t190 + t177 * t187 + t222, t202 + t312, t76 * t148 - t214 - g(2) * (t159 + t259) + t147 * t312 + (t307 * qJD(2) + t293) * pkin(1), t79, t116, -t297, 0, -t115, t80, t131 + t41 + (qJD(3) * t200 - t245) * t173 + (t191 - t292) * t177, t202 + t298, (t245 + (-t200 - t48) * qJD(3)) * t177 + t191 * t173 + t261, t26 * t91 + t48 * t60 - g(2) * t206 + (t207 * t252 + t293) * pkin(1) + t186 * t147 + t189, t6, t1, t32, t7, t31, 0, -t11 * t164 - t37 * t162 + t78 * t24 + t47 * t73 + t199, t10 * t164 + t38 * t162 - t78 * t23 - t47 * t75 + t198, -t10 * t73 + t11 * t75 + t37 * t23 - t38 * t24 + t210, t4 * t38 + t21 * t10 + t5 * t37 + t20 * t11 + t14 * t78 + t40 * t47 - g(1) * (-t175 * pkin(1) + t227) - g(2) * (t126 + t206) + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, t204 + t219, (-t246 + (-qJD(2) + t165) * t254) * pkin(1) + t299, 0, 0, t79, t297, t116, t80, t115, 0, t87 + (-pkin(2) * t251 - t277) * t173 + (-t213 + t229) * t177 + t237, (-t277 + (t241 - t290) * qJD(3)) * t177 + (t213 - t219) * t173 + t222, t188 + t312, -t76 * pkin(2) + pkin(7) * t312 - g(2) * t259 - t307 * t282 - t214, t79, t116, -t297, 0, -t115, t80, t41 + (t114 * t251 - t277) * t173 + (-t165 * t84 + t201 - t292) * t177 + t237, t188 + t298, (t277 + (-t114 * t165 - t241 - t48) * qJD(3)) * t177 + ((-t84 + t242) * t165 + t201) * t173 + t261, t26 * t114 + t48 * t84 - g(2) * t217 + (-t174 * t48 - t178 * t207) * t282 + t186 * pkin(7) + t189, t6, t1, t32, t7, t31, 0, -t51 * t162 - t285 * t164 + t215 * t73 + t92 * t24 + t199, t52 * t162 + t286 * t164 - t215 * t75 - t92 * t23 + t198, t51 * t23 - t52 * t24 + t285 * t75 - t286 * t73 + t210, t4 * t52 + t5 * t51 + t14 * t92 - g(1) * t227 - g(2) * (t126 + t217) + t215 * t40 + t286 * t21 + t285 * t20 + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t93, t137, t120, t264, qJDD(3), -t106 * t268 + t218, g(3) * t173 - t66 + (-t106 * t165 + t299) * t177, 0, 0, -t120, t137, t93, qJDD(3), -t264, t120, 0.2e1 * t275 - qJDD(4) + (-t173 * t48 + t177 * t85) * t165 + t218, -t209 * t163 + ((t81 - t168) * t173 + (t223 - t72) * t177) * t165, 0.2e1 * t166 + 0.2e1 * t167 + t66 + (t165 * t85 - g(3)) * t173 + (t165 * t48 - t299) * t177, t35 * qJ(4) + t81 * qJD(4) - t36 * pkin(3) - t48 * t85 - g(1) * (-pkin(3) * t270 + t125) - g(2) * (-pkin(3) * t272 + t123) - g(3) * t257 + t208 * t105, t288, t240, t304, -t288, t303, t162, -t107 * t162 - t283 * t164 - t56 * t73 - t310, t108 * t162 + t164 * t284 + t56 * t75 + t183, t107 * t23 - t108 * t24 - (-t21 - t283) * t75 + (t20 - t284) * t73, -g(1) * t125 - g(2) * t123 - g(3) * t236 + t5 * t107 + t4 * t108 + t283 * t20 + t284 * t21 + t299 * t300 - t40 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) - t120, t137, -t169 * t161 - t181, -t81 * qJD(3) + t48 * t268 + t235 + t36, 0, 0, 0, 0, 0, 0, -t176 * t162 - t172 * t221 - t73 * t268, t172 * t162 - t176 * t221 + t268 * t75, -t303 * t172 + t304 * t176, -t40 * t268 + (t5 - t280) * t176 + (t4 + t281) * t172 + t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t288, -t240, -t304, t288, -t303, -t162, -t280 + t310, -t183 - t281, 0, 0;];
tau_reg = t2;
