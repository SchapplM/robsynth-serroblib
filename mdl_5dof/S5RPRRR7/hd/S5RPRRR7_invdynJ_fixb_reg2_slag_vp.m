% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:15
% EndTime: 2019-12-31 19:04:23
% DurationCPUTime: 4.59s
% Computational Cost: add. (5227->440), mult. (11114->587), div. (0->0), fcn. (7248->14), ass. (0->221)
t180 = cos(qJ(3));
t261 = qJD(1) * t180;
t146 = -qJD(4) + t261;
t176 = sin(qJ(3));
t168 = qJ(1) + pkin(9);
t158 = sin(t168);
t159 = cos(t168);
t215 = g(1) * t159 + g(2) * t158;
t293 = g(3) * t180;
t194 = t215 * t176 - t293;
t172 = sin(pkin(9));
t150 = pkin(1) * t172 + pkin(6);
t131 = t150 * qJD(1);
t256 = qJD(3) * t180;
t126 = t150 * qJDD(1);
t312 = qJD(2) * qJD(3) + t126;
t222 = -t180 * qJDD(2) + t131 * t256 + t312 * t176;
t279 = qJDD(3) * pkin(3);
t43 = t222 - t279;
t191 = -t43 + t194;
t315 = -qJD(4) * pkin(7) * t146 - t191;
t175 = sin(qJ(4));
t179 = cos(qJ(4));
t249 = t179 * qJD(3);
t262 = qJD(1) * t176;
t114 = t175 * t262 - t249;
t258 = qJD(3) * t175;
t116 = t179 * t262 + t258;
t174 = sin(qJ(5));
t178 = cos(qJ(5));
t209 = t114 * t174 - t178 * t116;
t61 = t178 * t114 + t116 * t174;
t292 = t61 * t209;
t182 = -pkin(8) - pkin(7);
t238 = qJD(4) * t182;
t216 = pkin(3) * t176 - pkin(7) * t180;
t120 = t216 * qJD(1);
t113 = t176 * t131;
t91 = qJD(2) * t180 - t113;
t53 = t175 * t120 + t179 * t91;
t314 = -t53 + (pkin(8) * t261 + t238) * t175;
t265 = t179 * t180;
t205 = pkin(4) * t176 - pkin(8) * t265;
t52 = t179 * t120 - t175 * t91;
t313 = -qJD(1) * t205 + t179 * t238 - t52;
t248 = qJD(1) * qJD(3);
t230 = t180 * t248;
t245 = t176 * qJDD(1);
t311 = qJD(3) * qJD(4) + t230 + t245;
t310 = t209 ^ 2 - t61 ^ 2;
t141 = -qJD(5) + t146;
t253 = qJD(4) * t176;
t229 = qJD(1) * t253;
t221 = t311 * t175 + t179 * t229;
t202 = t179 * qJDD(3) - t221;
t250 = qJD(5) * t178;
t251 = qJD(5) * t174;
t56 = (-qJDD(3) + t229) * t175 - t311 * t179;
t15 = t114 * t250 + t116 * t251 - t174 * t202 + t178 * t56;
t309 = -t141 * t61 - t15;
t260 = qJD(2) * t176;
t92 = t131 * t180 + t260;
t81 = qJD(3) * pkin(7) + t92;
t217 = pkin(3) * t180 + pkin(7) * t176;
t206 = -pkin(2) - t217;
t173 = cos(pkin(9));
t300 = pkin(1) * t173;
t107 = t206 - t300;
t82 = t107 * qJD(1);
t37 = -t175 * t81 + t179 * t82;
t28 = -pkin(8) * t116 + t37;
t25 = -pkin(4) * t146 + t28;
t38 = t175 * t82 + t179 * t81;
t29 = -pkin(8) * t114 + t38;
t240 = -t176 * qJDD(2) - t312 * t180;
t257 = qJD(3) * t176;
t48 = -t131 * t257 - t240;
t42 = qJDD(3) * pkin(7) + t48;
t123 = t216 * qJD(3);
t58 = qJD(1) * t123 + qJDD(1) * t107;
t51 = t179 * t58;
t11 = -qJD(4) * t38 - t175 * t42 + t51;
t163 = t180 * qJDD(1);
t110 = t176 * t248 + qJDD(4) - t163;
t6 = pkin(4) * t110 + pkin(8) * t56 + t11;
t252 = qJD(4) * t179;
t254 = qJD(4) * t175;
t10 = t175 * t58 + t179 * t42 + t82 * t252 - t81 * t254;
t7 = pkin(8) * t202 + t10;
t1 = (qJD(5) * t25 + t7) * t178 + t174 * t6 - t29 * t251;
t171 = qJ(4) + qJ(5);
t166 = cos(t171);
t294 = g(3) * t176;
t80 = -qJD(3) * pkin(3) - t91;
t57 = pkin(4) * t114 + t80;
t165 = sin(t171);
t271 = t166 * t180;
t76 = -t158 * t271 + t159 * t165;
t78 = t158 * t165 + t159 * t271;
t308 = g(1) * t78 - g(2) * t76 + t166 * t294 + t57 * t61 - t1;
t281 = t178 * t29;
t9 = t174 * t25 + t281;
t2 = -qJD(5) * t9 - t174 * t7 + t178 * t6;
t272 = t165 * t180;
t75 = t158 * t272 + t159 * t166;
t77 = t158 * t166 - t159 * t272;
t307 = -g(1) * t77 + g(2) * t75 + t165 * t294 + t57 * t209 + t2;
t190 = qJD(5) * t209 + t174 * t56 + t178 * t202;
t306 = t141 * t209 + t190;
t305 = t146 * t37 + t10;
t280 = pkin(1) * qJDD(1);
t234 = t180 * t249;
t200 = -t175 * t253 + t234;
t157 = pkin(4) * t179 + pkin(3);
t208 = t157 * t180 - t176 * t182;
t268 = t175 * t180;
t86 = t158 * t268 + t159 * t179;
t88 = t158 * t179 - t159 * t268;
t304 = -g(1) * t88 + g(2) * t86;
t244 = qJD(4) + qJD(5);
t267 = t176 * t179;
t303 = t110 * t267 - t146 * t200;
t299 = pkin(4) * t175;
t298 = g(1) * t158;
t295 = g(2) * t159;
t235 = t175 * t256;
t118 = t174 * t179 + t175 * t178;
t68 = t244 * t118;
t32 = t174 * t235 + t176 * t68 - t178 * t234;
t270 = t174 * t175;
t117 = -t178 * t179 + t270;
t94 = t117 * t176;
t291 = -t190 * t94 + t32 * t61;
t106 = qJDD(5) + t110;
t269 = t175 * t176;
t33 = -t251 * t269 + (t244 * t267 + t235) * t178 + t200 * t174;
t93 = t118 * t176;
t290 = -t93 * t106 + t33 * t141;
t133 = t182 * t175;
t134 = t182 * t179;
t71 = t133 * t178 + t134 * t174;
t289 = qJD(5) * t71 + t313 * t174 + t314 * t178;
t72 = t133 * t174 - t134 * t178;
t288 = -qJD(5) * t72 - t314 * t174 + t313 * t178;
t197 = t202 * t179;
t287 = -t114 * t234 + t176 * t197;
t286 = -t117 * t261 - t178 * t252 - t179 * t250 + t244 * t270;
t285 = -t118 * t261 + t68;
t283 = t146 * t38;
t282 = t174 * t29;
t119 = t150 * t265;
t66 = t175 * t107 + t119;
t278 = t114 * t146;
t277 = t116 * t114;
t276 = t116 * t146;
t275 = t150 * t175;
t273 = t159 * t175;
t264 = t179 * t123 + t257 * t275;
t169 = t176 ^ 2;
t170 = t180 ^ 2;
t263 = t169 - t170;
t151 = -pkin(2) - t300;
t132 = qJD(1) * t151;
t259 = qJD(3) * t114;
t255 = qJD(4) * t114;
t127 = qJDD(1) * t151;
t184 = qJD(1) ^ 2;
t241 = t176 * t184 * t180;
t181 = cos(qJ(1));
t239 = t181 * pkin(1) + t159 * pkin(2) + t158 * pkin(6);
t237 = t116 * t256;
t236 = t146 * t258;
t232 = t176 * t252;
t177 = sin(qJ(1));
t228 = -pkin(1) * t177 + t159 * pkin(6);
t226 = t15 * t180 - t209 * t257;
t224 = t116 * t257 + t180 * t56;
t223 = -t56 + t255;
t220 = t116 * t232;
t219 = t176 * t230;
t218 = pkin(4) * t254 - t260 - (qJD(1) * t299 + t131) * t180;
t214 = g(1) * t177 - g(2) * t181;
t213 = -t15 * t93 - t209 * t33;
t212 = t106 * t94 - t141 * t32;
t96 = t179 * t107;
t44 = -pkin(8) * t267 + t96 + (-pkin(4) - t275) * t180;
t55 = -pkin(8) * t269 + t66;
t20 = -t174 * t55 + t178 * t44;
t21 = t174 * t44 + t178 * t55;
t211 = -t175 * t38 - t179 * t37;
t210 = t175 * t37 - t179 * t38;
t204 = -t180 * t190 - t61 * t257;
t203 = -t110 * t175 + t146 * t252;
t201 = -qJD(1) * t132 + t215;
t199 = t232 + t235;
t198 = -pkin(7) * t110 - t146 * t80;
t183 = qJD(3) ^ 2;
t196 = t150 * t183 + 0.2e1 * t127 + t295;
t195 = 0.2e1 * t132 * qJD(3) - qJDD(3) * t150;
t192 = -t215 * t180 - t294;
t30 = t175 * t123 + t107 * t252 + (-t176 * t249 - t180 * t254) * t150;
t189 = qJD(4) * t211 + t10 * t179 - t11 * t175;
t188 = t222 * t176 + t48 * t180 + (-t176 * t92 - t180 * t91) * qJD(3);
t135 = t176 * t298;
t125 = qJDD(3) * t180 - t176 * t183;
t124 = qJDD(3) * t176 + t180 * t183;
t99 = (t150 + t299) * t176;
t89 = t158 * t175 + t159 * t265;
t87 = -t158 * t265 + t273;
t69 = pkin(4) * t199 + t150 * t256;
t65 = -t150 * t268 + t96;
t31 = -t66 * qJD(4) + t264;
t24 = -pkin(8) * t199 + t30;
t23 = -pkin(4) * t202 + t43;
t22 = t205 * qJD(3) + (-t119 + (pkin(8) * t176 - t107) * t175) * qJD(4) + t264;
t13 = t178 * t28 - t282;
t12 = -t174 * t28 - t281;
t8 = t178 * t25 - t282;
t4 = -qJD(5) * t21 - t174 * t24 + t178 * t22;
t3 = qJD(5) * t20 + t174 * t22 + t178 * t24;
t5 = [0, 0, 0, 0, 0, qJDD(1), t214, g(1) * t181 + g(2) * t177, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t173 * t280 - t295 + t298, -0.2e1 * t172 * t280 + t215, 0, (t214 + (t172 ^ 2 + t173 ^ 2) * t280) * pkin(1), qJDD(1) * t169 + 0.2e1 * t219, 0.2e1 * t176 * t163 - 0.2e1 * t263 * t248, t124, qJDD(1) * t170 - 0.2e1 * t219, t125, 0, t195 * t176 + (-t196 + t298) * t180, t176 * t196 + t180 * t195 - t135, (t169 + t170) * t126 + t188 - t215, t127 * t151 - g(1) * (-pkin(2) * t158 + t228) - g(2) * t239 + t188 * t150, t116 * t200 - t56 * t267, -t220 + (-t237 + (t56 + t255) * t176) * t175 + t287, t224 + t303, t114 * t199 - t202 * t269, (-t202 + t236) * t180 + (t203 - t259) * t176, -t110 * t180 - t146 * t257, -g(1) * t87 - g(2) * t89 + t65 * t110 - t31 * t146 + (-t11 + (t114 * t150 + t175 * t80) * qJD(3)) * t180 + (t37 * qJD(3) - t150 * t202 + t43 * t175 + t80 * t252) * t176, -g(1) * t86 - g(2) * t88 - t110 * t66 + t146 * t30 + (t10 + (t116 * t150 + t179 * t80) * qJD(3)) * t180 + (-qJD(3) * t38 - t150 * t56 + t43 * t179 - t80 * t254) * t176, -t30 * t114 + t66 * t202 - t31 * t116 + t65 * t56 + t135 + t211 * t256 + (qJD(4) * t210 - t10 * t175 - t11 * t179 - t295) * t176, t10 * t66 + t38 * t30 + t11 * t65 + t37 * t31 - g(1) * t228 - g(2) * (t159 * t217 + t239) - t206 * t298 + (t176 * t43 + t80 * t256) * t150, t15 * t94 + t209 * t32, -t213 + t291, -t212 + t226, -t190 * t93 + t33 * t61, t204 + t290, -t106 * t180 - t141 * t257, -g(1) * t76 - g(2) * t78 + t106 * t20 - t141 * t4 - t180 * t2 - t190 * t99 + t23 * t93 + t8 * t257 + t33 * t57 + t61 * t69, -g(1) * t75 - g(2) * t77 + t1 * t180 - t106 * t21 + t141 * t3 - t15 * t99 - t209 * t69 - t23 * t94 - t9 * t257 - t32 * t57, -t1 * t93 + t15 * t20 - t176 * t295 + t190 * t21 + t2 * t94 + t209 * t4 - t3 * t61 + t32 * t8 - t33 * t9 + t135, t1 * t21 + t9 * t3 + t2 * t20 + t8 * t4 + t23 * t99 + t57 * t69 - g(1) * (pkin(4) * t273 + t228) - g(2) * (t159 * t208 + t239) + (-g(1) * (-pkin(2) - t208) - g(2) * t299) * t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t125, -t124, 0, t176 * t48 - t180 * t222 - g(3) + (-t176 * t91 + t180 * t92) * qJD(3), 0, 0, 0, 0, 0, 0, (t202 + t236) * t180 + (t203 + t259) * t176, t224 - t303, t220 + (t176 * t223 + t237) * t175 + t287, -g(3) + (-qJD(3) * t210 - t43) * t180 + (qJD(3) * t80 + t189) * t176, 0, 0, 0, 0, 0, 0, -t204 + t290, t212 + t226, t213 + t291, -t1 * t94 - t180 * t23 - t2 * t93 + t57 * t257 - t32 * t9 - t33 * t8 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, t263 * t184, t245, t241, t163, qJDD(3), qJD(3) * t92 + t176 * t201 - t222 - t293, t294 + (t91 + t113) * qJD(3) + t201 * t180 + t240, 0, 0, -t56 * t175 - t179 * t276, (-t56 + t278) * t179 + (t202 + t276) * t175, (-t116 * t176 + t146 * t265) * qJD(1) - t203, -t175 * t278 + t197, t146 * t254 + t110 * t179 + (t114 * t176 - t146 * t268) * qJD(1), t146 * t262, -pkin(3) * t221 + t52 * t146 - t37 * t262 - t92 * t114 + t198 * t175 + (t279 - t315) * t179, pkin(3) * t56 - t116 * t92 - t146 * t53 + t315 * t175 + t198 * t179 + t38 * t262, t53 * t114 + t52 * t116 + ((qJD(4) * t116 + t202) * pkin(7) + t305) * t179 + (pkin(7) * t223 - t11 + t283) * t175 + t192, -t37 * t52 - t38 * t53 - t80 * t92 + t191 * pkin(3) + (t189 + t192) * pkin(7), -t118 * t15 + t209 * t286, t117 * t15 + t118 * t190 + t209 * t285 + t286 * t61, t106 * t118 + t286 * t141 + t209 * t262, -t117 * t190 + t285 * t61, -t106 * t117 + t285 * t141 + t61 * t262, t141 * t262, t106 * t71 + t117 * t23 - t288 * t141 + t157 * t190 + t194 * t166 + t218 * t61 - t8 * t262 + t285 * t57, -t106 * t72 + t118 * t23 + t289 * t141 + t15 * t157 - t165 * t194 - t209 * t218 + t9 * t262 - t286 * t57, -t1 * t117 - t118 * t2 + t15 * t71 + t190 * t72 + t209 * t288 - t285 * t9 + t286 * t8 - t289 * t61 + t192, -g(3) * t208 + t1 * t72 - t23 * t157 + t2 * t71 + t218 * t57 + t288 * t8 + t289 * t9 + t215 * (t157 * t176 + t180 * t182); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t277, -t114 ^ 2 + t116 ^ 2, -t56 - t278, -t277, t202 - t276, t110, -t81 * t252 - t116 * t80 - t283 + t51 + (-qJD(4) * t82 + t294 - t42) * t175 + t304, g(1) * t89 - g(2) * t87 + g(3) * t267 + t114 * t80 - t305, 0, 0, -t292, t310, t309, t292, t306, t106, t12 * t141 + (t106 * t178 - t116 * t61 + t141 * t251) * pkin(4) + t307, -t13 * t141 + (-t106 * t174 + t116 * t209 + t141 * t250) * pkin(4) + t308, -t12 * t209 + t13 * t61 - t209 * t9 - t61 * t8 + (t15 * t178 + t190 * t174 + (-t174 * t209 - t178 * t61) * qJD(5)) * pkin(4), -t8 * t12 - t9 * t13 + (t1 * t174 + t2 * t178 - t57 * t116 + g(3) * t269 + (-t174 * t8 + t178 * t9) * qJD(5) + t304) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t292, t310, t309, t292, t306, t106, -t141 * t9 + t307, -t141 * t8 + t308, 0, 0;];
tau_reg = t5;
