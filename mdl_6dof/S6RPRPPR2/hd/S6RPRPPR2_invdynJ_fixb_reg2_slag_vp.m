% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:49
% EndTime: 2019-03-09 02:42:54
% DurationCPUTime: 4.28s
% Computational Cost: add. (5608->469), mult. (11904->556), div. (0->0), fcn. (8246->14), ass. (0->249)
t167 = sin(pkin(10));
t169 = cos(pkin(10));
t173 = sin(qJ(3));
t176 = cos(qJ(3));
t116 = t167 * t176 + t169 * t173;
t315 = t116 * qJD(1);
t324 = qJD(6) + t315;
t266 = t169 * t176;
t236 = qJD(1) * t266;
t259 = qJD(1) * t173;
t107 = t167 * t259 - t236;
t172 = sin(qJ(6));
t175 = cos(qJ(6));
t85 = qJD(3) * t172 - t175 * t107;
t232 = t324 * t85;
t255 = qJD(6) * t175;
t256 = qJD(6) * t172;
t109 = t116 * qJD(3);
t247 = t176 * qJDD(1);
t248 = t173 * qJDD(1);
t209 = t167 * t248 - t169 * t247;
t71 = qJD(1) * t109 + t209;
t28 = qJD(3) * t256 - t175 * qJDD(3) - t107 * t255 - t172 * t71;
t328 = t28 - t232;
t168 = sin(pkin(9));
t145 = pkin(1) * t168 + pkin(7);
t131 = t145 * qJD(1);
t226 = qJ(4) * qJD(1) + t131;
t254 = t173 * qJD(2);
t90 = t176 * t226 + t254;
t80 = t167 * t90;
t159 = t176 * qJD(2);
t89 = -t173 * t226 + t159;
t83 = qJD(3) * pkin(3) + t89;
t41 = t169 * t83 - t80;
t220 = qJD(5) - t41;
t297 = t315 * pkin(5);
t307 = pkin(4) + pkin(8);
t21 = -t307 * qJD(3) + t220 + t297;
t170 = cos(pkin(9));
t147 = -t170 * pkin(1) - pkin(2);
t160 = t176 * pkin(3);
t319 = t147 - t160;
t106 = qJD(1) * t319 + qJD(4);
t187 = -qJ(5) * t315 + t106;
t32 = t307 * t107 + t187;
t8 = -t172 * t32 + t175 * t21;
t332 = t8 * t324;
t251 = qJD(1) * qJD(3);
t235 = t173 * t251;
t191 = qJDD(1) * t116 - t167 * t235;
t234 = t176 * t251;
t72 = t169 * t234 + t191;
t285 = t72 * qJ(5);
t88 = pkin(3) * t235 + qJDD(1) * t319 + qJDD(4);
t181 = -qJD(5) * t315 - t285 + t88;
t10 = t307 * t71 + t181;
t158 = t176 * qJDD(2);
t129 = t145 * qJDD(1);
t224 = -qJD(2) * qJD(3) - t129;
t250 = qJD(1) * qJD(4);
t257 = qJD(3) * t176;
t39 = qJDD(3) * pkin(3) + t158 - t226 * t257 + (-qJ(4) * qJDD(1) + t224 - t250) * t173;
t240 = -t173 * qJDD(2) + t224 * t176;
t258 = qJD(3) * t173;
t63 = -t131 * t258 - t240;
t43 = t176 * t250 + (-t235 + t247) * qJ(4) + t63;
t294 = t167 * t43 - t169 * t39;
t241 = qJDD(5) + t294;
t6 = t72 * pkin(5) - t307 * qJDD(3) + t241;
t9 = t172 * t21 + t175 * t32;
t2 = -qJD(6) * t9 - t172 * t10 + t175 * t6;
t314 = t324 * t9 + t2;
t87 = qJD(3) * t175 + t107 * t172;
t29 = qJD(6) * t87 + t172 * qJDD(3) - t175 * t71;
t331 = t324 * t87 - t29;
t330 = 0.2e1 * qJD(3);
t321 = t175 * t324;
t329 = t87 * t321;
t229 = t172 * t324;
t70 = qJDD(6) + t72;
t62 = t175 * t70;
t201 = -t229 * t324 + t62;
t115 = t167 * t173 - t266;
t206 = -t109 * t315 - t115 * t72;
t112 = -t167 * t258 + t169 * t257;
t318 = -t112 * t107 - t116 * t71;
t327 = -t206 + t318;
t326 = t206 + t318;
t279 = t109 * t172;
t325 = t115 * t255 + t279;
t323 = t315 * t330 + t209;
t283 = pkin(1) * qJDD(1);
t46 = t169 * t89 - t80;
t263 = -qJD(5) + t46;
t162 = qJ(3) + pkin(10);
t152 = sin(t162);
t154 = cos(t162);
t320 = -t154 * pkin(4) - t152 * qJ(5);
t163 = qJ(1) + pkin(9);
t153 = sin(t163);
t155 = cos(t163);
t233 = -g(1) * t153 + g(2) * t155;
t300 = g(2) * t153;
t218 = g(1) * t155 + t300;
t105 = t315 ^ 2;
t308 = t107 ^ 2;
t317 = -t308 - t105;
t316 = -t308 + t105;
t146 = -pkin(3) * t169 - pkin(4);
t137 = -pkin(8) + t146;
t304 = pkin(5) * t107;
t287 = t169 * t90;
t42 = t167 * t83 + t287;
t36 = -qJD(3) * qJ(5) - t42;
t23 = -t36 - t304;
t313 = t137 * t70 + t23 * t324;
t276 = t115 * t172;
t312 = -t70 * t276 - t324 * t325;
t278 = t109 * t175;
t286 = t175 * t28;
t311 = -t115 * (t256 * t87 + t286) + t87 * t278;
t264 = qJ(4) + t145;
t113 = t264 * t173;
t114 = t264 * t176;
t68 = -t113 * t167 + t114 * t169;
t310 = -t68 * qJDD(3) + t152 * t233;
t189 = -g(3) * t152 - t218 * t154;
t309 = qJD(3) * (-t107 + t236) + t191;
t306 = t71 * pkin(4);
t305 = pkin(3) * t173;
t303 = pkin(8) * t154;
t142 = g(3) * t154;
t298 = g(3) * t176;
t174 = sin(qJ(1));
t296 = t174 * pkin(1);
t295 = t87 * t85;
t16 = t167 * t39 + t169 * t43;
t291 = t87 * t112 - t28 * t116;
t290 = t107 * t85;
t27 = t175 * t29;
t284 = t87 * t107;
t281 = qJDD(3) * pkin(4);
t280 = t107 * t315;
t275 = t131 * t173;
t274 = t131 * t176;
t273 = t152 * t155;
t272 = t153 * t154;
t271 = t153 * t172;
t270 = t153 * t175;
t269 = t154 * t155;
t268 = t155 * t172;
t267 = t155 * t175;
t45 = t167 * t89 + t287;
t265 = t45 * qJD(3);
t262 = t297 - t263;
t149 = t160 + pkin(2);
t177 = cos(qJ(1));
t161 = t177 * pkin(1);
t261 = t155 * t149 + t161;
t165 = t173 ^ 2;
t166 = t176 ^ 2;
t260 = t165 - t166;
t132 = qJD(1) * t147;
t130 = qJDD(1) * t147;
t246 = -t29 * t276 - t325 * t85;
t245 = qJDD(3) * qJ(5) + t16;
t151 = pkin(3) * t258;
t179 = qJD(1) ^ 2;
t242 = t173 * t179 * t176;
t239 = t160 - t320;
t238 = -g(1) * t273 - t152 * t300 + t142;
t230 = qJD(3) * t264;
t91 = t176 * qJD(4) - t173 * t230;
t92 = -t173 * qJD(4) - t176 * t230;
t50 = t167 * t91 - t169 * t92;
t231 = pkin(3) * t259 + t107 * qJ(5);
t67 = t169 * t113 + t114 * t167;
t223 = pkin(4) * t269 + qJ(5) * t273 + t261;
t222 = t173 * t234;
t221 = g(1) * t272 - g(2) * t269;
t219 = -pkin(4) * t152 - t305;
t216 = g(1) * t174 - g(2) * t177;
t215 = t172 * t9 + t175 * t8;
t214 = -t172 * t8 + t175 * t9;
t1 = qJD(6) * t8 + t175 * t10 + t172 * t6;
t213 = t1 - t332;
t212 = -t238 - t294;
t171 = -qJ(4) - pkin(7);
t210 = -t155 * t171 - t296;
t208 = -t112 * t85 - t116 * t29;
t51 = t167 * t92 + t169 * t91;
t197 = -t116 * qJ(5) + t319;
t47 = t307 * t115 + t197;
t53 = pkin(5) * t116 + t67;
t19 = t172 * t53 + t175 * t47;
t18 = -t172 * t47 + t175 * t53;
t207 = t107 * t109 + t115 * t71;
t205 = t112 * t315 + t116 * t72;
t100 = t254 + t274;
t75 = qJD(3) * t109 + qJDD(3) * t115;
t74 = qJD(3) * t112 + qJDD(3) * t116;
t13 = -qJD(3) * qJD(5) - t245;
t203 = -t149 + t320;
t199 = -t112 * qJ(5) - t116 * qJD(5) + t151;
t198 = t115 * t62 + (-t115 * t256 + t278) * t324;
t195 = -t67 * qJDD(3) + t221;
t194 = -qJD(1) * t132 + t218;
t193 = -t172 * t70 - t321 * t324;
t52 = t107 * pkin(4) + t187;
t192 = t315 * t52 + qJDD(5) - t212;
t190 = -qJDD(3) * t145 + t132 * t330;
t178 = qJD(3) ^ 2;
t188 = -t145 * t178 - 0.2e1 * t130 - t233;
t186 = qJD(6) * t214 + t1 * t172 + t2 * t175;
t64 = -t100 * qJD(3) - t173 * t129 + t158;
t99 = t159 - t275;
t185 = -t64 * t173 + t63 * t176 + (-t100 * t173 - t176 * t99) * qJD(3);
t184 = t233 + t88;
t7 = -pkin(5) * t71 - t13;
t183 = -qJD(6) * t137 * t324 + t189 + t7;
t182 = -t107 * t51 + t315 * t50 + t67 * t72 - t68 * t71 - t218;
t48 = (t107 + t236) * qJD(3) + t191;
t140 = pkin(3) * t167 + qJ(5);
t128 = qJDD(3) * t176 - t173 * t178;
t127 = qJDD(3) * t173 + t176 * t178;
t120 = qJ(5) * t269;
t118 = qJ(5) * t272;
t96 = -t152 * t271 + t267;
t95 = t152 * t270 + t268;
t94 = t152 * t268 + t270;
t93 = t152 * t267 - t271;
t61 = pkin(4) * t115 + t197;
t59 = pkin(4) * t315 + t231;
t54 = -pkin(5) * t115 + t68;
t44 = pkin(4) * t109 + t199;
t40 = t307 * t315 + t231;
t35 = -qJD(3) * pkin(4) + t220;
t31 = -pkin(5) * t109 + t51;
t30 = pkin(5) * t112 + t50;
t25 = t45 - t304;
t24 = t307 * t109 + t199;
t17 = t181 + t306;
t14 = t241 - t281;
t12 = t172 * t25 + t175 * t40;
t11 = -t172 * t40 + t175 * t25;
t4 = -qJD(6) * t19 - t172 * t24 + t175 * t30;
t3 = qJD(6) * t18 + t172 * t30 + t175 * t24;
t5 = [0, 0, 0, 0, 0, qJDD(1), t216, g(1) * t177 + g(2) * t174, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t170 * t283 - t233, -0.2e1 * t168 * t283 + t218, 0 (t216 + (t168 ^ 2 + t170 ^ 2) * t283) * pkin(1), qJDD(1) * t165 + 0.2e1 * t222, 0.2e1 * t173 * t247 - 0.2e1 * t251 * t260, t127, qJDD(1) * t166 - 0.2e1 * t222, t128, 0, t173 * t190 + t176 * t188, -t173 * t188 + t176 * t190 (t165 + t166) * t129 + t185 - t218, t130 * t147 - g(1) * (-pkin(2) * t153 + pkin(7) * t155 - t296) - g(2) * (pkin(2) * t155 + pkin(7) * t153 + t161) + t185 * t145, t205, t326, t74, t207, -t75, 0, t106 * t109 + t88 * t115 + t319 * t71 + (t107 * t305 - t50) * qJD(3) + t195, t106 * t112 + t88 * t116 + t319 * t72 + (t305 * t315 - t51) * qJD(3) + t310, -t109 * t42 - t112 * t41 - t115 * t16 + t116 * t294 + t182, t16 * t68 + t42 * t51 + t294 * t67 - t41 * t50 + t88 * t319 + t106 * t151 - g(1) * (-t153 * t149 + t210) - g(2) * (-t153 * t171 + t261) 0, -t74, t75, t205, t326, t207, t109 * t36 + t112 * t35 + t115 * t13 + t116 * t14 + t182, qJD(3) * t50 - t107 * t44 - t109 * t52 - t115 * t17 - t61 * t71 - t195, t51 * qJD(3) - t52 * t112 - t17 * t116 - t315 * t44 - t61 * t72 - t310, t17 * t61 + t52 * t44 - t13 * t68 - t36 * t51 + t14 * t67 + t35 * t50 - g(1) * t210 - g(2) * t223 + (-g(1) * t203 + g(2) * t171) * t153, t87 * t279 + (-t172 * t28 + t255 * t87) * t115, t246 + t311, t291 - t312, -t85 * t278 + (t256 * t85 - t27) * t115, t198 + t208, t112 * t324 + t116 * t70, -t23 * t278 - g(1) * t96 - g(2) * t94 + t4 * t324 + t8 * t112 + t2 * t116 + t18 * t70 + t54 * t29 + t31 * t85 + (-t7 * t175 + t23 * t256) * t115, t23 * t279 + g(1) * t95 - g(2) * t93 - t1 * t116 - t3 * t324 - t9 * t112 - t19 * t70 - t54 * t28 + t31 * t87 + (t7 * t172 + t23 * t255) * t115, t18 * t28 - t19 * t29 - t3 * t85 - t4 * t87 + t214 * t109 + (-qJD(6) * t215 + t1 * t175 - t2 * t172) * t115 + t221, t1 * t19 + t9 * t3 + t2 * t18 + t8 * t4 + t7 * t54 + t23 * t31 - g(1) * (t155 * pkin(5) + t210) - g(2) * (pkin(8) * t269 + t223) + (-g(1) * (t203 - t303) - g(2) * (pkin(5) - t171)) * t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t128, -t127, 0, t63 * t173 + t64 * t176 - g(3) + (t100 * t176 - t173 * t99) * qJD(3), 0, 0, 0, 0, 0, 0, -t75, -t74, t327, -t109 * t41 + t112 * t42 + t115 * t294 + t116 * t16 - g(3), 0, 0, 0, 0, 0, 0, t327, t75, t74, t109 * t35 - t112 * t36 + t115 * t14 - t116 * t13 - g(3), 0, 0, 0, 0, 0, 0, t198 - t208, t291 + t312, t246 - t311, t109 * t215 + t23 * t112 + t115 * t186 + t7 * t116 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t242, t260 * t179, t248, t242, t247, qJDD(3), -t298 + t158 + (t100 - t274) * qJD(3) + (t194 + t224) * t173, g(3) * t173 + (t99 + t275) * qJD(3) + t194 * t176 + t240, 0, 0, t280, t316, t48, -t280, -t209, qJDD(3), t265 - t106 * t315 + (qJDD(3) * t169 - t107 * t259) * pkin(3) + t212, t46 * qJD(3) + t106 * t107 + (-qJDD(3) * t167 - t259 * t315) * pkin(3) - t16 - t189 (t42 - t45) * t315 + (-t41 + t46) * t107 + (-t167 * t71 - t169 * t72) * pkin(3), t41 * t45 - t42 * t46 + (-t298 - t294 * t169 + t16 * t167 + (-qJD(1) * t106 + t218) * t173) * pkin(3), qJDD(3), -t48, t209, t280, t316, -t280, -t140 * t71 + t146 * t72 + (-t36 - t45) * t315 + (t35 + t263) * t107, -t265 + t59 * t107 + (-pkin(4) + t146) * qJDD(3) + t192, t140 * qJDD(3) - t52 * t107 + t59 * t315 + (0.2e1 * qJD(5) - t46) * qJD(3) + t189 + t245, -t13 * t140 + t14 * t146 - t52 * t59 - t35 * t45 - g(1) * (t155 * t219 + t120) - g(2) * (t153 * t219 + t118) - g(3) * t239 + t263 * t36, -t229 * t87 - t286, -t27 - t329 + (t28 + t232) * t172, t201 + t284, t29 * t172 + t321 * t85, t193 - t290, t324 * t107, t8 * t107 - t11 * t324 + t140 * t29 + t183 * t172 + t175 * t313 + t262 * t85, -t9 * t107 + t12 * t324 - t140 * t28 - t172 * t313 + t183 * t175 + t262 * t87, t11 * t87 + t12 * t85 + (-t315 * t9 + t137 * t28 - t2 + (-t137 * t85 - t9) * qJD(6)) * t175 + (t315 * t8 - t137 * t29 - t1 + (t137 * t87 + t8) * qJD(6)) * t172 - t238, t7 * t140 - t9 * t12 - t8 * t11 - g(1) * (-t155 * t305 + t120) - g(2) * (-t153 * t305 + t118) - g(3) * (t239 + t303) + t262 * t23 + t186 * t137 + t218 * t152 * t307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t323, t309, t317, t42 * t107 + t315 * t41 + t184, 0, 0, 0, 0, 0, 0, t317, -t323, -t309, t306 - t285 - t36 * t107 + (-qJD(5) - t35) * t315 + t184, 0, 0, 0, 0, 0, 0, t193 + t290, t284 - t201, -t328 * t172 - t27 + t329, t23 * t107 - t172 * t314 + t213 * t175 + t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, qJDD(3) - t280, -t105 - t178, qJD(3) * t36 + t192 - t281, 0, 0, 0, 0, 0, 0, -qJD(3) * t85 + t201, -qJD(3) * t87 + t193, t331 * t172 + t328 * t175, -t23 * qJD(3) + t213 * t172 + t175 * t314 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, -t85 ^ 2 + t87 ^ 2, -t328, -t295, t331, t70, -g(1) * t93 - g(2) * t95 + t175 * t142 - t23 * t87 + t314, g(1) * t94 - g(2) * t96 + t332 + t23 * t85 + (-qJD(6) * t21 - t10) * t175 + (qJD(6) * t32 - t142 - t6) * t172, 0, 0;];
tau_reg  = t5;
