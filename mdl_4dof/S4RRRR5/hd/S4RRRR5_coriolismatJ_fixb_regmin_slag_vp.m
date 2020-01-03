% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:15
% EndTime: 2019-12-31 17:28:21
% DurationCPUTime: 2.58s
% Computational Cost: add. (1921->265), mult. (4692->412), div. (0->0), fcn. (4559->6), ass. (0->234)
t253 = qJD(3) + qJD(4);
t180 = sin(qJ(2));
t179 = sin(qJ(3));
t309 = cos(qJ(4));
t238 = t309 * t179;
t178 = sin(qJ(4));
t181 = cos(qJ(3));
t286 = t178 * t181;
t202 = t238 + t286;
t115 = t202 * t180;
t314 = pkin(6) + pkin(7);
t146 = t314 * t179;
t147 = t314 * t181;
t102 = -t178 * t146 + t309 * t147;
t318 = t253 * t102;
t203 = -t309 * t146 - t178 * t147;
t317 = t253 * t203;
t312 = -t202 / 0.2e1;
t237 = t309 * t181;
t287 = t178 * t179;
t201 = t237 - t287;
t308 = pkin(3) * t179;
t239 = pkin(5) + t308;
t215 = t239 * t180;
t209 = -t215 / 0.2e1;
t303 = t181 * pkin(3);
t173 = -pkin(2) - t303;
t288 = t173 * t115;
t316 = -t201 * t209 - t288 / 0.2e1;
t174 = t179 ^ 2;
t176 = t181 ^ 2;
t158 = t176 - t174;
t257 = t180 * qJD(1);
t236 = t181 * t257;
t315 = t158 * qJD(2) - 0.2e1 * t179 * t236;
t313 = -t201 / 0.2e1;
t182 = cos(qJ(2));
t311 = -t182 / 0.2e1;
t310 = t182 / 0.2e1;
t307 = pkin(5) * t179;
t306 = t178 * pkin(3);
t305 = t180 * pkin(2);
t304 = t180 * pkin(3);
t302 = t182 * pkin(6);
t285 = t179 * t180;
t114 = t178 * t285 - t180 * t237;
t29 = t114 * t202 - t115 * t201;
t301 = t253 * t29;
t40 = t114 * t313 + t115 * t312;
t300 = t253 * t40;
t117 = t202 * t182;
t139 = t239 * t182;
t148 = -t302 + t305;
t138 = t181 * t148;
t164 = pkin(5) * t285;
t281 = t181 * t182;
t82 = -pkin(7) * t281 + t138 + t164 + t304;
t243 = t309 * t82;
t137 = t179 * t148;
t283 = t180 * t181;
t251 = pkin(5) * t283;
t284 = t179 * t182;
t96 = -pkin(7) * t284 + t137 - t251;
t295 = t178 * t96;
t214 = -t182 * pkin(2) - t180 * pkin(6);
t141 = -pkin(1) + t214;
t130 = t181 * t141;
t217 = -pkin(7) * t283 + t130;
t77 = (-pkin(3) - t307) * t182 + t217;
t244 = t309 * t77;
t250 = pkin(5) * t281;
t110 = t179 * t141 + t250;
t93 = -pkin(7) * t285 + t110;
t296 = t178 * t93;
t36 = -t244 + t296;
t1 = (t243 - t295) * t182 + t36 * t180 - t139 * t115 - t117 * t215;
t299 = t1 * qJD(1);
t298 = t178 * t82;
t252 = pkin(5) * t284;
t92 = t217 - t252;
t297 = t178 * t92;
t118 = t201 * t182;
t240 = t309 * t96;
t241 = t309 * t93;
t37 = t178 * t77 + t241;
t2 = (t240 + t298) * t182 - t139 * t114 + (t239 * t118 - t37) * t180;
t294 = t2 * qJD(1);
t293 = t203 * t182;
t292 = t102 * t182;
t218 = pkin(3) * t310 - t77 / 0.2e1;
t15 = (t92 / 0.2e1 + t218) * t178;
t291 = t15 * qJD(1);
t222 = t309 * pkin(3) / 0.2e1;
t187 = -t244 / 0.2e1 + t182 * t222;
t242 = t309 * t92;
t17 = t242 / 0.2e1 + t187;
t290 = t17 * qJD(1);
t289 = t173 * t114;
t175 = t180 ^ 2;
t282 = t181 * t175;
t41 = -t241 - t297;
t88 = t114 * t215;
t19 = -pkin(3) * t115 * t283 + t41 * t182 + t88;
t280 = t19 * qJD(1);
t42 = t242 - t296;
t20 = t42 * t182 + (-t114 * t303 - t239 * t115) * t180;
t279 = t20 * qJD(1);
t25 = -t115 * t215 - t36 * t182;
t278 = t25 * qJD(1);
t26 = -t37 * t182 + t88;
t277 = t26 * qJD(1);
t32 = t114 * t117 - t118 * t115;
t276 = t32 * qJD(1);
t109 = -t130 + t252;
t45 = t109 * t180 + (-t164 + t138) * t182;
t275 = t45 * qJD(1);
t46 = t137 * t182 + (-t110 + t250) * t180;
t274 = t46 * qJD(1);
t48 = t115 * t180 - t117 * t182;
t273 = t48 * qJD(1);
t49 = t114 * t180 + t118 * t182;
t272 = t49 * qJD(1);
t193 = -t286 / 0.2e1 - t238 / 0.2e1;
t59 = (t312 + t193) * t182;
t53 = t59 * qJD(1);
t192 = t237 / 0.2e1 - t287 / 0.2e1;
t60 = (t201 / 0.2e1 + t192) * t182;
t54 = t60 * qJD(1);
t78 = -t109 * t182 - t175 * t307;
t271 = t78 * qJD(1);
t79 = -pkin(5) * t282 - t110 * t182;
t270 = t79 * qJD(1);
t177 = t182 ^ 2;
t159 = t177 - t175;
t269 = qJD(1) * t114;
t268 = qJD(2) * t202;
t267 = qJD(2) * t173;
t266 = qJD(2) * t179;
t265 = qJD(2) * t181;
t264 = qJD(3) * t179;
t263 = qJD(3) * t181;
t262 = qJD(3) * t182;
t261 = qJD(4) * t173;
t135 = t159 * t179;
t260 = t135 * qJD(1);
t136 = t181 * t177 - t282;
t259 = t136 * qJD(1);
t258 = t159 * qJD(1);
t256 = t180 * qJD(2);
t255 = t182 * qJD(1);
t254 = t182 * qJD(2);
t249 = pkin(1) * t257;
t248 = pkin(1) * t255;
t247 = t308 / 0.2e1;
t246 = -t303 / 0.2e1;
t245 = t303 / 0.2e1;
t235 = t179 * t262;
t234 = t181 * t262;
t233 = t179 * t263;
t232 = t179 * t265;
t231 = t180 * t254;
t230 = t180 * t255;
t229 = t181 * t256;
t228 = -t292 / 0.2e1;
t227 = t309 * qJD(3);
t226 = t309 * qJD(4);
t225 = t253 * t115;
t224 = t253 * t202;
t223 = -qJD(3) + t255;
t220 = t179 * t229;
t219 = -t240 / 0.2e1;
t216 = t255 - qJD(3) / 0.2e1;
t183 = t114 * t247 + t203 * t311;
t197 = -t298 / 0.2e1 + t219;
t185 = t288 / 0.2e1 + t197;
t3 = (-t306 / 0.2e1 - t202 * t245 + t239 * t313) * t180 + t183 + t185;
t56 = t173 * t201 + t202 * t308;
t213 = t3 * qJD(1) - t56 * qJD(2);
t184 = -t115 * t247 + t228;
t196 = -t295 / 0.2e1 + t243 / 0.2e1;
t186 = t289 / 0.2e1 + t196;
t4 = (-t201 * t246 + t239 * t312 + t222) * t180 + t184 + t186;
t55 = t173 * t202 - t201 * t308;
t212 = t4 * qJD(1) - t55 * qJD(2);
t211 = -qJD(4) + t223;
t47 = -t114 ^ 2 + t115 ^ 2;
t7 = t47 * qJD(1) + t29 * qJD(2);
t50 = t201 ^ 2 - t202 ^ 2;
t13 = t29 * qJD(1) + t50 * qJD(2);
t210 = t223 * t180;
t208 = t215 / 0.2e1;
t207 = t302 / 0.2e1 - t305 / 0.2e1;
t194 = t207 * t181;
t91 = -t138 / 0.2e1 + t194;
t206 = pkin(2) * t266 - t91 * qJD(1);
t195 = t207 * t179;
t90 = t137 / 0.2e1 - t195;
t205 = pkin(2) * t265 - t90 * qJD(1);
t9 = -t202 * t208 + t186 + t228;
t204 = t9 * qJD(1) - t202 * t267;
t10 = -t293 / 0.2e1 - t201 * t208 + t185;
t200 = t10 * qJD(1) - t201 * t267;
t22 = -t40 * qJD(2) - t115 * t269;
t30 = t40 * qJD(1) + t201 * t268;
t199 = t181 * t210;
t122 = (t174 / 0.2e1 - t176 / 0.2e1) * t180;
t198 = -t122 * qJD(1) + t232;
t191 = t179 * qJD(1) * t282 + t122 * qJD(2);
t134 = t158 * t175;
t190 = t134 * qJD(1) + 0.2e1 * t220;
t188 = -t202 * t209 + t196 - t289 / 0.2e1;
t168 = -t257 / 0.2e1;
t167 = t257 / 0.2e1;
t166 = t256 / 0.2e1;
t129 = t216 * t180;
t121 = t122 * qJD(3);
t113 = (-qJD(4) / 0.2e1 + t216) * t180;
t64 = t164 + t138 / 0.2e1 + t194;
t63 = t251 - t137 / 0.2e1 - t195;
t62 = t193 * t182 - t202 * t311;
t61 = t192 * t182 - t201 * t310;
t44 = t60 * qJD(2) - t115 * t255;
t43 = t59 * qJD(2) + t114 * t255;
t35 = -t224 - t53;
t34 = t201 * t253 - t54;
t24 = t61 * qJD(2) + t211 * t115;
t23 = t62 * qJD(2) - t211 * t114;
t18 = t296 - t242 / 0.2e1 + t187;
t16 = -t241 - t297 / 0.2e1 + t218 * t178;
t12 = t292 / 0.2e1 + t188;
t11 = t293 / 0.2e1 + t197 + t316;
t6 = -t246 * t115 + t219 + (-t304 / 0.2e1 - t82 / 0.2e1) * t178 - t183 + t316;
t5 = -t184 + t188 + (-t201 * t245 + t222) * t180;
t8 = [0, 0, 0, t231, t159 * qJD(2), 0, 0, 0, -pkin(1) * t256, -pkin(1) * t254, -t175 * t233 + t176 * t231, -t134 * qJD(3) - 0.2e1 * t182 * t220, -t136 * qJD(2) + t180 * t235, t135 * qJD(2) + t180 * t234, -t231, -t45 * qJD(2) - t79 * qJD(3), t46 * qJD(2) + t78 * qJD(3), (-qJD(2) * t118 + t225) * t114, t32 * qJD(2) + t253 * t47, -t49 * qJD(2) + t182 * t225, -t114 * t182 * t253 - t48 * qJD(2), -t231, -t1 * qJD(2) - t19 * qJD(3) - t26 * qJD(4), t2 * qJD(2) + t20 * qJD(3) + t25 * qJD(4); 0, 0, 0, t230, t258, t254, -t256, 0, -pkin(5) * t254 - t249, pkin(5) * t256 - t248, -t121 + (t176 * t257 + t232) * t182, -0.2e1 * t180 * t233 + t315 * t182, t179 * t256 - t259, t229 + t260, -t129, -t275 + (t179 * t214 - t250) * qJD(2) + t64 * qJD(3), t274 + (t181 * t214 + t252) * qJD(2) + t63 * qJD(3), (t268 - t269) * t118 + t300, t276 + (-t117 * t202 + t118 * t201) * qJD(2) + t301, t202 * t256 + t253 * t61 - t272, t201 * t256 + t253 * t62 - t273, -t113, -t299 + (t173 * t117 - t139 * t201 + t180 * t203) * qJD(2) + t5 * qJD(3) + t12 * qJD(4), t294 + (-t102 * t180 + t173 * t118 + t139 * t202) * qJD(2) + t6 * qJD(3) + t11 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, -t190, t179 * t210, t199, t166, t64 * qJD(2) - t110 * qJD(3) - t270, t63 * qJD(2) + t109 * qJD(3) + t271, -t22, t7, t24, t23, t166, t5 * qJD(2) + t41 * qJD(3) + t16 * qJD(4) - t280, t6 * qJD(2) - t42 * qJD(3) + t18 * qJD(4) + t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t7, t24, t23, t166, t12 * qJD(2) + t16 * qJD(3) - t37 * qJD(4) - t277, t11 * qJD(2) + t18 * qJD(3) + t36 * qJD(4) + t278; 0, 0, 0, -t230, -t258, 0, 0, 0, t249, t248, -t176 * t230 - t121, 0.2e1 * t179 * t199, -t234 + t259, t235 - t260, t129, t91 * qJD(3) + t275, t90 * qJD(3) - t274, t118 * t269 + t300, -t276 + t301, -t253 * t60 + t272, -t253 * t59 + t273, t113, -t4 * qJD(3) - t9 * qJD(4) + t299, -t3 * qJD(3) - t10 * qJD(4) - t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, t158 * qJD(3), 0, 0, 0, -pkin(2) * t264, -pkin(2) * t263, t201 * t224, t253 * t50, 0, 0, 0, t55 * qJD(3) + t202 * t261, t56 * qJD(3) + t201 * t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t315, -t223 * t181, t223 * t179, t168, -pkin(6) * t263 - t206, pkin(6) * t264 - t205, t30, t13, t34, t35, t168, -t212 - t318, -t213 - t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t13, t34, t35, t168, -t204 - t318, -t200 - t317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, t190, (-t179 * t257 + t265) * t182, (-t236 - t266) * t182, t166, -t91 * qJD(2) + t270, -t90 * qJD(2) - t271, t22, -t7, t44, t43, t166, t4 * qJD(2) + t15 * qJD(4) + t280, t3 * qJD(2) + t17 * qJD(4) - t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, -t315, t181 * t255, -t179 * t255, t167, t206, t205, -t30, -t13, t54, t53, t167, t212, t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t306, -pkin(3) * t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t253 * t306 + t291, t290 + (-t227 - t226) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t7, t44, t43, t166, t9 * qJD(2) - t15 * qJD(3) + t277, t10 * qJD(2) - t17 * qJD(3) - t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t13, t54, t53, t167, t204, t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t306 - t291, pkin(3) * t227 - t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t8;
