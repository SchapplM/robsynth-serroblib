% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:27
% EndTime: 2019-12-31 21:49:35
% DurationCPUTime: 3.31s
% Computational Cost: add. (2296->321), mult. (4726->379), div. (0->0), fcn. (3494->6), ass. (0->268)
t178 = sin(qJ(3));
t309 = cos(qJ(2));
t239 = t309 * t178;
t179 = sin(qJ(2));
t181 = cos(qJ(3));
t276 = t179 * t181;
t137 = (t239 + t276) * pkin(1);
t315 = -t137 / 0.2e1;
t177 = sin(qJ(4));
t180 = cos(qJ(4));
t209 = -t180 * pkin(4) - t177 * qJ(5);
t145 = -pkin(3) + t209;
t238 = t309 * t181;
t277 = t178 * t179;
t139 = (t238 - t277) * pkin(1);
t175 = t177 ^ 2;
t176 = t180 ^ 2;
t227 = t176 / 0.2e1 + t175 / 0.2e1;
t320 = t227 * t139;
t323 = -t320 * pkin(8) + t145 * t315;
t164 = t176 - t175;
t253 = -qJD(2) - qJD(3);
t224 = qJD(1) - t253;
t321 = t224 * t164;
t304 = pkin(2) * qJD(2);
t251 = t178 * t304;
t252 = t309 * pkin(1);
t215 = t252 + pkin(2);
t201 = t178 * t215;
t128 = pkin(1) * t276 + t201;
t81 = t178 * pkin(2);
t214 = t81 / 0.2e1 + t128 / 0.2e1;
t199 = t315 + t214;
t46 = t199 * t177;
t293 = t46 * qJD(1) + t177 * t251;
t47 = t199 * t180;
t292 = t47 * qJD(1) + t180 * t251;
t265 = t175 + t176;
t319 = pkin(3) / 0.2e1;
t150 = t181 * t215;
t127 = pkin(1) * t277 - t150;
t76 = t127 + t145;
t318 = t76 / 0.2e1;
t124 = -pkin(3) + t127;
t317 = -t124 / 0.2e1;
t305 = t181 * pkin(2);
t136 = t145 - t305;
t316 = t136 / 0.2e1;
t314 = t145 / 0.2e1;
t168 = -pkin(3) - t305;
t313 = -t168 / 0.2e1;
t310 = t177 / 0.2e1;
t308 = pkin(1) * t179;
t307 = pkin(3) * t180;
t306 = t177 * pkin(4);
t303 = pkin(2) * qJD(3);
t125 = pkin(8) + t128;
t55 = t265 * t127;
t3 = -t125 * t55 + t76 * t128;
t302 = t3 * qJD(1);
t58 = t265 * t139;
t4 = t125 * t58 + t76 * t137;
t301 = t4 * qJD(1);
t275 = t180 * qJ(5);
t146 = -t275 + t306;
t300 = t76 * t146;
t299 = t76 * t180;
t259 = qJD(1) * t180;
t102 = t128 * t259;
t298 = t47 * qJD(2) + t102;
t123 = t136 * t177;
t107 = t123 / 0.2e1;
t68 = t76 * t177;
t63 = t68 / 0.2e1;
t297 = t107 + t63;
t108 = -t123 / 0.2e1;
t64 = -t68 / 0.2e1;
t296 = t108 + t64;
t120 = t137 * t259;
t295 = -t47 * qJD(3) + t120;
t282 = t146 * t177;
t50 = t282 + t299;
t291 = qJD(1) * t50;
t281 = t146 * t180;
t51 = -t68 + t281;
t290 = qJD(1) * t51;
t289 = qJD(1) * t55;
t288 = qJD(1) * t58;
t287 = t136 * t180;
t285 = t145 * t180;
t284 = t146 * t136;
t283 = t146 * t145;
t280 = t168 * t180;
t279 = t177 * t127;
t278 = t177 * t139;
t274 = t180 * t127;
t273 = t180 * t139;
t270 = t81 * qJD(1);
t83 = t150 / 0.2e1 + (-t252 / 0.2e1 + pkin(2) / 0.2e1) * t181;
t269 = t83 * qJD(1);
t122 = t128 * qJD(3);
t100 = t177 * t122;
t129 = t137 * qJD(2);
t118 = t177 * t129;
t268 = t100 + t118;
t267 = -t129 - t122;
t248 = t178 * t303;
t163 = t177 * t248;
t173 = t175 * qJD(5);
t266 = t173 - t163;
t264 = qJD(1) * t127;
t263 = qJD(1) * t128;
t262 = qJD(1) * t137;
t261 = qJD(1) * t139;
t260 = qJD(1) * t177;
t258 = qJD(2) * t177;
t257 = qJD(3) * t177;
t256 = qJD(4) * qJ(5);
t255 = t177 * qJD(4);
t174 = t180 * qJD(4);
t254 = t180 * qJD(5);
t250 = pkin(8) * t255;
t249 = pkin(8) * t174;
t247 = -t305 / 0.2e1;
t246 = t305 / 0.2e1;
t245 = t319 + t317;
t244 = t319 + t313;
t243 = t76 * t260;
t242 = qJD(1) * t300;
t241 = t316 + t318;
t240 = t314 + t318;
t237 = t124 * t260;
t236 = t124 * t259;
t235 = t125 * t255;
t234 = t125 * t174;
t167 = pkin(8) + t81;
t233 = t167 * t255;
t232 = t167 * t174;
t229 = t314 + t316;
t228 = t313 + t317;
t226 = t309 * qJD(1);
t225 = t309 * qJD(2);
t223 = pkin(2) * t253;
t222 = t265 * t181;
t101 = t128 * t260;
t221 = -qJD(2) * t46 - t101;
t119 = t137 * t260;
t220 = qJD(3) * t46 - t119;
t140 = t145 * t177;
t132 = t140 / 0.2e1;
t219 = t132 - t281;
t133 = -t140 / 0.2e1;
t218 = t133 + t281;
t217 = -t285 / 0.2e1 - t282;
t216 = t180 * t248;
t45 = t137 * t310 + t214 * t177;
t213 = qJD(2) * t45 + t100 + t101;
t212 = qJD(3) * t45 + t118 + t119;
t211 = t227 * t127;
t210 = t227 * t181;
t182 = (t125 * t210 + t178 * t318) * pkin(2) - t167 * t211 + t128 * t316;
t2 = t182 + t323;
t49 = (t136 * t178 + t167 * t222) * pkin(2);
t208 = t2 * qJD(1) + t49 * qJD(2);
t112 = -t278 / 0.2e1;
t11 = t112 + t281 + t296;
t62 = -t123 + t281;
t207 = qJD(1) * t11 + qJD(2) * t62;
t114 = t273 / 0.2e1;
t12 = t180 * t241 + t114 + t282;
t61 = t282 + t287;
t206 = qJD(1) * t12 + qJD(2) * t61;
t10 = t265 * (t246 - t127 / 0.2e1 - t139 / 0.2e1);
t138 = pkin(2) * t222;
t205 = -qJD(1) * t10 - qJD(2) * t138;
t204 = qJD(4) * t146 - qJD(5) * t177;
t48 = (-t214 + t315) * t180;
t203 = t48 * qJD(3) - t129 * t180 - t120;
t202 = t48 * qJD(2) - t122 * t180 - t102;
t200 = t275 / 0.2e1 - t306 / 0.2e1;
t190 = t200 * t139;
t7 = -t146 * t241 + t190;
t197 = -t7 * qJD(1) + qJD(2) * t284;
t111 = t278 / 0.2e1;
t19 = t111 + t297;
t196 = qJD(1) * t19 + t136 * t258;
t29 = t177 * t228 + t112;
t195 = qJD(1) * t29 - t168 * t258;
t115 = -t273 / 0.2e1;
t30 = t180 * t228 + t115;
t194 = qJD(1) * t30 - qJD(2) * t280;
t193 = -t216 - t292;
t192 = t267 * t180;
t191 = t200 * t127;
t94 = t279 / 0.2e1;
t15 = t218 + t64 + t94;
t153 = t177 * t247;
t25 = t108 + t153 + t218;
t72 = -t140 + t281;
t189 = qJD(1) * t15 + qJD(2) * t25 + qJD(3) * t72;
t96 = -t274 / 0.2e1;
t16 = t180 * t240 + t282 + t96;
t158 = t180 * t246;
t26 = t180 * t229 + t158 + t282;
t71 = t282 + t285;
t188 = qJD(1) * t16 + qJD(2) * t26 + qJD(3) * t71;
t187 = t200 * t305;
t41 = t177 * t245 + t94;
t80 = t177 * t244 + t153;
t186 = pkin(3) * t257 + qJD(1) * t41 + qJD(2) * t80;
t97 = t274 / 0.2e1;
t42 = t180 * t245 + t97;
t159 = t180 * t247;
t82 = t180 * t244 + t159;
t185 = qJD(1) * t42 + qJD(2) * t82 + qJD(3) * t307;
t23 = -t146 * t229 + t187;
t5 = -t146 * t240 - t191;
t184 = -t5 * qJD(1) - t23 * qJD(2) + qJD(3) * t283;
t93 = -t279 / 0.2e1;
t21 = t93 + t132 + t63;
t152 = t177 * t246;
t53 = t152 + t132 + t107;
t183 = qJD(1) * t21 + qJD(2) * t53 + t145 * t257;
t126 = qJD(4) * t209 + t254;
t172 = -t307 / 0.2e1;
t171 = -pkin(3) * t177 / 0.2e1;
t166 = t177 * t174;
t165 = t177 * t254;
t148 = t164 * qJD(4);
t144 = t280 / 0.2e1;
t143 = t168 * t310;
t135 = t224 * t175;
t131 = t139 * qJD(2);
t130 = t138 * qJD(3);
t121 = t127 * qJD(3);
t109 = -t287 / 0.2e1;
t91 = t283 / 0.2e1;
t90 = t124 * t180 / 0.2e1;
t89 = t124 * t310;
t86 = t224 * t180 * t177;
t85 = t172 + t144 + t159;
t84 = t171 + t143 + t153;
t77 = t284 / 0.2e1;
t67 = t247 - t150 / 0.2e1 + (t277 - t238 / 0.2e1) * pkin(1);
t66 = -t81 / 0.2e1 - t201 / 0.2e1 + (-t276 - t239 / 0.2e1) * pkin(1);
t65 = -t299 / 0.2e1;
t57 = t300 / 0.2e1;
t56 = t58 * qJD(2);
t54 = t133 + t108 + t152;
t52 = t55 * qJD(3);
t44 = t172 + t90 + t97;
t43 = t171 + t89 + t94;
t32 = t144 + t90 + t115;
t31 = t143 + t89 + t112;
t28 = t107 + t153 + t219;
t27 = t109 + t158 + t217;
t24 = t91 + t77 + t187;
t22 = t133 + t64 + t93;
t20 = t111 + t296;
t18 = t219 + t63 + t94;
t17 = t217 + t65 + t96;
t14 = t112 - t281 + t297;
t13 = t109 + t114 + t65 - t282;
t9 = pkin(2) * t210 - t211 + t320;
t8 = t77 + t57 + t190;
t6 = t91 + t57 - t191;
t1 = t182 - t323;
t33 = [0, 0, 0, 0, -qJD(2) * t308, -pkin(1) * t225, 0, t267, -t131 + t121, t166, t148, 0, 0, 0, t124 * t255 + t192, t124 * t174 + t268, -t51 * qJD(4) + t165 + t192, t56 - t52, -qJD(4) * t50 + t173 - t268, t4 * qJD(2) + t3 * qJD(3) + t204 * t76; 0, 0, 0, 0, (-qJD(1) - qJD(2)) * t308, (-t226 - t225) * pkin(1), 0, qJD(3) * t66 - t129 - t262, qJD(3) * t67 - t131 - t261, t166, t148, 0, 0, 0, qJD(4) * t31 + t203, qJD(4) * t32 + t212, qJD(4) * t14 + t165 + t203, qJD(3) * t9 + t288 + t56, qJD(4) * t13 + t173 - t212, t301 + (t137 * t136 + t167 * t58) * qJD(2) + t1 * qJD(3) + t8 * qJD(4) + t20 * qJD(5); 0, 0, 0, 0, 0, 0, 0, qJD(2) * t66 - t122 - t263, qJD(2) * t67 + t121 + t264, t166, t148, 0, 0, 0, qJD(4) * t43 + t202, qJD(4) * t44 + t213, qJD(4) * t18 + t165 + t202, qJD(2) * t9 - t289 - t52, qJD(4) * t17 + t173 - t213, t302 + t1 * qJD(2) + (-pkin(8) * t55 + t128 * t145) * qJD(3) + t6 * qJD(4) + t22 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t321, t174, -t255, 0, qJD(2) * t31 + qJD(3) * t43 - t234 + t237, qJD(2) * t32 + qJD(3) * t44 + t235 + t236, qJD(2) * t14 + qJD(3) * t18 - t234 - t290, t126, qJD(2) * t13 + qJD(3) * t17 - t235 - t291, t8 * qJD(2) + t6 * qJD(3) + t125 * t126 + t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t174, t135, qJD(2) * t20 + qJD(3) * t22 + t234 - t243; 0, 0, 0, 0, qJD(1) * t308, pkin(1) * t226, 0, -qJD(3) * t81 + t262, -qJD(3) * t83 + t261, t166, t148, 0, 0, 0, -qJD(4) * t29 + t295, -qJD(4) * t30 + t220, -qJD(4) * t11 + t165 + t295, qJD(3) * t10 - t288, -qJD(4) * t12 + t173 - t220, qJD(3) * t2 - qJD(4) * t7 - qJD(5) * t19 - t301; 0, 0, 0, 0, 0, 0, 0, -t248, -t181 * t303, t166, t148, 0, 0, 0, t168 * t255 - t216, t168 * t174 + t163, -qJD(4) * t62 + t165 - t216, t130, -qJD(4) * t61 + t266, t49 * qJD(3) + t136 * t204; 0, 0, 0, 0, 0, 0, 0, t178 * t223 - t270, t181 * t223 - t269, t166, t148, 0, 0, 0, qJD(4) * t84 + t193, qJD(4) * t85 + t163 + t293, qJD(4) * t28 + t165 + t193, t130 - t205, qJD(4) * t27 + t266 - t293, t24 * qJD(4) + t54 * qJD(5) + (pkin(8) * t222 + t145 * t178) * t303 + t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t321, t174, -t255, 0, qJD(3) * t84 - t195 - t232, qJD(3) * t85 - t194 + t233, qJD(3) * t28 - t207 - t232, t126, qJD(3) * t27 - t206 - t233, t24 * qJD(3) + t126 * t167 + t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t174, t135, qJD(3) * t54 - t196 + t232; 0, 0, 0, 0, 0, 0, 0, qJD(2) * t81 + t263, qJD(2) * t83 - t264, t166, t148, 0, 0, 0, -qJD(4) * t41 + t298, -qJD(4) * t42 + t221, -qJD(4) * t15 + t165 + t298, -qJD(2) * t10 + t289, -qJD(4) * t16 + t173 - t221, -qJD(2) * t2 - qJD(4) * t5 - qJD(5) * t21 - t302; 0, 0, 0, 0, 0, 0, 0, t251 + t270, t181 * t304 + t269, t166, t148, 0, 0, 0, -qJD(4) * t80 + t292, -qJD(4) * t82 - t293, -qJD(4) * t25 + t165 + t292, t205, -qJD(4) * t26 + t173 + t293, -qJD(4) * t23 - qJD(5) * t53 - t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t148, 0, 0, 0, -pkin(3) * t255, -pkin(3) * t174, -qJD(4) * t72 + t165, 0, -qJD(4) * t71 + t173, t204 * t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t321, t174, -t255, 0, -t186 - t249, -t185 + t250, -t189 - t249, t126, -t188 - t250, pkin(8) * t126 + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t174, t135, -t183 + t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t321, 0, 0, 0, qJD(2) * t29 + qJD(3) * t41 - t237, qJD(2) * t30 + qJD(3) * t42 - t236, qJD(2) * t11 + qJD(3) * t15 + t290, 0, qJD(2) * t12 + qJD(3) * t16 + t291, qJD(2) * t7 + qJD(3) * t5 - t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t321, 0, 0, 0, qJD(3) * t80 + t195, qJD(3) * t82 + t194, qJD(3) * t25 + t207, 0, qJD(3) * t26 + t206, qJD(3) * t23 - t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t321, 0, 0, 0, t186, t185, t189, 0, t188, -t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, 0, -t135, qJD(2) * t19 + qJD(3) * t21 + t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, 0, -t135, qJD(3) * t53 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, 0, -t135, t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t33;
