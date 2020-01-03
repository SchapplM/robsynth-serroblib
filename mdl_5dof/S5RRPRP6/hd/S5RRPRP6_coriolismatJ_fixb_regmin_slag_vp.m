% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:29
% EndTime: 2019-12-31 19:58:36
% DurationCPUTime: 2.38s
% Computational Cost: add. (3711->231), mult. (7199->354), div. (0->0), fcn. (7622->6), ass. (0->212)
t160 = sin(qJ(2));
t275 = -qJ(3) - pkin(6);
t145 = t275 * t160;
t162 = cos(qJ(2));
t146 = t275 * t162;
t158 = sin(pkin(8));
t259 = cos(pkin(8));
t103 = -t259 * t145 - t158 * t146;
t161 = cos(qJ(4));
t290 = t103 * t161;
t291 = t290 / 0.2e1;
t159 = sin(qJ(4));
t143 = t259 * t146;
t248 = t158 * t145;
t286 = -t143 + t248;
t258 = t286 * t159;
t257 = t286 * t161;
t190 = t259 * t160;
t247 = t158 * t162;
t140 = t190 + t247;
t81 = t159 * t140;
t189 = 0.2e1 * t161 * t81;
t138 = t158 * t160 - t259 * t162;
t278 = t138 * pkin(4);
t250 = t140 * t161;
t200 = -t162 * pkin(2) - pkin(1);
t76 = t138 * pkin(3) - t140 * pkin(7) + t200;
t50 = -t161 * t76 + t258;
t45 = -qJ(5) * t250 - t50;
t36 = t45 + t278;
t289 = -t36 + t45;
t153 = t158 * pkin(2) + pkin(7);
t226 = qJ(5) + t153;
t133 = t226 * t161;
t203 = t45 / 0.2e1 - t36 / 0.2e1;
t288 = t203 * t133;
t136 = t138 ^ 2;
t137 = t140 ^ 2;
t287 = -t137 - t136;
t132 = t190 / 0.2e1 + t247 / 0.2e1;
t246 = t103 * t159;
t277 = t160 * pkin(2);
t84 = t140 * pkin(3) + t138 * pkin(7) + t277;
t69 = t161 * t84;
t83 = t161 * t138;
t37 = t140 * pkin(4) + qJ(5) * t83 + t246 + t69;
t284 = t37 / 0.2e1;
t282 = -t140 / 0.2e1;
t157 = t161 ^ 2;
t281 = -t157 / 0.2e1;
t280 = t159 / 0.2e1;
t279 = pkin(4) * t159;
t276 = t161 * pkin(4);
t51 = t159 * t76 + t257;
t46 = -qJ(5) * t81 + t51;
t266 = t46 * t161;
t269 = t36 * t159;
t274 = -t269 / 0.2e1 + t266 / 0.2e1;
t68 = t159 * t84;
t273 = -t68 + t290;
t272 = pkin(4) * qJD(4);
t271 = qJD(2) * pkin(2);
t78 = t159 * t138;
t47 = qJ(5) * t78 - t273;
t58 = -pkin(4) * t78 + t286;
t59 = pkin(4) * t81 + t103;
t3 = t36 * t37 + t46 * t47 + t59 * t58;
t270 = t3 * qJD(1);
t268 = t37 * t161;
t204 = t278 / 0.2e1;
t179 = t204 - t203;
t4 = t179 * t161;
t267 = t4 * qJD(1);
t265 = t47 * t159;
t264 = t59 * t159;
t184 = t46 * t159 + t36 * t161;
t6 = (t265 + t268) * t140 - t184 * t138;
t263 = t6 * qJD(1);
t7 = t289 * t81;
t262 = t7 * qJD(1);
t208 = pkin(4) * t250;
t8 = t59 * t208 + t289 * t46;
t261 = t8 * qJD(1);
t9 = t179 * t159;
t260 = t9 * qJD(1);
t131 = t226 * t159;
t255 = t131 * t159;
t254 = t131 * t161;
t253 = t133 * t159;
t252 = t133 * t161;
t154 = -t259 * pkin(2) - pkin(3);
t144 = t154 - t276;
t191 = t140 * t144 / 0.2e1;
t164 = (-t252 / 0.2e1 - t255 / 0.2e1) * t138 + t191;
t178 = -t268 / 0.2e1 - t265 / 0.2e1;
t14 = t164 + t178;
t251 = t14 * qJD(1);
t15 = t59 * t140 + (-t266 + t269) * t138;
t249 = t15 * qJD(1);
t16 = (-t50 + t258) * t140 + t69 * t138;
t245 = t16 * qJD(1);
t17 = (-t51 + t257) * t140 + (-t290 + t273) * t138;
t244 = t17 * qJD(1);
t18 = t184 * t140;
t243 = t18 * qJD(1);
t24 = -t103 * t81 + t50 * t138;
t241 = t24 * qJD(1);
t25 = t103 * t250 - t51 * t138;
t240 = t25 * qJD(1);
t26 = t200 * t277;
t239 = t26 * qJD(1);
t48 = t103 * t140 - t138 * t286;
t238 = t48 * qJD(1);
t156 = t159 ^ 2;
t186 = t140 * t281 + t156 * t282;
t54 = t186 - t132;
t237 = t54 * qJD(1);
t209 = t137 - t136;
t55 = t209 * t159;
t236 = t55 * qJD(1);
t56 = t287 * t159;
t235 = t56 * qJD(1);
t57 = t209 * t161;
t234 = t57 * qJD(1);
t168 = -t158 * t138 / 0.2e1 + t259 * t282;
t62 = (-t160 / 0.2e1 + t168) * pkin(2);
t233 = t62 * qJD(1);
t232 = t78 * qJD(1);
t231 = t81 * qJD(1);
t230 = t83 * qJD(1);
t129 = t156 * t138;
t130 = t157 * t138;
t85 = t129 + t130;
t229 = t85 * qJD(1);
t149 = t156 + t157;
t86 = t149 * t137;
t228 = t86 * qJD(1);
t88 = t287 * t161;
t227 = t88 * qJD(1);
t150 = t157 - t156;
t224 = qJD(1) * t138;
t223 = qJD(1) * t140;
t222 = qJD(1) * t161;
t221 = qJD(1) * t162;
t220 = qJD(2) * t159;
t219 = qJD(2) * t161;
t218 = qJD(3) * t161;
t217 = qJD(4) * t159;
t216 = qJD(4) * t161;
t215 = t287 * qJD(1);
t214 = t132 * qJD(1);
t213 = t149 * qJD(2);
t151 = -t160 ^ 2 + t162 ^ 2;
t212 = t151 * qJD(1);
t211 = t160 * qJD(2);
t210 = t162 * qJD(2);
t207 = pkin(1) * t160 * qJD(1);
t206 = pkin(1) * t221;
t205 = pkin(4) * t217;
t202 = -t68 / 0.2e1 + t291;
t199 = t140 * t217;
t198 = t140 * t216;
t197 = t138 * t223;
t196 = t138 * t140 * qJD(2);
t195 = t159 * t216;
t194 = t159 * t219;
t193 = t160 * t221;
t192 = t140 * t222;
t188 = t159 * t204;
t187 = -qJD(4) - t224;
t185 = qJD(2) * t189;
t1 = -t288 + (t284 - t144 * t250 / 0.2e1 - t264 / 0.2e1) * pkin(4);
t49 = t144 * t279;
t183 = -t1 * qJD(1) + t49 * qJD(2);
t165 = (-t253 / 0.2e1 + t254 / 0.2e1) * t140 + t274;
t166 = t143 / 0.2e1 - t248 / 0.2e1 + t188;
t12 = t165 + t166;
t75 = t252 + t255;
t182 = t12 * qJD(1) + t75 * qJD(2);
t181 = -t138 * t154 - t140 * t153;
t180 = t187 * t161;
t177 = t153 * t138 / 0.2e1 + t154 * t282;
t167 = t177 * t161;
t22 = -t69 / 0.2e1 - t167;
t176 = -t22 * qJD(1) - t154 * t220;
t163 = t177 * t159 + t291;
t20 = t163 - t202;
t175 = -t20 * qJD(1) - t154 * t219;
t77 = (t156 / 0.2e1 + t281) * t140;
t174 = -t77 * qJD(1) + t194;
t173 = t140 * t180;
t172 = t132 * qJD(4) + t197;
t171 = t137 * t159 * t222 + t77 * qJD(2);
t87 = t150 * t137;
t170 = t87 * qJD(1) + t185;
t169 = qJD(1) * t189 - t150 * qJD(2);
t128 = t132 * qJD(2);
t126 = t140 * t219;
t107 = (t192 + t220) * pkin(4);
t71 = t78 * qJD(4);
t70 = t77 * qJD(4);
t61 = t277 / 0.2e1 + t168 * pkin(2);
t60 = -t217 - t232;
t53 = t186 + t132;
t23 = t103 * t280 + t246 / 0.2e1 + t69 / 0.2e1 - t167;
t21 = t163 + t202;
t13 = t164 - t178;
t11 = t165 - t166;
t10 = t45 * t280 - t266 / 0.2e1 + t188 + t274;
t5 = (t203 + t204) * t161;
t2 = t191 * t276 + t288 + (t264 / 0.2e1 + t284) * pkin(4);
t19 = [0, 0, 0, t160 * t210, t151 * qJD(2), 0, 0, 0, -pkin(1) * t211, -pkin(1) * t210, -t287 * qJD(3), t26 * qJD(2) + t48 * qJD(3), -t137 * t195 - t157 * t196, -t87 * qJD(4) + t138 * t185, t57 * qJD(2) - t138 * t199, -t55 * qJD(2) - t138 * t198, t196, t16 * qJD(2) - t56 * qJD(3) + t25 * qJD(4), t17 * qJD(2) - t88 * qJD(3) + t24 * qJD(4), -t6 * qJD(2) - t7 * qJD(4) + t86 * qJD(5), t3 * qJD(2) + t15 * qJD(3) + t8 * qJD(4) - t18 * qJD(5); 0, 0, 0, t193, t212, t210, -t211, 0, -pkin(6) * t210 - t207, pkin(6) * t211 - t206, (t259 * t138 - t140 * t158) * t271, t239 + (-t103 * t158 - t259 * t286) * t271 + t61 * qJD(3), -t70 + (-t157 * t223 - t194) * t138, (t129 - t130) * qJD(2) + (-qJD(4) + t224) * t189, t140 * t220 + t234, t126 - t236, t172, t245 + (t159 * t181 - t257) * qJD(2) + t23 * qJD(4), t244 + (t161 * t181 + t258) * qJD(2) + t21 * qJD(4), -t263 + (-t37 * t159 + t47 * t161 + (t253 - t254) * t138) * qJD(2) + t5 * qJD(4), t270 + (-t37 * t131 + t47 * t133 + t58 * t144) * qJD(2) + t13 * qJD(3) + t2 * qJD(4) + t11 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, t61 * qJD(2) + t238, 0, 0, 0, 0, 0, -t235, -t227, 0, t13 * qJD(2) + t10 * qJD(4) + t53 * qJD(5) + t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171, -t170, t187 * t81, t173, t128, t23 * qJD(2) - t51 * qJD(4) + t240, t21 * qJD(2) + t50 * qJD(4) + t241, pkin(4) * t199 + t5 * qJD(2) - t262, t2 * qJD(2) + t10 * qJD(3) - t46 * t272 + t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t11 * qJD(2) + t53 * qJD(3) - t243; 0, 0, 0, -t193, -t212, 0, 0, 0, t207, t206, 0, t62 * qJD(3) - t239, t157 * t197 - t70, 0.2e1 * t159 * t173, t83 * qJD(4) - t234, -t71 + t236, -t172, t22 * qJD(4) - t140 * t218 - t245, t81 * qJD(3) + t20 * qJD(4) - t244, -t85 * qJD(3) - t4 * qJD(4) + t263, t14 * qJD(3) - t1 * qJD(4) + t12 * qJD(5) - t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, t150 * qJD(4), 0, 0, 0, t154 * t217, t154 * t216, t149 * qJD(5), t49 * qJD(4) + t75 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, 0, 0, 0, 0, 0, -t192, t231, -t229, t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, -t169, t216 + t230, t60, -t214, -t153 * t216 - t176, t153 * t217 - t175, -pkin(4) * t216 - t267, -t133 * t272 + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, -t62 * qJD(2) - t238, 0, 0, 0, 0, 0, t126 - t71 + t235, -t81 * qJD(2) - t138 * t216 + t227, t85 * qJD(2), -t14 * qJD(2) - t9 * qJD(4) + t54 * qJD(5) - t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t233, 0, 0, 0, 0, 0, t192, -t231, t229, -t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t180, 0, -t205 - t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, t170, -t83 * qJD(2) + t159 * t197, t78 * qJD(2) + t138 * t192, t128, -t22 * qJD(2) + t78 * qJD(3) - t240, -t20 * qJD(2) + t138 * t218 - t241, t4 * qJD(2) + t262, t1 * qJD(2) + t9 * qJD(3) - qJD(5) * t208 - t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, t169, -t230, t232, t214, t176, t175, t267, -qJD(5) * t279 - t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, t138 * t222, 0, t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t228, pkin(4) * t198 - t12 * qJD(2) - t54 * qJD(3) + t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213, -t182 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t19;
