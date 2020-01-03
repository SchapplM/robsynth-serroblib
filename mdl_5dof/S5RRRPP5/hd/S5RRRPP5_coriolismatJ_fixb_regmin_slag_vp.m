% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:31
% EndTime: 2019-12-31 20:58:37
% DurationCPUTime: 2.18s
% Computational Cost: add. (2423->253), mult. (4558->280), div. (0->0), fcn. (4516->4), ass. (0->186)
t159 = qJD(2) + qJD(3);
t161 = sin(qJ(3));
t162 = sin(qJ(2));
t227 = t161 * t162;
t268 = -pkin(7) - pkin(6);
t138 = t268 * t227;
t163 = cos(qJ(2));
t142 = t268 * t163;
t256 = cos(qJ(3));
t188 = t256 * t142;
t276 = -t188 + t138;
t284 = t276 * pkin(3);
t194 = t256 * pkin(2);
t153 = -t194 - pkin(3);
t283 = t276 * t153;
t134 = -t256 * t163 + t227;
t123 = t134 * qJ(5);
t282 = t123 + t276;
t148 = t159 * qJ(4);
t266 = t123 / 0.2e1;
t66 = t266 - t123 / 0.2e1;
t281 = -t66 * qJD(1) - t148;
t186 = t256 * t162;
t179 = -t161 * t142 - t268 * t186;
t280 = t159 * t179;
t136 = t161 * t163 + t186;
t51 = -t136 * qJ(5) + t179;
t279 = t159 * t51;
t278 = -t188 / 0.2e1;
t69 = t159 * t134;
t164 = -pkin(3) - pkin(4);
t226 = t164 * t134;
t270 = t136 ^ 2;
t274 = qJD(4) * t270;
t273 = t270 * qJD(1);
t195 = t256 / 0.2e1;
t253 = t161 * pkin(2);
t151 = qJ(4) + t253;
t230 = t151 * t136;
t263 = -t136 / 0.2e1;
t272 = -t230 / 0.2e1 - (t134 * t195 + t161 * t263) * pkin(2);
t271 = t134 ^ 2;
t269 = -pkin(4) / 0.2e1;
t267 = -qJ(4) / 0.2e1;
t124 = qJ(4) * t134;
t265 = -t124 / 0.2e1;
t126 = t136 * qJ(4);
t264 = -t126 / 0.2e1;
t261 = t188 / 0.2e1;
t150 = -pkin(4) + t153;
t260 = t150 / 0.2e1;
t259 = -t151 / 0.2e1;
t258 = t161 / 0.2e1;
t257 = -t164 / 0.2e1;
t255 = t134 * pkin(3);
t254 = t136 * pkin(4);
t252 = t162 * pkin(2);
t196 = t138 + t278;
t57 = t278 + t196;
t251 = -qJD(2) * t276 - t57 * qJD(3);
t250 = -t57 * qJD(2) - qJD(3) * t276;
t154 = -t163 * pkin(2) - pkin(1);
t177 = t126 - t154;
t36 = t177 + t226;
t127 = t136 * pkin(3);
t71 = t124 + t127;
t42 = -t71 - t254;
t38 = t42 - t252;
t3 = t36 * t38;
t248 = t3 * qJD(1);
t247 = t36 * t134;
t246 = t36 * t136;
t4 = t36 * t42;
t245 = t4 * qJD(1);
t62 = -t177 + t255;
t242 = t62 * t134;
t241 = t62 * t136;
t64 = t71 + t252;
t9 = t62 * t64;
t240 = t9 * qJD(1);
t10 = t62 * t71;
t237 = t10 * qJD(1);
t13 = -t38 * t134 - t246;
t235 = t13 * qJD(1);
t14 = t38 * t136 - t247;
t234 = t14 * qJD(1);
t15 = -t42 * t134 - t246;
t233 = t15 * qJD(1);
t232 = t150 * t134;
t231 = t151 * t134;
t229 = t153 * t134;
t16 = t42 * t136 - t247;
t228 = t16 * qJD(1);
t176 = t151 / 0.2e1 - t253 / 0.2e1;
t18 = t264 + t176 * t136 + (t194 / 0.2e1 + t260 + t257) * t134;
t225 = t18 * qJD(1);
t19 = t134 * t282 - t51 * t136;
t224 = t19 * qJD(1);
t21 = (qJ(4) / 0.2e1 - t176) * t136 + (-t194 / 0.2e1 - t153 / 0.2e1 - pkin(3) / 0.2e1) * t134;
t223 = t21 * qJD(1);
t22 = t64 * t134 + t241;
t222 = t22 * qJD(1);
t23 = -t64 * t136 + t242;
t221 = t23 * qJD(1);
t180 = t265 - t127 / 0.2e1;
t168 = -t252 / 0.2e1 + t180;
t24 = -t231 / 0.2e1 + (t269 + t260) * t136 + t168;
t220 = t24 * qJD(1);
t26 = t71 * t134 + t241;
t219 = t26 * qJD(1);
t27 = -t71 * t136 + t242;
t218 = t27 * qJD(1);
t28 = t265 + (t269 + t164 / 0.2e1) * t136 + t180;
t217 = t28 * qJD(1);
t87 = t278 + t261;
t31 = t66 + t87;
t216 = t31 * qJD(1);
t55 = -t270 + t271;
t215 = t55 * qJD(1);
t59 = t134 * t252 + t154 * t136;
t214 = t59 * qJD(1);
t60 = -t154 * t134 + t136 * t252;
t213 = t60 * qJD(1);
t84 = t270 + t271;
t211 = t84 * qJD(1);
t210 = t87 * qJD(1);
t76 = t87 * qJD(2);
t199 = t151 * qJD(2);
t209 = t151 * qJD(3) + t199;
t208 = qJD(1) * t163;
t207 = qJD(3) * t154;
t206 = qJD(4) * t134;
t203 = t134 * qJD(1);
t202 = t136 * qJD(1);
t201 = t136 * qJD(4);
t149 = -t162 ^ 2 + t163 ^ 2;
t200 = t149 * qJD(1);
t198 = t162 * qJD(2);
t197 = t163 * qJD(2);
t157 = qJD(3) * t194;
t147 = t157 + qJD(4);
t193 = pkin(1) * t162 * qJD(1);
t192 = pkin(1) * t208;
t191 = qJD(3) * t253;
t190 = -t254 / 0.2e1;
t189 = t62 * t202;
t187 = t256 * t151;
t85 = t134 * t202;
t185 = t134 * t201;
t184 = t154 * t203;
t183 = t154 * t202;
t182 = t162 * t208;
t178 = -qJ(4) * qJD(3) - t199;
t70 = t159 * t136;
t166 = (t195 * t282 + t51 * t258) * pkin(2) + t51 * t259 + t282 * t260;
t169 = t257 * t282 - t267 * t51;
t2 = t166 + t169;
t92 = (t150 * t161 + t187) * pkin(2);
t175 = t2 * qJD(1) + t92 * qJD(2);
t165 = (t179 * t258 + t195 * t276) * pkin(2) + t179 * t259 + t283 / 0.2e1;
t170 = t284 / 0.2e1 - t179 * t267;
t8 = t165 + t170;
t95 = (t153 * t161 + t187) * pkin(2);
t174 = t8 * qJD(1) + t95 * qJD(2);
t173 = 0.2e1 * t266 + t138;
t160 = qJ(4) * qJD(4);
t158 = qJ(4) * t194;
t156 = qJD(2) * t194;
t143 = t151 * qJD(4);
t119 = t136 * qJD(5);
t115 = t134 * qJD(5);
t78 = t87 * qJD(3);
t61 = qJD(2) * t253 + t210;
t58 = qJD(4) - t156;
t49 = -t159 * t253 - t210;
t44 = t147 + t156;
t34 = -t123 + t261 - t196;
t32 = -t188 + t173;
t30 = 0.2e1 * t278 + t173;
t29 = t124 / 0.2e1 + t136 * t257 + t190 + t180;
t25 = t231 / 0.2e1 + t150 * t263 + t190 + t168;
t20 = -t229 / 0.2e1 + t264 + t255 / 0.2e1 + t272;
t17 = t232 / 0.2e1 + t126 / 0.2e1 + t226 / 0.2e1 - t272;
t7 = t165 - t170;
t1 = t166 - t169;
t5 = [0, 0, 0, t162 * t197, t149 * qJD(2), 0, 0, 0, -pkin(1) * t198, -pkin(1) * t197, -t134 * t70, t159 * t55, 0, 0, 0, t59 * qJD(2) + t136 * t207, t60 * qJD(2) - t134 * t207, t22 * qJD(2) + t26 * qJD(3) - t185, 0, t23 * qJD(2) + t27 * qJD(3) + t274, t9 * qJD(2) + t10 * qJD(3) - t201 * t62, t13 * qJD(2) + t15 * qJD(3) - t185, t14 * qJD(2) + t16 * qJD(3) + t274, t84 * qJD(5), t3 * qJD(2) + t4 * qJD(3) + t19 * qJD(5) + t201 * t36; 0, 0, 0, t182, t200, t197, -t198, 0, -pkin(6) * t197 - t193, pkin(6) * t198 - t192, -t85, t215, -t69, -t70, 0, t214 + t251, t280 + t213, t222 + t251, (-t229 - t230) * qJD(2) + t20 * qJD(3) - t206, -t280 + t221, t240 + (-t151 * t179 + t283) * qJD(2) + t7 * qJD(3) + t57 * qJD(4), -qJD(2) * t282 + t34 * qJD(3) + t235, t234 - t279, (t230 + t232) * qJD(2) + t17 * qJD(3) + t206, t248 + (t150 * t282 - t51 * t151) * qJD(2) + t1 * qJD(3) + t30 * qJD(4) + t25 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t215, -t69, -t70, 0, t183 + t250, t280 - t184, t219 + t250, t20 * qJD(2) + (-t126 + t255) * qJD(3) - t206, -t280 + t218, t237 + t7 * qJD(2) + (-qJ(4) * t179 - t284) * qJD(3) + t276 * qJD(4), t34 * qJD(2) - qJD(3) * t282 + t233, t228 - t279, t17 * qJD(2) + (t126 + t226) * qJD(3) + t206, t245 + t1 * qJD(2) + (-t51 * qJ(4) + t164 * t282) * qJD(3) + t32 * qJD(4) + t29 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t69, t273, -t189 - t250, -t85, t273, t69, t30 * qJD(2) + t32 * qJD(3) + t202 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t25 * qJD(2) + t29 * qJD(3) + t224; 0, 0, 0, -t182, -t200, 0, 0, 0, t193, t192, t85, -t215, 0, 0, 0, -t78 - t214, -t213, -t78 - t222, t21 * qJD(3), -t221, t8 * qJD(3) + t87 * qJD(4) - t240, t119 - t78 - t235, t115 - t234, t18 * qJD(3), t2 * qJD(3) + t31 * qJD(4) - t24 * qJD(5) - t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, -t157, -t191, 0, t147, t95 * qJD(3) + t143, -t191, t147, 0, t92 * qJD(3) + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t157 - t156, t49, t223, t44, (-pkin(3) * t253 + t158) * qJD(3) + t143 + t174, t49, t44, t225, (t164 * t253 + t158) * qJD(3) + t143 + t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t210 + t209, 0, t159, 0, t209 + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, t203, 0, -t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t215, 0, 0, 0, t76 - t183, t184, t76 - t219, -t21 * qJD(2), -t218, -t8 * qJD(2) - t237, t119 + t76 - t233, t115 - t228, -t18 * qJD(2), -t2 * qJD(2) + t66 * qJD(4) - t28 * qJD(5) - t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t156, t61, -t223, t58, t160 - t174, t61, t58, -t225, t160 - t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t160, 0, qJD(4), 0, t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t148, 0, t159, 0, -t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, t203, 0, -t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, -t273, t189 - t76, t85, -t273, 0, -t31 * qJD(2) - t66 * qJD(3) + (-qJD(1) * t36 - qJD(5)) * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159, t178 - t210, 0, -t159, 0, t178 - t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159, -t148, 0, -t159, 0, t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t69, -t211, t24 * qJD(2) + t28 * qJD(3) + t201 - t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, -t203, 0, t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, -t203, 0, t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
