% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRPR5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:47
% EndTime: 2021-01-15 16:04:59
% DurationCPUTime: 3.04s
% Computational Cost: add. (2135->217), mult. (5348->385), div. (0->0), fcn. (6097->10), ass. (0->190)
t301 = qJ(4) + pkin(7);
t167 = sin(qJ(5));
t171 = cos(qJ(3));
t252 = t301 * t171;
t165 = sin(pkin(10));
t168 = sin(qJ(3));
t269 = t165 * t168;
t277 = cos(pkin(10));
t298 = t277 * t252 - t269 * t301;
t275 = t298 * t167;
t170 = cos(qJ(5));
t274 = t298 * t170;
t211 = t277 * t171;
t187 = t211 - t269;
t135 = t187 ^ 2;
t212 = t277 * t168;
t268 = t165 * t171;
t139 = t212 + t268;
t136 = t139 ^ 2;
t299 = -t136 - t135;
t107 = t165 * t252 + t212 * t301;
t163 = t167 ^ 2;
t164 = t170 ^ 2;
t153 = t164 - t163;
t264 = t167 * t170;
t210 = 0.2e1 * t139 * t264;
t184 = qJD(2) * t210 - t153 * qJD(3);
t294 = t187 / 0.2e1;
t293 = -t139 / 0.2e1;
t292 = t139 / 0.2e1;
t289 = -t170 / 0.2e1;
t288 = t170 / 0.2e1;
t287 = t168 * pkin(3);
t285 = qJD(3) * pkin(3);
t92 = t139 * pkin(4) - pkin(8) * t187 + t287;
t284 = t167 * t92;
t283 = t170 * t92;
t166 = sin(pkin(5));
t169 = sin(qJ(2));
t267 = t166 * t169;
t278 = cos(pkin(5));
t132 = t168 * t267 - t278 * t171;
t133 = t278 * t168 + t171 * t267;
t188 = -t165 * t132 + t277 * t133;
t172 = cos(qJ(2));
t266 = t166 * t172;
t56 = t167 * t188 + t170 * t266;
t282 = t56 * t187;
t57 = -t167 * t266 + t170 * t188;
t281 = t57 * t187;
t73 = t277 * t132 + t165 * t133;
t280 = t73 * t139;
t279 = t188 * t187;
t110 = t139 * t266;
t273 = t110 * t167;
t272 = t110 * t170;
t271 = t139 * t170;
t111 = t187 * t266;
t265 = t167 * t111;
t89 = t167 * t139;
t263 = t170 * t111;
t21 = -t166 ^ 2 * t172 * t169 + t73 * t110 + t111 * t188;
t262 = t21 * qJD(1);
t235 = t136 - t135;
t58 = t235 * t167;
t261 = t58 * qJD(2);
t59 = t299 * t167;
t260 = t59 * qJD(2);
t60 = t235 * t170;
t259 = t60 * qJD(2);
t228 = -t277 / 0.2e1;
t181 = t139 * t228 + t165 * t294;
t71 = (-t168 / 0.2e1 + t181) * pkin(3);
t258 = t71 * qJD(2);
t86 = t167 * t187;
t257 = t86 * qJD(2);
t256 = t89 * qJD(2);
t91 = t170 * t187;
t255 = t91 * qJD(2);
t94 = t299 * t170;
t254 = t94 * qJD(2);
t251 = qJD(2) * t169;
t250 = qJD(2) * t171;
t249 = qJD(3) * t167;
t248 = qJD(3) * t170;
t247 = qJD(5) * t167;
t246 = qJD(5) * t170;
t245 = t299 * qJD(2);
t134 = t212 / 0.2e1 + t268 / 0.2e1;
t244 = t134 * qJD(2);
t243 = t187 * qJD(2);
t242 = t187 * qJD(4);
t241 = t139 * qJD(2);
t240 = t139 * qJD(3);
t239 = t139 * qJD(4);
t154 = -t168 ^ 2 + t171 ^ 2;
t238 = t154 * qJD(2);
t237 = t168 * qJD(3);
t236 = t171 * qJD(3);
t234 = pkin(2) * t168 * qJD(2);
t233 = pkin(2) * t250;
t231 = t167 * t267;
t230 = -t280 / 0.2e1;
t229 = t73 * t288;
t159 = -t171 * pkin(3) - pkin(2);
t227 = t187 * t246;
t226 = t187 * t241;
t225 = t187 * t240;
t224 = t166 * t251;
t223 = qJD(2) * t266;
t222 = t167 * t246;
t221 = t167 * t248;
t220 = t168 * t250;
t219 = t170 * t241;
t218 = t267 / 0.2e1;
t217 = -t266 / 0.2e1;
t216 = t266 / 0.2e1;
t215 = -t265 / 0.2e1;
t214 = -t263 / 0.2e1;
t209 = -qJD(5) + t243;
t205 = qJD(3) * t210;
t178 = t188 * t292;
t174 = t178 * t167 + t56 * t293;
t3 = t272 / 0.2e1 + t174;
t183 = -pkin(4) * t187 - t139 * pkin(8) + t159;
t32 = -t170 * t183 + t275;
t7 = (-t32 + t275) * t139 - t283 * t187;
t204 = t3 * qJD(1) + t7 * qJD(2);
t173 = t178 * t170 + t57 * t293;
t6 = -t273 / 0.2e1 + t173;
t33 = t167 * t183 + t274;
t8 = (-t33 + t274) * t139 + t284 * t187;
t203 = t6 * qJD(1) + t8 * qJD(2);
t182 = t111 * t165 / 0.2e1 + t110 * t228;
t1 = (t168 * t216 + t182) * pkin(3);
t24 = t159 * t287;
t202 = -t1 * qJD(1) + t24 * qJD(2);
t193 = t218 + t230;
t12 = t215 - t281 / 0.2e1 + t193 * t170;
t23 = t107 * t271 + t187 * t33;
t200 = t12 * qJD(1) - t23 * qJD(2);
t13 = t214 + t282 / 0.2e1 - t193 * t167;
t22 = -t107 * t89 - t187 * t32;
t199 = -t13 * qJD(1) + t22 * qJD(2);
t25 = -t279 / 0.2e1 + t193;
t31 = t107 * t139 + t187 * t298;
t198 = -t25 * qJD(1) + t31 * qJD(2);
t61 = (t292 - t134) * t266;
t79 = t159 * t139 - t187 * t287;
t197 = -t61 * qJD(1) + t79 * qJD(2);
t180 = -t211 / 0.2e1 + t269 / 0.2e1;
t62 = (t294 + t180) * t266;
t80 = t139 * t287 + t159 * t187;
t196 = -t62 * qJD(1) + t80 * qJD(2);
t157 = t165 * pkin(3) + pkin(8);
t158 = -t277 * pkin(3) - pkin(4);
t195 = -t139 * t157 + t158 * t187;
t194 = t209 * t170;
t192 = -t157 * t187 / 0.2e1 + t158 * t293;
t85 = (t163 / 0.2e1 - t164 / 0.2e1) * t139;
t191 = -t85 * qJD(2) + t221;
t190 = t139 * t194;
t189 = t134 * qJD(5) - t226;
t186 = t136 * qJD(2) * t264 + t85 * qJD(3);
t93 = t153 * t136;
t185 = t93 * qJD(2) + t205;
t179 = t92 / 0.2e1 + t192;
t19 = t179 * t170;
t177 = t19 * qJD(2) - t158 * t249;
t17 = t179 * t167;
t176 = -t17 * qJD(2) - t158 * t248;
t131 = t134 * qJD(3);
t130 = t170 * t240;
t82 = t86 * qJD(5);
t81 = t85 * qJD(5);
t70 = t287 / 0.2e1 + t181 * pkin(3);
t69 = -t247 + t257;
t64 = -t134 * t266 + t139 * t217;
t63 = t180 * t266 - t187 * t216;
t30 = -t289 * t73 + t229;
t29 = t73 * t167;
t26 = t279 / 0.2e1 + t280 / 0.2e1 + t218;
t20 = t283 / 0.2e1 - t192 * t170 + t107 * t167;
t18 = -t284 / 0.2e1 + t192 * t167 + (t288 - t289) * t107;
t15 = t281 / 0.2e1 + t139 * t229 + t215 + t170 * t218;
t14 = -t282 / 0.2e1 + t167 * t230 + t214 - t231 / 0.2e1;
t5 = t273 / 0.2e1 + t173;
t4 = -t272 / 0.2e1 + t174;
t2 = t182 * pkin(3) + t217 * t287;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t224, -t223, 0, 0, 0, 0, 0, (-t169 * t250 - t172 * t237) * t166, (t168 * t251 - t172 * t236) * t166, t64 * qJD(3) - t187 * t224, t63 * qJD(3) + t139 * t224, (t110 * t139 + t111 * t187) * qJD(2), t262 + (t110 * t107 + t111 * t298 + t159 * t267) * qJD(2) + t2 * qJD(3) + t26 * qJD(4), 0, 0, 0, 0, 0, (-(t170 * t267 - t265) * t187 + t110 * t89) * qJD(2) + t4 * qJD(3) + t15 * qJD(5), ((t231 + t263) * t187 + t110 * t271) * qJD(2) + t5 * qJD(3) + t14 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133 * qJD(3) - t168 * t223, t132 * qJD(3) - t171 * t223, t64 * qJD(2) - qJD(3) * t188, t63 * qJD(2) + qJD(3) * t73, 0, t2 * qJD(2) + (-t165 * t73 - t188 * t277) * t285, 0, 0, 0, 0, 0, t4 * qJD(2) + t29 * qJD(5) - t188 * t248, t5 * qJD(2) + t30 * qJD(5) + t188 * t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * qJD(2) + t29 * qJD(3) - t57 * qJD(5), t14 * qJD(2) + t30 * qJD(3) + t56 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 * qJD(3), -t62 * qJD(3), 0, -t1 * qJD(3) - t25 * qJD(4) - t262, 0, 0, 0, 0, 0, t3 * qJD(3) - t12 * qJD(5), t6 * qJD(3) - t13 * qJD(5); 0, 0, 0, 0, t168 * t236, t154 * qJD(3), 0, 0, 0, -pkin(2) * t237, -pkin(2) * t236, t79 * qJD(3), t80 * qJD(3), -qJD(4) * t299, t24 * qJD(3) + t31 * qJD(4), -t136 * t222 + t164 * t225, -t93 * qJD(5) - t187 * t205, t139 * t187 * t247 + t60 * qJD(3), -t58 * qJD(3) + t139 * t227, -t225, t7 * qJD(3) - t59 * qJD(4) + t23 * qJD(5), t8 * qJD(3) - t94 * qJD(4) + t22 * qJD(5); 0, 0, 0, 0, t220, t238, t236, -t237, 0, -pkin(7) * t236 - t234, pkin(7) * t237 - t233, -qJD(3) * t298 + t197, qJD(3) * t107 + t196, (-t139 * t165 - t187 * t277) * t285, (-t107 * t165 - t277 * t298) * t285 + t70 * qJD(4) + t202, -t81 - (-t164 * t241 - t221) * t187, -0.2e1 * t139 * t222 - t184 * t187, t167 * t240 + t259, t130 - t261, t189, (t167 * t195 - t274) * qJD(3) + t20 * qJD(5) + t204, (t170 * t195 + t275) * qJD(3) + t18 * qJD(5) + t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, t70 * qJD(3) + t198, 0, 0, 0, 0, 0, -t260, -t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t186, -t185, t209 * t89, t190, t131, t20 * qJD(3) - t33 * qJD(5) - t200, t18 * qJD(3) + t32 * qJD(5) + t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61 * qJD(2), t62 * qJD(2), 0, t1 * qJD(2), 0, 0, 0, 0, 0, -t3 * qJD(2), -t6 * qJD(2); 0, 0, 0, 0, -t220, -t238, 0, 0, 0, t234, t233, -t197 - t239, -t196 - t242, 0, t71 * qJD(4) - t202, -t164 * t226 - t81, 0.2e1 * t167 * t190, -t91 * qJD(5) - t259, t82 + t261, -t189, -t19 * qJD(5) - t170 * t239 - t204, t89 * qJD(4) + t17 * qJD(5) - t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, t153 * qJD(5), 0, 0, 0, t158 * t247, t158 * t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, -t243, 0, t258, 0, 0, 0, 0, 0, -t219, t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, -t184, t246 - t255, t69, -t244, -t157 * t246 - t177, t157 * t247 - t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 * qJD(2), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, t187 * qJD(3), t245, -t71 * qJD(3) - t198, 0, 0, 0, 0, 0, t130 + t82 + t260, -t89 * qJD(3) + t227 + t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, t243, 0, -t258, 0, 0, 0, 0, 0, t219, -t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * qJD(2), t13 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, t185, t91 * qJD(3) - t167 * t226, -t86 * qJD(3) - t187 * t219, t131, t19 * qJD(3) - t86 * qJD(4) + t200, -t17 * qJD(3) - t170 * t242 - t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, t184, t255, -t257, t244, t177, t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t257, -t170 * t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
