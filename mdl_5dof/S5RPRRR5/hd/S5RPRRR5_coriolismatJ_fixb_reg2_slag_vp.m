% Calculate inertial parameters regressor of coriolis matrix for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRR5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:49:09
% EndTime: 2022-01-20 09:49:13
% DurationCPUTime: 2.58s
% Computational Cost: add. (5123->208), mult. (9660->260), div. (0->0), fcn. (9384->8), ass. (0->162)
t245 = qJD(4) + qJD(5);
t202 = sin(qJ(4));
t288 = -pkin(8) - pkin(7);
t189 = t288 * t202;
t204 = cos(qJ(4));
t190 = t288 * t204;
t201 = sin(qJ(5));
t282 = cos(qJ(5));
t223 = -t282 * t189 - t201 * t190;
t303 = t245 * t223;
t203 = sin(qJ(3));
t222 = cos(pkin(9)) * pkin(1) + pkin(2);
t281 = sin(pkin(9)) * pkin(1);
t283 = cos(qJ(3));
t174 = t203 * t222 + t283 * t281;
t168 = pkin(7) + t174;
t277 = pkin(8) + t168;
t145 = t277 * t202;
t146 = t277 * t204;
t224 = t282 * t145 + t201 * t146;
t302 = t245 * t224;
t82 = -t201 * t145 + t282 * t146;
t301 = t245 * t82;
t150 = t201 * t189 - t282 * t190;
t300 = t245 * t150;
t246 = qJD(1) + qJD(3);
t299 = t204 * t246;
t238 = t282 * t204;
t265 = t201 * t202;
t183 = -t238 + t265;
t298 = t245 * t183;
t239 = t282 * t202;
t264 = t201 * t204;
t184 = t239 + t264;
t114 = t183 ^ 2 - t184 ^ 2;
t297 = t246 * t114;
t198 = t202 ^ 2;
t199 = t204 ^ 2;
t191 = t199 - t198;
t296 = t246 * t191;
t295 = t246 * t183 * t184;
t292 = -pkin(3) / 0.2e1;
t291 = t224 / 0.2e1;
t290 = -t82 / 0.2e1;
t173 = t203 * t281 - t283 * t222;
t99 = t183 * t173;
t289 = t99 / 0.2e1;
t285 = -t150 / 0.2e1;
t228 = t285 + t150 / 0.2e1;
t286 = t223 / 0.2e1;
t229 = -t223 / 0.2e1 + t286;
t241 = t290 + t82 / 0.2e1;
t242 = -t224 / 0.2e1 + t291;
t7 = (t228 + t241) * t184 + (-t229 - t242) * t183;
t287 = t7 * qJD(4);
t284 = t204 / 0.2e1;
t280 = pkin(4) * t201;
t279 = pkin(4) * t202;
t278 = t204 * pkin(4);
t23 = t228 * t183 + t229 * t184;
t276 = t23 * qJD(4);
t275 = pkin(3) * qJD(3);
t274 = qJD(4) * pkin(4);
t167 = -pkin(3) + t173;
t157 = t167 - t278;
t270 = t157 * t183;
t269 = t157 * t184;
t194 = -pkin(3) - t278;
t267 = t194 * t183;
t266 = t194 * t184;
t263 = t202 * t173;
t98 = t184 * t173;
t29 = t184 * t289 - t98 * t183 / 0.2e1;
t262 = t29 * qJD(1);
t30 = -t99 * t183 - t98 * t184;
t261 = t30 * qJD(1);
t101 = (-t198 - t199) * t173;
t35 = t168 * t101 + t167 * t174;
t260 = t35 * qJD(1);
t175 = t183 * t279;
t89 = t175 + t269;
t259 = t89 * qJD(1);
t244 = t184 * t279;
t90 = t244 - t270;
t258 = t90 * qJD(1);
t257 = -t270 / 0.2e1 - t267 / 0.2e1;
t230 = t173 * t284;
t256 = -t201 * t263 / 0.2e1 + t282 * t230;
t255 = qJD(1) * t157;
t254 = qJD(1) * t167;
t253 = qJD(3) * t194;
t252 = qJD(5) * t157;
t251 = qJD(5) * t194;
t250 = t101 * qJD(1);
t249 = t173 * qJD(1);
t248 = t174 * qJD(1);
t166 = t174 * qJD(3);
t247 = t202 * qJD(4);
t197 = t204 * qJD(4);
t240 = (t282 * t183 - t184 * t201) * t274;
t237 = t183 * t255;
t236 = t184 * t255;
t235 = t202 * t254;
t234 = t204 * t254;
t233 = t183 * t248;
t232 = t184 * t248;
t231 = t202 * t248;
t227 = t194 / 0.2e1 + t157 / 0.2e1;
t226 = t282 * qJD(4);
t225 = t282 * qJD(5);
t138 = t245 * t184;
t221 = t173 / 0.2e1 + pkin(3) / 0.2e1 - t167 / 0.2e1;
t220 = t7 * qJD(3);
t219 = t7 * qJD(1);
t11 = t241 * t183 + t242 * t184;
t14 = t157 * t279;
t218 = t14 * qJD(1) + t11 * qJD(2);
t15 = t157 * t174 - t224 * t98 + t82 * t99;
t217 = t15 * qJD(1) + t29 * qJD(2);
t216 = t11 * qJD(1) + t23 * qJD(3);
t123 = t175 + t266;
t206 = (t264 / 0.2e1 + t239 / 0.2e1) * t173;
t37 = -t227 * t184 + t206;
t31 = -t175 + t37;
t215 = t31 * qJD(1) - t123 * qJD(3);
t124 = t244 - t267;
t213 = t244 + t257;
t33 = (-t238 / 0.2e1 + t265 / 0.2e1) * t173 + t213;
t214 = -t33 * qJD(1) - t124 * qJD(3);
t94 = t221 * t202;
t212 = t94 * qJD(1) + t202 * t275;
t95 = t221 * t204;
t211 = t95 * qJD(1) + t204 * t275;
t38 = t227 * t183 + t256;
t210 = t38 * qJD(1) + t183 * t253;
t209 = t37 * qJD(1) - t184 * t253;
t208 = t201 * t289 + t98 * t282 / 0.2e1;
t205 = -t150 * t291 - t223 * t290 - t224 * t285 - t286 * t82;
t2 = (-t227 * t202 + t208) * pkin(4) + t205;
t36 = t194 * t279;
t207 = -t2 * qJD(1) + t23 * qJD(2) + t36 * qJD(3);
t40 = t266 / 0.2e1 + t269 / 0.2e1 + t206;
t192 = t202 * t197;
t188 = t191 * qJD(4);
t176 = t202 * t299;
t165 = t173 * qJD(3);
t160 = t202 * t166;
t129 = t184 * t166;
t128 = t183 * t166;
t100 = t101 * qJD(3);
t97 = t167 * t284 + t204 * t292 + t230;
t96 = t263 / 0.2e1 + (t292 + t167 / 0.2e1) * t202;
t85 = t183 * t138;
t83 = t184 * t298;
t52 = t245 * t114;
t39 = t256 + t257;
t34 = t213 + t256;
t32 = t175 + t40;
t28 = t30 * qJD(3);
t3 = t208 * pkin(4) - t205 + (t157 + t194) * t279 / 0.2e1;
t1 = t29 * qJD(3) + t11 * qJD(4);
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, t165, 0, 0, t192, t188, 0, -t192, 0, 0, -t204 * t166 + t167 * t247, t167 * t197 + t160, t100, t35 * qJD(3), -t85, t52, 0, t83, 0, 0, t89 * qJD(4) + t184 * t252 + t128, t90 * qJD(4) - t183 * t252 + t129, t28, t15 * qJD(3) + t14 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166 - t248, t165 + t249, 0, 0, t192, t188, 0, -t192, 0, 0, t96 * qJD(4) - t174 * t299, t97 * qJD(4) + t160 + t231, t100 + t250, t260 + (-t174 * pkin(3) + pkin(7) * t101) * qJD(3), -t85, t52, 0, t83, 0, 0, t32 * qJD(4) + t40 * qJD(5) + t128 + t233, t34 * qJD(4) + t39 * qJD(5) + t129 + t232, t28 + t261 + t287, (t99 * t150 + t174 * t194 - t223 * t98) * qJD(3) + t3 * qJD(4) + t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t296, t197, -t176, -t247, 0, t96 * qJD(3) - t168 * t197 + t235, t97 * qJD(3) + t168 * t247 + t234, 0, 0, -t295, t297, -t298, t295, -t138, 0, t32 * qJD(3) + t259 - t301, t34 * qJD(3) + t258 + t302, t220 + t240, t3 * qJD(3) + (-t201 * t224 - t282 * t82) * t274 + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t295, t297, -t298, t295, -t138, 0, t40 * qJD(3) + t236 - t301, t39 * qJD(3) - t237 + t302, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262 + t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t247, -t197, 0, 0, 0, 0, 0, 0, 0, 0, -t138, t298, 0, (-t183 * t201 - t282 * t184) * t274 + t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, t298, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, -t249, 0, 0, t192, t188, 0, -t192, 0, 0, -t94 * qJD(4) + t204 * t248, -t95 * qJD(4) - t231, -t250, -t260, -t85, t52, 0, t83, 0, 0, -t31 * qJD(4) - t37 * qJD(5) - t233, t33 * qJD(4) - t38 * qJD(5) - t232, -t261 + t287, -t2 * qJD(4) - t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262 + t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, t188, 0, -t192, 0, 0, -pkin(3) * t247, -pkin(3) * t197, 0, 0, -t85, t52, 0, t83, 0, 0, t123 * qJD(4) + t184 * t251, t124 * qJD(4) - t183 * t251, 0, t36 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t296, t197, -t176, -t247, 0, -pkin(7) * t197 - t212, pkin(7) * t247 - t211, 0, 0, -t295, t297, -t298, t295, -t138, 0, -t215 - t300, -t214 + t303, t219 + t240, (-t150 * t282 - t201 * t223) * t274 + t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t295, t297, -t298, t295, -t138, 0, -t209 - t300, -t210 + t303, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, -t296, 0, t176, 0, 0, t94 * qJD(3) - t235, t95 * qJD(3) - t234, 0, 0, t295, -t297, 0, -t295, 0, 0, t31 * qJD(3) - t259, -t33 * qJD(3) - t258, -t220, t2 * qJD(3) - t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, -t296, 0, t176, 0, 0, t212, t211, 0, 0, t295, -t297, 0, -t295, 0, 0, t215, t214, -t219, -t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t280, -pkin(4) * t225, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245 * t280, (-t226 - t225) * pkin(4), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, -t297, 0, -t295, 0, 0, t37 * qJD(3) - t236, t38 * qJD(3) + t237, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, -t297, 0, -t295, 0, 0, t209, t210, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201 * t274, pkin(4) * t226, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t4;
