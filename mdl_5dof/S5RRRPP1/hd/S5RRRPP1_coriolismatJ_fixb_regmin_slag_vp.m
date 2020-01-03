% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:59
% EndTime: 2019-12-31 20:50:05
% DurationCPUTime: 2.20s
% Computational Cost: add. (3603->237), mult. (6935->303), div. (0->0), fcn. (6762->6), ass. (0->193)
t196 = sin(qJ(3));
t291 = cos(pkin(8));
t249 = t291 * t196;
t195 = sin(pkin(8));
t198 = cos(qJ(3));
t279 = t195 * t198;
t160 = t249 + t279;
t199 = cos(qJ(2));
t305 = t199 * pkin(1);
t128 = t160 * t305;
t248 = t291 * t198;
t280 = t195 * t196;
t211 = t248 - t280;
t129 = t211 * t305;
t219 = t128 * t160 + t129 * t211;
t320 = t219 * qJD(1);
t157 = t160 ^ 2;
t332 = t211 ^ 2 + t157;
t342 = t332 * qJD(4);
t349 = qJD(2) * t219 + t342;
t350 = t349 + t320;
t348 = t342 - t320;
t261 = qJD(1) + qJD(2);
t347 = t261 * t332;
t197 = sin(qJ(2));
t306 = t197 * pkin(1);
t258 = t306 / 0.2e1;
t193 = t198 * qJ(4);
t174 = pkin(7) * t198 + t193;
t163 = t291 * t174;
t304 = -qJ(4) - pkin(7);
t240 = t304 * t280;
t331 = t163 + t240;
t282 = t331 * t211;
t117 = t174 * t195 - t249 * t304;
t283 = t117 * t160;
t185 = pkin(7) + t306;
t155 = t185 * t198 + t193;
t139 = t291 * t155;
t245 = (-qJ(4) - t185) * t196;
t239 = t195 * t245;
t333 = t139 + t239;
t287 = t333 * t211;
t100 = t155 * t195 - t245 * t291;
t290 = t100 * t160;
t203 = t282 / 0.2e1 + t287 / 0.2e1 + t283 / 0.2e1 + t290 / 0.2e1 + t258;
t222 = t287 + t290;
t319 = t222 * qJD(1);
t346 = qJD(2) * t203 + t319;
t201 = (-t117 / 0.2e1 - t100 / 0.2e1) * t160 - t211 * (t331 / 0.2e1 + t333 / 0.2e1) + t258;
t343 = qJD(2) * t201 - t319;
t184 = -pkin(3) * t291 - pkin(4);
t307 = t196 * pkin(3);
t188 = t307 / 0.2e1;
t182 = pkin(3) * t195 + qJ(5);
t312 = t182 / 0.2e1;
t41 = t188 + (pkin(4) / 0.2e1 - t184 / 0.2e1) * t160 - (qJ(5) / 0.2e1 + t312) * t211;
t341 = t261 * t41;
t255 = -t291 / 0.2e1;
t313 = t211 / 0.2e1;
t209 = t160 * t255 + t195 * t313;
t309 = -t196 / 0.2e1;
t85 = (t309 + t209) * pkin(3);
t340 = t261 * t85;
t187 = -pkin(3) * t198 - pkin(2);
t103 = -pkin(4) * t211 - qJ(5) * t160 + t187;
t96 = t103 - t305;
t254 = t103 / 0.2e1 + t96 / 0.2e1;
t339 = t254 * t160;
t336 = t261 * t211;
t335 = t261 * t160;
t180 = -t196 ^ 2 + t198 ^ 2;
t334 = t261 * t180;
t325 = qJD(4) * t201;
t220 = t282 + t283;
t324 = qJD(4) * t220;
t323 = qJD(4) * t222;
t321 = t203 * qJD(4);
t316 = -qJD(1) * t201 + qJD(2) * t220;
t237 = t139 / 0.2e1;
t236 = t163 / 0.2e1;
t311 = t184 / 0.2e1;
t186 = -pkin(2) - t305;
t310 = t186 / 0.2e1;
t308 = t160 * pkin(4);
t303 = pkin(1) * qJD(1);
t302 = pkin(1) * qJD(2);
t301 = pkin(2) * qJD(2);
t300 = qJD(3) * pkin(3);
t281 = t211 * qJ(5);
t104 = -t281 + t307 + t308;
t298 = t96 * t104;
t299 = t298 * qJD(1);
t297 = t96 * t211;
t296 = t96 * t160;
t150 = t160 * qJD(5);
t293 = qJD(3) * t41 - t150;
t292 = (-t160 * t182 + t184 * t211) * qJD(3) + qJD(5) * t211;
t286 = t103 * t211;
t285 = t103 * t160;
t284 = t104 * t103;
t173 = t187 - t305;
t23 = t173 * t307;
t278 = t23 * qJD(1);
t223 = t100 * t128 + t129 * t333;
t24 = t306 * t96 + t223;
t277 = t24 * qJD(1);
t25 = t173 * t306 + t223;
t276 = t25 * qJD(1);
t76 = t104 * t211;
t29 = -t76 + t296;
t275 = t29 * qJD(1);
t77 = t104 * t160;
t30 = -t77 - t297;
t274 = t30 * qJD(1);
t113 = t211 * t150;
t260 = t197 * t302;
t269 = -t211 * t260 + t113;
t149 = t157 * qJD(5);
t268 = -t160 * t260 + t149;
t228 = t249 * t305;
t267 = t279 * t305 / 0.2e1 + t228 / 0.2e1;
t257 = -t305 / 0.2e1;
t241 = t198 * t257;
t266 = t195 * t241 - t228 / 0.2e1;
t265 = qJD(1) * t186;
t264 = t211 * qJD(3);
t263 = t160 * qJD(4);
t262 = t196 * qJD(3);
t192 = t198 * qJD(3);
t259 = t197 * t303;
t256 = qJD(1) * t296;
t253 = t196 * t265;
t252 = t198 * t265;
t247 = t117 * t128 + t129 * t331;
t246 = pkin(1) * t261;
t244 = t211 * t259;
t243 = t160 * t259;
t242 = t196 * t259;
t238 = t197 * t246;
t213 = t128 * t311 + t129 * t312;
t1 = -t104 * t254 + t213;
t235 = -qJD(1) * t1 + qJD(2) * t284;
t28 = t187 * t307;
t210 = t129 * t195 / 0.2e1 + t128 * t255;
t7 = ((-t187 / 0.2e1 - t173 / 0.2e1) * t196 + t210) * pkin(3);
t230 = -t7 * qJD(1) + t28 * qJD(2);
t17 = t266 + t76 - t339;
t33 = -t76 + t285;
t225 = qJD(1) * t17 - qJD(2) * t33;
t204 = (t248 / 0.2e1 - t280 / 0.2e1) * t305;
t18 = t211 * t254 + t204 + t77;
t34 = -t77 - t286;
t224 = qJD(1) * t18 - qJD(2) * t34;
t218 = t257 + pkin(2) / 0.2e1 - t186 / 0.2e1;
t216 = -t285 / 0.2e1 - t296 / 0.2e1;
t122 = t218 * t196;
t215 = qJD(1) * t122 + t196 * t301;
t123 = t218 * t198;
t214 = qJD(1) * t123 + t198 * t301;
t26 = t267 + t339;
t212 = -qJD(1) * t26 - qJD(2) * t285;
t115 = t236 - t163 / 0.2e1;
t98 = t237 - t139 / 0.2e1;
t207 = qJD(1) * t98 + qJD(2) * t115 + qJD(3) * t182;
t181 = t196 * t192;
t179 = t196 * t260;
t175 = t180 * qJD(3);
t148 = t160 * qJD(3);
t145 = t211 * qJD(4);
t143 = t261 * t198 * t196;
t125 = t241 + (-pkin(2) / 0.2e1 + t310) * t198;
t124 = pkin(2) * t309 + (t257 + t310) * t196;
t106 = t261 * t157;
t97 = (-t160 * t195 - t211 * t291) * t300;
t95 = 0.2e1 * t236 + t240;
t84 = pkin(3) * t209 + t188;
t82 = t85 * qJD(3);
t81 = t85 * qJD(4);
t80 = t84 * qJD(3);
t79 = t84 * qJD(4);
t60 = 0.2e1 * t237 + t239;
t59 = t211 * t335;
t42 = t182 * t313 + t160 * t311 + t188 - t281 / 0.2e1 + t308 / 0.2e1;
t40 = t41 * qJD(4);
t39 = t42 * qJD(3);
t38 = t42 * qJD(4);
t27 = t216 + t267;
t20 = -t77 - t286 / 0.2e1 - t297 / 0.2e1 + t204;
t19 = -t76 - t216 + t266;
t8 = t210 * pkin(3) + (t173 + t187) * t188;
t2 = t284 / 0.2e1 + t298 / 0.2e1 + t213;
t3 = [0, 0, 0, 0, -t260, -t199 * t302, t181, t175, 0, 0, 0, t186 * t262 - t198 * t260, t186 * t192 + t179, t349, qJD(2) * t25 + qJD(3) * t23 + t323, qJD(3) * t29 + t269, t349, qJD(3) * t30 + t268, qJD(2) * t24 + qJD(3) * t298 - t150 * t96 + t323; 0, 0, 0, 0, -t238, -t199 * t246, t181, t175, 0, 0, 0, t124 * qJD(3) - t198 * t238, qJD(3) * t125 + t179 + t242, t350, t276 + (t187 * t306 + t247) * qJD(2) + t8 * qJD(3) + t321, qJD(3) * t19 - t244 + t269, t350, qJD(3) * t20 - t243 + t268, t277 + (t103 * t306 + t247) * qJD(2) + t2 * qJD(3) + t321 + t27 * qJD(5); 0, 0, 0, 0, 0, 0, t143, t334, t192, -t262, 0, qJD(2) * t124 - t185 * t192 + t253, qJD(2) * t125 + t185 * t262 + t252, t97, t278 + t8 * qJD(2) + (-t100 * t195 - t291 * t333) * t300 + t79, qJD(2) * t19 - qJD(3) * t333 + t275, t292, qJD(2) * t20 - qJD(3) * t100 + t274, t299 + t2 * qJD(2) + (-t100 * t182 + t184 * t333) * qJD(3) + t38 + t60 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t347, t80 + t346, 0, t347, 0, t39 + t346; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t264, t106, qJD(2) * t27 + qJD(3) * t60 - t256; 0, 0, 0, 0, t259, t199 * t303, t181, t175, 0, 0, 0, -qJD(3) * t122 + t198 * t259, -qJD(3) * t123 - t242, t348, -qJD(3) * t7 - t276 - t325, -qJD(3) * t17 + t113 + t244, t348, -qJD(3) * t18 + t149 + t243, -qJD(3) * t1 - qJD(5) * t26 - t277 - t325; 0, 0, 0, 0, 0, 0, t181, t175, 0, 0, 0, -pkin(2) * t262, -pkin(2) * t192, t342, qJD(3) * t28 + t324, qJD(3) * t33 + t113, t342, qJD(3) * t34 + t149, qJD(3) * t284 - t103 * t150 + t324; 0, 0, 0, 0, 0, 0, t143, t334, t192, -t262, 0, -pkin(7) * t192 - t215, pkin(7) * t262 - t214, t97, (-t117 * t195 - t291 * t331) * t300 + t79 + t230, -qJD(3) * t331 - t225, t292, -qJD(3) * t117 - t224, (-t117 * t182 + t184 * t331) * qJD(3) + t38 + t95 * qJD(5) + t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t347, t316 + t80, 0, t347, 0, t316 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t264, t106, qJD(3) * t95 + t212; 0, 0, 0, 0, 0, 0, -t143, -t334, 0, 0, 0, qJD(2) * t122 - t253, qJD(2) * t123 - t252, 0, qJD(2) * t7 - t278 + t81, qJD(2) * t17 - t263 - t275, 0, qJD(2) * t18 + t145 - t274, qJD(2) * t1 + qJD(5) * t98 - t299 - t40; 0, 0, 0, 0, 0, 0, -t143, -t334, 0, 0, 0, t215, t214, 0, -t230 + t81, -t263 + t225, 0, t145 + t224, qJD(5) * t115 - t235 - t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t182 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, -t335, 0, t336, -t341; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t347, -t82 + t343, t148, -t347, -t264, t293 + t343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t347, -t316 - t82, t148, -t347, -t264, -t316 + t293; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340, t335, 0, -t336, t341; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, -t106, qJD(2) * t26 - qJD(3) * t98 + t256 + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, -t106, -qJD(3) * t115 - t212 + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
