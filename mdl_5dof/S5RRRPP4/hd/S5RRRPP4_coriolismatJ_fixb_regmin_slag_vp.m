% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPP4
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
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:48
% EndTime: 2019-12-31 20:55:56
% DurationCPUTime: 3.25s
% Computational Cost: add. (5033->226), mult. (9660->308), div. (0->0), fcn. (10680->6), ass. (0->188)
t176 = qJD(2) + qJD(3);
t177 = sin(pkin(8));
t178 = sin(qJ(3));
t179 = sin(qJ(2));
t180 = cos(qJ(2));
t299 = cos(qJ(3));
t159 = -t178 * t180 - t299 * t179;
t311 = pkin(6) + pkin(7);
t161 = t311 * t179;
t162 = t311 * t180;
t199 = -t299 * t161 - t178 * t162;
t321 = t159 * qJ(4) + t199;
t329 = t177 * t321;
t103 = t178 * t161 - t299 * t162;
t157 = t178 * t179 - t299 * t180;
t100 = -t157 * qJ(4) - t103;
t279 = cos(pkin(8));
t91 = t279 * t100;
t335 = t329 / 0.2e1 + t91 / 0.2e1;
t342 = 0.2e1 * t335;
t343 = t342 * qJD(5);
t340 = -t177 * t100 + t279 * t321;
t341 = t176 * t340;
t336 = t91 + t329;
t322 = -t177 * t157 - t279 * t159;
t123 = t322 ^ 2;
t332 = t322 / 0.2e1;
t244 = qJD(5) * t322;
t327 = t322 * qJD(1);
t116 = t322 * qJD(4);
t326 = t176 * t199;
t198 = -t279 * t157 + t177 * t159;
t325 = t176 * t198;
t320 = t176 * t103;
t208 = t198 * t336 - t322 * t340;
t319 = qJD(1) * t208;
t318 = qJD(4) * t208;
t205 = t198 ^ 2 + t123;
t315 = t205 * qJD(1);
t314 = qJD(4) * t205;
t313 = pkin(4) / 0.2e1;
t310 = -qJ(5) / 0.2e1;
t214 = t279 * t178;
t230 = t299 * pkin(2);
t172 = t230 + pkin(3);
t256 = t177 * t172;
t144 = pkin(2) * t214 + t256;
t137 = qJ(5) + t144;
t309 = -t137 / 0.2e1;
t308 = t137 / 0.2e1;
t298 = pkin(2) * t178;
t165 = t177 * t298;
t143 = t279 * t172 - t165;
t138 = -pkin(4) - t143;
t307 = t138 / 0.2e1;
t306 = t143 / 0.2e1;
t305 = -t144 / 0.2e1;
t295 = t177 * pkin(3);
t166 = qJ(5) + t295;
t304 = -t166 / 0.2e1;
t303 = t166 / 0.2e1;
t167 = -t279 * pkin(3) - pkin(4);
t302 = -t167 / 0.2e1;
t175 = t179 * pkin(2);
t301 = t175 / 0.2e1;
t300 = t177 / 0.2e1;
t297 = t322 * pkin(4);
t296 = t159 * pkin(3);
t291 = pkin(3) * qJD(3);
t173 = -pkin(2) * t180 - pkin(1);
t136 = t157 * pkin(3) + t173;
t52 = -pkin(4) * t198 - qJ(5) * t322 + t136;
t290 = t322 * t52;
t289 = t198 * t52;
t268 = t198 * qJ(5);
t60 = -t268 - t296 + t297;
t53 = t175 + t60;
t5 = t52 * t53;
t286 = t5 * qJD(1);
t225 = t299 * t177;
t150 = (t225 + t214) * pkin(2);
t285 = t340 * t150;
t151 = t279 * t230 - t165;
t284 = t336 * t151;
t6 = t52 * t60;
t283 = t6 * qJD(1);
t17 = -t198 * t53 + t290;
t278 = qJD(1) * t17;
t18 = -t322 * t53 - t289;
t277 = qJD(1) * t18;
t21 = -t198 * t60 + t290;
t274 = qJD(1) * t21;
t22 = -t322 * t60 - t289;
t273 = qJD(1) * t22;
t11 = t136 * (t175 - t296);
t271 = t11 * qJD(1);
t12 = t136 * t296;
t270 = t12 * qJD(1);
t269 = t322 * t177;
t267 = t137 * t322;
t266 = t138 * t198;
t220 = t150 * t332;
t261 = t151 * t198;
t76 = t261 / 0.2e1;
t213 = t76 + t220;
t14 = (t307 + t302) * t198 + (t309 + t303) * t322 + t213;
t265 = t14 * qJD(1);
t264 = t143 * t198;
t263 = t144 * t322;
t188 = -t264 / 0.2e1 - t263 / 0.2e1 + t220;
t215 = t279 * t198;
t191 = (-t269 / 0.2e1 - t215 / 0.2e1) * pkin(3);
t15 = -t261 / 0.2e1 + t191 - t188;
t262 = t15 * qJD(1);
t260 = t166 * t322;
t259 = t167 * t198;
t149 = -t296 / 0.2e1;
t231 = t301 + t149;
t23 = (t310 + t309) * t198 + (t313 - t138 / 0.2e1) * t322 + t231;
t255 = t23 * qJD(1);
t33 = t149 + (t310 + t304) * t198 + (t313 + t302) * t322;
t254 = t33 * qJD(1);
t200 = t198 * t305 + t306 * t322;
t37 = t200 + t231;
t251 = t37 * qJD(1);
t226 = -t279 / 0.2e1;
t190 = (t198 * t300 + t226 * t322) * pkin(3);
t50 = t296 / 0.2e1 + t190;
t250 = t50 * qJD(1);
t247 = qJD(1) * t173;
t246 = qJD(1) * t180;
t245 = qJD(3) * t173;
t101 = t157 ^ 2 - t159 ^ 2;
t243 = t101 * qJD(1);
t111 = t157 * t175 - t159 * t173;
t240 = t111 * qJD(1);
t112 = -t157 * t173 - t159 * t175;
t239 = t112 * qJD(1);
t238 = t123 * qJD(1);
t236 = t198 * qJD(1);
t164 = -t179 ^ 2 + t180 ^ 2;
t235 = t164 * qJD(1);
t234 = t179 * qJD(2);
t233 = t180 * qJD(2);
t232 = t151 * qJD(3) + qJD(5);
t229 = pkin(1) * t179 * qJD(1);
t228 = pkin(1) * t246;
t227 = t52 * t327;
t224 = t198 * t327;
t223 = t157 * t247;
t222 = t159 * t247;
t221 = t179 * t246;
t219 = t299 * qJD(2);
t218 = t299 * qJD(3);
t131 = t176 * t159;
t45 = t284 / 0.2e1;
t185 = t45 + t336 * t307 - t285 / 0.2e1 + t340 * t308;
t201 = t302 * t336 + t304 * t340;
t2 = t185 + t201;
t61 = t137 * t151 + t138 * t150;
t207 = t2 * qJD(1) + t61 * qJD(2);
t189 = t336 * t306 + t285 / 0.2e1 + t340 * t305;
t193 = (t226 * t336 + t300 * t340) * pkin(3);
t3 = -t284 / 0.2e1 + t193 + t189;
t62 = -t143 * t150 + t144 * t151;
t206 = t3 * qJD(1) - t62 * qJD(2);
t204 = qJD(2) * t137;
t203 = qJD(2) * t150;
t202 = qJD(2) * t151;
t197 = t149 + t297 / 0.2e1 - t268 / 0.2e1;
t135 = -qJ(5) + (t230 / 0.2e1 - pkin(3) / 0.2e1 - t172 / 0.2e1) * t177;
t194 = qJD(2) * t135 - qJD(3) * t166;
t145 = t150 * qJD(3);
t132 = t159 * t157 * qJD(1);
t130 = t176 * t157;
t129 = t295 / 0.2e1 + qJ(5) + t256 / 0.2e1 + (t214 + t225 / 0.2e1) * pkin(2);
t118 = t198 * qJD(4);
t67 = t198 * qJD(5);
t51 = t149 + t190;
t38 = -t200 + t231;
t34 = t167 * t332 + t198 * t303 + t197;
t27 = -0.2e1 * t335;
t24 = t198 * t308 + t307 * t322 + t197 + t301;
t16 = t76 + t191 + t188;
t13 = -t267 / 0.2e1 + t266 / 0.2e1 - t260 / 0.2e1 + t259 / 0.2e1 + t213;
t4 = t45 + t193 - t189;
t1 = t185 - t201;
t7 = [0, 0, 0, t179 * t233, t164 * qJD(2), 0, 0, 0, -pkin(1) * t234, -pkin(1) * t233, t157 * t131, t176 * t101, 0, 0, 0, qJD(2) * t111 - t159 * t245, qJD(2) * t112 - t157 * t245, t314, qJD(2) * t11 - qJD(3) * t12 + t318, qJD(2) * t17 + qJD(3) * t21 + t198 * t244, t314, qJD(2) * t18 + qJD(3) * t22 + qJD(5) * t123, qJD(2) * t5 + qJD(3) * t6 - t244 * t52 + t318; 0, 0, 0, t221, t235, t233, -t234, 0, -pkin(6) * t233 - t229, pkin(6) * t234 - t228, t132, t243, -t130, t131, 0, t240 + t320, t239 - t326, (-t263 - t264) * qJD(2) + t16 * qJD(3), t271 + (-t143 * t336 + t144 * t340) * qJD(2) + t4 * qJD(3) + t38 * qJD(4), -qJD(2) * t336 + qJD(3) * t27 + t278, (t266 - t267) * qJD(2) + t13 * qJD(3) + t67, t277 + t341, t286 + (t137 * t340 + t138 * t336) * qJD(2) + t1 * qJD(3) + t24 * qJD(4) + t343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t243, -t130, t131, 0, -t222 + t320, -t223 - t326, t16 * qJD(2) + (-t215 - t269) * t291, -t270 + t4 * qJD(2) + (t177 * t340 - t279 * t336) * t291 + t51 * qJD(4), qJD(2) * t27 - qJD(3) * t336 + t274, t13 * qJD(2) + (t259 - t260) * qJD(3) + t67, t273 + t341, t283 + t1 * qJD(2) + (t166 * t340 + t167 * t336) * qJD(3) + t34 * qJD(4) + t343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t315, qJD(2) * t38 + qJD(3) * t51 + t319, 0, t315, 0, qJD(2) * t24 + qJD(3) * t34 + t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, t325, t238, t176 * t342 - t227; 0, 0, 0, -t221, -t235, 0, 0, 0, t229, t228, -t132, -t243, 0, 0, 0, -t240, -t239, -qJD(3) * t15, -qJD(3) * t3 - qJD(4) * t37 - t271, -t116 - t278, qJD(3) * t14, t118 - t277, qJD(3) * t2 - qJD(4) * t23 - t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t298, -pkin(2) * t218, 0, t62 * qJD(3), -t145, 0, t232, qJD(3) * t61 + qJD(5) * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176 * t298, (-t219 - t218) * pkin(2), -t262, (-t150 * t279 + t151 * t177) * t291 - t206, -t145 - t203, t265, t202 + t232, (t150 * t167 + t151 * t166) * qJD(3) + t129 * qJD(5) + t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t251, -t327, 0, t236, -t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, qJD(3) * t129 + t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t243, 0, 0, 0, t222, t223, qJD(2) * t15, qJD(2) * t3 + qJD(4) * t50 + t270, -t116 - t274, -qJD(2) * t14, t118 - t273, -qJD(2) * t2 - qJD(4) * t33 - t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t298, pkin(2) * t219, t262, t206, t203, -t265, qJD(5) - t202, -qJD(5) * t135 - t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t166 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t250, -t327, 0, t236, -t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, -t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t315, qJD(2) * t37 - qJD(3) * t50 - t319, t176 * t322, -t315, -t325, qJD(2) * t23 + qJD(3) * t33 - t244 - t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t327, 0, -t236, t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t250, t327, 0, -t236, t254; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224, 0, -t238, t116 + t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, qJD(3) * t135 - t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
