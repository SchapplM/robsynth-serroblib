% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x28]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:50
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:49:02
% EndTime: 2021-01-15 22:49:22
% DurationCPUTime: 5.99s
% Computational Cost: add. (6690->234), mult. (12999->319), div. (0->0), fcn. (15145->8), ass. (0->208)
t281 = qJD(2) + qJD(3);
t449 = qJD(5) + t281;
t229 = sin(qJ(5));
t232 = cos(qJ(5));
t230 = sin(qJ(3));
t231 = sin(qJ(2));
t233 = cos(qJ(2));
t372 = cos(qJ(3));
t205 = t230 * t231 - t372 * t233;
t207 = -t230 * t233 - t372 * t231;
t227 = sin(pkin(9));
t228 = cos(pkin(9));
t168 = -t228 * t205 + t227 * t207;
t380 = pkin(6) + pkin(7);
t213 = t380 * t231;
t214 = t380 * t233;
t140 = t372 * t213 + t230 * t214;
t388 = t207 * qJ(4) - t140;
t392 = t227 * t388;
t209 = t372 * t214;
t324 = t230 * t213;
t387 = -t209 + t324;
t132 = t205 * qJ(4) + t387;
t399 = t228 * t132;
t403 = t399 - t392;
t409 = -t168 * pkin(8) + t403;
t250 = -t227 * t205 - t228 * t207;
t404 = t227 * t132 + t228 * t388;
t411 = -t250 * pkin(8) + t404;
t446 = t449 * (-t229 * t409 - t232 * t411);
t445 = t449 * (-t229 * t411 + t232 * t409);
t419 = t229 * t168 + t232 * t250;
t293 = t419 * qJD(5);
t385 = t281 * t419 + t293;
t426 = t419 * qJD(1);
t91 = t419 * qJD(4);
t390 = t229 * t250;
t395 = t232 * t168;
t99 = t395 - t390;
t408 = -t419 ^ 2 + t99 ^ 2;
t418 = t408 * qJD(1);
t416 = t281 * t404;
t413 = t99 * t426;
t412 = t403 * t228;
t379 = t399 / 0.2e1;
t371 = pkin(2) * t230;
t218 = t227 * t371;
t279 = t372 * pkin(2);
t199 = t228 * t279 - t218;
t374 = t199 / 0.2e1;
t378 = t395 / 0.2e1;
t393 = -t250 / 0.2e1;
t386 = t281 * t140;
t377 = t250 / 0.2e1;
t223 = t279 + pkin(3);
t331 = t228 * t230;
t191 = pkin(2) * t331 + t227 * t223;
t376 = -t191 / 0.2e1;
t198 = (t372 * t227 + t331) * pkin(2);
t375 = -t198 / 0.2e1;
t261 = -t209 / 0.2e1;
t373 = t227 / 0.2e1;
t369 = t207 * pkin(3);
t368 = t227 * pkin(3);
t226 = t231 * pkin(2);
t363 = pkin(3) * qJD(3);
t224 = -t233 * pkin(2) - pkin(1);
t179 = t205 * pkin(3) + t224;
t102 = -pkin(4) * t168 + t179;
t362 = t102 * t419;
t361 = t102 * t99;
t346 = t250 * t227;
t345 = t168 * t228;
t186 = t226 - t369;
t17 = t179 * t186;
t344 = t17 * qJD(1);
t18 = t179 * t369;
t343 = t18 * qJD(1);
t189 = t228 * t223 - t218;
t342 = t189 * t168;
t341 = t191 * t250;
t327 = t229 * t198;
t326 = t229 * t199;
t321 = t232 * t199;
t234 = t342 / 0.2e1 + t341 / 0.2e1 + t250 * t375 - t168 * t374;
t239 = (-t346 / 0.2e1 - t345 / 0.2e1) * pkin(3);
t24 = t239 + t234;
t320 = t24 * qJD(1);
t112 = pkin(4) * t250 - t369;
t103 = t112 + t226;
t26 = -t103 * t99 + t362;
t319 = t26 * qJD(1);
t27 = t103 * t419 + t361;
t318 = t27 * qJD(1);
t28 = -t112 * t99 + t362;
t317 = t28 * qJD(1);
t29 = t112 * t419 + t361;
t316 = t29 * qJD(1);
t30 = -t168 * t403 - t250 * t404;
t315 = t30 * qJD(1);
t42 = -t390 + 0.2e1 * t378;
t313 = t42 * qJD(1);
t43 = t378 - t395 / 0.2e1 + (t393 + t377) * t229;
t312 = t43 * qJD(1);
t47 = t168 ^ 2 + t250 ^ 2;
t308 = t47 * qJD(1);
t243 = t168 * t376 + t189 * t377;
t197 = -t369 / 0.2e1;
t280 = t197 + t226 / 0.2e1;
t60 = t243 + t280;
t307 = t60 * qJD(1);
t104 = t179 * t250;
t68 = -t168 * t186 + t104;
t306 = t68 * qJD(1);
t105 = t179 * t168;
t69 = t186 * t250 + t105;
t305 = t69 * qJD(1);
t74 = -t168 * t369 - t104;
t304 = t74 * qJD(1);
t75 = t250 * t369 - t105;
t303 = t75 * qJD(1);
t242 = t168 * t373 + t228 * t393;
t82 = (t207 / 0.2e1 + t242) * pkin(3);
t302 = t82 * qJD(1);
t300 = t99 * qJD(5);
t299 = t99 * qJD(1);
t297 = qJD(1) * t102;
t296 = qJD(1) * t224;
t295 = qJD(1) * t233;
t294 = qJD(3) * t224;
t137 = t205 ^ 2 - t207 ^ 2;
t292 = t137 * qJD(1);
t156 = t205 * t226 - t224 * t207;
t289 = t156 * qJD(1);
t157 = -t224 * t205 - t207 * t226;
t288 = t157 * qJD(1);
t287 = t250 * qJD(1);
t286 = t168 * qJD(1);
t174 = t261 + t209 / 0.2e1;
t285 = t174 * qJD(1);
t217 = -t231 ^ 2 + t233 ^ 2;
t284 = t217 * qJD(1);
t283 = t231 * qJD(2);
t282 = t233 * qJD(2);
t278 = pkin(1) * t231 * qJD(1);
t277 = pkin(1) * t295;
t273 = t99 * t297;
t272 = t419 * t297;
t271 = t205 * t296;
t270 = t207 * t296;
t269 = t231 * t295;
t188 = pkin(4) + t189;
t219 = t228 * pkin(3) + pkin(4);
t267 = -t219 / 0.2e1 - t188 / 0.2e1;
t266 = t372 * qJD(2);
t265 = t372 * qJD(3);
t262 = t281 * t99;
t171 = t281 * t207;
t260 = t368 / 0.2e1 + t191 / 0.2e1;
t257 = t374 + t267;
t235 = (t376 + t198 / 0.2e1) * t404 + (t374 - t189 / 0.2e1) * t403;
t240 = (t404 * t373 + t412 / 0.2e1) * pkin(3);
t1 = t240 + t235;
t101 = -t189 * t198 + t191 * t199;
t254 = t1 * qJD(1) - t101 * qJD(2);
t128 = -t232 * t188 + t229 * t191;
t253 = t128 * qJD(2);
t129 = t229 * t188 + t232 * t191;
t252 = t129 * qJD(2);
t187 = t232 * t198;
t154 = t187 + t326;
t251 = t154 * qJD(2);
t155 = t321 - t327;
t249 = t155 * qJD(2);
t32 = t379 - t399 / 0.2e1;
t248 = -t32 * qJD(1) + t198 * qJD(2);
t247 = t199 * qJD(2);
t246 = t260 * t232;
t190 = -t232 * t219 + t229 * t368;
t76 = t257 * t232 + (t375 + t260) * t229;
t238 = -t76 * qJD(2) - t190 * qJD(3);
t192 = t229 * t219 + t232 * t368;
t78 = t187 / 0.2e1 - t246 + t257 * t229;
t237 = -t78 * qJD(2) + t192 * qJD(3);
t194 = t199 * qJD(3);
t193 = t198 * qJD(3);
t185 = t192 * qJD(5);
t184 = t190 * qJD(5);
t172 = t207 * t205 * qJD(1);
t170 = t281 * t205;
t164 = t168 * qJD(4);
t161 = t250 * qJD(4);
t141 = 0.2e1 * t261 + t324;
t139 = t155 * qJD(3);
t138 = t154 * qJD(3);
t122 = t129 * qJD(5);
t121 = t128 * qJD(5);
t94 = t99 * qJD(4);
t83 = pkin(3) * t242 + t197;
t79 = -t326 / 0.2e1 - t187 / 0.2e1 - t246 + t267 * t229;
t77 = -t321 / 0.2e1 + t327 / 0.2e1 + t267 * t232 + t260 * t229;
t61 = -t243 + t280;
t36 = t43 * qJD(5);
t33 = -t392 + 0.2e1 * t379;
t25 = t239 - t234;
t20 = t42 * qJD(5) + t262;
t2 = t240 - t235;
t3 = [0, 0, 0, t231 * t282, t217 * qJD(2), 0, 0, 0, -pkin(1) * t283, -pkin(1) * t282, t205 * t171, t281 * t137, 0, 0, 0, t156 * qJD(2) - t207 * t294, t157 * qJD(2) - t205 * t294, qJD(2) * t68 - qJD(3) * t74, qJD(2) * t69 - qJD(3) * t75, qJD(4) * t47, qJD(2) * t17 - qJD(3) * t18 + qJD(4) * t30, (t262 + t300) * t419, t449 * t408, 0, 0, 0, t26 * qJD(2) + t28 * qJD(3) + t102 * t293, t27 * qJD(2) + t29 * qJD(3) + t102 * t300; 0, 0, 0, t269, t284, t282, -t283, 0, -pkin(6) * t282 - t278, pkin(6) * t283 - t277, t172, t292, -t170, t171, 0, qJD(2) * t387 + t141 * qJD(3) + t289, t288 + t386, qJD(2) * t403 + t33 * qJD(3) + t306, t305 - t416, (-t341 - t342) * qJD(2) + t25 * qJD(3), t344 + (t189 * t403 + t191 * t404) * qJD(2) + t2 * qJD(3) + t61 * qJD(4), t413, t418, t20, -t385, 0, t319 + t445, t318 + t446; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t292, -t170, t171, 0, t141 * qJD(2) + qJD(3) * t387 - t270, -t271 + t386, t33 * qJD(2) + qJD(3) * t403 - t304, -t303 - t416, t25 * qJD(2) + (-t345 - t346) * t363, -t343 + t2 * qJD(2) + t83 * qJD(4) + (t227 * t404 + t412) * t363, t413, t418, t20, -t385, 0, t317 + t445, t316 + t446; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t308, t61 * qJD(2) + t83 * qJD(3) + t315, 0, 0, 0, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t413, t418, t281 * t42 + t300, -t385, 0, t272 + t445, t43 * qJD(4) + t273 + t446; 0, 0, 0, -t269, -t284, 0, 0, 0, t278, t277, -t172, -t292, 0, 0, 0, t174 * qJD(3) - t289, -t288, t32 * qJD(3) - t161 - t306, -t164 - t305, -qJD(3) * t24, -qJD(3) * t1 - qJD(4) * t60 - t344, -t413, -t418, t36, 0, 0, -t319 - t91, -t318 - t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t371, -pkin(2) * t265, -t193, -t194, 0, t101 * qJD(3), 0, 0, 0, 0, 0, -t138 - t122, -t139 + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t281 * t371 + t285, (-t266 - t265) * pkin(2), -t193 - t248, -t194 - t247, -t320, (-t198 * t228 + t199 * t227) * t363 - t254, 0, 0, 0, 0, 0, t79 * qJD(5) - t138 - t251, t77 * qJD(5) - t139 - t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, -t286, 0, -t307, 0, 0, 0, 0, 0, -t426, -t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312, 0, 0, t79 * qJD(3) - t122 - t252, t77 * qJD(3) + t121 + t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t172, -t292, 0, 0, 0, -t174 * qJD(2) + t270, t271, -t32 * qJD(2) - t161 + t304, -t164 + t303, qJD(2) * t24, t1 * qJD(2) + t82 * qJD(4) + t343, -t413, -t418, t36, 0, 0, -t317 - t91, -t316 - t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t371 - t285, pkin(2) * t266, t248, t247, t320, t254, 0, 0, 0, 0, 0, t78 * qJD(5) + t251, t76 * qJD(5) + t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185, t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, -t286, 0, t302, 0, 0, 0, 0, 0, -t426, -t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312, 0, 0, -t185 - t237, t184 - t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t281 * t250, t281 * t168, -t308, t60 * qJD(2) - t82 * qJD(3) - t315, 0, 0, 0, 0, 0, t385, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t287, t286, 0, t307, 0, 0, 0, 0, 0, t426, t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t287, t286, 0, -t302, 0, 0, 0, 0, 0, t426, t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t426, t313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t413, -t418, -t281 * t43, 0, 0, -t272 - t91, -t42 * qJD(4) - t273; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t312, 0, 0, -t78 * qJD(3) + t252, -t76 * qJD(3) - t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t312, 0, 0, t237, t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t426, -t313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
