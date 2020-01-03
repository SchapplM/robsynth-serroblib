% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x26]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRR8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:16
% EndTime: 2019-12-31 20:18:27
% DurationCPUTime: 5.15s
% Computational Cost: add. (5498->312), mult. (11117->452), div. (0->0), fcn. (13056->8), ass. (0->257)
t224 = cos(qJ(5));
t219 = sin(pkin(9));
t220 = cos(pkin(9));
t223 = sin(qJ(2));
t225 = cos(qJ(2));
t196 = -t219 * t223 + t220 * t225;
t197 = t219 * t225 + t220 * t223;
t222 = sin(qJ(4));
t379 = cos(qJ(4));
t166 = -t196 * t379 + t197 * t222;
t356 = t166 ^ 2;
t195 = t379 * t197;
t346 = t222 * t196;
t395 = t195 + t346;
t408 = t395 ^ 2;
t409 = t408 - t356;
t420 = t409 * t224;
t423 = t420 * qJD(1);
t405 = t224 * t166;
t414 = t405 / 0.2e1;
t415 = -t405 / 0.2e1;
t417 = t415 + t414;
t424 = qJD(5) * t417;
t427 = t423 + t424;
t292 = qJD(2) + qJD(4);
t221 = sin(qJ(5));
t372 = -qJ(3) - pkin(6);
t204 = t372 * t223;
t205 = t372 * t225;
t174 = t204 * t219 - t205 * t220;
t126 = -pkin(7) * t196 - t174;
t411 = t379 * t126;
t394 = t204 * t220 + t205 * t219;
t403 = -pkin(7) * t197 + t394;
t412 = t222 * t403;
t418 = -t411 + t412;
t426 = t418 * t221;
t425 = t418 * t224;
t421 = t409 * t221;
t422 = t421 * qJD(1);
t410 = t379 * t403;
t389 = -t410 / 0.2e1;
t259 = t411 / 0.2e1;
t398 = t395 * qJD(1);
t277 = t166 * t398;
t258 = t195 / 0.2e1;
t393 = t258 + t346 / 0.2e1;
t419 = qJD(5) * t393 + t277;
t413 = t222 * t126;
t56 = t410 + t413;
t416 = t409 * qJD(1);
t406 = t221 * t395;
t271 = t406 / 0.2e1;
t374 = t395 * pkin(4);
t375 = t166 * pkin(8);
t93 = t374 + t375;
t407 = -t166 / 0.2e1;
t386 = t166 / 0.2e1;
t388 = t395 / 0.2e1;
t376 = pkin(8) * t395;
t301 = t166 * qJD(4);
t404 = -qJD(2) * t166 - t301;
t402 = qJD(1) * t166;
t401 = qJD(3) * t166;
t399 = t393 * qJD(1);
t397 = 0.2e1 * t395;
t217 = t221 ^ 2;
t218 = t224 ^ 2;
t208 = t218 - t217;
t396 = t208 * t292;
t87 = 0.2e1 * t415;
t284 = pkin(2) * t220 + pkin(3);
t378 = pkin(2) * t219;
t193 = t222 * t284 + t378 * t379;
t187 = pkin(8) + t193;
t192 = t222 * t378 - t284 * t379;
t383 = t192 / 0.2e1;
t186 = -pkin(4) + t192;
t384 = -t186 / 0.2e1;
t265 = t383 + t384;
t227 = (-t187 / 0.2e1 + t193 / 0.2e1) * t395 + t265 * t166;
t392 = -t376 / 0.2e1 + t227;
t351 = t221 * t224;
t279 = qJD(1) * t351;
t264 = t217 / 0.2e1 - t218 / 0.2e1;
t70 = t264 * t395;
t391 = t279 * t408 + t292 * t70;
t382 = -t221 / 0.2e1;
t381 = -t224 / 0.2e1;
t380 = t224 / 0.2e1;
t377 = pkin(4) * t221;
t216 = t223 * pkin(2);
t371 = qJD(2) * pkin(2);
t285 = -pkin(2) * t225 - pkin(1);
t180 = -pkin(3) * t196 + t285;
t256 = pkin(4) * t166 - t376;
t228 = t180 + t256;
t19 = -t224 * t228 + t426;
t233 = t395 * t418;
t181 = pkin(3) * t197 + t216;
t58 = t181 + t93;
t1 = -t19 * t395 + t221 * t233 + t405 * t58;
t370 = t1 * qJD(1);
t20 = t221 * t228 + t425;
t75 = t221 * t166;
t2 = -t20 * t395 + t224 * t233 - t58 * t75;
t369 = t2 * qJD(1);
t3 = (-t19 + t426) * t395 + t93 * t405;
t366 = t3 * qJD(1);
t4 = (-t20 + t425) * t395 - t93 * t75;
t365 = t4 * qJD(1);
t226 = t376 / 0.2e1 + pkin(4) * t407 + t227;
t5 = t221 * t226;
t364 = t5 * qJD(1);
t361 = t221 * t56;
t360 = t224 * t56;
t303 = t395 * qJD(2);
t266 = 0.2e1 * t388;
t82 = t266 * t224;
t359 = qJD(4) * t82 + t224 * t303;
t13 = t166 * t19 + t406 * t56;
t358 = t13 * qJD(1);
t344 = t224 * t395;
t14 = -t166 * t20 - t344 * t56;
t357 = t14 * qJD(1);
t354 = t186 * t221;
t250 = t408 + t356;
t29 = t250 * t221;
t340 = t29 * qJD(1);
t31 = t250 * t224;
t338 = t31 * qJD(1);
t33 = t216 * t285;
t336 = t33 * qJD(1);
t34 = t166 * t181 + t180 * t395;
t335 = t34 * qJD(1);
t35 = -t166 * t180 + t181 * t395;
t334 = t35 * qJD(1);
t61 = t174 * t196 - t197 * t394;
t330 = t61 * qJD(1);
t63 = (t407 + t386) * t351;
t329 = t63 * qJD(1);
t328 = t70 * qJD(1);
t327 = t406 * qJD(1);
t267 = t388 - t395 / 0.2e1;
t73 = t267 * t221;
t326 = t73 * qJD(1);
t74 = t266 * t221;
t325 = t74 * qJD(1);
t324 = t75 * qJD(1);
t77 = 0.2e1 * t407 * t221;
t67 = t77 * qJD(1);
t81 = t267 * t224;
t323 = t81 * qJD(1);
t322 = t82 * qJD(1);
t321 = t405 * qJD(1);
t86 = 0.2e1 * t414;
t320 = t86 * qJD(1);
t319 = t87 * qJD(1);
t92 = t208 * t408;
t318 = t92 * qJD(1);
t316 = qJD(1) * t180;
t315 = qJD(1) * t225;
t314 = qJD(2) * t224;
t313 = qJD(4) * t224;
t312 = qJD(5) * t221;
t215 = qJD(5) * t224;
t109 = 0.2e1 * t258 + t346;
t310 = t109 * qJD(1);
t232 = (t219 * t196 / 0.2e1 - t220 * t197 / 0.2e1) * pkin(2);
t124 = -t216 / 0.2e1 + t232;
t309 = t124 * qJD(1);
t163 = t258 - t195 / 0.2e1;
t306 = t163 * qJD(1);
t305 = t163 * qJD(4);
t302 = t395 * qJD(3);
t298 = t395 * qJD(4);
t170 = t196 ^ 2 + t197 ^ 2;
t297 = t170 * qJD(1);
t296 = t193 * qJD(2);
t183 = t193 * qJD(4);
t209 = -t223 ^ 2 + t225 ^ 2;
t295 = t209 * qJD(1);
t294 = t223 * qJD(2);
t293 = t225 * qJD(2);
t290 = pkin(1) * t223 * qJD(1);
t289 = pkin(1) * t315;
t280 = t218 * t398;
t278 = qJD(5) * t166 * t395;
t276 = t395 * t402;
t275 = t166 * t316;
t274 = t395 * t316;
t210 = t221 * t215;
t273 = t223 * t315;
t272 = t224 * t398;
t270 = t166 * t382;
t262 = t292 * t224;
t260 = t395 * t279;
t257 = pkin(4) / 0.2e1 + t265;
t254 = -0.2e1 * t260;
t253 = 0.2e1 * t260;
t252 = t221 * t262;
t251 = t292 * t351;
t248 = -t186 * t166 - t187 * t395;
t22 = t389 + t410 / 0.2e1;
t247 = -qJD(1) * t22 - qJD(2) * t192;
t21 = t259 - t411 / 0.2e1;
t246 = qJD(1) * t21 + t296;
t245 = t395 * (-qJD(5) - t402);
t244 = qJD(4) * t109 + t303;
t242 = t375 / 0.2e1 + t374 / 0.2e1;
t239 = t187 * t386 + t384 * t395;
t231 = t58 / 0.2e1 + t239;
t9 = t221 * t231;
t241 = -qJD(1) * t9 - t186 * t314;
t8 = t224 * t226;
t240 = -qJD(1) * t8 - t221 * t296;
t11 = t224 * t231;
t238 = qJD(1) * t11 - qJD(2) * t354;
t235 = t404 * t395;
t234 = t93 / 0.2e1 + t242;
t110 = t257 * t221;
t17 = t234 * t224;
t230 = qJD(1) * t17 + qJD(2) * t110 + qJD(4) * t377;
t111 = t257 * t224;
t15 = t234 * t221;
t229 = pkin(4) * t313 - qJD(1) * t15 + qJD(2) * t111;
t207 = t208 * qJD(5);
t182 = t192 * qJD(4);
t179 = t221 * t183;
t125 = -0.2e1 * t395 * t210;
t123 = t216 / 0.2e1 + t232;
t113 = pkin(4) * t381 + (t186 + t192) * t380;
t112 = -t377 / 0.2e1 + t354 / 0.2e1 + t221 * t383;
t95 = t254 + t396;
t94 = t253 - t396;
t91 = t292 * t393;
t80 = t271 - t406 / 0.2e1;
t79 = 0.2e1 * t271;
t78 = -t75 / 0.2e1 - t270;
t69 = t81 * qJD(4);
t66 = t78 * qJD(5);
t65 = t77 * qJD(5);
t64 = t70 * qJD(5);
t62 = t87 * t221;
t60 = -t312 + t67;
t53 = t251 - t328;
t52 = -t252 + t328;
t27 = t166 * t264 + t217 * t386 + t218 * t407;
t24 = 0.2e1 * t259 - t412;
t23 = -t413 + 0.2e1 * t389;
t18 = -t224 * t242 + t380 * t93 - t361;
t16 = t221 * t242 + t382 * t93 - t360;
t12 = -t361 / 0.2e1 + t56 * t382 + t58 * t380 - t239 * t224;
t10 = -t360 / 0.2e1 + t56 * t381 + t58 * t382 + t239 * t221;
t7 = pkin(4) * t414 + t392 * t224 + t426;
t6 = -pkin(4) * t270 + t392 * t221 - t425;
t25 = [0, 0, 0, t223 * t293, t209 * qJD(2), 0, 0, 0, -pkin(1) * t294, -pkin(1) * t293, t170 * qJD(3), qJD(2) * t33 + qJD(3) * t61, t235, -t292 * t409, 0, 0, 0, qJD(2) * t34 + t180 * t298, qJD(2) * t35 - t180 * t301, -t210 * t408 + t218 * t235, -0.2e1 * t221 * t344 * t404 - t92 * qJD(5), -t221 * t278 + t292 * t420, -t224 * t278 - t292 * t421, (t298 + t303) * t166, qJD(2) * t1 + qJD(3) * t29 + qJD(4) * t3 + qJD(5) * t14, qJD(2) * t2 + qJD(3) * t31 + qJD(4) * t4 + qJD(5) * t13; 0, 0, 0, t273, t295, t293, -t294, 0, -pkin(6) * t293 - t290, pkin(6) * t294 - t289, (-t196 * t220 - t197 * t219) * t371, t336 + t123 * qJD(3) + (-t174 * t220 + t219 * t394) * t371, -t276, -t416, t404, -t244, 0, -qJD(2) * t418 + qJD(4) * t24 + t335, -qJD(2) * t56 + qJD(4) * t23 + t334, t62 * qJD(4) - t64 - (t221 * t314 + t280) * t166, t27 * qJD(4) + t125 - (qJD(2) * t208 + t254) * t166, qJD(4) * t79 + t221 * t303 + t427, t66 - t422 + t359, t419, t370 + (t221 * t248 - t425) * qJD(2) + t6 * qJD(4) + t12 * qJD(5), t369 + (t224 * t248 + t426) * qJD(2) + t7 * qJD(4) + t10 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297, qJD(2) * t123 + t330, 0, 0, 0, 0, 0, t305, 0, 0, 0, 0, 0, 0, t66 - t69 + t340, qJD(4) * t80 + t338 + t424; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t277, -t416, t404, -qJD(2) * t109 - t298, 0, qJD(2) * t24 + qJD(3) * t163 - qJD(4) * t418 + t274, qJD(2) * t23 - qJD(4) * t56 - t275, t62 * qJD(2) - t64 + (-t221 * t313 - t280) * t166, t27 * qJD(2) + t125 + (-qJD(4) * t208 + t253) * t166, qJD(2) * t79 + t221 * t298 + t427, qJD(2) * t82 + t224 * t298 - t422, t419, t366 + t6 * qJD(2) - t81 * qJD(3) + (t221 * t256 - t425) * qJD(4) + t18 * qJD(5), t365 + t7 * qJD(2) + t80 * qJD(3) + (t224 * t256 + t426) * qJD(4) + t16 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t391, -t252 * t397 - t318, t221 * t245 + t292 * t417, t78 * qJD(2) + t224 * t245, t91, qJD(2) * t12 + qJD(3) * t78 + qJD(4) * t18 - qJD(5) * t20 + t357, t10 * qJD(2) + qJD(3) * t417 + t16 * qJD(4) + t19 * qJD(5) + t358; 0, 0, 0, -t273, -t295, 0, 0, 0, t290, t289, 0, qJD(3) * t124 - t336, t276, t416, 0, -t305, 0, -qJD(4) * t21 - t302 - t335, qJD(4) * t22 - t334 + t401, qJD(4) * t63 + t218 * t276 - t64, -t166 * t253 + t125, -qJD(4) * t73 + qJD(5) * t86 - t423, t65 - t69 + t422, -t419, qJD(4) * t5 - qJD(5) * t11 - t224 * t302 - t370, qJD(3) * t406 + qJD(4) * t8 + qJD(5) * t9 - t369; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183, t182, t210, t207, 0, 0, 0, -t183 * t224 + t186 * t312, t186 * t215 + t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t309, 0, 0, 0, 0, 0, -t398, t402, 0, 0, 0, 0, 0, -t272, t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t306, 0, -t183 - t246, t182 - t247, t210 + t329, t207, -t326, -t323, 0, t112 * qJD(5) - t193 * t262 + t364, qJD(5) * t113 + t179 - t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t95, t215 + t320, t60, -t399, qJD(4) * t112 - t187 * t215 - t238, qJD(4) * t113 + t187 * t312 - t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t297, -qJD(2) * t124 - t330, 0, 0, 0, 0, 0, t244, t404, 0, 0, 0, 0, 0, t65 - t340 + t359, -qJD(2) * t406 - qJD(4) * t74 + qJD(5) * t87 - t338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t309, 0, 0, 0, 0, 0, t398, -t402, 0, 0, 0, 0, 0, t272, -t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, -t402, 0, 0, 0, 0, 0, t322, -t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t215 + t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t277, t416, 0, t163 * qJD(2), 0, qJD(2) * t21 - qJD(3) * t109 - t274, -qJD(2) * t22 + t275 + t401, -qJD(2) * t63 + t218 * t277 - t64, t166 * t254 + t125, qJD(2) * t73 + qJD(5) * t405 - t423, qJD(2) * t81 - qJD(5) * t75 + t422, -t419, -qJD(2) * t5 - qJD(3) * t82 - qJD(5) * t17 - t366, -qJD(2) * t8 + qJD(3) * t74 + qJD(5) * t15 - t365; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t306, 0, t246, t247, t210 - t329, t207, t326, t323, 0, -qJD(5) * t110 + t224 * t296 - t364, -qJD(5) * t111 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t310, t402, 0, 0, 0, 0, 0, -t322, t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, t207, 0, 0, 0, -pkin(4) * t312, -pkin(4) * t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t95, t215 + t321, -t312 - t324, -t399, -pkin(8) * t215 - t230, pkin(8) * t312 - t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t391, t251 * t397 + t318, -qJD(2) * t86 - qJD(4) * t405 + t221 * t277, -qJD(2) * t77 + qJD(4) * t75 + t224 * t277, t91, qJD(2) * t11 - qJD(3) * t77 + qJD(4) * t17 - t357, -qJD(2) * t9 - qJD(3) * t87 - qJD(4) * t15 - t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t94, -t320, -t67, t399, qJD(4) * t110 + t238, qJD(4) * t111 + t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t319; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t94, -t321, t324, t399, t230, t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t25;
