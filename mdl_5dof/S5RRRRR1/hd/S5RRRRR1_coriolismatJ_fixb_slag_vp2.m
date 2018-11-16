% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:51:15
% EndTime: 2018-11-16 14:51:24
% DurationCPUTime: 5.29s
% Computational Cost: add. (9863->400), mult. (21102->538), div. (0->0), fcn. (24558->8), ass. (0->257)
t234 = cos(qJ(5));
t400 = t234 / 0.2e1;
t232 = sin(qJ(5));
t386 = Ifges(6,4) * t234;
t217 = Ifges(6,1) * t232 + t386;
t343 = t234 * t217;
t387 = Ifges(6,4) * t232;
t216 = Ifges(6,2) * t234 + t387;
t346 = t232 * t216;
t273 = t346 / 0.2e1 - t343 / 0.2e1;
t289 = -Ifges(6,2) * t232 + t386;
t290 = Ifges(6,1) * t234 - t387;
t401 = -t232 / 0.2e1;
t299 = t289 * t400 - t290 * t401 - t273;
t215 = Ifges(6,5) * t232 + Ifges(6,6) * t234;
t395 = sin(qJ(3));
t396 = sin(qJ(2));
t398 = cos(qJ(3));
t399 = cos(qJ(2));
t212 = t395 * t399 + t396 * t398;
t233 = sin(qJ(4));
t253 = t395 * t396 - t398 * t399;
t397 = cos(qJ(4));
t245 = -t212 * t397 + t233 * t253;
t439 = t245 * t215;
t441 = Ifges(5,6) * t245;
t173 = t212 * t233 + t397 * t253;
t450 = Ifges(6,6) * t245 + t173 * t289;
t461 = t234 * t450;
t449 = Ifges(6,5) * t245 + t173 * t290;
t462 = t232 * t449;
t468 = t439 / 0.2e1 - t441 + t461 / 0.2e1 + t462 / 0.2e1;
t414 = -t173 / 0.4e1;
t457 = t173 / 0.4e1;
t467 = t343 * t457 + t346 * t414 + t461 / 0.4e1 + t462 / 0.4e1;
t466 = -t449 / 0.4e1;
t465 = -t450 / 0.4e1;
t456 = Ifges(5,5) * t173;
t464 = -t456 / 0.2e1;
t372 = t232 * mrSges(6,1);
t391 = mrSges(6,2) * t234;
t291 = t372 + t391;
t455 = t291 * t173;
t463 = pkin(4) * t455;
t460 = -t439 / 0.4e1 + t441 / 0.2e1 + t464;
t415 = -t173 / 0.2e1;
t228 = Ifges(6,5) * t234;
t381 = Ifges(6,6) * t232;
t437 = t228 - t381;
t447 = t245 / 0.2e1;
t458 = t400 * t449 + t401 * t450 + t437 * t447 - Ifges(5,4) * t245 + (Ifges(5,2) + Ifges(6,3)) * t415;
t384 = Ifges(6,5) * t173;
t330 = -t384 / 0.2e1;
t321 = t395 * t233;
t204 = (t397 * t398 - t321) * pkin(2);
t230 = t232 ^ 2;
t231 = t234 ^ 2;
t342 = t230 + t231;
t303 = t342 * t204;
t358 = t173 * t234;
t444 = mrSges(6,1) * t245;
t113 = -mrSges(6,3) * t358 + t444;
t362 = t113 * t232;
t347 = t232 * t173;
t443 = mrSges(6,2) * t245;
t110 = -mrSges(6,3) * t347 - t443;
t363 = t110 * t234;
t454 = t363 / 0.2e1 - t362 / 0.2e1;
t120 = pkin(4) * t245 - t173 * pkin(6);
t226 = pkin(2) * t398 + pkin(3);
t198 = -pkin(2) * t321 + t226 * t397;
t193 = -pkin(4) - t198;
t411 = t193 / 0.2e1;
t390 = mrSges(6,3) * t232;
t111 = mrSges(6,2) * t173 - t245 * t390;
t344 = t234 * t111;
t389 = mrSges(6,3) * t234;
t114 = -mrSges(6,1) * t173 - t245 * t389;
t348 = t232 * t114;
t436 = t344 / 0.2e1 - t348 / 0.2e1;
t109 = -t173 * t390 - t443;
t345 = t234 * t109;
t112 = -t173 * t389 + t444;
t349 = t232 * t112;
t435 = -t349 / 0.2e1 + t345 / 0.2e1;
t367 = t234 * mrSges(6,1);
t371 = t232 * mrSges(6,2);
t433 = -t371 / 0.2e1 + t367 / 0.2e1;
t382 = Ifges(6,6) * t173;
t49 = t245 * t289 - t382;
t370 = t232 * t49;
t52 = t245 * t290 - t384;
t431 = t52 * t400 - t370 / 0.2e1 + Ifges(5,1) * t447 + Ifges(5,4) * t173;
t429 = (-t395 * mrSges(4,1) - t398 * mrSges(4,2)) * pkin(2);
t428 = mrSges(6,3) * t303;
t427 = (-t228 / 0.2e1 + t381 / 0.2e1) * t173;
t426 = t463 / 0.2e1 + t460;
t425 = -pkin(6) / 0.2e1;
t424 = m(6) * pkin(3);
t423 = -mrSges(6,1) / 0.2e1;
t422 = mrSges(6,2) / 0.2e1;
t421 = Ifges(6,3) / 0.2e1;
t420 = -t49 / 0.4e1;
t419 = t52 / 0.4e1;
t418 = -t111 / 0.2e1;
t417 = t111 / 0.2e1;
t413 = t173 / 0.2e1;
t297 = t395 * t397;
t199 = pkin(2) * t297 + t233 * t226;
t194 = pkin(6) + t199;
t410 = -t194 / 0.2e1;
t409 = t198 / 0.2e1;
t203 = (t233 * t398 + t297) * pkin(2);
t408 = t203 / 0.2e1;
t407 = -t204 / 0.2e1;
t406 = -t216 / 0.4e1;
t405 = -t217 / 0.4e1;
t394 = pkin(3) * t233;
t224 = pkin(6) + t394;
t404 = -t224 / 0.2e1;
t403 = t224 / 0.2e1;
t339 = t397 * pkin(3);
t225 = -t339 - pkin(4);
t402 = t225 / 0.2e1;
t392 = t212 * pkin(3);
t227 = t399 * pkin(2) + pkin(1);
t380 = Ifges(6,3) * t245;
t116 = mrSges(5,1) * t245 + mrSges(5,2) * t173;
t117 = -mrSges(5,1) * t173 + mrSges(5,2) * t245;
t187 = -pkin(3) * t253 + t227;
t338 = t396 * pkin(2);
t192 = -t338 - t392;
t333 = t421 + Ifges(5,2) / 0.2e1;
t236 = t227 * (-t212 * mrSges(4,1) + mrSges(4,2) * t253) + t253 ^ 2 * Ifges(4,4) + (t415 * t437 + t431) * t173 + (Ifges(5,1) * t413 - t173 * t333 + t458) * t245;
t239 = -Ifges(4,4) * t212 + (-Ifges(4,1) + Ifges(4,2)) * t253;
t350 = t232 * t111;
t281 = t234 * t114 + t350;
t282 = t232 * t110 + t234 * t113;
t306 = m(6) * t342;
t78 = -pkin(4) * t173 - pkin(6) * t245 + t187;
t80 = t120 - t392;
t79 = -t338 + t80;
t1 = t192 * t117 + t281 * t79 + (m(5) * t192 + t116) * t187 + (t306 * t79 + t282) * t78 + (mrSges(4,2) * t338 + t239) * t212 + t236 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t399) * t399 + (-pkin(1) * mrSges(3,1) + (-m(4) * t227 + mrSges(4,1) * t253) * pkin(2) + (-Ifges(3,2) + Ifges(3,1)) * t399 - Ifges(3,4) * t396) * t396;
t379 = t1 * qJD(1);
t378 = t198 * mrSges(5,2);
t377 = t199 * mrSges(5,1);
t2 = t187 * t116 + t281 * t80 + (t306 * t80 + t282) * t78 + ((-m(5) * t187 - t117) * pkin(3) + t239) * t212 + t236;
t376 = t2 * qJD(1);
t375 = t203 * mrSges(5,1);
t374 = t204 * mrSges(5,2);
t3 = (t232 * t109 + t234 * t112) * t78 + (t187 * mrSges(5,1) + t458) * t245 - (-t187 * mrSges(5,2) - t427 + (-Ifges(5,1) / 0.2e1 + t333) * t245 - t431) * t173 + (t306 * t78 + t281) * t120;
t364 = t3 * qJD(1);
t16 = (t232 * t52 / 0.2e1 + t49 * t400 + t215 * t415 - t273 * t245) * t245 + (t348 - t344) * t78;
t361 = t16 * qJD(1);
t359 = t245 * t233;
t356 = t193 * t455;
t355 = t198 * t173;
t354 = t199 * t245;
t292 = t367 - t371;
t353 = t199 * t292;
t352 = t203 * t292;
t351 = t225 * t455;
t341 = t292 * t394;
t340 = -t397 / 0.2e1;
t337 = -t394 / 0.2e1;
t105 = t292 * t245;
t336 = -pkin(4) * t105 / 0.2e1;
t335 = t114 * t425;
t334 = Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.4e1;
t322 = t397 * t173;
t318 = t228 * t414;
t317 = t105 * t411;
t108 = t291 * t245;
t316 = t108 * t408;
t315 = t105 * t402;
t314 = t455 * t402;
t313 = -t350 / 0.2e1;
t310 = t233 * t108 / 0.2e1;
t307 = t403 + t410;
t304 = t342 * t198;
t302 = t342 * t224;
t301 = mrSges(6,3) * t339;
t300 = -t339 / 0.2e1;
t296 = t397 * t417;
t295 = t232 * t340;
t294 = mrSges(6,3) * (-t231 / 0.2e1 - t230 / 0.2e1);
t276 = pkin(4) * t291;
t287 = -t276 / 0.2e1 + t299;
t249 = mrSges(6,3) * t304 - t353 - t377 - t378;
t19 = m(6) * (t193 * t199 + t194 * t304) + t249;
t254 = t454 * pkin(6) - t426 + t467;
t274 = -t455 * t411 - t199 * t108 / 0.2e1;
t5 = (t109 * t410 + t198 * t418 + t217 * t414 + t465) * t234 + (t466 + t194 * t112 / 0.2e1 + t114 * t409 - t173 * t406) * t232 + (Ifges(5,6) / 0.2e1 - t215 / 0.4e1) * t245 + t254 + t464 + t274;
t286 = -t5 * qJD(1) + t19 * qJD(2);
t268 = -t234 * t334 + t406;
t248 = (t268 - t387) * t245 + t419;
t244 = t114 * t410 + t248;
t267 = t232 * t334 + t405;
t246 = (t457 + t413) * Ifges(6,6) + t267 * t245 + t420;
t275 = -t380 / 0.2e1 + t318;
t279 = t194 * t294;
t10 = t317 + t245 * t279 + (t111 * t410 + t422 * t79 + t246) * t232 + (t423 * t79 + t244 + t330) * t234 + t275;
t272 = t193 * t291;
t90 = t272 + t299;
t285 = t10 * qJD(1) + t90 * qJD(2);
t247 = -t354 / 0.2e1 + t245 * t408 - t173 * t407 - t355 / 0.2e1;
t15 = t316 + (t411 - t225 / 0.2e1) * t455 + (-t110 * t307 + t204 * t417) * t234 + (t113 * t307 + t114 * t407) * t232 + ((t359 / 0.2e1 + t322 / 0.2e1) * pkin(3) + t247) * mrSges(5,3);
t20 = (-mrSges(5,1) - t292) * t203 + t429 + (mrSges(6,3) * t342 - mrSges(5,2)) * t204 + m(6) * (t193 * t203 + t194 * t303) + m(5) * (-t198 * t203 + t199 * t204);
t284 = t15 * qJD(1) + t20 * qJD(2);
t283 = -t362 + t363;
t280 = pkin(6) * t294;
t278 = t224 * t294;
t277 = t372 / 0.2e1 + t391 / 0.2e1;
t271 = t225 * t291;
t266 = -t460 + t467;
t265 = t276 / 0.2e1;
t264 = t342 * t397;
t263 = t277 * t198;
t262 = t277 * t204;
t261 = -t272 / 0.2e1;
t235 = m(6) * (t225 * t199 + t198 * t302 + (t193 * t233 + t194 * t264) * pkin(3)) / 0.2e1 - t378 / 0.2e1 - t377 / 0.2e1 - t353 / 0.2e1 + mrSges(5,1) * t337 - t341 / 0.2e1 + mrSges(5,2) * t300 + t342 * (mrSges(6,3) * t409 + t301 / 0.2e1);
t238 = -m(6) * (-pkin(4) * t203 + pkin(6) * t303) / 0.2e1 + t375 / 0.2e1 + t352 / 0.2e1 + t374 / 0.2e1 - t428 / 0.2e1;
t18 = t235 + t238;
t7 = t314 + pkin(3) * t310 + (pkin(3) * t296 + t109 * t403 + t110 * t425 + t173 * t405 + t465) * t234 + (t114 * t300 + t466 + t216 * t457 + t112 * t404 + pkin(6) * t113 / 0.2e1) * t232 + t266 + t426;
t241 = -mrSges(5,1) * t394 - mrSges(5,2) * t339 + t301 * t342 - t341;
t91 = (t224 * t264 + t225 * t233) * t424 + t241;
t259 = t7 * qJD(1) + t18 * qJD(2) + t91 * qJD(3);
t115 = t271 + t299;
t243 = t114 * t404 + t248;
t12 = t315 + t245 * t278 + (t111 * t404 + t422 * t80 + t246) * t232 + (t423 * t80 + t243 + t330) * t234 + t275;
t237 = -t271 / 0.2e1 - t299;
t31 = t261 - t262 + t237;
t258 = t12 * qJD(1) - t31 * qJD(2) + t115 * qJD(3);
t257 = t267 * t232;
t256 = Ifges(4,5) * t253 + Ifges(4,6) * t212 - t173 * t273 + t456 + t468;
t255 = t437 * t414 - t370 / 0.4e1;
t121 = t276 - t299;
t14 = t336 + t318 + (-Ifges(6,3) / 0.2e1 + t280) * t245 + (t120 * t423 + t268 * t245 + t330 + t335 + t419) * t234 + (0.3e1 / 0.4e1 * t382 + pkin(6) * t418 + t420 + t120 * t422 + (t267 - t386) * t245) * t232;
t33 = t265 + t261 - t263 - t299;
t251 = (mrSges(6,1) * t295 + t340 * t391) * pkin(3);
t57 = t265 + t251 + t237;
t252 = t14 * qJD(1) - t33 * qJD(2) - t57 * qJD(3) - t121 * qJD(4);
t250 = -Ifges(6,6) * t347 / 0.2e1 + Ifges(6,5) * t358 / 0.2e1 + t380 / 0.2e1 + t255;
t242 = t254 + t266;
t197 = t271 / 0.2e1;
t176 = t272 / 0.2e1;
t58 = t197 + t251 + t287;
t34 = t176 - t263 + t287;
t32 = t197 + t176 - t262 + t299;
t17 = t235 - t238;
t13 = t336 + pkin(6) * t313 + (t335 + t248) * t234 + t255 + (t421 + t280 + t257) * t245 - t427 + t433 * t120;
t11 = t315 + t224 * t313 + (t278 + t257) * t245 + t243 * t234 + t250 + t433 * t80;
t9 = t194 * t313 + t317 + (t279 + t257) * t245 + t244 * t234 + t250 + t433 * t79;
t8 = t256 + t351 / 0.2e1 + t316 + t356 / 0.2e1 + t436 * t204 + (t173 * t300 + t245 * t337 + t247) * mrSges(5,3) + (t224 + t194) * t454;
t6 = t435 * t224 + (t114 * t295 + t234 * t296 + t310) * pkin(3) + t314 + t242;
t4 = t194 * t435 + t198 * t436 + t242 - t274;
t21 = [qJD(2) * t1 + qJD(3) * t2 + qJD(4) * t3 - qJD(5) * t16, t379 + (t256 + (t212 * t395 - t253 * t398) * pkin(2) * mrSges(4,3) + (-t354 - t355) * mrSges(5,3) + t283 * t194 - Ifges(3,5) * t399 + Ifges(3,6) * t396 + t356) * qJD(2) + t8 * qJD(3) + t4 * qJD(4) + t9 * qJD(5), t376 + t8 * qJD(2) + (t351 + t283 * t224 + (-t322 - t359) * pkin(3) * mrSges(5,3) + t256) * qJD(3) + t6 * qJD(4) + t11 * qJD(5), t4 * qJD(2) + t6 * qJD(3) + t13 * qJD(5) + t364 + (-t463 - (-Ifges(5,5) + t273) * t173 + (t345 - t349) * pkin(6) + t468) * qJD(4), -t361 + t9 * qJD(2) + t11 * qJD(3) + t13 * qJD(4) + (-t291 * t78 - t439) * qJD(5); qJD(3) * t15 - qJD(4) * t5 + qJD(5) * t10 - t379, qJD(3) * t20 + qJD(4) * t19 + qJD(5) * t90 (m(6) * (t225 * t203 + t204 * t302) - t352 + m(5) * (-t203 * t397 + t204 * t233) * pkin(3) - t374 - t375 + t428 + t429) * qJD(3) + t17 * qJD(4) + t32 * qJD(5) + t284, t17 * qJD(3) + (m(6) * (-pkin(4) * t199 + pkin(6) * t304) + t249) * qJD(4) + t34 * qJD(5) + t286, t32 * qJD(3) + t34 * qJD(4) + (-t194 * t292 + t437) * qJD(5) + t285; -qJD(2) * t15 + qJD(4) * t7 + qJD(5) * t12 - t376, qJD(4) * t18 - qJD(5) * t31 - t284, qJD(4) * t91 + qJD(5) * t115 ((-pkin(4) * t233 + pkin(6) * t264) * t424 + t241) * qJD(4) + t58 * qJD(5) + t259, t58 * qJD(4) + (-t224 * t292 + t437) * qJD(5) + t258; qJD(2) * t5 - qJD(3) * t7 + qJD(5) * t14 - t364, -qJD(3) * t18 - qJD(5) * t33 - t286, -qJD(5) * t57 - t259, -t121 * qJD(5) (-pkin(6) * t292 + t437) * qJD(5) + t252; -qJD(2) * t10 - qJD(3) * t12 - qJD(4) * t14 + t361, qJD(3) * t31 + qJD(4) * t33 - t285, qJD(4) * t57 - t258, -t252, 0;];
Cq  = t21;
