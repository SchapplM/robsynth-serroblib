% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:04:14
% EndTime: 2019-03-08 21:04:23
% DurationCPUTime: 4.30s
% Computational Cost: add. (8268->396), mult. (18761->557), div. (0->0), fcn. (20496->10), ass. (0->238)
t236 = sin(pkin(11));
t240 = sin(qJ(3));
t243 = cos(qJ(3));
t344 = cos(pkin(11));
t294 = t344 * t243;
t206 = t236 * t240 - t294;
t208 = t236 * t243 + t240 * t344;
t228 = -pkin(3) * t243 - pkin(2);
t273 = -qJ(5) * t208 + t228;
t146 = pkin(4) * t206 + t273;
t158 = -mrSges(6,2) * t206 - mrSges(6,3) * t208;
t428 = m(6) * t146 + t158;
t242 = cos(qJ(6));
t234 = t242 ^ 2;
t357 = t234 * mrSges(7,3);
t239 = sin(qJ(6));
t232 = t239 ^ 2;
t358 = t232 * mrSges(7,3);
t427 = -t357 / 0.2e1 - t358 / 0.2e1;
t395 = m(6) / 0.2e1;
t426 = 0.2e1 * t395;
t425 = m(6) + m(5);
t424 = m(6) + m(7);
t371 = mrSges(5,3) + mrSges(6,1);
t422 = t206 * t371;
t421 = t208 * t371;
t391 = m(5) * pkin(3);
t310 = t391 / 0.2e1;
t420 = t236 * t310;
t237 = sin(pkin(6));
t244 = cos(qJ(2));
t325 = t237 * t244;
t419 = t325 * t395;
t370 = -qJ(4) - pkin(8);
t215 = t370 * t240;
t217 = t370 * t243;
t169 = -t344 * t215 - t217 * t236;
t116 = pkin(5) * t208 + t169;
t350 = t242 * mrSges(7,1);
t355 = t239 * mrSges(7,2);
t214 = -t350 + t355;
t142 = t214 * t206;
t143 = t214 * t208;
t230 = t240 * pkin(3);
t292 = qJ(5) * t206 + t230;
t147 = pkin(4) * t208 + t292;
t361 = t208 * mrSges(7,3);
t149 = -mrSges(7,1) * t206 - t239 * t361;
t151 = mrSges(7,2) * t206 + t242 * t361;
t238 = cos(pkin(6));
t241 = sin(qJ(2));
t326 = t237 * t241;
t200 = t238 * t243 - t240 * t326;
t201 = t238 * t240 + t243 * t326;
t262 = t236 * t200 + t201 * t344;
t295 = t344 * t200;
t328 = t236 * t201;
t263 = -t295 + t328;
t332 = t206 * t239;
t148 = t208 * mrSges(7,1) - mrSges(7,3) * t332;
t320 = t242 * t148;
t331 = t206 * t242;
t150 = -t208 * mrSges(7,2) + mrSges(7,3) * t331;
t323 = t239 * t150;
t264 = t320 / 0.2e1 + t323 / 0.2e1;
t386 = pkin(4) + pkin(9);
t99 = t206 * t386 + t273;
t55 = t116 * t242 - t239 * t99;
t56 = t116 * t239 + t242 * t99;
t282 = t239 * t56 + t242 * t55;
t306 = t240 * t325;
t384 = t262 / 0.2e1;
t92 = t239 * t325 + t242 * t263;
t387 = t92 / 0.2e1;
t393 = m(7) / 0.2e1;
t396 = m(5) / 0.2e1;
t407 = t236 * t215 - t344 * t217;
t415 = -t206 * pkin(5) + t407;
t100 = t208 * t386 + t292;
t57 = -t100 * t239 + t242 * t415;
t58 = t100 * t242 + t239 * t415;
t93 = -t239 * t263 + t242 * t325;
t418 = -pkin(3) * t306 * t396 - t147 * t419 + (-t415 * t263 + t57 * t92 - t58 * t93 + (-t116 + t282) * t262) * t393 - t263 * t142 / 0.2e1 + t143 * t384 + t149 * t387 - t93 * t151 / 0.2e1 + t264 * t262;
t348 = t242 * t57;
t354 = t239 * t58;
t281 = t348 + t354;
t417 = -mrSges(4,1) * t243 + mrSges(4,2) * t240;
t416 = t395 + t396;
t414 = t236 * t391 - mrSges(5,2);
t303 = t344 * pkin(3);
t227 = -t303 - pkin(4);
t413 = t227 * t395 - t344 * t310;
t412 = m(6) * t227 - t344 * t391 - mrSges(5,1) + mrSges(6,2);
t349 = t242 * mrSges(7,2);
t356 = t239 * mrSges(7,1);
t216 = t349 + t356;
t410 = mrSges(6,3) + t216;
t155 = -t208 * mrSges(6,2) + t206 * mrSges(6,3);
t156 = t208 * mrSges(5,1) - t206 * mrSges(5,2);
t408 = t155 + t156;
t406 = t427 * t206;
t176 = t208 * t325;
t140 = t176 * t242 - t239 * t326;
t141 = t176 * t239 + t242 * t326;
t277 = t140 * t242 + t141 * t239;
t218 = t240 * mrSges(4,1) + t243 * mrSges(4,2);
t364 = Ifges(7,6) * t242;
t365 = Ifges(7,5) * t239;
t270 = t365 / 0.2e1 + t364 / 0.2e1;
t405 = Ifges(5,4) + Ifges(6,6) - t270;
t402 = (t419 + (t239 * t92 + t242 * t93) * t393) * t208;
t401 = -t142 + t422;
t397 = 0.2e1 * t208;
t394 = m(6) / 0.4e1;
t392 = m(7) / 0.4e1;
t390 = mrSges(7,1) / 0.2e1;
t389 = -mrSges(7,2) / 0.2e1;
t385 = t415 / 0.2e1;
t383 = t176 / 0.2e1;
t381 = -t206 / 0.2e1;
t380 = -t208 / 0.2e1;
t379 = t208 / 0.2e1;
t223 = -pkin(9) + t227;
t378 = -t223 / 0.2e1;
t377 = t223 / 0.2e1;
t376 = -t239 / 0.2e1;
t375 = t239 / 0.2e1;
t374 = t242 / 0.2e1;
t373 = pkin(3) * t236;
t367 = Ifges(7,4) * t239;
t366 = Ifges(7,4) * t242;
t363 = t206 * mrSges(6,1);
t362 = t206 * mrSges(5,3);
t360 = t208 * Ifges(7,5);
t359 = t208 * Ifges(7,6);
t353 = t239 * t93;
t197 = Ifges(7,4) * t331;
t97 = Ifges(7,1) * t332 + t197 + t360;
t352 = t239 * t97;
t347 = t242 * t92;
t286 = Ifges(7,2) * t242 + t367;
t95 = t206 * t286 + t359;
t346 = t242 * t95;
t280 = t347 - t353;
t11 = m(7) * (-t263 + t280) * t262;
t343 = t11 * qJD(1);
t177 = -t236 * t306 + t294 * t325;
t77 = t262 * t177;
t341 = t262 * t214;
t224 = qJ(5) + t373;
t312 = t232 + t234;
t293 = t312 * t223;
t248 = t216 * t381 + (t293 * t393 + t413 + t427) * t208 + (-t420 + (-t395 - t393) * t224) * t206;
t249 = t147 * t395 + (-t239 * t57 + t242 * t58) * t393 + t149 * t376 + t151 * t374 + t240 * t310;
t13 = t248 - t249 - t408;
t340 = t13 * qJD(2);
t329 = t237 ^ 2 * t241;
t14 = m(7) * (t140 * t92 - t141 * t93 + t77) + m(4) * (-t329 + (-t200 * t240 + t201 * t243) * t237) * t244 + t425 * (t176 * t263 - t244 * t329 + t77);
t339 = t14 * qJD(1);
t335 = t177 * t224;
t334 = t206 * t262;
t259 = (t356 / 0.2e1 + t349 / 0.2e1) * t208;
t318 = t242 * t150;
t324 = t239 * t148;
t266 = t324 / 0.2e1 - t318 / 0.2e1;
t296 = t232 / 0.2e1 + t234 / 0.2e1;
t288 = mrSges(7,3) * t296;
t23 = t206 * t288 + t259 + t266;
t330 = t23 * qJD(2);
t322 = t239 * t151;
t220 = Ifges(7,1) * t242 - t367;
t321 = t239 * t220;
t319 = t242 * t149;
t219 = -Ifges(7,2) * t239 + t366;
t317 = t242 * t219;
t268 = t355 / 0.2e1 - t350 / 0.2e1;
t257 = t268 * t208;
t25 = -t257 + t264;
t316 = t25 * qJD(2);
t49 = (m(7) * t296 + t395) * t397;
t315 = t49 * qJD(2);
t311 = t240 ^ 2 + t243 ^ 2;
t309 = t392 + t394;
t308 = -Ifges(7,2) / 0.4e1 + Ifges(7,1) / 0.4e1;
t299 = -t325 / 0.2e1;
t287 = Ifges(7,1) * t239 + t366;
t285 = -t364 - t365;
t246 = t335 * t395 + (t223 * t277 + t335) * t393 + mrSges(6,2) * t383 - t277 * mrSges(7,3) / 0.2e1 + t218 * t299 + (-mrSges(5,1) / 0.2e1 + t413) * t176 + (-mrSges(5,2) / 0.2e1 + t420 + t410 / 0.2e1) * t177;
t2 = -t246 + (-t362 / 0.2e1 - t363 / 0.2e1) * t263 + (t218 + t408) * t299 + t371 * (-t263 * t381 + (t379 + t380) * t262) + t418;
t157 = mrSges(5,1) * t206 + mrSges(5,2) * t208;
t96 = -Ifges(7,6) * t206 + t208 * t286;
t98 = -Ifges(7,5) * t206 + t208 * t287;
t3 = -pkin(2) * t218 - t116 * t142 + t415 * t143 + t146 * t155 + t57 * t148 + t55 * t149 + t58 * t150 + t56 * t151 + (-Ifges(4,4) * t240 + pkin(3) * t157) * t240 + m(7) * (-t116 * t415 + t55 * t57 + t56 * t58) + (t405 * t206 + t96 * t374 + t98 * t375) * t206 + (Ifges(4,4) * t243 + (Ifges(4,1) - Ifges(4,2)) * t240) * t243 + (t346 / 0.2e1 + t352 / 0.2e1 - t405 * t208 + (-Ifges(7,3) - Ifges(5,1) + Ifges(6,3) + Ifges(5,2) - Ifges(6,2)) * t206) * t208 + (m(5) * t230 + t156) * t228 + t428 * t147;
t284 = t2 * qJD(1) + t3 * qJD(2);
t144 = t206 * t216;
t196 = Ifges(7,5) * t331;
t4 = t55 * t150 - t56 * t148 + t196 * t379 + t415 * t144 + ((-t55 * mrSges(7,3) + t97 / 0.2e1 + t197 / 0.2e1) * t242 + (-t56 * mrSges(7,3) - t359 / 0.2e1 - t95 / 0.2e1 + (-t367 / 0.2e1 + (Ifges(7,1) / 0.2e1 - Ifges(7,2) / 0.2e1) * t242) * t206) * t239) * t206;
t250 = (t353 / 0.2e1 - t347 / 0.2e1) * t206 * mrSges(7,3) + t144 * t384 + t150 * t387 + t93 * t148 / 0.2e1;
t269 = t140 * t390 + t141 * t389;
t7 = t250 - t269;
t283 = t7 * qJD(1) + t4 * qJD(2);
t10 = t401 * t206 + (t320 + t323 + t421) * t208 + m(7) * (-t206 * t415 + t208 * t282) + t425 * (t169 * t208 - t206 * t407);
t247 = (t208 * t280 - t334) * t393 + t416 * (t208 * t263 - t334);
t252 = t416 * t326 + (-t140 * t239 + t141 * t242) * t393;
t16 = -t247 + t252;
t279 = qJD(1) * t16 - qJD(2) * t10;
t19 = (m(7) * (t239 * t55 - t242 * t56) - t318 + t324 - t428) * t208;
t253 = m(6) * t383 + t277 * t393;
t30 = -t253 + t402;
t278 = qJD(1) * t30 + qJD(2) * t19;
t275 = t169 * t176 + t177 * t407;
t272 = t312 * t392 + t394;
t271 = t389 * t58 + t390 * t57;
t267 = t214 * t385 - t224 * t144 / 0.2e1;
t265 = -t322 / 0.2e1 - t319 / 0.2e1;
t261 = t272 * t262;
t258 = t268 * t262;
t21 = t341 / 0.2e1 - t258;
t5 = (-Ifges(7,3) / 0.2e1 + t223 * t288) * t206 + (0.3e1 / 0.4e1 * t360 + t97 / 0.4e1 + t197 / 0.4e1 + t148 * t377 + (t219 / 0.4e1 + t308 * t239) * t206) * t239 + (0.3e1 / 0.4e1 * t359 + t95 / 0.4e1 + t150 * t378 + (0.3e1 / 0.4e1 * t367 - t220 / 0.4e1 - t308 * t242) * t206) * t242 + t267 + t271;
t75 = -t224 * t214 - t321 / 0.2e1 - t242 * t287 / 0.2e1 - t317 / 0.2e1 + t286 * t375;
t256 = t21 * qJD(1) + t5 * qJD(2) - t75 * qJD(3);
t152 = 0.4e1 * t224 * t309 + t410;
t20 = t268 * t206 + 0.2e1 * (t415 / 0.4e1 - t354 / 0.4e1 - t348 / 0.4e1) * m(7) + t265;
t32 = -0.2e1 * t262 * t309 + 0.2e1 * t261;
t254 = qJD(1) * t32 - qJD(2) * t20 - qJD(3) * t152;
t50 = t272 * t397 + (m(7) * t312 + m(6)) * t380;
t33 = t384 * t424 + 0.2e1 * t261;
t29 = t402 + t253;
t26 = -t257 - t264;
t24 = t259 - t266 + t406;
t22 = -t341 / 0.2e1 - t258;
t18 = t281 * t393 + m(7) * t385 + (-mrSges(6,1) + t268) * t206 - t265 + t407 * t426;
t17 = t247 + t252;
t15 = t248 + t249;
t8 = t250 + t269;
t6 = -t352 / 0.4e1 - t346 / 0.4e1 - t239 * (-Ifges(7,2) * t332 + t197) / 0.4e1 + t318 * t377 + t324 * t378 + Ifges(7,3) * t381 - t267 + t271 + (t220 / 0.2e1 - t286 / 0.4e1) * t331 + t406 * t223 + (t285 / 0.4e1 + t270) * t208 - (t287 + t219) * t332 / 0.4e1;
t1 = (-t218 / 0.2e1 - t156 / 0.2e1 - t155 / 0.2e1) * t325 + t246 + (t263 / 0.2e1 + t295 / 0.2e1 - t328 / 0.2e1) * t422 + t418;
t9 = [t14 * qJD(2) + t11 * qJD(3), t1 * qJD(3) + t17 * qJD(4) + t29 * qJD(5) + t8 * qJD(6) + t339 + (t140 * t148 + t141 * t150 - t401 * t177 + (t146 * t326 + t275) * t426 + 0.2e1 * (t140 * t55 + t141 * t56 + t177 * t415) * t393 + 0.2e1 * (t228 * t326 + t275) * t396 + t176 * t421 + ((mrSges(4,3) * t311 - mrSges(3,2)) * t244 + (-mrSges(3,1) + t157 + t158 + t417) * t241 + m(4) * (pkin(8) * t244 * t311 - pkin(2) * t241)) * t237) * qJD(2), t1 * qJD(2) + t33 * qJD(5) + t22 * qJD(6) + t343 + (-t201 * mrSges(4,1) - t200 * mrSges(4,2) + (m(7) * t293 - t357 - t358 + t412) * t262 + (-t224 * t424 - t410 - t414) * t263) * qJD(3), qJD(2) * t17, qJD(2) * t29 + qJD(3) * t33, t8 * qJD(2) + t22 * qJD(3) + (mrSges(7,1) * t93 - mrSges(7,2) * t92) * qJD(6); qJD(3) * t2 - qJD(4) * t16 + qJD(5) * t30 + qJD(6) * t7 - t339, qJD(3) * t3 + qJD(4) * t10 + qJD(5) * t19 + qJD(6) * t4, t15 * qJD(4) + t18 * qJD(5) + t6 * qJD(6) + t284 + ((-t224 * mrSges(6,1) - mrSges(5,3) * t373 + Ifges(6,5) - Ifges(5,6)) * t208 + Ifges(4,5) * t243 - Ifges(4,6) * t240 + t224 * t143 + t303 * t362 + t98 * t374 + t96 * t376 + (Ifges(7,5) * t242 - Ifges(7,6) * t239) * t381 - t227 * t363 + (t321 + t317) * t379 + (-m(7) * t224 - t216) * t116 + (m(7) * t281 + t319 + t322) * t223 + (-Ifges(5,5) + Ifges(6,4)) * t206 + (-m(6) * t224 - mrSges(6,3) - t414) * t169 + t412 * t407 + t417 * pkin(8) - t281 * mrSges(7,3)) * qJD(3), qJD(3) * t15 + qJD(5) * t50 + qJD(6) * t26 - t279, qJD(3) * t18 + qJD(4) * t50 + qJD(6) * t24 + t278, t6 * qJD(3) + t26 * qJD(4) + t24 * qJD(5) + (-mrSges(7,1) * t56 - mrSges(7,2) * t55 - Ifges(7,6) * t332 + t196) * qJD(6) + t283; -qJD(2) * t2 - qJD(5) * t32 - qJD(6) * t21 - t343, qJD(4) * t13 + qJD(5) * t20 - qJD(6) * t5 - t284, qJD(5) * t152 + qJD(6) * t75, t340, -t254 (-t216 * t223 + t285) * qJD(6) - t256; qJD(2) * t16, -qJD(3) * t13 - qJD(5) * t49 - qJD(6) * t25 + t279, -t340, 0, -t315, t214 * qJD(6) - t316; -qJD(2) * t30 + qJD(3) * t32, -qJD(3) * t20 + qJD(4) * t49 - qJD(6) * t23 - t278, t254, t315, 0, -t216 * qJD(6) - t330; -t7 * qJD(2) + t21 * qJD(3), qJD(3) * t5 + qJD(4) * t25 + qJD(5) * t23 - t283, t256, t316, t330, 0;];
Cq  = t9;
