% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:08:00
% EndTime: 2019-03-09 10:08:27
% DurationCPUTime: 15.92s
% Computational Cost: add. (16750->650), mult. (44058->894), div. (0->0), fcn. (34230->10), ass. (0->316)
t288 = sin(pkin(11));
t290 = cos(pkin(11));
t291 = sin(qJ(6));
t294 = cos(qJ(6));
t259 = t288 * t294 + t290 * t291;
t248 = t259 * qJD(6);
t289 = sin(pkin(10));
t293 = sin(qJ(2));
t295 = cos(qJ(2));
t355 = cos(pkin(10));
t256 = -t289 * t293 + t295 * t355;
t303 = qJD(1) * t256;
t389 = cos(qJ(4));
t229 = t389 * t303;
t320 = t355 * t293;
t340 = qJD(1) * t295;
t244 = -qJD(1) * t320 - t289 * t340;
t292 = sin(qJ(4));
t194 = t244 * t292 + t229;
t448 = t259 * t194;
t425 = t448 - t248;
t309 = t288 * t291 - t290 * t294;
t247 = t309 * qJD(6);
t447 = t309 * t194;
t424 = t447 - t247;
t189 = qJD(6) - t194;
t287 = qJD(2) + qJD(4);
t299 = t292 * t303;
t297 = -t244 * t389 + t299;
t173 = t287 * t288 + t290 * t297;
t317 = t290 * t287 - t288 * t297;
t117 = t173 * t294 + t291 * t317;
t373 = Ifges(7,4) * t117;
t445 = -t173 * t291 + t294 * t317;
t59 = Ifges(7,2) * t445 + Ifges(7,6) * t189 + t373;
t408 = t59 / 0.2e1;
t111 = Ifges(7,4) * t445;
t60 = Ifges(7,1) * t117 + Ifges(7,5) * t189 + t111;
t406 = t60 / 0.2e1;
t383 = -qJ(3) - pkin(7);
t272 = t383 * t295;
t264 = qJD(1) * t272;
t249 = t289 * t264;
t271 = t383 * t293;
t263 = qJD(1) * t271;
t255 = qJD(2) * pkin(2) + t263;
t197 = t355 * t255 + t249;
t387 = pkin(8) * t244;
t168 = qJD(2) * pkin(3) + t197 + t387;
t321 = t355 * t264;
t198 = t289 * t255 - t321;
t301 = pkin(8) * t303;
t170 = t301 + t198;
t107 = t292 * t168 + t170 * t389;
t101 = t287 * qJ(5) + t107;
t203 = -qJD(1) * pkin(1) - pkin(2) * t340 - pkin(3) * t303 + qJD(3);
t112 = -pkin(4) * t194 - qJ(5) * t297 + t203;
t68 = -t101 * t288 + t290 * t112;
t33 = -pkin(5) * t194 - pkin(9) * t173 + t68;
t69 = t290 * t101 + t288 * t112;
t43 = pkin(9) * t317 + t69;
t12 = -t291 * t43 + t294 * t33;
t13 = t291 * t33 + t294 * t43;
t258 = t289 * t295 + t320;
t245 = t258 * qJD(2);
t232 = qJD(1) * t245;
t246 = t256 * qJD(2);
t233 = qJD(1) * t246;
t338 = qJD(4) * t292;
t148 = qJD(4) * t229 - t292 * t232 + t233 * t389 + t244 * t338;
t328 = qJD(4) * t389;
t149 = qJD(4) * t299 + t232 * t389 + t292 * t233 - t244 * t328;
t322 = qJD(2) * t383;
t241 = qJD(3) * t295 + t293 * t322;
t220 = t241 * qJD(1);
t242 = -t293 * qJD(3) + t295 * t322;
t221 = t242 * qJD(1);
t176 = t355 * t220 + t289 * t221;
t160 = -pkin(8) * t232 + t176;
t175 = -t289 * t220 + t221 * t355;
t302 = -t233 * pkin(8) + t175;
t55 = t389 * t160 + t168 * t328 - t170 * t338 + t292 * t302;
t51 = qJD(5) * t287 + t55;
t337 = qJD(1) * qJD(2);
t327 = t293 * t337;
t279 = pkin(2) * t327;
t204 = pkin(3) * t232 + t279;
t67 = pkin(4) * t149 - qJ(5) * t148 - qJD(5) * t297 + t204;
t23 = t288 * t67 + t290 * t51;
t353 = t148 * t288;
t14 = -pkin(9) * t353 + t23;
t22 = -t288 * t51 + t290 * t67;
t352 = t148 * t290;
t6 = pkin(5) * t149 - pkin(9) * t352 + t22;
t2 = qJD(6) * t12 + t14 * t294 + t291 * t6;
t3 = -qJD(6) * t13 - t14 * t291 + t294 * t6;
t56 = qJD(4) * t107 + t292 * t160 - t389 * t302;
t32 = pkin(5) * t353 + t56;
t329 = t352 / 0.2e1;
t330 = -t353 / 0.2e1;
t359 = t23 * t290;
t374 = Ifges(6,4) * t290;
t375 = Ifges(6,4) * t288;
t390 = t290 / 0.2e1;
t396 = t189 / 0.2e1;
t399 = t149 / 0.2e1;
t402 = t117 / 0.2e1;
t404 = t445 / 0.2e1;
t315 = Ifges(6,1) * t290 - t375;
t440 = t315 / 0.2e1;
t410 = Ifges(6,5) * t399 + t148 * t440;
t42 = -qJD(6) * t117 - t148 * t259;
t411 = t42 / 0.2e1;
t41 = qJD(6) * t445 - t148 * t309;
t412 = t41 / 0.2e1;
t413 = Ifges(7,1) * t412 + Ifges(7,4) * t411 + Ifges(7,5) * t399;
t414 = Ifges(7,4) * t412 + Ifges(7,2) * t411 + Ifges(7,6) * t399;
t314 = -Ifges(6,2) * t288 + t374;
t52 = t149 * Ifges(6,6) + t148 * t314;
t106 = t168 * t389 - t292 * t170;
t100 = -t287 * pkin(4) + qJD(5) - t106;
t88 = -pkin(5) * t317 + t100;
t458 = mrSges(6,3) * t359 + t32 * (mrSges(7,1) * t309 + mrSges(7,2) * t259) + (Ifges(7,4) * t259 - Ifges(7,2) * t309) * t411 + (Ifges(7,1) * t259 - Ifges(7,4) * t309) * t412 - t309 * t414 + (Ifges(6,5) * t288 + Ifges(7,5) * t259 + Ifges(6,6) * t290 - Ifges(7,6) * t309) * t399 + (-mrSges(6,1) * t290 + mrSges(6,2) * t288 - mrSges(5,1)) * t56 + Ifges(5,5) * t148 - Ifges(5,6) * t149 - t55 * mrSges(5,2) + (Ifges(6,1) * t288 + t374) * t329 + (Ifges(6,2) * t290 + t375) * t330 + (-Ifges(7,5) * t247 - Ifges(7,6) * t248) * t396 + (-Ifges(7,1) * t247 - Ifges(7,4) * t248) * t402 + (-Ifges(7,4) * t247 - Ifges(7,2) * t248) * t404 + t52 * t390 + t288 * t410 + t259 * t413 + (-t425 * mrSges(7,1) + mrSges(7,2) * t424) * t88 + t425 * t408 + t424 * t406 + (-t12 * t424 + t425 * t13 - t2 * t309 - t3 * t259) * mrSges(7,3);
t449 = t194 * t288;
t457 = pkin(5) * t449;
t456 = pkin(9) * t449;
t453 = -Ifges(7,1) * t447 - Ifges(7,4) * t448;
t452 = -Ifges(7,4) * t447 - Ifges(7,2) * t448;
t451 = -Ifges(7,5) * t447 - Ifges(7,6) * t448;
t393 = -t194 / 0.2e1;
t188 = Ifges(5,4) * t194;
t357 = t287 * Ifges(5,5);
t433 = t297 * Ifges(5,1);
t141 = t188 + t357 + t433;
t316 = mrSges(6,1) * t288 + mrSges(6,2) * t290;
t306 = t100 * t316;
t370 = t106 * mrSges(5,3);
t96 = t173 * Ifges(6,4) + Ifges(6,2) * t317 - Ifges(6,6) * t194;
t97 = t173 * Ifges(6,1) + Ifges(6,4) * t317 - Ifges(6,5) * t194;
t446 = t306 + t97 * t390 - t288 * t96 / 0.2e1 + t141 / 0.2e1 + t203 * mrSges(5,2) - t370 + t173 * t440 + t317 * t314 / 0.2e1 + t357 / 0.2e1;
t150 = pkin(4) * t297 - qJ(5) * t194;
t356 = t287 * Ifges(5,6);
t363 = t173 * Ifges(6,5);
t364 = t317 * Ifges(6,6);
t369 = t107 * mrSges(5,3);
t435 = t13 * mrSges(7,2);
t436 = t12 * mrSges(7,1);
t444 = t356 / 0.2e1 - t363 / 0.2e1 - t364 / 0.2e1 - t203 * mrSges(5,1) - t68 * mrSges(6,1) + t69 * mrSges(6,2) + t369 + t435 - t436;
t443 = t246 / 0.2e1;
t442 = -t297 / 0.2e1;
t441 = -t303 / 0.2e1;
t437 = pkin(5) * t297;
t362 = t189 * Ifges(7,3);
t367 = t117 * Ifges(7,5);
t368 = t445 * Ifges(7,6);
t58 = t362 + t367 + t368;
t95 = -Ifges(6,3) * t194 + t363 + t364;
t434 = t95 + t58;
t376 = Ifges(5,4) * t297;
t331 = t355 * pkin(2);
t280 = t331 + pkin(3);
t388 = pkin(2) * t289;
t240 = t292 * t280 + t389 * t388;
t236 = qJ(5) + t240;
t207 = (-pkin(9) - t236) * t288;
t286 = t290 * pkin(9);
t344 = t236 * t290;
t208 = t286 + t344;
t162 = t207 * t294 - t208 * t291;
t276 = t292 * t388;
t227 = -qJD(4) * t276 + t280 * t328;
t219 = qJD(5) + t227;
t347 = t194 * t290;
t201 = -t289 * t263 + t321;
t177 = t201 - t301;
t202 = t355 * t263 + t249;
t178 = t202 + t387;
t120 = t292 * t177 + t178 * t389;
t341 = qJD(1) * t293;
t284 = pkin(2) * t341;
t213 = -pkin(3) * t244 + t284;
t121 = t150 + t213;
t73 = -t120 * t288 + t290 * t121;
t36 = -pkin(9) * t347 + t437 + t73;
t74 = t290 * t120 + t288 * t121;
t57 = t74 - t456;
t432 = qJD(6) * t162 - t219 * t309 - t291 * t36 - t294 * t57;
t163 = t207 * t291 + t208 * t294;
t431 = -qJD(6) * t163 - t219 * t259 + t291 * t57 - t294 * t36;
t267 = (-pkin(9) - qJ(5)) * t288;
t354 = qJ(5) * t290;
t268 = t286 + t354;
t205 = t267 * t294 - t268 * t291;
t81 = -t106 * t288 + t290 * t150;
t44 = -t194 * t286 + t437 + t81;
t82 = t290 * t106 + t288 * t150;
t62 = t82 - t456;
t430 = -qJD(5) * t309 + qJD(6) * t205 - t291 * t44 - t294 * t62;
t206 = t267 * t291 + t268 * t294;
t429 = -qJD(5) * t259 - qJD(6) * t206 + t291 * t62 - t294 * t44;
t426 = -t120 + t227;
t423 = mrSges(5,1) * t287 + mrSges(6,1) * t317 - mrSges(6,2) * t173 - mrSges(5,3) * t297;
t209 = t355 * t271 + t272 * t289;
t186 = -pkin(8) * t258 + t209;
t210 = t289 * t271 - t355 * t272;
t187 = pkin(8) * t256 + t210;
t422 = t389 * t186 - t292 * t187;
t421 = (t295 * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t341) + t293 * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t340)) * pkin(7);
t420 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t419 = -m(5) * t106 + m(6) * t100 - t423;
t72 = -mrSges(7,1) * t445 + mrSges(7,2) * t117;
t418 = -m(7) * t88 - t419 - t72;
t405 = -t445 / 0.2e1;
t403 = -t117 / 0.2e1;
t140 = Ifges(5,2) * t194 + t356 + t376;
t401 = t140 / 0.2e1;
t397 = -t189 / 0.2e1;
t394 = t194 / 0.2e1;
t391 = -t244 / 0.2e1;
t386 = t290 * pkin(5);
t184 = -t241 * t289 + t355 * t242;
t165 = -pkin(8) * t246 + t184;
t185 = t355 * t241 + t289 * t242;
t166 = -pkin(8) * t245 + t185;
t75 = qJD(4) * t422 + t292 * t165 + t389 * t166;
t307 = t256 * t389 - t292 * t258;
t155 = qJD(4) * t307 - t292 * t245 + t246 * t389;
t200 = t292 * t256 + t258 * t389;
t156 = qJD(4) * t200 + t245 * t389 + t292 * t246;
t339 = qJD(2) * t293;
t285 = pkin(2) * t339;
t214 = pkin(3) * t245 + t285;
t79 = pkin(4) * t156 - qJ(5) * t155 - qJD(5) * t200 + t214;
t28 = t288 * t79 + t290 * t75;
t382 = mrSges(4,3) * t232;
t381 = mrSges(4,3) * t233;
t380 = mrSges(5,3) * t149;
t379 = Ifges(3,4) * t293;
t378 = Ifges(3,4) * t295;
t377 = Ifges(4,4) * t244;
t372 = Ifges(6,5) * t290;
t371 = Ifges(6,6) * t288;
t366 = t422 * t56;
t365 = t148 * mrSges(5,3);
t361 = t200 * t56;
t360 = t22 * t288;
t358 = t244 * mrSges(4,3);
t351 = t155 * t288;
t350 = t155 * t290;
t346 = t200 * t288;
t345 = t200 * t290;
t282 = -t295 * pkin(2) - pkin(1);
t222 = -t256 * pkin(3) + t282;
t135 = -pkin(4) * t307 - t200 * qJ(5) + t222;
t139 = t292 * t186 + t187 * t389;
t84 = t288 * t135 + t290 * t139;
t85 = mrSges(6,1) * t353 + mrSges(6,2) * t352;
t342 = qJD(1) * t282;
t334 = Ifges(7,5) * t41 + Ifges(7,6) * t42 + Ifges(7,3) * t149;
t15 = -t42 * mrSges(7,1) + t41 * mrSges(7,2);
t326 = t295 * t337;
t27 = -t288 * t75 + t290 * t79;
t319 = t232 * mrSges(4,1) + t233 * mrSges(4,2);
t318 = t149 * mrSges(5,1) + t148 * mrSges(5,2);
t83 = t290 * t135 - t139 * t288;
t119 = -t389 * t177 + t178 * t292;
t239 = t280 * t389 - t276;
t313 = Ifges(3,5) * t295 - Ifges(3,6) * t293;
t312 = -t371 + t372;
t311 = t359 - t360;
t310 = -t288 * t68 + t290 * t69;
t61 = -pkin(5) * t307 - pkin(9) * t345 + t83;
t70 = -pkin(9) * t346 + t84;
t25 = -t291 * t70 + t294 * t61;
t26 = t291 * t61 + t294 * t70;
t237 = -pkin(4) - t239;
t308 = pkin(1) * (mrSges(3,1) * t293 + mrSges(3,2) * t295);
t305 = t293 * (Ifges(3,1) * t295 - t379);
t304 = (Ifges(3,2) * t295 + t379) * qJD(1);
t300 = mrSges(4,3) * t303;
t76 = qJD(4) * t139 - t389 * t165 + t292 * t166;
t283 = Ifges(3,4) * t340;
t281 = -pkin(4) - t386;
t266 = qJD(3) + t342;
t254 = Ifges(3,1) * t341 + Ifges(3,5) * qJD(2) + t283;
t253 = Ifges(3,6) * qJD(2) + t304;
t238 = Ifges(4,4) * t303;
t218 = qJD(2) * mrSges(4,1) + t358;
t217 = -qJD(2) * mrSges(4,2) + t300;
t215 = t237 - t386;
t196 = -mrSges(4,1) * t303 - t244 * mrSges(4,2);
t191 = -t244 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t238;
t190 = Ifges(4,2) * t303 + Ifges(4,6) * qJD(2) - t377;
t179 = -mrSges(5,2) * t287 + mrSges(5,3) * t194;
t153 = t309 * t200;
t152 = t259 * t200;
t151 = -mrSges(5,1) * t194 + mrSges(5,2) * t297;
t129 = -mrSges(6,1) * t194 - mrSges(6,3) * t173;
t128 = mrSges(6,2) * t194 + mrSges(6,3) * t317;
t98 = pkin(5) * t346 - t422;
t92 = t119 + t457;
t91 = mrSges(7,1) * t189 - mrSges(7,3) * t117;
t90 = -mrSges(7,2) * t189 + mrSges(7,3) * t445;
t89 = t107 + t457;
t87 = mrSges(6,1) * t149 - mrSges(6,3) * t352;
t86 = -mrSges(6,2) * t149 - mrSges(6,3) * t353;
t66 = -t155 * t259 + t200 * t247;
t65 = -t155 * t309 - t200 * t248;
t40 = pkin(5) * t351 + t76;
t31 = -mrSges(7,2) * t149 + mrSges(7,3) * t42;
t30 = mrSges(7,1) * t149 - mrSges(7,3) * t41;
t24 = -pkin(9) * t351 + t28;
t16 = pkin(5) * t156 - pkin(9) * t350 + t27;
t5 = -qJD(6) * t26 + t16 * t294 - t24 * t291;
t4 = qJD(6) * t25 + t16 * t291 + t24 * t294;
t1 = [(t305 - 0.2e1 * t308) * t337 - t421 * qJD(2) - (t365 + t85) * t422 + (Ifges(7,4) * t65 + Ifges(7,2) * t66) * t404 + (t233 * t258 + t246 * t391) * Ifges(4,1) - (t304 + t253) * t339 / 0.2e1 + t191 * t443 + (-Ifges(7,4) * t153 - Ifges(7,2) * t152) * t411 + (-t12 * t65 + t13 * t66 - t152 * t2 + t153 * t3) * mrSges(7,3) + (-Ifges(7,1) * t153 - Ifges(7,4) * t152) * t412 + (-Ifges(7,5) * t153 - Ifges(7,6) * t152) * t399 + t32 * (mrSges(7,1) * t152 - mrSges(7,2) * t153) + m(7) * (t12 * t5 + t13 * t4 + t2 * t26 + t25 * t3 + t32 * t98 + t40 * t88) + t316 * t361 + (-t256 * t232 + t245 * t441) * Ifges(4,2) + (-t258 * t232 + t256 * t233 - t245 * t391 + t303 * t443) * Ifges(4,4) + m(5) * (t107 * t75 + t139 * t55 + t203 * t214 + t204 * t222 - t366) + m(6) * (t22 * t83 + t23 * t84 + t27 * t68 + t28 * t69 - t366) + (Ifges(7,1) * t65 + Ifges(7,4) * t66) * t402 + (t204 * mrSges(5,2) + Ifges(5,1) * t148 + t312 * t399 + t314 * t330 + t315 * t329) * t200 + t214 * t151 + t185 * t217 + t184 * t218 + (Ifges(7,5) * t65 + Ifges(7,6) * t66) * t396 + (t433 / 0.2e1 + t312 * t393 + t446) * t155 + (-mrSges(4,1) * t256 + mrSges(4,2) * t258) * t279 + ((t254 + qJD(1) * (Ifges(3,1) * t293 + t378)) * t295 + t313 * qJD(2) + Ifges(4,5) * t246 - Ifges(4,6) * t245) * qJD(2) / 0.2e1 - (Ifges(6,3) * t149 + t148 * t312 + t334) * t307 / 0.2e1 + (t148 * t307 - t149 * t200 + t155 * t394) * Ifges(5,4) - (t204 * mrSges(5,1) + t22 * mrSges(6,1) - t23 * mrSges(6,2) + Ifges(6,5) * t329 + Ifges(7,5) * t412 + Ifges(5,2) * t149 + Ifges(6,6) * t330 + Ifges(7,6) * t411 + (Ifges(6,3) + Ifges(7,3)) * t399 + t420) * t307 + (t307 * t55 + t361) * mrSges(5,3) + (-t22 * t345 - t23 * t346 - t350 * t68 - t351 * t69) * mrSges(6,3) + t75 * t179 + t28 * t128 + t27 * t129 + t98 * t15 + t4 * t90 + t5 * t91 + t84 * t86 + t83 * t87 + t88 * (-mrSges(7,1) * t66 + mrSges(7,2) * t65) + t40 * t72 + t419 * t76 + t25 * t30 + t26 * t31 - t139 * t380 - t209 * t381 + (-Ifges(3,2) * t293 + t378) * t326 - t210 * t382 + (-t175 * t258 + t176 * t256 - t197 * t246 - t198 * t245) * mrSges(4,3) + t266 * (mrSges(4,1) * t245 + mrSges(4,2) * t246) - t245 * t190 / 0.2e1 + t196 * t285 + t65 * t406 + t66 * t408 + t345 * t410 - t153 * t413 - t152 * t414 + t222 * t318 + t282 * t319 + (t442 * Ifges(5,4) + t434 / 0.2e1 - Ifges(5,2) * t394 - t140 / 0.2e1 + Ifges(6,3) * t393 + Ifges(7,3) * t396 + Ifges(7,5) * t402 + Ifges(7,6) * t404 - t444) * t156 + m(4) * (t175 * t209 + t176 * t210 + t184 * t197 + t185 * t198 + (t266 + t342) * t285) - t52 * t346 / 0.2e1; (Ifges(7,3) * t297 + t451) * t397 + (Ifges(7,6) * t297 + t452) * t405 + (Ifges(7,5) * t297 + t453) * t403 - t69 * (-mrSges(6,2) * t297 - mrSges(6,3) * t449) + t96 * t449 / 0.2e1 - t317 * (Ifges(6,6) * t297 + t194 * t314) / 0.2e1 + t297 * t369 - t173 * (Ifges(6,5) * t297 + t194 * t315) / 0.2e1 - t203 * (mrSges(5,1) * t297 + mrSges(5,2) * t194) - t287 * (Ifges(5,5) * t194 - Ifges(5,6) * t297) / 0.2e1 + (Ifges(6,3) * t297 + t194 * t312) * t394 + t297 * t401 - t68 * (mrSges(6,1) * t297 - mrSges(6,3) * t347) + (-Ifges(5,2) * t297 + t141 + t188) * t393 + (Ifges(5,1) * t194 - t376 + t434) * t442 - (-Ifges(3,2) * t341 + t254 + t283) * t340 / 0.2e1 + t194 * t370 + t197 * t300 + t297 * t435 + (Ifges(4,2) * t244 + t191 + t238) * t441 + (-t219 * t288 - t73) * t129 + (-t100 * t119 + t219 * t310 + t236 * t311 + t237 * t56 - t68 * t73 - t69 * t74) * m(6) - t194 * t306 - t266 * (-t244 * mrSges(4,1) + mrSges(4,2) * t303) - qJD(2) * (Ifges(4,5) * t303 + Ifges(4,6) * t244) / 0.2e1 + (-mrSges(3,1) * t326 + mrSges(3,2) * t327) * pkin(7) + (t219 * t290 - t74) * t128 + t86 * t344 - t213 * t151 + t215 * t15 - t202 * t217 - t201 * t218 + t458 + (t106 * t119 + t107 * t426 - t203 * t213 - t239 * t56) * m(5) + (m(5) * t55 - qJD(4) * t418 - t380) * t240 + t175 * mrSges(4,1) - t176 * mrSges(4,2) + t162 * t30 + t163 * t31 - t297 * t436 + t431 * t91 + t432 * t90 + (t12 * t431 + t13 * t432 + t162 * t3 + t163 * t2 + t215 * t32 - t88 * t92) * m(7) - t92 * t72 + t426 * t179 + t423 * t119 - t288 * t236 * t87 + ((t175 * t355 + t176 * t289) * pkin(2) - t197 * t201 - t198 * t202 - t266 * t284) * m(4) - t331 * t381 + t244 * (Ifges(4,1) * t303 + t377) / 0.2e1 - t239 * t365 - mrSges(6,3) * t360 - t198 * t358 - Ifges(4,6) * t232 + Ifges(4,5) * t233 + t237 * t85 + Ifges(3,5) * t326 - t382 * t388 + t190 * t391 + (t421 + (-t305 / 0.2e1 + t308) * qJD(1)) * qJD(1) - Ifges(3,6) * t327 - t196 * t284 - t313 * t337 / 0.2e1 + t253 * t341 / 0.2e1 - t97 * t347 / 0.2e1; -t217 * t303 + t129 * t449 - t244 * t218 - t309 * t30 + t259 * t31 + t288 * t86 + t290 * t87 + t318 + t319 + m(6) * (t22 * t290 + t23 * t288) + m(5) * t204 - t128 * t347 + t425 * t91 + t424 * t90 + (t12 * t425 + t13 * t424 + t2 * t259 - t3 * t309) * m(7) + (-t197 * t244 - t198 * t303 + t279) * m(4) + (-m(5) * t107 - m(6) * t310 - t179) * t194 + t418 * t297; t451 * t397 + t452 * t405 + t453 * t403 + (-pkin(4) * t56 + qJ(5) * t311 + qJD(5) * t310 - t100 * t107 - t68 * t81 - t69 * t82) * m(6) + (qJD(5) * t290 - t82) * t128 + t205 * t30 + t206 * t31 + t86 * t354 + (-t22 * mrSges(6,3) - qJ(5) * t87 - qJD(5) * t129) * t288 + t458 - (t188 / 0.2e1 - (t372 / 0.2e1 - t371 / 0.2e1) * t194 + (-t288 * t69 - t290 * t68) * mrSges(6,3) + (Ifges(5,1) / 0.2e1 - Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1) * t297 + t446) * t194 - t106 * t179 - t81 * t129 + t429 * t91 + t430 * t90 + (t12 * t429 + t13 * t430 + t2 * t206 + t205 * t3 + t281 * t32 - t88 * t89) * m(7) - t89 * t72 + t423 * t107 - pkin(4) * t85 + t281 * t15 + (t401 - t362 / 0.2e1 - t368 / 0.2e1 - t367 / 0.2e1 - t95 / 0.2e1 - t58 / 0.2e1 + t376 / 0.2e1 + t444) * t297; t117 * t91 - t445 * t90 - t317 * t128 + t173 * t129 + t15 + t85 + (t117 * t12 - t13 * t445 + t32) * m(7) + (t173 * t68 - t317 * t69 + t56) * m(6); -t88 * (mrSges(7,1) * t117 + mrSges(7,2) * t445) + (Ifges(7,1) * t445 - t373) * t403 + t59 * t402 + (Ifges(7,5) * t445 - Ifges(7,6) * t117) * t397 - t12 * t90 + t13 * t91 + (t117 * t13 + t12 * t445) * mrSges(7,3) + t334 + (-Ifges(7,2) * t117 + t111 + t60) * t405 + t420;];
tauc  = t1(:);
