% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:22:53
% EndTime: 2019-03-09 09:23:17
% DurationCPUTime: 12.99s
% Computational Cost: add. (9109->669), mult. (22470->925), div. (0->0), fcn. (15521->8), ass. (0->302)
t263 = sin(pkin(10));
t264 = cos(pkin(10));
t267 = sin(qJ(5));
t270 = cos(qJ(5));
t219 = t263 * t270 - t264 * t267;
t271 = cos(qJ(2));
t282 = t271 * t219;
t182 = qJD(1) * t282;
t200 = t219 * qJD(5);
t333 = t182 + t200;
t286 = t263 * t267 + t264 * t270;
t281 = t271 * t286;
t183 = qJD(1) * t281;
t199 = t286 * qJD(5);
t332 = t183 + t199;
t268 = sin(qJ(2));
t311 = -pkin(7) * t263 - pkin(3);
t336 = t264 * t271;
t316 = pkin(8) * t336;
t277 = -t316 + (-pkin(4) + t311) * t268;
t292 = pkin(2) * t268 - qJ(3) * t271;
t224 = t292 * qJD(1);
t338 = t264 * t224;
t111 = qJD(1) * t277 - t338;
t201 = t263 * t224;
t328 = qJD(1) * t268;
t251 = qJ(4) * t328;
t337 = t264 * t268;
t340 = t263 * t271;
t283 = -pkin(7) * t337 + pkin(8) * t340;
t129 = qJD(1) * t283 + t201 + t251;
t366 = -pkin(8) + qJ(3);
t230 = t366 * t263;
t231 = t366 * t264;
t167 = t267 * t230 + t270 * t231;
t324 = qJD(3) * t270;
t325 = qJD(3) * t267;
t430 = -qJD(5) * t167 - t270 * t111 + t129 * t267 + t263 * t324 - t264 * t325;
t320 = qJD(5) * t270;
t321 = qJD(5) * t267;
t429 = -t267 * t111 - t270 * t129 + t230 * t320 - t231 * t321 + t263 * t325 + t264 * t324;
t440 = Ifges(4,1) + Ifges(5,1);
t442 = pkin(5) * t328 + pkin(9) * t332 + t430;
t441 = pkin(9) * t333 - t429;
t259 = t271 * qJD(1);
t250 = t259 + qJD(5);
t241 = qJD(6) + t250;
t433 = Ifges(5,4) + Ifges(4,5);
t307 = t263 * t259;
t255 = pkin(7) * t259;
t309 = qJ(4) * t259;
t330 = t264 * t309 - t255;
t389 = pkin(3) + pkin(4);
t151 = -t307 * t389 + t330;
t439 = qJD(4) * t263 - t151;
t266 = sin(qJ(6));
t308 = t263 * t328;
t318 = t264 * qJD(2);
t213 = t308 - t318;
t306 = t264 * t328;
t319 = t263 * qJD(2);
t214 = t306 + t319;
t287 = t213 * t267 + t214 * t270;
t228 = -pkin(2) * t271 - t268 * qJ(3) - pkin(1);
t207 = t228 * qJD(1);
t237 = qJD(2) * qJ(3) + t255;
t149 = t264 * t207 - t263 * t237;
t128 = pkin(3) * t259 + qJD(4) - t149;
t88 = pkin(4) * t259 - t214 * pkin(8) + t128;
t150 = t263 * t207 + t264 * t237;
t130 = -t309 + t150;
t94 = t213 * pkin(8) + t130;
t45 = -t267 * t94 + t270 * t88;
t37 = -pkin(9) * t287 + t45;
t33 = pkin(5) * t250 + t37;
t269 = cos(qJ(6));
t143 = t213 * t270 - t214 * t267;
t46 = t267 * t88 + t270 * t94;
t38 = pkin(9) * t143 + t46;
t348 = t269 * t38;
t10 = t266 * t33 + t348;
t302 = t269 * t143 - t287 * t266;
t425 = t143 * t266 + t269 * t287;
t65 = Ifges(7,4) * t302;
t36 = Ifges(7,1) * t425 + Ifges(7,5) * t241 + t65;
t373 = Ifges(7,4) * t425;
t380 = -t241 / 0.2e1;
t393 = -t425 / 0.2e1;
t395 = -t302 / 0.2e1;
t254 = pkin(7) * t328;
t227 = -qJD(2) * pkin(2) + qJD(3) + t254;
t121 = t213 * pkin(3) - t214 * qJ(4) + t227;
t91 = -pkin(4) * t213 - t121;
t54 = -pkin(5) * t143 + t91;
t349 = t266 * t38;
t9 = t269 * t33 - t349;
t437 = (Ifges(7,5) * t302 - Ifges(7,6) * t425) * t380 + (t10 * t425 + t302 * t9) * mrSges(7,3) + (-Ifges(7,2) * t425 + t36 + t65) * t395 - t54 * (mrSges(7,1) * t425 + mrSges(7,2) * t302) + (Ifges(7,1) * t302 - t373) * t393;
t436 = -qJD(2) / 0.2e1;
t435 = qJD(2) / 0.2e1;
t434 = -Ifges(4,4) + Ifges(5,5);
t166 = t270 * t230 - t231 * t267;
t116 = -pkin(9) * t219 + t166;
t117 = -pkin(9) * t286 + t167;
t53 = t116 * t266 + t117 * t269;
t432 = -qJD(6) * t53 + t441 * t266 + t442 * t269;
t52 = t116 * t269 - t117 * t266;
t431 = qJD(6) * t52 + t442 * t266 - t441 * t269;
t428 = t333 * pkin(5) + t439;
t357 = Ifges(5,5) * t263;
t362 = Ifges(4,4) * t263;
t427 = t440 * t264 + t357 - t362;
t198 = qJD(2) * t292 - t268 * qJD(3);
t184 = t198 * qJD(1);
t225 = (qJD(3) - t254) * qJD(2);
t125 = t264 * t184 - t263 * t225;
t317 = qJD(1) * qJD(2);
t81 = (-t268 * t389 - t316) * t317 - t125;
t126 = t263 * t184 + t264 * t225;
t304 = t268 * t317;
t312 = qJ(4) * t304 + t126;
t82 = (pkin(8) * t319 - qJD(4)) * t259 + t312;
t16 = -qJD(5) * t46 - t267 * t82 + t270 * t81;
t278 = qJD(2) * t281;
t89 = qJD(1) * t278 + qJD(5) * t143;
t11 = -pkin(5) * t304 - pkin(9) * t89 + t16;
t15 = t267 * t81 + t270 * t82 + t88 * t320 - t321 * t94;
t279 = qJD(2) * t282;
t90 = qJD(1) * t279 - qJD(5) * t287;
t14 = pkin(9) * t90 + t15;
t2 = qJD(6) * t9 + t11 * t266 + t14 * t269;
t29 = qJD(6) * t302 + t266 * t90 + t269 * t89;
t3 = -qJD(6) * t10 + t11 * t269 - t14 * t266;
t30 = -qJD(6) * t425 - t266 * t89 + t269 * t90;
t426 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t29 + Ifges(7,6) * t30;
t399 = t29 / 0.2e1;
t398 = t30 / 0.2e1;
t35 = Ifges(7,2) * t302 + Ifges(7,6) * t241 + t373;
t421 = t35 / 0.2e1;
t391 = t89 / 0.2e1;
t390 = t90 / 0.2e1;
t420 = t214 / 0.2e1;
t418 = -t304 / 0.2e1;
t417 = t317 / 0.2e1;
t305 = Ifges(3,6) * t436;
t416 = Ifges(3,5) * t435;
t248 = pkin(7) * t340;
t262 = t271 * pkin(3);
t372 = pkin(8) * t268;
t132 = pkin(4) * t271 + t248 + t262 + (-t228 - t372) * t264;
t177 = pkin(7) * t336 + t263 * t228;
t343 = qJ(4) * t271;
t169 = t177 - t343;
t148 = t263 * t372 + t169;
t64 = t267 * t132 + t270 * t148;
t285 = t266 * t267 - t269 * t270;
t410 = t241 * t285;
t223 = t266 * t270 + t267 * t269;
t409 = t241 * t223;
t327 = qJD(2) * t268;
t408 = qJ(4) * t327 - qJD(4) * t271;
t405 = t263 * t389 + pkin(7);
t404 = t16 * mrSges(6,1) - t15 * mrSges(6,2) + Ifges(6,5) * t89 + Ifges(6,6) * t90 + t426;
t363 = Ifges(3,4) * t268;
t403 = 0.2e1 * (Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * t213 + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t214 + t10 * mrSges(7,2) + t130 * mrSges(5,3) + t149 * mrSges(4,1) + t46 * mrSges(6,2) + t305 - (t271 * Ifges(3,2) + t363) * qJD(1) / 0.2e1 - t241 * Ifges(7,3) - t425 * Ifges(7,5) - t302 * Ifges(7,6) - t250 * Ifges(6,3) - t287 * Ifges(6,5) - t143 * Ifges(6,6) - t128 * mrSges(5,1) - t150 * mrSges(4,2) - t45 * mrSges(6,1) - t9 * mrSges(7,1) + t433 * t420 - (Ifges(4,3) + Ifges(5,2)) * t259 / 0.2e1;
t253 = Ifges(3,4) * t259;
t356 = Ifges(5,5) * t264;
t293 = Ifges(5,3) * t263 + t356;
t361 = Ifges(4,4) * t264;
t294 = -Ifges(4,2) * t263 + t361;
t364 = mrSges(5,3) * t264;
t374 = t264 / 0.2e1;
t375 = t263 / 0.2e1;
t376 = -t263 / 0.2e1;
t402 = (t128 * t264 - t130 * t263) * mrSges(5,2) - (t149 * t264 + t150 * t263) * mrSges(4,3) + (m(4) * t227 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t213 + mrSges(4,2) * t214 + mrSges(3,3) * t328) * pkin(7) + t121 * (mrSges(5,1) * t263 - t364) + t227 * (mrSges(4,1) * t263 + mrSges(4,2) * t264) + Ifges(3,1) * t328 / 0.2e1 + t253 / 0.2e1 + t416 + (Ifges(5,5) * t214 - Ifges(5,6) * t259) * t375 + (Ifges(4,4) * t214 - Ifges(4,6) * t259) * t376 + t427 * t420 + (t214 * t440 - t433 * t259) * t374 + (Ifges(5,3) * t375 - Ifges(4,2) * t376 + t434 * t374 - t294 / 0.2e1 + t293 / 0.2e1) * t213;
t401 = Ifges(7,4) * t399 + Ifges(7,2) * t398 + Ifges(7,6) * t418;
t400 = Ifges(7,1) * t399 + Ifges(7,4) * t398 + Ifges(7,5) * t418;
t397 = Ifges(6,4) * t391 + Ifges(6,2) * t390 + Ifges(6,6) * t418;
t396 = Ifges(6,1) * t391 + Ifges(6,4) * t390 + Ifges(6,5) * t418;
t394 = t302 / 0.2e1;
t392 = t425 / 0.2e1;
t388 = pkin(1) * mrSges(3,1);
t387 = pkin(1) * mrSges(3,2);
t385 = -t143 / 0.2e1;
t384 = t143 / 0.2e1;
t383 = -t287 / 0.2e1;
t382 = t287 / 0.2e1;
t379 = t241 / 0.2e1;
t378 = -t250 / 0.2e1;
t377 = t250 / 0.2e1;
t367 = -Ifges(6,3) - Ifges(7,3);
t147 = t219 * t269 - t266 * t286;
t73 = -qJD(6) * t147 + t199 * t266 - t200 * t269;
t99 = t182 * t269 - t183 * t266;
t365 = t73 - t99;
t360 = Ifges(5,4) * t264;
t359 = Ifges(6,4) * t287;
t358 = Ifges(4,5) * t264;
t355 = Ifges(4,6) * t263;
t354 = Ifges(5,6) * t263;
t100 = t182 * t266 + t183 * t269;
t146 = -t219 * t266 - t269 * t286;
t72 = qJD(6) * t146 - t199 * t269 - t200 * t266;
t347 = t100 - t72;
t154 = mrSges(5,1) * t213 - mrSges(5,3) * t214;
t74 = -mrSges(6,1) * t143 + mrSges(6,2) * t287;
t346 = -t154 + t74;
t342 = qJD(2) * mrSges(3,2);
t339 = t264 * t198;
t178 = -t213 * mrSges(5,2) - mrSges(5,3) * t259;
t179 = mrSges(4,2) * t259 - t213 * mrSges(4,3);
t335 = t179 + t178;
t180 = -mrSges(4,1) * t259 - t214 * mrSges(4,3);
t181 = mrSges(5,1) * t259 + t214 * mrSges(5,2);
t334 = t180 - t181;
t303 = t271 * t317;
t300 = t264 * t303;
t331 = -qJ(4) * t300 - t214 * qJD(4);
t301 = t263 * t303;
t186 = mrSges(4,1) * t301 + mrSges(4,2) * t300;
t329 = -qJD(4) * t337 - t318 * t343;
t326 = qJD(2) * t271;
t226 = -t264 * pkin(3) - t263 * qJ(4) - pkin(2);
t315 = pkin(7) * t327;
t310 = pkin(3) * t263 + pkin(7);
t63 = t270 * t132 - t267 * t148;
t176 = t264 * t228 - t248;
t202 = t264 * pkin(4) - t226;
t297 = t311 * t268;
t44 = -t90 * mrSges(6,1) + t89 * mrSges(6,2);
t8 = -t30 * mrSges(7,1) + t29 * mrSges(7,2);
t193 = t286 * t268;
t50 = pkin(5) * t271 - t193 * pkin(9) + t63;
t192 = t219 * t268;
t51 = pkin(9) * t192 + t64;
t25 = -t266 * t51 + t269 * t50;
t26 = t266 * t50 + t269 * t51;
t109 = -mrSges(6,2) * t250 + mrSges(6,3) * t143;
t110 = mrSges(6,1) * t250 - mrSges(6,3) * t287;
t290 = t109 * t270 - t110 * t267;
t112 = t192 * t269 - t193 * t266;
t113 = t192 * t266 + t193 * t269;
t172 = -pkin(7) * t306 + t201;
t189 = t263 * t198;
t161 = -t264 * t315 + t189;
t284 = t310 * t326;
t101 = qJD(2) * t277 - t339;
t102 = qJD(2) * t283 + t189 + t408;
t31 = t267 * t101 + t270 * t102 + t132 * t320 - t148 * t321;
t246 = qJ(4) * t337;
t168 = -t268 * t405 + t246;
t127 = -t326 * t405 - t329;
t114 = qJD(1) * t284 + t331;
t93 = -qJD(4) * t259 + t312;
t276 = t114 * mrSges(5,1) + (Ifges(5,6) * t268 + t271 * t293) * t417 - (Ifges(4,6) * t268 + t271 * t294) * t317 / 0.2e1 - t93 * mrSges(5,2) - t126 * mrSges(4,3);
t105 = -pkin(3) * t304 - t125;
t275 = t105 * mrSges(5,2) - t125 * mrSges(4,3) - t114 * mrSges(5,3) + (t268 * t433 + t271 * t427) * t417;
t95 = t303 * t405 + t331;
t32 = -qJD(5) * t64 + t270 * t101 - t102 * t267;
t274 = t358 / 0.2e1 - t355 / 0.2e1 + t360 / 0.2e1 + t354 / 0.2e1;
t238 = mrSges(3,3) * t259 - t342;
t234 = mrSges(5,2) * t300;
t232 = mrSges(5,1) * t301;
t197 = (-mrSges(5,2) * t340 + mrSges(5,3) * t268) * t317;
t196 = -mrSges(5,1) * t304 + t234;
t195 = (mrSges(4,1) * t268 - mrSges(4,3) * t336) * t317;
t194 = (-mrSges(4,2) * t268 - mrSges(4,3) * t340) * t317;
t190 = t268 * t310 - t246;
t185 = -mrSges(5,3) * t300 + t232;
t174 = pkin(3) * t307 - t330;
t171 = pkin(7) * t308 + t338;
t170 = -t176 + t262;
t160 = t263 * t315 + t339;
t157 = pkin(5) * t286 + t202;
t156 = qJD(1) * t297 - t338;
t153 = t172 + t251;
t152 = t284 + t329;
t141 = qJD(2) * t297 - t339;
t140 = Ifges(6,4) * t143;
t124 = t161 + t408;
t120 = qJD(5) * t192 + t278;
t119 = -t199 * t268 + t279;
t106 = -pkin(5) * t192 + t168;
t76 = mrSges(6,2) * t304 + mrSges(6,3) * t90;
t75 = -mrSges(6,1) * t304 - mrSges(6,3) * t89;
t62 = Ifges(6,1) * t287 + t250 * Ifges(6,5) + t140;
t61 = t143 * Ifges(6,2) + t250 * Ifges(6,6) + t359;
t59 = -t119 * pkin(5) + t127;
t58 = mrSges(7,1) * t241 - mrSges(7,3) * t425;
t57 = -mrSges(7,2) * t241 + mrSges(7,3) * t302;
t48 = t90 * pkin(5) + t95;
t41 = -qJD(6) * t113 + t119 * t269 - t120 * t266;
t40 = qJD(6) * t112 + t119 * t266 + t120 * t269;
t39 = -mrSges(7,1) * t302 + mrSges(7,2) * t425;
t24 = mrSges(7,2) * t304 + mrSges(7,3) * t30;
t23 = -mrSges(7,1) * t304 - mrSges(7,3) * t29;
t20 = pkin(9) * t119 + t31;
t19 = -pkin(5) * t327 - pkin(9) * t120 + t32;
t13 = t269 * t37 - t349;
t12 = -t266 * t37 - t348;
t5 = -qJD(6) * t26 + t19 * t269 - t20 * t266;
t4 = qJD(6) * t25 + t19 * t266 + t20 * t269;
t1 = [m(6) * (t127 * t91 + t15 * t64 + t16 * t63 - t168 * t95 + t31 * t46 + t32 * t45) + m(7) * (t10 * t4 - t106 * t48 + t2 * t26 + t25 * t3 + t5 * t9 + t54 * t59) + m(5) * (t105 * t170 + t114 * t190 + t121 * t152 + t124 * t130 + t128 * t141 + t169 * t93) + (t10 * t41 + t112 * t2 - t113 * t3 - t40 * t9) * mrSges(7,3) + (t119 * t46 - t120 * t45 + t15 * t192 - t16 * t193) * mrSges(6,3) + (-t93 * mrSges(5,3) + t126 * mrSges(4,2) + t105 * mrSges(5,1) - t125 * mrSges(4,1) + (t416 + (-0.2e1 * t387 + (-0.3e1 / 0.2e1 * t360 - 0.3e1 / 0.2e1 * t354 - 0.3e1 / 0.2e1 * t358 + 0.3e1 / 0.2e1 * t355 + 0.3e1 / 0.2e1 * Ifges(3,4)) * t271 + (m(4) * pkin(7) ^ 2 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(5,2) + (pkin(7) * mrSges(4,2) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t264) * t264 + (pkin(7) * mrSges(4,1) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t263 + t434 * t264) * t263 + t367) * t268) * qJD(1) + t402) * qJD(2) + t404) * t271 + m(4) * (t125 * t176 + t126 * t177 + t149 * t160 + t150 * t161) + (pkin(7) * t186 + t275 * t264 + t276 * t263 + (((-0.3e1 / 0.2e1 * Ifges(3,4) + t274) * t268 - Ifges(7,5) * t113 / 0.2e1 - Ifges(7,6) * t112 / 0.2e1 - Ifges(6,5) * t193 / 0.2e1 - Ifges(6,6) * t192 / 0.2e1 - 0.2e1 * t388) * qJD(1) - pkin(7) * t238 + t305 + t403) * qJD(2)) * t268 + t25 * t23 + t26 * t24 + (Ifges(7,4) * t40 + Ifges(7,2) * t41) * t394 + t193 * t396 + t192 * t397 + (Ifges(7,4) * t113 + Ifges(7,2) * t112) * t398 + (Ifges(7,1) * t113 + Ifges(7,4) * t112) * t399 + t113 * t400 + t112 * t401 + (Ifges(6,4) * t120 + Ifges(6,2) * t119) * t384 + (Ifges(6,4) * t193 + Ifges(6,2) * t192) * t390 + (Ifges(6,1) * t193 + Ifges(6,4) * t192) * t391 + (Ifges(7,1) * t40 + Ifges(7,4) * t41) * t392 + t40 * t36 / 0.2e1 + t54 * (-mrSges(7,1) * t41 + mrSges(7,2) * t40) + t4 * t57 + t5 * t58 + t59 * t39 + t63 * t75 + t64 * t76 + t41 * t421 + (Ifges(6,5) * t120 + Ifges(6,6) * t119) * t377 + (Ifges(7,5) * t40 + Ifges(7,6) * t41) * t379 + (Ifges(6,1) * t120 + Ifges(6,4) * t119) * t382 + t106 * t8 + t31 * t109 + t32 * t110 - t48 * (-mrSges(7,1) * t112 + mrSges(7,2) * t113) + t119 * t61 / 0.2e1 + t120 * t62 / 0.2e1 + t91 * (-mrSges(6,1) * t119 + mrSges(6,2) * t120) + t127 * t74 + t152 * t154 + t168 * t44 + t124 * t178 + t161 * t179 + t160 * t180 + t141 * t181 + t190 * t185 - t95 * (-mrSges(6,1) * t192 + mrSges(6,2) * t193) + t177 * t194 + t176 * t195 + t170 * t196 + t169 * t197; (((t388 + t363 / 0.2e1) * qJD(1) + (t238 + t342) * pkin(7) + t305 + (Ifges(6,5) * t219 + Ifges(7,5) * t147 - Ifges(6,6) * t286 + Ifges(7,6) * t146) * t436 + ((Ifges(4,6) - Ifges(5,6)) * t264 + t433 * t263) * t435 - t403) * t268 + ((-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t328 + (t271 * t274 + t387) * qJD(1) - t253 / 0.2e1 + ((Ifges(4,2) * t264 + t362) * t376 + (-Ifges(5,3) * t264 + t357) * t375 + Ifges(3,5) / 0.2e1 + (-m(4) * pkin(2) - mrSges(4,1) * t264 + mrSges(4,2) * t263 - mrSges(3,1)) * pkin(7) + (t263 * t440 - t356 + t361) * t374) * qJD(2) - t402) * t271) * qJD(1) + (t72 / 0.2e1 - t100 / 0.2e1) * t36 + (-t15 * t286 - t16 * t219 + t332 * t45 - t333 * t46) * mrSges(6,3) + (Ifges(6,4) * t219 - Ifges(6,2) * t286) * t390 + (Ifges(6,1) * t219 - Ifges(6,4) * t286) * t391 - t95 * (mrSges(6,1) * t286 + mrSges(6,2) * t219) - t286 * t397 + (-Ifges(6,1) * t199 - Ifges(6,4) * t200) * t382 + (t73 / 0.2e1 - t99 / 0.2e1) * t35 - m(4) * (t149 * t171 + t150 * t172) + t431 * t57 + t432 * t58 + (t10 * t431 - t157 * t48 + t2 * t53 + t3 * t52 + t428 * t54 + t432 * t9) * m(7) + m(4) * ((-t149 * t263 + t150 * t264) * qJD(3) + (-t125 * t263 + t126 * t264) * qJ(3)) + (-t183 / 0.2e1 - t199 / 0.2e1) * t62 + (-t182 / 0.2e1 - t200 / 0.2e1) * t61 + (-Ifges(6,4) * t199 - Ifges(6,2) * t200) * t384 + (-Ifges(6,5) * t199 - Ifges(6,6) * t200) * t377 + (t15 * t167 + t16 * t166 - t202 * t95 + t429 * t46 + t430 * t45 + t439 * t91) * m(6) + t428 * t39 + t429 * t109 + t430 * t110 + (-mrSges(7,1) * t365 - mrSges(7,2) * t347) * t54 + (t10 * t365 + t146 * t2 - t147 * t3 + t347 * t9) * mrSges(7,3) + (t346 * qJD(4) - t334 * qJD(3) + (-t195 + t196) * qJ(3) + t275) * t263 + (mrSges(6,1) * t333 - mrSges(6,2) * t332) * t91 + (t335 * qJD(3) + (t194 + t197) * qJ(3) - t276) * t264 + (t114 * t226 + (qJ(3) * t93 + qJD(3) * t130) * t264 + (qJ(3) * t105 + qJD(3) * t128 - qJD(4) * t121) * t263 - t121 * t174 - t128 * t156 - t130 * t153) * m(5) + (Ifges(7,4) * t72 + Ifges(7,2) * t73) * t394 + (Ifges(7,4) * t100 + Ifges(7,2) * t99) * t395 + t219 * t396 + (Ifges(7,4) * t147 + Ifges(7,2) * t146) * t398 + (Ifges(7,1) * t147 + Ifges(7,4) * t146) * t399 + t147 * t400 + t146 * t401 + (Ifges(6,1) * t183 + Ifges(6,4) * t182) * t383 + (Ifges(6,4) * t183 + Ifges(6,2) * t182) * t385 + (Ifges(7,1) * t72 + Ifges(7,4) * t73) * t392 + (Ifges(7,1) * t100 + Ifges(7,4) * t99) * t393 + t52 * t23 + t53 * t24 + (Ifges(6,5) * t183 + Ifges(6,6) * t182) * t378 + (Ifges(7,5) * t72 + Ifges(7,6) * t73) * t379 + (Ifges(7,5) * t100 + Ifges(7,6) * t99) * t380 - t48 * (-mrSges(7,1) * t146 + mrSges(7,2) * t147) - t151 * t74 + t157 * t8 + t166 * t75 + t167 * t76 - t174 * t154 - t153 * t178 - t172 * t179 - t171 * t180 - t156 * t181 - pkin(2) * t186 + t202 * t44 + t226 * t185; t143 * t109 - t287 * t110 + t302 * t57 - t425 * t58 + t232 + t334 * t214 + t335 * t213 + (m(4) * pkin(7) - t364) * t303 - m(4) * (-t149 * t214 - t150 * t213) - t8 - t44 + t186 + (t10 * t302 - t425 * t9 + t48) * m(7) + (t143 * t46 - t287 * t45 + t95) * m(6) + (-t128 * t214 + t130 * t213 + t114) * m(5); -t285 * t23 + t223 * t24 + t267 * t76 + t270 * t75 + t234 - t409 * t58 - t410 * t57 + t290 * qJD(5) + (-t39 - t346) * t214 + (-mrSges(5,1) * t327 + (t178 + t290) * t271) * qJD(1) + (-t10 * t410 + t2 * t223 - t214 * t54 - t285 * t3 - t409 * t9) * m(7) + (t15 * t267 + t16 * t270 - t91 * t214 - t250 * (t267 * t45 - t270 * t46)) * m(6) + (t121 * t214 + t130 * t259 + t105) * m(5); (-Ifges(6,2) * t287 + t140 + t62) * t385 + t404 + (-t287 * t39 + t269 * t23 + t266 * t24 + (-t266 * t58 + t269 * t57) * qJD(6) + (-t287 * t54 + t2 * t266 + t269 * t3 + (t10 * t269 - t266 * t9) * qJD(6)) * m(7)) * pkin(5) + (t143 * t45 + t287 * t46) * mrSges(6,3) + (Ifges(6,5) * t143 - Ifges(6,6) * t287) * t378 - t91 * (mrSges(6,1) * t287 + mrSges(6,2) * t143) + t425 * t421 + t367 * t304 + (Ifges(6,1) * t143 - t359) * t383 + t61 * t382 - m(7) * (t10 * t13 + t12 * t9) - t13 * t57 - t12 * t58 - t45 * t109 + t46 * t110 + t437; -Ifges(7,3) * t304 + t10 * t58 + t35 * t392 - t9 * t57 + t426 + t437;];
tauc  = t1(:);
