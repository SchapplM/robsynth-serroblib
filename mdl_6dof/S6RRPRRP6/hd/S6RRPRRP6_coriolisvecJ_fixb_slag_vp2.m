% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:14:19
% EndTime: 2018-11-23 17:14:47
% DurationCPUTime: 28.69s
% Computational Cost: add. (15713->789), mult. (47241->1061), div. (0->0), fcn. (37355->10), ass. (0->368)
t262 = sin(pkin(6));
t261 = sin(pkin(11));
t266 = sin(qJ(2));
t269 = cos(qJ(2));
t389 = cos(pkin(11));
t328 = t389 * t269;
t282 = -t261 * t266 + t328;
t274 = t262 * t282;
t220 = qJD(1) * t274;
t481 = t220 - qJD(4);
t518 = Ifges(6,1) + Ifges(7,1);
t516 = Ifges(6,5) + Ifges(7,4);
t263 = cos(pkin(6));
t424 = pkin(1) * t263;
t257 = t269 * t424;
t252 = qJD(1) * t257;
t414 = pkin(8) + qJ(3);
t339 = t414 * t266;
t322 = t262 * t339;
t203 = -qJD(1) * t322 + t252;
t256 = t266 * t424;
t380 = t262 * t269;
t477 = t380 * t414 + t256;
t204 = t477 * qJD(1);
t330 = t389 * t204;
t154 = t203 * t261 + t330;
t265 = sin(qJ(4));
t268 = cos(qJ(4));
t542 = -t154 - t481 * (pkin(4) * t265 - pkin(10) * t268);
t329 = t389 * t266;
t369 = qJD(1) * t262;
t348 = t269 * t369;
t221 = -t261 * t348 - t329 * t369;
t254 = qJD(1) * t263 + qJD(2);
t183 = -t221 * t268 + t254 * t265;
t223 = qJD(2) * t274;
t215 = qJD(1) * t223;
t137 = qJD(4) * t183 + t215 * t265;
t445 = t137 / 0.2e1;
t227 = (t261 * t269 + t329) * t262;
t222 = qJD(2) * t227;
t214 = qJD(1) * t222;
t264 = sin(qJ(5));
t267 = cos(qJ(5));
t363 = qJD(5) * t264;
t182 = t221 * t265 + t254 * t268;
t136 = qJD(4) * t182 + t215 * t268;
t272 = -t282 * t369 + qJD(4);
t496 = qJD(5) * t272 + t136;
t70 = -t183 * t363 + t264 * t214 + t267 * t496;
t453 = t70 / 0.2e1;
t541 = t516 * t445 + t518 * t453;
t431 = t214 / 0.2e1;
t540 = 0.2e1 * t431;
t361 = qJD(5) * t267;
t71 = t183 * t361 - t267 * t214 + t264 * t496;
t452 = -t71 / 0.2e1;
t451 = t71 / 0.2e1;
t446 = -t137 / 0.2e1;
t517 = -Ifges(6,4) + Ifges(7,5);
t539 = Ifges(6,6) - Ifges(7,6);
t514 = Ifges(6,3) + Ifges(7,2);
t447 = t136 / 0.2e1;
t538 = Ifges(5,4) * t447;
t181 = qJD(5) - t182;
t138 = t264 * t183 - t267 * t272;
t139 = t267 * t183 + t264 * t272;
t511 = -t539 * t138 + t516 * t139 + t514 * t181;
t537 = t511 / 0.2e1;
t536 = t517 * t451 + t541;
t184 = pkin(2) * t254 + t203;
t126 = t261 * t184 + t330;
t119 = pkin(9) * t254 + t126;
t242 = (-pkin(2) * t269 - pkin(1)) * t262;
t370 = qJD(1) * t242;
t237 = qJD(3) + t370;
t148 = -t220 * pkin(3) + t221 * pkin(9) + t237;
t79 = t268 * t119 + t265 * t148;
t64 = pkin(10) * t272 + t79;
t194 = t261 * t204;
t125 = t184 * t389 - t194;
t118 = -t254 * pkin(3) - t125;
t75 = -t182 * pkin(4) - t183 * pkin(10) + t118;
t23 = t264 * t75 + t267 * t64;
t18 = qJ(6) * t181 + t23;
t535 = t18 * mrSges(7,2);
t132 = Ifges(6,4) * t138;
t398 = Ifges(7,5) * t138;
t510 = t139 * t518 + t516 * t181 - t132 + t398;
t534 = Ifges(5,4) * t446;
t350 = t389 * pkin(2);
t260 = -t350 - pkin(3);
t241 = -t268 * pkin(4) - t265 * pkin(10) + t260;
t155 = t203 * t389 - t194;
t349 = t266 * t369;
t326 = pkin(2) * t349;
t167 = -pkin(3) * t221 - pkin(9) * t220 + t326;
t94 = t268 * t155 + t265 * t167;
t81 = -pkin(10) * t221 + t94;
t533 = t241 * t361 + t264 * t542 - t267 * t81;
t423 = pkin(2) * t261;
t259 = pkin(9) + t423;
t376 = t267 * t268;
t494 = t264 * t241 + t259 * t376;
t532 = -qJD(5) * t494 + t264 * t81 + t267 * t542;
t173 = t220 * t376 - t221 * t264;
t360 = t268 * qJD(4);
t362 = qJD(5) * t265;
t276 = -t264 * t362 + t267 * t360;
t531 = t173 - t276;
t22 = -t264 * t64 + t267 * t75;
t249 = qJD(2) * t252;
t273 = (-qJD(2) * t339 + qJD(3) * t269) * t262;
t177 = qJD(1) * t273 + t249;
t381 = t262 * t266;
t188 = -qJD(2) * t477 - qJD(3) * t381;
t270 = qJD(1) * t188;
t112 = t177 * t389 + t261 * t270;
t368 = qJD(2) * t262;
t343 = t266 * t368;
t321 = qJD(1) * t343;
t289 = pkin(2) * t321;
t149 = pkin(3) * t214 - pkin(9) * t215 + t289;
t366 = qJD(4) * t265;
t27 = t268 * t112 - t119 * t366 + t148 * t360 + t265 * t149;
t24 = pkin(10) * t214 + t27;
t111 = t177 * t261 - t389 * t270;
t54 = pkin(4) * t137 - pkin(10) * t136 + t111;
t3 = t267 * t24 + t264 * t54 + t75 * t361 - t363 * t64;
t4 = -qJD(5) * t23 - t24 * t264 + t267 * t54;
t312 = -t264 * t4 + t267 * t3;
t530 = -t22 * t361 - t23 * t363 + t312;
t497 = qJD(6) - t22;
t17 = -pkin(5) * t181 + t497;
t1 = qJ(6) * t137 + qJD(6) * t181 + t3;
t2 = -pkin(5) * t137 - t4;
t313 = t1 * t267 + t2 * t264;
t529 = t17 * t361 - t18 * t363 + t313;
t46 = mrSges(6,1) * t137 - mrSges(6,3) * t70;
t47 = -t137 * mrSges(7,1) + t70 * mrSges(7,2);
t48 = -mrSges(6,2) * t137 - mrSges(6,3) * t71;
t49 = -mrSges(7,2) * t71 + mrSges(7,3) * t137;
t528 = (t48 + t49) * t267 + (-t46 + t47) * t264;
t409 = mrSges(4,3) * t221;
t495 = -mrSges(4,1) * t254 - mrSges(5,1) * t182 + mrSges(5,2) * t183 - t409;
t527 = m(5) * t118 + t495;
t526 = mrSges(7,2) * t2 - mrSges(6,3) * t4 + t536;
t525 = -t70 * Ifges(7,5) / 0.2e1 + Ifges(7,6) * t446 + Ifges(6,4) * t453 + Ifges(6,6) * t445 + (Ifges(7,3) + Ifges(6,2)) * t452;
t468 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t524 = 0.2e1 * Ifges(5,2) * t446 + Ifges(5,6) * t540 - Ifges(6,6) * t452 - Ifges(7,6) * t451 - t514 * t445 - t516 * t453 - t468 + t538;
t523 = t22 * mrSges(6,1) - t23 * mrSges(6,2) + t18 * mrSges(7,3);
t78 = -t265 * t119 + t268 * t148;
t63 = -pkin(4) * t272 - t78;
t522 = m(6) * t63;
t520 = t78 * mrSges(5,1);
t519 = t79 * mrSges(5,2);
t513 = t514 * t137 + t516 * t70 - t539 * t71;
t509 = Ifges(4,5) * t215;
t508 = Ifges(4,6) * t214;
t507 = t272 * Ifges(5,5);
t506 = t272 * Ifges(5,6);
t107 = mrSges(5,1) * t214 - mrSges(5,3) * t136;
t20 = mrSges(6,1) * t71 + mrSges(6,2) * t70;
t505 = -t107 + t20;
t394 = t183 * mrSges(5,3);
t147 = mrSges(5,1) * t272 - t394;
t90 = mrSges(6,1) * t138 + mrSges(6,2) * t139;
t504 = t147 - t90;
t386 = t220 * t265;
t503 = -qJ(6) * t386 + (-t259 * t363 - qJD(6)) * t268 + (-t259 * t267 + qJ(6)) * t366 + t533;
t331 = t259 * t264 + pkin(5);
t502 = pkin(5) * t386 - t331 * t366 - t532;
t346 = t259 * t366;
t501 = t264 * t346 + t532;
t365 = qJD(4) * t267;
t500 = (-t265 * t365 - t268 * t363) * t259 + t533;
t378 = t264 * t268;
t172 = t220 * t378 + t267 * t221;
t293 = pkin(5) * t264 - qJ(6) * t267;
t288 = t259 + t293;
t294 = pkin(5) * t267 + qJ(6) * t264;
t93 = -t265 * t155 + t167 * t268;
t80 = pkin(4) * t221 - t93;
t499 = -pkin(5) * t172 + qJ(6) * t173 + (qJD(5) * t294 - qJD(6) * t267) * t265 + t288 * t360 - t80;
t498 = -qJD(6) * t264 + t181 * t293 - t79;
t26 = t138 * pkin(5) - t139 * qJ(6) + t63;
t310 = mrSges(7,1) * t264 - mrSges(7,3) * t267;
t311 = mrSges(6,1) * t264 + mrSges(6,2) * t267;
t493 = t26 * t310 + t63 * t311;
t492 = t516 * t264 + t539 * t267;
t491 = -t539 * t264 + t516 * t267;
t396 = Ifges(7,5) * t267;
t399 = Ifges(6,4) * t267;
t490 = t518 * t264 - t396 + t399;
t397 = Ifges(7,5) * t264;
t400 = Ifges(6,4) * t264;
t489 = t518 * t267 + t397 - t400;
t486 = -t366 + t386;
t359 = qJD(1) * qJD(2);
t338 = t269 * t359;
t320 = t262 * t338;
t485 = Ifges(3,5) * t320 - Ifges(3,6) * t321 - t508 + t509;
t28 = -t265 * t112 - t119 * t360 - t148 * t366 + t149 * t268;
t484 = -t265 * t28 + t268 * t27;
t131 = Ifges(7,5) * t139;
t57 = Ifges(7,6) * t181 + Ifges(7,3) * t138 + t131;
t401 = Ifges(6,4) * t139;
t60 = -Ifges(6,2) * t138 + Ifges(6,6) * t181 + t401;
t483 = t60 / 0.2e1 - t57 / 0.2e1;
t455 = -t60 / 0.2e1;
t482 = t455 + t57 / 0.2e1;
t253 = qJD(2) * t257;
t187 = t253 + t273;
t122 = t187 * t389 + t261 * t188;
t200 = pkin(2) * t263 + t257 - t322;
t371 = pkin(8) * t380 + t256;
t218 = qJ(3) * t380 + t371;
t163 = t261 * t200 + t389 * t218;
t151 = pkin(9) * t263 + t163;
t325 = pkin(2) * t343;
t168 = pkin(3) * t222 - pkin(9) * t223 + t325;
t226 = t261 * t381 - t262 * t328;
t171 = t226 * pkin(3) - t227 * pkin(9) + t242;
t42 = t268 * t122 - t151 * t366 + t265 * t168 + t171 * t360;
t34 = pkin(10) * t222 + t42;
t97 = t268 * t151 + t265 * t171;
t86 = pkin(10) * t226 + t97;
t162 = t200 * t389 - t261 * t218;
t150 = -t263 * pkin(3) - t162;
t197 = t227 * t265 - t263 * t268;
t198 = t227 * t268 + t263 * t265;
t95 = t197 * pkin(4) - t198 * pkin(10) + t150;
t413 = t264 * t95 + t267 * t86;
t121 = t187 * t261 - t389 * t188;
t160 = -qJD(4) * t197 + t223 * t268;
t161 = qJD(4) * t198 + t223 * t265;
t56 = pkin(4) * t161 - pkin(10) * t160 + t121;
t8 = -qJD(5) * t413 - t264 * t34 + t267 * t56;
t25 = -pkin(4) * t214 - t28;
t9 = pkin(5) * t71 - qJ(6) * t70 - qJD(6) * t139 + t25;
t474 = -t25 * mrSges(6,1) - t9 * mrSges(7,1);
t473 = mrSges(6,2) * t25 - mrSges(7,3) * t9;
t472 = t28 * mrSges(5,1) - t27 * mrSges(5,2);
t471 = 0.2e1 * Ifges(5,1) * t447 + Ifges(5,5) * t540 + t534;
t470 = -mrSges(7,2) * t1 - mrSges(6,3) * t3 - t525;
t224 = -pkin(8) * t321 + t249;
t234 = t371 * qJD(2);
t225 = qJD(1) * t234;
t469 = -t225 * mrSges(3,1) - t224 * mrSges(3,2) - t112 * mrSges(4,2);
t467 = -t63 * mrSges(6,1) - t26 * mrSges(7,1) + t23 * mrSges(6,3) + t535;
t437 = -t181 / 0.2e1;
t442 = -t139 / 0.2e1;
t443 = t138 / 0.2e1;
t444 = -t138 / 0.2e1;
t466 = Ifges(6,6) * t443 + Ifges(7,6) * t444 + t437 * t514 + t442 * t516 - t523;
t464 = t262 ^ 2;
t462 = m(5) / 0.2e1;
t461 = m(6) / 0.2e1;
t460 = m(7) / 0.2e1;
t404 = Ifges(5,4) * t183;
t105 = Ifges(5,2) * t182 + t404 + t506;
t448 = -t105 / 0.2e1;
t441 = t139 / 0.2e1;
t436 = t181 / 0.2e1;
t435 = -t182 / 0.2e1;
t434 = -t183 / 0.2e1;
t433 = t183 / 0.2e1;
t430 = -t220 / 0.2e1;
t429 = -t221 / 0.2e1;
t412 = mrSges(4,3) * t214;
t411 = mrSges(4,3) * t215;
t410 = mrSges(4,3) * t220;
t408 = mrSges(6,3) * t138;
t407 = mrSges(6,3) * t139;
t406 = Ifges(3,4) * t266;
t405 = Ifges(3,4) * t269;
t403 = Ifges(5,4) * t265;
t402 = Ifges(5,4) * t268;
t395 = t182 * mrSges(5,3);
t393 = t221 * Ifges(4,4);
t392 = t25 * t265;
t123 = pkin(4) * t183 - pkin(10) * t182;
t51 = t264 * t123 + t267 * t78;
t388 = t182 * t264;
t387 = t182 * t267;
t385 = t220 * t268;
t384 = t241 * t267;
t100 = -mrSges(7,2) * t138 + mrSges(7,3) * t181;
t101 = -mrSges(6,2) * t181 - t408;
t375 = t100 + t101;
t102 = mrSges(6,1) * t181 - t407;
t103 = -mrSges(7,1) * t181 + mrSges(7,2) * t139;
t374 = -t103 + t102;
t367 = qJD(4) * t264;
t357 = pkin(10) * t363;
t356 = pkin(10) * t361;
t351 = Ifges(5,5) * t136 - Ifges(5,6) * t137 + Ifges(5,3) * t214;
t347 = t269 * t368;
t345 = t259 * t360;
t344 = t264 * t360;
t334 = -t362 / 0.2e1;
t333 = t361 / 0.2e1;
t332 = t360 / 0.2e1;
t327 = t214 * mrSges(4,1) + t215 * mrSges(4,2);
t96 = -t265 * t151 + t171 * t268;
t324 = mrSges(3,3) * t349;
t323 = mrSges(3,3) * t348;
t318 = t347 / 0.2e1;
t315 = -t343 / 0.2e1;
t309 = Ifges(5,1) * t268 - t403;
t304 = -Ifges(5,2) * t265 + t402;
t303 = -Ifges(6,2) * t264 + t399;
t302 = Ifges(6,2) * t267 + t400;
t299 = Ifges(5,5) * t268 - Ifges(5,6) * t265;
t296 = Ifges(7,3) * t264 + t396;
t295 = -Ifges(7,3) * t267 + t397;
t36 = -t264 * t86 + t267 * t95;
t50 = t123 * t267 - t264 * t78;
t170 = t198 * t267 + t226 * t264;
t169 = t198 * t264 - t267 * t226;
t43 = -t265 * t122 - t151 * t360 + t168 * t268 - t171 * t366;
t85 = -pkin(4) * t226 - t96;
t7 = t264 * t56 + t267 * t34 + t95 * t361 - t363 * t86;
t283 = (Ifges(3,2) * t269 + t406) * t262;
t279 = pkin(1) * t464 * (mrSges(3,1) * t266 + mrSges(3,2) * t269);
t278 = t254 * t262 * (Ifges(3,5) * t269 - Ifges(3,6) * t266);
t277 = t266 * t464 * (Ifges(3,1) * t269 - t406);
t275 = t265 * t361 + t344;
t35 = -pkin(4) * t222 - t43;
t251 = Ifges(3,4) * t348;
t247 = -pkin(4) - t294;
t239 = -pkin(8) * t381 + t257;
t233 = -pkin(8) * t343 + t253;
t232 = t371 * qJD(1);
t231 = -pkin(8) * t349 + t252;
t230 = -t254 * mrSges(3,2) + t323;
t229 = mrSges(3,1) * t254 - t324;
t217 = Ifges(4,4) * t220;
t216 = t288 * t265;
t202 = Ifges(3,1) * t349 + t254 * Ifges(3,5) + t251;
t201 = Ifges(3,6) * t254 + qJD(1) * t283;
t192 = -t259 * t378 + t384;
t189 = -mrSges(4,2) * t254 + t410;
t186 = t268 * t331 - t384;
t185 = -qJ(6) * t268 + t494;
t179 = Ifges(5,4) * t182;
t175 = -mrSges(4,1) * t220 - mrSges(4,2) * t221;
t166 = -t221 * Ifges(4,1) + t254 * Ifges(4,5) + t217;
t165 = t220 * Ifges(4,2) + t254 * Ifges(4,6) - t393;
t146 = -mrSges(5,2) * t272 + t395;
t108 = -mrSges(5,2) * t214 - mrSges(5,3) * t137;
t106 = Ifges(5,1) * t183 + t179 + t507;
t104 = Ifges(5,5) * t183 + Ifges(5,6) * t182 + t272 * Ifges(5,3);
t89 = mrSges(7,1) * t138 - mrSges(7,3) * t139;
t88 = pkin(5) * t139 + qJ(6) * t138;
t87 = mrSges(5,1) * t137 + mrSges(5,2) * t136;
t84 = qJD(5) * t170 + t160 * t264 - t267 * t222;
t83 = -qJD(5) * t169 + t160 * t267 + t222 * t264;
t45 = pkin(5) * t169 - qJ(6) * t170 + t85;
t39 = -pkin(5) * t183 - t50;
t38 = qJ(6) * t183 + t51;
t30 = -pkin(5) * t197 - t36;
t29 = qJ(6) * t197 + t413;
t19 = mrSges(7,1) * t71 - mrSges(7,3) * t70;
t10 = pkin(5) * t84 - qJ(6) * t83 - qJD(6) * t170 + t35;
t6 = -pkin(5) * t161 - t8;
t5 = qJ(6) * t161 + qJD(6) * t197 + t7;
t11 = [(t224 * t380 + t225 * t381 - t231 * t347 - t232 * t343 - t239 * t320 - t321 * t371) * mrSges(3,3) + m(3) * (t224 * t371 - t225 * t239 - t231 * t234 + t232 * t233) + (t111 * mrSges(5,2) + t471 + t534) * t198 + (-m(4) * t125 + t527) * t121 + (t161 * t514 + t516 * t83 - t539 * t84) * t436 + (-Ifges(6,2) * t452 + Ifges(7,3) * t451 - t445 * t539 + t453 * t517 + t470 - t474) * t169 + (mrSges(4,2) * t289 + mrSges(4,3) * t111 + Ifges(4,1) * t215 - Ifges(4,4) * t214) * t227 + m(4) * (-t111 * t162 + t112 * t163 + t122 * t126 + (t237 + t370) * t325) + m(5) * (t111 * t150 + t27 * t97 + t28 * t96 + t42 * t79 + t43 * t78) + (Ifges(7,5) * t83 + Ifges(7,6) * t161) * t443 + (Ifges(6,4) * t83 + Ifges(6,6) * t161) * t444 + ((Ifges(3,1) * t266 + t405) * t262 * t318 + t283 * t315) * qJD(1) + (-t125 * t223 - t126 * t222) * mrSges(4,3) + (-t160 * t78 - t161 * t79 - t198 * t28) * mrSges(5,3) + (-Ifges(6,2) * t444 + Ifges(7,3) * t443 - t467 + t482) * t84 + m(6) * (t22 * t8 + t23 * t7 + t25 * t85 + t3 * t413 + t35 * t63 + t36 * t4) + t413 * t48 + (-0.2e1 * t279 + t277) * t359 + (Ifges(5,6) * t446 + Ifges(5,5) * t447 + mrSges(4,1) * t289 + t351 / 0.2e1 + Ifges(5,3) * t431 - Ifges(4,4) * t215 + Ifges(4,2) * t214 - t112 * mrSges(4,3) + t472) * t226 + (t161 * t18 - t26 * t83) * mrSges(7,3) + (Ifges(6,4) * t452 + Ifges(7,5) * t451 + t473 + t526 + t541) * t170 + (-t27 * mrSges(5,3) + t513 / 0.2e1 + t111 * mrSges(5,1) - t538 - t524) * t197 - t162 * t411 - t163 * t412 + t175 * t325 + t464 * (-Ifges(3,2) * t266 + t405) * t338 + (-t161 * t23 + t63 * t83) * mrSges(6,2) + qJD(2) * t278 / 0.2e1 + m(7) * (t1 * t29 + t10 * t26 + t17 * t6 + t18 * t5 + t2 * t30 + t45 * t9) + t242 * t327 + t30 * t47 + t29 * t49 + t45 * t19 + t36 * t46 + t272 * (Ifges(5,5) * t160 - Ifges(5,6) * t161 + Ifges(5,3) * t222) / 0.2e1 + (-t111 * mrSges(4,1) - t508 / 0.2e1 + t509 / 0.2e1 + (Ifges(3,5) * t318 + Ifges(3,6) * t315) * qJD(1) + t469 + t485 / 0.2e1) * t263 + (Ifges(4,1) * t223 - Ifges(4,4) * t222) * t429 + (Ifges(5,1) * t160 - Ifges(5,4) * t161 + Ifges(5,5) * t222) * t433 + t222 * t520 + t161 * t448 + t201 * t315 + t510 * t83 / 0.2e1 + (t161 * t516 + t517 * t84 + t518 * t83) * t441 - t222 * t519 + t85 * t20 + t10 * t89 + t35 * t90 + t5 * t100 + t7 * t101 + t8 * t102 + t6 * t103 + t96 * t107 + t97 * t108 + t161 * t537 + t42 * t146 + t43 * t147 + t150 * t87 + t202 * t318 + t160 * t106 / 0.2e1 + t22 * (mrSges(6,1) * t161 - mrSges(6,3) * t83) + t17 * (-mrSges(7,1) * t161 + mrSges(7,2) * t83) + t118 * (mrSges(5,1) * t161 + mrSges(5,2) * t160) + t122 * t189 + t222 * t104 / 0.2e1 + t182 * (Ifges(5,4) * t160 - Ifges(5,2) * t161 + Ifges(5,6) * t222) / 0.2e1 - t222 * t165 / 0.2e1 + t223 * t166 / 0.2e1 + t220 * (Ifges(4,4) * t223 - Ifges(4,2) * t222) / 0.2e1 + t233 * t230 - t234 * t229 + t237 * (mrSges(4,1) * t222 + mrSges(4,2) * t223) + t254 * (Ifges(4,5) * t223 - Ifges(4,6) * t222) / 0.2e1; (t448 + t537 + t523) * t366 + (t265 * t470 + t332 * t57) * t264 + t510 * (t264 * t334 - t173 / 0.2e1) + (t265 * t526 + t332 * t510 + t334 * t60) * t267 + (-t63 * t80 + t192 * t4 + t494 * t3 + (t360 * t63 + t392) * t259 + t500 * t23 + t501 * t22) * m(6) + t494 * t48 - t275 * t535 + (t22 * t531 - t23 * t275) * mrSges(6,3) + (t486 * mrSges(7,1) - mrSges(7,2) * t531) * t17 + (-Ifges(6,2) * t443 + Ifges(7,3) * t444 - t437 * t539 + t442 * t517 + t467 + t483) * t172 + (t259 * t108 + (t296 * t443 + t303 * t444) * qJD(4) + t524) * t268 + (-t385 / 0.2e1 + t332) * t106 + (-t63 * mrSges(6,2) + t26 * mrSges(7,3) + Ifges(6,4) * t443 + Ifges(7,5) * t444 + t437 * t516 + t442 * t518) * t173 + (-t278 / 0.2e1 + (-t277 / 0.2e1 + t279) * qJD(1)) * qJD(1) + (-t345 - t93) * t147 - t495 * t154 + t499 * t89 + t500 * t101 + t501 * t102 + t502 * t103 + t503 * t100 + (t1 * t185 + t17 * t502 + t18 * t503 + t186 * t2 + t216 * t9 + t26 * t499) * m(7) + (t324 + t229) * t232 + (t486 * t79 + (-t360 + t385) * t78 + t484) * mrSges(5,3) + (t296 * t451 + t303 * t452 + t310 * t9 + t333 * t57 + t489 * t453 + t491 * t445 + (Ifges(6,6) * t444 + Ifges(7,6) * t443) * qJD(4) + t471) * t265 + (-t118 * t154 - t78 * t93 - t79 * t94 + t111 * t260 + ((-t265 * t79 - t268 * t78) * qJD(4) + t484) * t259) * m(5) - t481 * t118 * (mrSges(5,1) * t265 + mrSges(5,2) * t268) + t505 * t259 * t265 + (Ifges(4,1) * t220 + t104 + t393) * t221 / 0.2e1 + (t345 - t80) * t90 + (t323 - t230) * t231 + t403 * t446 + t485 + t201 * t349 / 0.2e1 - t175 * t326 - t412 * t423 + (-t295 * t443 - t302 * t444) * t362 - t126 * t409 - t350 * t411 + t402 * t447 + (t105 / 0.2e1 + t466 - t511 / 0.2e1) * t386 + t469 - (-Ifges(3,2) * t349 + t202 + t251) * t348 / 0.2e1 + (t182 * t304 + t183 * t309 + t272 * t299) * qJD(4) / 0.2e1 + (-mrSges(5,1) * t268 + mrSges(5,2) * t265 - mrSges(4,1)) * t111 - t272 * (-Ifges(5,3) * t221 + t220 * t299) / 0.2e1 + t63 * (mrSges(6,1) * t275 + mrSges(6,2) * t276) + t26 * (mrSges(7,1) * t275 - mrSges(7,3) * t276) + (t125 * t154 - t126 * t155 - t237 * t326 + (-t111 * t389 + t112 * t261) * pkin(2)) * m(4) + t165 * t429 + (-Ifges(5,5) * t221 + t220 * t309) * t434 + (-Ifges(5,6) * t221 + t220 * t304) * t435 + t125 * t410 + t221 * t520 + t311 * t392 + t344 * t455 - t513 * t268 / 0.2e1 + (-t492 * t362 + (t265 * t514 + t268 * t491) * qJD(4)) * t436 + (-t490 * t362 + (t265 * t516 + t268 * t489) * qJD(4)) * t441 - t221 * t519 + (-t346 - t94) * t146 + (Ifges(4,2) * t221 + t166 + t217) * t430 + t185 * t49 + t186 * t47 - t155 * t189 + t192 * t46 + t216 * t19 - t237 * (-mrSges(4,1) * t221 + mrSges(4,2) * t220) - t254 * (Ifges(4,5) * t220 + Ifges(4,6) * t221) / 0.2e1 + t260 * t87; -t220 * t189 - t375 * t173 + t374 * t172 + (-t220 * t146 - t19 + (-t264 * t374 + t267 * t375 + t146) * qJD(4) - t505) * t268 + (t108 + (-t264 * t375 - t267 * t374) * qJD(5) + t481 * (-t89 + t504) + t528) * t265 - m(7) * (t17 * t172 + t173 * t18) - m(6) * (-t172 * t22 + t173 * t23) + 0.2e1 * ((t17 * t367 + t18 * t365 - t9) * t460 + (-t22 * t367 + t23 * t365 - t25) * t461 + (qJD(4) * t79 + t28) * t462 + m(5) * t79 * t430) * t268 + 0.2e1 * ((qJD(4) * t26 + t529) * t460 + (qJD(4) * t63 + t530) * t461 + (-qJD(4) * t78 + t27) * t462 + (t78 * t462 - m(7) * t26 / 0.2e1 - t522 / 0.2e1) * t220) * t265 + t327 + t527 * t221 + (-t125 * t221 - t126 * t220 + t289) * m(4); (t473 + t536) * t264 + (t139 * t489 + t181 * t491) * qJD(5) / 0.2e1 + (-t17 * t387 + t18 * t388 + t529) * mrSges(7,2) + (t22 * t387 + t23 * t388 + t530) * mrSges(6,3) + (((t17 * t267 - t18 * t264) * qJD(5) + t313) * m(7) + ((-t22 * t267 - t23 * t264) * qJD(5) + t312) * m(6) + t528) * pkin(10) + (t474 + t525) * t267 + (t506 / 0.2e1 - Ifges(5,2) * t435 + t17 * mrSges(7,1) - t118 * mrSges(5,1) + t466) * t183 + (t296 * t444 - t507 / 0.2e1 + t303 * t443 + Ifges(5,1) * t434 - t118 * mrSges(5,2) + t489 * t442 + t491 * t437 - t493) * t182 + t498 * t89 + t490 * t453 + t492 * t445 + t482 * t363 + t483 * t388 + t472 + (t493 + (-t303 / 0.2e1 + t296 / 0.2e1) * t138) * qJD(5) + (t333 - t387 / 0.2e1) * t510 + (t179 + t106) * t435 + (t395 - t146) * t78 + (-t17 * t39 - t18 * t38 + t247 * t9 + t498 * t26) * m(7) + (-pkin(4) * t25 - t22 * t50 - t23 * t51) * m(6) + t351 - pkin(4) * t20 + (t356 - t39) * t103 + (-t356 - t50) * t102 + (-t357 - t51) * t101 + (-t357 - t38) * t100 + t105 * t433 + t295 * t451 + t302 * t452 + (-t404 + t511) * t434 + (t394 + t504 - t522) * t79 + t247 * t19; (-pkin(5) * t2 + qJ(6) * t1 - t17 * t23 + t18 * t497 - t26 * t88) * m(7) + (Ifges(7,3) * t139 - t398) * t444 + (-t138 * t518 + t131 - t401 + t57) * t442 + t468 + (t138 * t17 + t139 * t18) * mrSges(7,2) + (-t138 * t516 - t139 * t539) * t437 + (t374 + t407) * t23 + (-t375 - t408) * t22 + (-Ifges(6,2) * t139 - t132 + t510) * t443 - pkin(5) * t47 + qJ(6) * t49 + t60 * t441 - t88 * t89 + qJD(6) * t100 - t26 * (mrSges(7,1) * t139 + mrSges(7,3) * t138) - t63 * (mrSges(6,1) * t139 - mrSges(6,2) * t138) + t513; -t181 * t100 + t139 * t89 + 0.2e1 * (t2 / 0.2e1 + t26 * t441 + t18 * t437) * m(7) + t47;];
tauc  = t11(:);
