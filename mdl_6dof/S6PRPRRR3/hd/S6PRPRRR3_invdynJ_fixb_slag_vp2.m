% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:48
% EndTime: 2019-03-08 20:32:17
% DurationCPUTime: 17.96s
% Computational Cost: add. (11839->635), mult. (28641->881), div. (0->0), fcn. (23422->18), ass. (0->302)
t251 = sin(pkin(12));
t400 = pkin(8) + qJ(3);
t222 = t400 * t251;
t254 = cos(pkin(12));
t223 = t400 * t254;
t259 = sin(qJ(4));
t263 = cos(qJ(4));
t170 = -t259 * t222 + t263 * t223;
t215 = t251 * t263 + t254 * t259;
t253 = sin(pkin(6));
t264 = cos(qJ(2));
t358 = t253 * t264;
t273 = t215 * t358;
t441 = qJD(1) * t273 - t215 * qJD(3) - qJD(4) * t170;
t357 = t254 * t263;
t214 = -t251 * t259 + t357;
t272 = t214 * t358;
t346 = qJD(4) * t263;
t440 = -qJD(1) * t272 - t222 * t346 + qJD(3) * t357 + (-qJD(3) * t251 - qJD(4) * t223) * t259;
t207 = t215 * qJD(4);
t474 = pkin(9) * t207 - t440;
t206 = t214 * qJD(4);
t473 = -pkin(9) * t206 + t441;
t472 = mrSges(6,2) - mrSges(7,3);
t257 = sin(qJ(6));
t261 = cos(qJ(6));
t464 = t261 * mrSges(7,1) - t257 * mrSges(7,2);
t471 = -mrSges(6,1) - t464;
t258 = sin(qJ(5));
t262 = cos(qJ(5));
t169 = -t263 * t222 - t223 * t259;
t123 = -pkin(9) * t215 + t169;
t124 = pkin(9) * t214 + t170;
t295 = t262 * t123 - t124 * t258;
t470 = -qJD(5) * t295 - t473 * t258 + t262 * t474;
t260 = sin(qJ(2));
t350 = qJD(1) * t253;
t321 = t260 * t350;
t407 = pkin(4) * t207;
t290 = t262 * t214 - t215 * t258;
t93 = qJD(5) * t290 + t206 * t262 - t207 * t258;
t160 = t214 * t258 + t215 * t262;
t94 = qJD(5) * t160 + t206 * t258 + t262 * t207;
t469 = pkin(5) * t94 - pkin(10) * t93 - t321 + t407;
t468 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t79 = t123 * t258 + t124 * t262;
t446 = -qJD(5) * t79 + t258 * t474 + t473 * t262;
t250 = qJD(4) + qJD(5);
t204 = t214 * qJD(2);
t205 = t215 * qJD(2);
t291 = t204 * t258 + t262 * t205;
t110 = t250 * t261 - t257 * t291;
t111 = t250 * t257 + t261 * t291;
t374 = mrSges(6,1) * t250 + mrSges(7,1) * t110 - mrSges(7,2) * t111 - mrSges(6,3) * t291;
t467 = t251 ^ 2 + t254 ^ 2;
t255 = cos(pkin(6));
t373 = cos(pkin(11));
t311 = t373 * t264;
t252 = sin(pkin(11));
t361 = t252 * t260;
t203 = -t255 * t361 + t311;
t249 = pkin(12) + qJ(4);
t242 = sin(t249);
t243 = cos(t249);
t362 = t252 * t253;
t466 = -t203 * t242 + t243 * t362;
t359 = t253 * t260;
t465 = -t242 * t359 + t243 * t255;
t305 = mrSges(7,1) * t257 + mrSges(7,2) * t261;
t428 = -m(4) * qJ(3) - m(5) * t400 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - t305;
t244 = qJ(5) + t249;
t236 = sin(t244);
t237 = cos(t244);
t238 = pkin(3) * t254 + pkin(2);
t308 = -mrSges(4,1) * t254 + mrSges(4,2) * t251;
t423 = m(7) * pkin(10);
t424 = m(7) * pkin(5);
t463 = -m(4) * pkin(2) - m(5) * t238 - t243 * mrSges(5,1) + t242 * mrSges(5,2) - mrSges(3,1) + t308 + (-t424 + t471) * t237 + (-t423 + t472) * t236;
t301 = -t264 * t350 + qJD(3);
t193 = -qJD(2) * t238 + t301;
t153 = -pkin(4) * t204 + t193;
t385 = t250 * Ifges(6,6);
t395 = Ifges(6,4) * t291;
t221 = qJD(2) * qJ(3) + t321;
t349 = qJD(1) * t255;
t232 = t254 * t349;
t392 = pkin(8) * qJD(2);
t163 = t232 + (-t221 - t392) * t251;
t180 = t254 * t221 + t251 * t349;
t164 = t254 * t392 + t180;
t100 = t163 * t259 + t164 * t263;
t88 = pkin(9) * t204 + t100;
t375 = t262 * t88;
t368 = t164 * t259;
t99 = t263 * t163 - t368;
t87 = -pkin(9) * t205 + t99;
t86 = qJD(4) * pkin(4) + t87;
t49 = t258 * t86 + t375;
t402 = t49 * mrSges(6,3);
t43 = pkin(10) * t250 + t49;
t310 = t262 * t204 - t205 * t258;
t68 = -pkin(5) * t310 - pkin(10) * t291 + t153;
t20 = t257 * t68 + t261 * t43;
t452 = t20 * mrSges(7,2);
t19 = -t257 * t43 + t261 * t68;
t453 = t19 * mrSges(7,1);
t126 = qJD(6) - t310;
t387 = t126 * Ifges(7,3);
t388 = t111 * Ifges(7,5);
t390 = t110 * Ifges(7,6);
t57 = t387 + t388 + t390;
t448 = t310 * Ifges(6,2);
t81 = t385 + t395 + t448;
t462 = -t153 * mrSges(6,1) + t402 + t81 / 0.2e1 - t57 / 0.2e1 + t452 - t453 + t395 / 0.2e1 + t385 / 0.2e1;
t125 = Ifges(6,4) * t310;
t380 = t258 * t88;
t48 = t262 * t86 - t380;
t42 = -pkin(5) * t250 - t48;
t281 = t42 * t305;
t109 = Ifges(7,4) * t110;
t59 = t111 * Ifges(7,1) + t126 * Ifges(7,5) + t109;
t377 = t261 * t59;
t386 = t250 * Ifges(6,5);
t403 = t48 * mrSges(6,3);
t409 = t257 / 0.2e1;
t389 = t111 * Ifges(7,4);
t58 = t110 * Ifges(7,2) + t126 * Ifges(7,6) + t389;
t449 = t291 * Ifges(6,1);
t82 = t125 + t386 + t449;
t461 = -t153 * mrSges(6,2) - t281 - t377 / 0.2e1 + t58 * t409 + t403 - t82 / 0.2e1 - t386 / 0.2e1 - t125 / 0.2e1;
t460 = -m(7) - m(6);
t186 = -pkin(4) * t214 - t238;
t80 = -pkin(5) * t290 - pkin(10) * t160 + t186;
t41 = t257 * t80 + t261 * t79;
t455 = -qJD(6) * t41 + t257 * t470 + t261 * t469;
t40 = -t257 * t79 + t261 * t80;
t454 = qJD(6) * t40 + t257 * t469 - t261 * t470;
t245 = qJDD(4) + qJDD(5);
t157 = qJD(2) * t206 + qJDD(2) * t215;
t158 = -qJD(2) * t207 + qJDD(2) * t214;
t69 = qJD(5) * t310 + t157 * t262 + t158 * t258;
t46 = qJD(6) * t110 + t245 * t257 + t261 * t69;
t47 = -qJD(6) * t111 + t245 * t261 - t257 * t69;
t16 = -mrSges(7,1) * t47 + mrSges(7,2) * t46;
t61 = mrSges(6,1) * t245 - mrSges(6,3) * t69;
t451 = t16 - t61;
t450 = mrSges(4,3) * t467;
t292 = t180 * t254 - (-t221 * t251 + t232) * t251;
t443 = t264 * t292;
t112 = -mrSges(6,2) * t250 + mrSges(6,3) * t310;
t75 = -mrSges(7,2) * t126 + mrSges(7,3) * t110;
t76 = mrSges(7,1) * t126 - mrSges(7,3) * t111;
t298 = -t257 * t76 + t261 * t75;
t442 = -t112 - t298;
t181 = -t236 * t359 + t237 * t255;
t182 = t236 * t255 + t237 * t359;
t439 = t471 * t181 + t182 * t472;
t151 = -t203 * t236 + t237 * t362;
t152 = t203 * t237 + t236 * t362;
t438 = t471 * t151 + t152 * t472;
t312 = t373 * t260;
t360 = t252 * t264;
t201 = t255 * t312 + t360;
t313 = t253 * t373;
t149 = -t201 * t236 - t237 * t313;
t150 = t201 * t237 - t236 * t313;
t437 = t471 * t149 + t150 * t472;
t342 = qJD(6) * t261;
t284 = t160 * t342 + t257 * t93;
t316 = qJD(2) * t350;
t228 = t264 * t316;
t340 = qJDD(1) * t253;
t196 = t260 * t340 + t228;
t174 = t196 + t468;
t339 = qJDD(1) * t255;
t230 = t254 * t339;
t147 = -t174 * t251 + t230;
t148 = t254 * t174 + t251 * t339;
t436 = -t147 * t251 + t148 * t254;
t70 = -qJD(5) * t291 - t157 * t258 + t158 * t262;
t67 = qJDD(6) - t70;
t21 = mrSges(7,1) * t67 - mrSges(7,3) * t46;
t22 = -mrSges(7,2) * t67 + mrSges(7,3) * t47;
t435 = -t257 * t21 + t261 * t22;
t336 = qJDD(2) * t254;
t337 = qJDD(2) * t251;
t213 = -mrSges(4,1) * t336 + mrSges(4,2) * t337;
t28 = -t70 * mrSges(6,1) + t69 * mrSges(6,2);
t95 = -t158 * mrSges(5,1) + t157 * mrSges(5,2);
t434 = -t213 - t28 - t95;
t89 = -mrSges(6,1) * t310 + mrSges(6,2) * t291;
t431 = -mrSges(5,1) * t204 + mrSges(5,2) * t205 + t308 * qJD(2) + t89;
t121 = t230 + (-pkin(8) * qJDD(2) - t174) * t251;
t122 = pkin(8) * t336 + t148;
t54 = -qJD(4) * t100 + t263 * t121 - t122 * t259;
t35 = qJDD(4) * pkin(4) - pkin(9) * t157 + t54;
t53 = -qJD(4) * t368 + t259 * t121 + t263 * t122 + t163 * t346;
t36 = pkin(9) * t158 + t53;
t11 = -qJD(5) * t49 - t258 * t36 + t262 * t35;
t227 = t260 * t316;
t195 = t264 * t340 - t227;
t282 = qJDD(3) - t195;
t167 = -qJDD(2) * t238 + t282;
t102 = -pkin(4) * t158 + t167;
t23 = -pkin(5) * t70 - pkin(10) * t69 + t102;
t344 = qJD(5) * t262;
t345 = qJD(5) * t258;
t10 = t258 * t35 + t262 * t36 + t86 * t344 - t345 * t88;
t7 = pkin(10) * t245 + t10;
t2 = qJD(6) * t19 + t23 * t257 + t261 * t7;
t3 = -qJD(6) * t20 + t23 * t261 - t257 * t7;
t429 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t90 = pkin(5) * t291 - pkin(10) * t310;
t422 = t46 / 0.2e1;
t421 = t47 / 0.2e1;
t418 = t67 / 0.2e1;
t415 = -t110 / 0.2e1;
t414 = -t111 / 0.2e1;
t413 = t111 / 0.2e1;
t412 = -t126 / 0.2e1;
t410 = t205 / 0.2e1;
t408 = pkin(4) * t205;
t406 = pkin(4) * t258;
t405 = pkin(4) * t262;
t399 = mrSges(5,3) * t204;
t398 = mrSges(7,3) * t257;
t397 = mrSges(7,3) * t261;
t396 = Ifges(5,4) * t205;
t394 = Ifges(7,4) * t257;
t393 = Ifges(7,4) * t261;
t391 = t100 * mrSges(5,3);
t370 = t160 * t257;
t369 = t160 * t261;
t348 = qJD(2) * t260;
t343 = qJD(6) * t257;
t335 = Ifges(7,5) * t46 + Ifges(7,6) * t47 + Ifges(7,3) * t67;
t332 = m(4) + m(5) - t460;
t328 = t257 * t358;
t327 = t261 * t358;
t324 = t377 / 0.2e1;
t320 = t253 * t348;
t315 = -t343 / 0.2e1;
t314 = t151 * pkin(5) + pkin(10) * t152;
t309 = t466 * pkin(4);
t304 = Ifges(7,1) * t261 - t394;
t303 = -Ifges(7,2) * t257 + t393;
t302 = Ifges(7,5) * t261 - Ifges(7,6) * t257;
t300 = t19 * t261 + t20 * t257;
t299 = -t19 * t257 + t20 * t261;
t197 = -t251 * t359 + t254 * t255;
t198 = t251 * t255 + t254 * t359;
t129 = t197 * t263 - t198 * t259;
t130 = t197 * t259 + t198 * t263;
t294 = t262 * t129 - t130 * t258;
t84 = t129 * t258 + t130 * t262;
t289 = t465 * pkin(4);
t72 = -t257 * t84 - t327;
t285 = -t261 * t84 + t328;
t283 = t160 * t343 - t261 * t93;
t280 = t110 * t303;
t279 = t111 * t304;
t278 = t126 * t302;
t274 = -t201 * t242 - t243 * t313;
t271 = t274 * pkin(4);
t270 = -qJD(6) * t300 - t3 * t257;
t269 = t2 * t261 + t270;
t14 = t46 * Ifges(7,4) + t47 * Ifges(7,2) + t67 * Ifges(7,6);
t15 = t46 * Ifges(7,1) + t47 * Ifges(7,4) + t67 * Ifges(7,5);
t8 = -pkin(5) * t245 - t11;
t267 = t11 * mrSges(6,1) - t10 * mrSges(6,2) + t15 * t409 + t2 * t397 - t8 * t464 + t261 * t14 / 0.2e1 + Ifges(6,3) * t245 + (Ifges(7,1) * t257 + t393) * t422 + (Ifges(7,2) * t261 + t394) * t421 + t58 * t315 + (Ifges(7,5) * t257 + Ifges(7,6) * t261) * t418 + Ifges(6,6) * t70 + Ifges(6,5) * t69 + (t281 + t324) * qJD(6) + (t280 + t279 + t278) * qJD(6) / 0.2e1;
t265 = qJD(2) ^ 2;
t246 = -pkin(9) - t400;
t240 = -pkin(5) - t405;
t219 = pkin(4) * t243 + t238;
t217 = -qJD(2) * pkin(2) + t301;
t202 = t255 * t360 + t312;
t200 = -t255 * t311 + t361;
t194 = Ifges(5,4) * t204;
t185 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t205;
t184 = -qJD(4) * mrSges(5,2) + t399;
t183 = -qJDD(2) * pkin(2) + t282;
t178 = t181 * pkin(5);
t145 = t149 * pkin(5);
t134 = t205 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t194;
t133 = t204 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t396;
t132 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t158;
t131 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t157;
t98 = -qJD(2) * t273 - qJD(4) * t130;
t97 = qJD(2) * t272 + qJD(4) * t129;
t74 = t408 + t90;
t62 = -mrSges(6,2) * t245 + mrSges(6,3) * t70;
t51 = t262 * t87 - t380;
t50 = t258 * t87 + t375;
t30 = qJD(5) * t84 + t258 * t97 - t262 * t98;
t29 = qJD(5) * t294 + t258 * t98 + t262 * t97;
t27 = t257 * t90 + t261 * t48;
t26 = -t257 * t48 + t261 * t90;
t25 = t257 * t74 + t261 * t51;
t24 = -t257 * t51 + t261 * t74;
t18 = qJD(6) * t285 - t257 * t29 + t261 * t320;
t17 = qJD(6) * t72 + t257 * t320 + t261 * t29;
t1 = [t29 * t112 + t129 * t131 + t130 * t132 + t17 * t75 + t18 * t76 + t97 * t184 + t98 * t185 + t72 * t21 - t285 * t22 + t84 * t62 - t451 * t294 - t374 * t30 + (-t197 * t251 + t198 * t254) * qJDD(2) * mrSges(4,3) + (-m(2) - m(3) - t332) * g(3) + m(6) * (t10 * t84 + t11 * t294 + t29 * t49 - t30 * t48) + m(5) * (t100 * t97 + t129 * t54 + t130 * t53 + t98 * t99) + m(4) * (t147 * t197 + t148 * t198) + m(7) * (t17 * t20 + t18 * t19 - t2 * t285 - t294 * t8 + t3 * t72 + t30 * t42) + ((-t265 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + qJD(2) * t431) * t260 + (qJDD(2) * mrSges(3,1) + (-mrSges(3,2) + t450) * t265 + t434) * t264 + m(6) * (-t102 * t264 + t153 * t348) + m(5) * (-t167 * t264 + t193 * t348) + m(4) * (qJD(2) * t443 - t183 * t264 + t217 * t348) + m(3) * (t195 * t264 + t196 * t260)) * t253 + (m(3) * t255 ^ 2 + m(2)) * qJDD(1); (t19 * t455 + t2 * t41 + t20 * t454 - t295 * t8 + t3 * t40 - t42 * t446) * m(7) - t451 * t295 - (Ifges(7,6) * t421 + Ifges(7,5) * t422 + Ifges(7,3) * t418 - Ifges(6,4) * t69 - Ifges(6,2) * t70 - Ifges(6,6) * t245 + t102 * mrSges(6,1) + t335 / 0.2e1 - t10 * mrSges(6,3) + t429) * t290 + t454 * t75 + t455 * t76 - t94 * t452 + (mrSges(5,2) * t167 - mrSges(5,3) * t54 + Ifges(5,1) * t157 + Ifges(5,4) * t158 + Ifges(5,5) * qJDD(4)) * t215 + (Ifges(4,4) * t251 + Ifges(4,2) * t254) * t336 + (Ifges(4,1) * t251 + Ifges(4,4) * t254) * t337 - t470 * t112 + (t10 * t79 + t102 * t186 + t11 * t295 + t153 * t407 + t446 * t48 - t470 * t49) * m(6) + (-mrSges(5,1) * t167 + mrSges(5,3) * t53 + Ifges(5,4) * t157 + Ifges(5,2) * t158 + Ifges(5,6) * qJDD(4)) * t214 + (-m(5) * t193 - m(6) * t153 - t431) * t321 + (t460 * (-t202 * t219 - t203 * t246) + t428 * t203 - t463 * t202) * g(1) + (t460 * (-t200 * t219 - t201 * t246) + t428 * t201 - t463 * t200) * g(2) + (t102 * mrSges(6,2) - t11 * mrSges(6,3) + Ifges(6,1) * t69 + Ifges(6,4) * t70 + Ifges(6,5) * t245 + t302 * t418 + t303 * t421 + t304 * t422 + t305 * t8 + t315 * t59) * t160 + (-t196 + t228) * mrSges(3,2) + (t436 + t467 * (-t228 + t468)) * mrSges(4,3) + (-pkin(2) * t183 + t292 * qJD(3) + t436 * qJ(3) - (t217 * t260 + t443) * t350) * m(4) - t284 * t58 / 0.2e1 + t440 * t184 + t441 * t185 + (t100 * t440 - t167 * t238 + t169 * t54 + t170 * t53 + t441 * t99) * m(5) + t93 * t324 + (t19 * t283 - t2 * t370 - t20 * t284 - t3 * t369) * mrSges(7,3) + (-Ifges(7,1) * t283 - Ifges(7,4) * t284 + Ifges(7,5) * t94) * t413 - t207 * t391 + Ifges(3,3) * qJDD(2) + t250 * (Ifges(6,5) * t93 - Ifges(6,6) * t94) / 0.2e1 - t238 * t95 - pkin(2) * t213 + t206 * t134 / 0.2e1 - t207 * t133 / 0.2e1 + t186 * t28 + t169 * t131 + t170 * t132 + t89 * t407 - t93 * t403 + t153 * (mrSges(6,1) * t94 + mrSges(6,2) * t93) - t94 * t402 + t93 * t82 / 0.2e1 + t94 * t57 / 0.2e1 - t94 * t81 / 0.2e1 + t79 * t62 - t14 * t370 / 0.2e1 + t15 * t369 / 0.2e1 + t40 * t21 + t41 * t22 + t94 * t453 + (t460 * t219 * t358 + (t463 * t264 + (-t246 * t460 + t428) * t260) * t253) * g(3) + t446 * t374 - t99 * t206 * mrSges(5,3) + t42 * (mrSges(7,1) * t284 - mrSges(7,2) * t283) + t110 * (-Ifges(7,4) * t283 - Ifges(7,2) * t284 + Ifges(7,6) * t94) / 0.2e1 + t126 * (-Ifges(7,5) * t283 - Ifges(7,6) * t284 + Ifges(7,3) * t94) / 0.2e1 + t183 * t308 + t310 * (Ifges(6,4) * t93 - Ifges(6,2) * t94) / 0.2e1 + t291 * (Ifges(6,1) * t93 - Ifges(6,4) * t94) / 0.2e1 + (Ifges(5,1) * t206 - Ifges(5,4) * t207) * t410 + qJD(4) * (Ifges(5,5) * t206 - Ifges(5,6) * t207) / 0.2e1 + t204 * (Ifges(5,4) * t206 - Ifges(5,2) * t207) / 0.2e1 + t193 * (mrSges(5,1) * t207 + mrSges(5,2) * t206) + (t195 + t227) * mrSges(3,1); t374 * t291 - t265 * t450 + t261 * t21 + t257 * t22 + t205 * t185 - t204 * t184 + t442 * t310 + t298 * qJD(6) + (t126 * t299 + t2 * t257 + t261 * t3 - t291 * t42) * m(7) + (t291 * t48 - t310 * t49 + t102) * m(6) + (-t100 * t204 + t205 * t99 + t167) * m(5) + (-qJD(2) * t292 + t183) * m(4) + (-g(1) * t202 - g(2) * t200 + g(3) * t358) * t332 - t434; (t399 - t184) * t99 + (-t19 * t24 - t20 * t25 - t42 * t50 + t240 * t8 + (t258 * t42 + t262 * t299) * qJD(5) * pkin(4)) * m(7) + (m(7) * t269 - t342 * t76 - t343 * t75 + t435) * (pkin(10) + t406) - t442 * pkin(4) * t344 + (t302 * t412 + t304 * t414 + t303 * t415 - t449 / 0.2e1 + t19 * t397 + t20 * t398 + t461) * t310 + (Ifges(7,3) * t412 + Ifges(7,5) * t414 + Ifges(7,6) * t415 + t448 / 0.2e1 + t462) * t291 + t267 + (-m(6) * t271 - t274 * mrSges(5,1) - (-t201 * t243 + t242 * t313) * mrSges(5,2) - m(7) * (t150 * pkin(10) + t145 + t271) + t437) * g(2) + t61 * t405 + t62 * t406 + t133 * t410 + t205 * t391 + ((t10 * t258 + t11 * t262 + (-t258 * t48 + t262 * t49) * qJD(5)) * pkin(4) - t153 * t408 + t48 * t50 - t49 * t51) * m(6) + (-t19 * t342 - t20 * t343) * mrSges(7,3) + Ifges(5,3) * qJDD(4) + t240 * t16 - t193 * (mrSges(5,1) * t205 + mrSges(5,2) * t204) - qJD(4) * (Ifges(5,5) * t204 - Ifges(5,6) * t205) / 0.2e1 + t100 * t185 - t89 * t408 + Ifges(5,5) * t157 + Ifges(5,6) * t158 - t3 * t398 - t205 * (Ifges(5,1) * t204 - t396) / 0.2e1 - t51 * t112 - t24 * t76 - t25 * t75 - t53 * mrSges(5,2) + t54 * mrSges(5,1) + (-t465 * mrSges(5,1) - (-t242 * t255 - t243 * t359) * mrSges(5,2) - m(6) * t289 - m(7) * (pkin(10) * t182 + t178 + t289) + t439) * g(3) + (-t466 * mrSges(5,1) - (-t203 * t243 - t242 * t362) * mrSges(5,2) - m(6) * t309 - m(7) * (t309 + t314) + t438) * g(1) + t374 * (-pkin(4) * t345 + t50) - (-Ifges(5,2) * t205 + t134 + t194) * t204 / 0.2e1; ((-t257 * t75 - t261 * t76) * qJD(6) + (-g(2) * t150 - g(3) * t182) * m(7) + t435) * pkin(10) + (-m(7) * t145 + t437) * g(2) + (-m(7) * t314 + t438) * g(1) + (-m(7) * t178 + t439) * g(3) - m(7) * (t19 * t26 + t20 * t27 + t42 * t49) + t374 * t49 + t269 * t423 - t8 * t424 + (-t278 / 0.2e1 - t280 / 0.2e1 - t279 / 0.2e1 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t291 + t300 * mrSges(7,3) + t461) * t310 + t267 + (-t387 / 0.2e1 - t390 / 0.2e1 - t388 / 0.2e1 + t462) * t291 - t48 * t112 - t26 * t76 - t27 * t75 - pkin(5) * t16 + t270 * mrSges(7,3); -t42 * (mrSges(7,1) * t111 + mrSges(7,2) * t110) + (Ifges(7,1) * t110 - t389) * t414 + t58 * t413 + (Ifges(7,5) * t110 - Ifges(7,6) * t111) * t412 - t19 * t75 + t20 * t76 - g(1) * ((-t152 * t257 + t202 * t261) * mrSges(7,1) + (-t152 * t261 - t202 * t257) * mrSges(7,2)) - g(2) * ((-t150 * t257 + t200 * t261) * mrSges(7,1) + (-t150 * t261 - t200 * t257) * mrSges(7,2)) - g(3) * ((-t182 * t257 - t327) * mrSges(7,1) + (-t182 * t261 + t328) * mrSges(7,2)) + (t110 * t19 + t111 * t20) * mrSges(7,3) + t335 + (-Ifges(7,2) * t111 + t109 + t59) * t415 + t429;];
tau  = t1;
