% Calculate vector of inverse dynamics joint torques for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:55:07
% EndTime: 2019-12-05 18:55:34
% DurationCPUTime: 12.86s
% Computational Cost: add. (8055->606), mult. (17896->833), div. (0->0), fcn. (13566->14), ass. (0->300)
t238 = qJ(4) + qJ(5);
t231 = sin(t238);
t241 = sin(qJ(4));
t369 = mrSges(5,2) * t241;
t464 = -mrSges(6,2) * t231 - t369;
t461 = -mrSges(6,3) - mrSges(5,3);
t246 = cos(qJ(4));
t225 = pkin(3) * t246 + pkin(2);
t418 = m(5) * pkin(2);
t233 = cos(t238);
t373 = mrSges(5,1) * t246;
t459 = -mrSges(6,1) * t233 - t373;
t463 = m(6) * t225 + t418 - t459;
t462 = m(5) + m(6);
t239 = qJ(2) + qJ(3);
t232 = sin(t239);
t460 = t464 * t232;
t242 = sin(qJ(3));
t243 = sin(qJ(2));
t247 = cos(qJ(3));
t248 = cos(qJ(2));
t193 = t242 * t248 + t243 * t247;
t181 = t193 * qJD(1);
t237 = qJD(2) + qJD(3);
t152 = t181 * t246 + t237 * t241;
t311 = qJD(1) * t248;
t312 = qJD(1) * t243;
t180 = -t242 * t312 + t247 * t311;
t177 = qJD(4) - t180;
t170 = qJD(5) + t177;
t151 = -t181 * t241 + t237 * t246;
t240 = sin(qJ(5));
t245 = cos(qJ(5));
t283 = t245 * t151 - t152 * t240;
t444 = Ifges(5,3) * t177;
t445 = Ifges(5,6) * t151;
t99 = t151 * t240 + t152 * t245;
t457 = Ifges(5,5) * t152 + Ifges(6,5) * t99 + Ifges(6,6) * t283 + Ifges(6,3) * t170 + t444 + t445;
t297 = pkin(1) * t311;
t128 = -pkin(2) * t180 - pkin(5) * t181 - t297;
t356 = pkin(1) * qJD(2);
t229 = t242 * t356;
t199 = pkin(5) * t237 + t229;
t100 = t246 * t128 - t199 * t241;
t456 = t100 * mrSges(5,1);
t101 = t128 * t241 + t199 * t246;
t455 = t101 * mrSges(5,2);
t333 = t240 * t241;
t266 = -t245 * t246 + t333;
t429 = qJD(4) + qJD(5);
t454 = t429 * t266;
t453 = t461 * t232;
t235 = t248 * pkin(1);
t371 = mrSges(4,2) * t232;
t452 = -m(4) * t235 - mrSges(3,1) * t248 + mrSges(3,2) * t243 + t371;
t191 = t242 * t243 - t247 * t248;
t260 = t191 * qJD(3);
t146 = -qJD(2) * t191 - t260;
t308 = qJD(4) * t246;
t263 = t146 * t241 + t193 * t308;
t355 = pkin(1) * qJD(3);
t291 = qJD(2) * t355;
t348 = pkin(1) * qJDD(2);
t187 = t242 * t348 + t247 * t291;
t236 = qJDD(2) + qJDD(3);
t169 = pkin(5) * t236 + t187;
t307 = qJD(1) * qJD(2);
t195 = qJDD(1) * t248 - t243 * t307;
t196 = qJDD(1) * t243 + t248 * t307;
t114 = -qJD(1) * t260 + t195 * t242 + t196 * t247;
t115 = -qJD(3) * t181 + t195 * t247 - t242 * t196;
t383 = pkin(1) * t195;
t57 = -pkin(2) * t115 - pkin(5) * t114 - t383;
t26 = qJD(4) * t100 + t169 * t246 + t241 * t57;
t27 = -qJD(4) * t101 - t169 * t241 + t246 * t57;
t451 = -t27 * t241 + t246 * t26;
t113 = qJDD(4) - t115;
t109 = qJDD(5) + t113;
t400 = t109 / 0.2e1;
t72 = qJD(4) * t151 + t114 * t246 + t236 * t241;
t73 = -qJD(4) * t152 - t114 * t241 + t236 * t246;
t20 = -qJD(5) * t99 - t240 * t72 + t245 * t73;
t415 = t20 / 0.2e1;
t19 = qJD(5) * t283 + t240 * t73 + t245 * t72;
t416 = t19 / 0.2e1;
t419 = Ifges(6,1) * t416 + Ifges(6,4) * t415 + Ifges(6,5) * t400;
t420 = Ifges(6,4) * t416 + Ifges(6,2) * t415 + Ifges(6,6) * t400;
t409 = t72 / 0.2e1;
t408 = t73 / 0.2e1;
t399 = t113 / 0.2e1;
t139 = pkin(2) * t181 - pkin(5) * t180;
t298 = pkin(1) * t312;
t127 = t139 + t298;
t330 = t241 * t245;
t192 = t240 * t246 + t330;
t145 = t429 * t192;
t382 = pkin(1) * t242;
t224 = pkin(5) + t382;
t293 = t247 * t355;
t176 = t181 * pkin(3);
t95 = t127 * t246 + t176;
t450 = -t127 * t330 - t145 * t224 - t240 * t95 - t266 * t293;
t449 = t127 * t333 - t192 * t293 + t224 * t454 - t245 * t95;
t175 = Ifges(4,4) * t180;
t448 = Ifges(4,5) * t237;
t447 = Ifges(4,2) * t180;
t446 = Ifges(4,6) * t237;
t443 = t248 * Ifges(3,2);
t417 = m(6) * pkin(3);
t442 = mrSges(5,1) + t417;
t295 = t247 * t356;
t119 = t139 * t241 + t246 * t295;
t118 = t246 * t139 - t241 * t295;
t90 = t118 + t176;
t441 = -pkin(5) * t145 - t119 * t245 - t240 * t90;
t440 = pkin(5) * t454 + t119 * t240 - t245 * t90;
t200 = -pkin(2) * t237 - t295;
t275 = mrSges(5,1) * t241 + mrSges(5,2) * t246;
t439 = t200 * t275;
t367 = mrSges(4,3) * t181;
t438 = -mrSges(4,1) * t237 - mrSges(5,1) * t151 + mrSges(5,2) * t152 + t367;
t309 = qJD(4) * t241;
t292 = pkin(3) * t309;
t294 = t242 * t355;
t341 = t180 * t241;
t303 = pkin(3) * t341;
t437 = t292 + t294 - t303;
t381 = pkin(1) * t243;
t436 = t232 * t463 + t462 * t381;
t249 = cos(qJ(1));
t234 = cos(t239);
t334 = t234 * t249;
t435 = t249 * t460 + t334 * t461;
t244 = sin(qJ(1));
t335 = t234 * t244;
t434 = t244 * t460 + t335 * t461;
t45 = mrSges(5,1) * t113 - mrSges(5,3) * t72;
t46 = -mrSges(5,2) * t113 + mrSges(5,3) * t73;
t433 = -t241 * t45 + t246 * t46;
t430 = Ifges(3,6) * qJDD(2);
t428 = -mrSges(4,3) - mrSges(3,3) + mrSges(2,2);
t220 = t234 * mrSges(4,1);
t427 = -t220 - mrSges(2,1) + t452;
t16 = pkin(3) * t113 + t27;
t347 = t101 * t240;
t82 = pkin(3) * t177 + t100;
t42 = t245 * t82 - t347;
t3 = qJD(5) * t42 + t16 * t240 + t245 * t26;
t346 = t101 * t245;
t43 = t240 * t82 + t346;
t4 = -qJD(5) * t43 + t16 * t245 - t240 * t26;
t426 = t4 * mrSges(6,1) - t3 * mrSges(6,2);
t364 = mrSges(5,3) * t151;
t116 = -mrSges(5,2) * t177 + t364;
t363 = mrSges(5,3) * t152;
t117 = mrSges(5,1) * t177 - t363;
t268 = t100 * t246 + t101 * t241;
t425 = m(5) * t268 + t241 * t116 + t246 * t117;
t424 = t27 * mrSges(5,1) - t26 * mrSges(5,2);
t121 = t192 * t180;
t122 = t266 * t180;
t360 = Ifges(4,4) * t181;
t129 = t360 + t446 + t447;
t130 = Ifges(4,1) * t181 + t175 + t448;
t132 = -pkin(3) * t151 + t200;
t186 = -t242 * t291 + t247 * t348;
t168 = -pkin(2) * t236 - t186;
t23 = Ifges(5,4) * t72 + Ifges(5,2) * t73 + Ifges(5,6) * t113;
t269 = Ifges(5,5) * t246 - Ifges(5,6) * t241;
t357 = Ifges(5,4) * t246;
t271 = -Ifges(5,2) * t241 + t357;
t358 = Ifges(5,4) * t241;
t273 = Ifges(5,1) * t246 - t358;
t280 = mrSges(4,3) * t295;
t285 = -t309 / 0.2e1;
t387 = t246 / 0.2e1;
t150 = Ifges(5,4) * t151;
t86 = Ifges(5,1) * t152 + Ifges(5,5) * t177 + t150;
t290 = t86 * t387;
t340 = t180 * t246;
t385 = mrSges(6,3) * t43;
t389 = t181 / 0.2e1;
t392 = -t177 / 0.2e1;
t393 = t170 / 0.2e1;
t394 = -t170 / 0.2e1;
t396 = -t152 / 0.2e1;
t397 = -t151 / 0.2e1;
t403 = t99 / 0.2e1;
t404 = -t99 / 0.2e1;
t405 = t283 / 0.2e1;
t406 = -t283 / 0.2e1;
t93 = Ifges(6,4) * t283;
t41 = Ifges(6,1) * t99 + Ifges(6,5) * t170 + t93;
t410 = t41 / 0.2e1;
t411 = -t41 / 0.2e1;
t384 = Ifges(6,4) * t99;
t40 = Ifges(6,2) * t283 + Ifges(6,6) * t170 + t384;
t412 = t40 / 0.2e1;
t413 = -t40 / 0.2e1;
t414 = Ifges(5,1) * t409 + Ifges(5,4) * t408 + Ifges(5,5) * t399;
t60 = -pkin(3) * t73 + t168;
t359 = Ifges(5,4) * t152;
t85 = Ifges(5,2) * t151 + Ifges(5,6) * t177 + t359;
t423 = (t100 * t340 + t101 * t341 + t451) * mrSges(5,3) + (mrSges(6,1) * t60 - mrSges(6,3) * t3 - 0.2e1 * t420) * t266 + (mrSges(6,2) * t60 + 0.2e1 * t419) * t192 + (t285 + t341 / 0.2e1) * t85 - (Ifges(6,4) * t403 + Ifges(6,2) * t405 + Ifges(6,6) * t393 + t385 + t412) * t145 - t43 * (-mrSges(6,2) * t181 - mrSges(6,3) * t121) + (-Ifges(6,1) * t122 - Ifges(6,4) * t121 + Ifges(6,5) * t181) * t404 + (-Ifges(6,4) * t122 - Ifges(6,2) * t121 + Ifges(6,6) * t181) * t406 + (-Ifges(6,5) * t122 - Ifges(6,6) * t121 + Ifges(6,3) * t181) * t394 - t42 * (mrSges(6,1) * t181 + mrSges(6,3) * t122) + (t151 * t271 + t152 * t273 + t177 * t269) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t181 + t130 + t175) * t180 / 0.2e1 + t168 * (t369 - t373) + (Ifges(5,2) * t246 + t358) * t408 + (Ifges(5,1) * t241 + t357) * t409 - t122 * t411 - t121 * t413 + t241 * t414 + t23 * t387 + t129 * t389 + (Ifges(5,3) * t181 + t180 * t269) * t392 + (Ifges(5,5) * t181 + t180 * t273) * t396 + (Ifges(5,6) * t181 + t180 * t271) * t397 + (Ifges(5,5) * t241 + Ifges(5,6) * t246) * t399 + Ifges(4,5) * t114 + Ifges(4,6) * t115 + (t439 + t290) * qJD(4) - (Ifges(6,1) * t403 + Ifges(6,4) * t405 + Ifges(6,5) * t393 + t410) * t454 + (-t192 * t4 + t42 * t454) * mrSges(6,3) + ((t122 - t454) * mrSges(6,2) + (-t121 + t145) * mrSges(6,1)) * t132 + (mrSges(4,1) * t181 + mrSges(4,2) * t180) * t297 + t180 * t280 - t180 * t439 - t181 * t456 + t186 * mrSges(4,1) - t187 * mrSges(4,2) - (Ifges(4,1) * t180 - t360 + t457) * t181 / 0.2e1 + Ifges(4,3) * t236 - t237 * (Ifges(4,5) * t180 - Ifges(4,6) * t181) / 0.2e1 + t181 * t455 - t86 * t340 / 0.2e1;
t422 = -t220 + t453 + (t459 - t464) * t234;
t53 = -mrSges(6,1) * t283 + mrSges(6,2) * t99;
t401 = pkin(3) * t53;
t395 = t152 / 0.2e1;
t386 = t248 / 0.2e1;
t380 = pkin(1) * t247;
t379 = pkin(3) * t152;
t376 = g(3) * t232;
t221 = t232 * pkin(5);
t222 = t234 * pkin(2);
t370 = mrSges(4,2) * t234;
t366 = mrSges(5,3) * t100;
t365 = mrSges(5,3) * t101;
t362 = Ifges(3,4) * t243;
t361 = Ifges(3,4) * t248;
t142 = pkin(2) * t191 - pkin(5) * t193 - t235;
t345 = t142 * t241;
t344 = t142 * t246;
t339 = t193 * t241;
t338 = t193 * t246;
t336 = t232 * t249;
t201 = t234 * t225;
t331 = t241 * t244;
t329 = t241 * t249;
t326 = t244 * t246;
t324 = t246 * t249;
t164 = t231 * t335 + t233 * t249;
t165 = t231 * t249 - t233 * t335;
t323 = -t164 * mrSges(6,1) + t165 * mrSges(6,2);
t166 = -t231 * t334 + t233 * t244;
t167 = t231 * t244 + t233 * t334;
t322 = t166 * mrSges(6,1) - t167 * mrSges(6,2);
t317 = t201 + t221;
t314 = pkin(5) * t336 + t235 * t249;
t313 = t222 + t221;
t310 = qJD(2) * t243;
t305 = t241 * t417;
t304 = m(4) * t248 * pkin(1) ^ 2;
t302 = Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t109;
t301 = Ifges(5,5) * t72 + Ifges(5,6) * t73 + Ifges(5,3) * t113;
t296 = pkin(1) * t310;
t289 = t142 * t309;
t288 = t142 * t308;
t284 = t307 / 0.2e1;
t282 = t243 * t304;
t281 = mrSges(4,3) * t229;
t279 = -t235 - t221;
t147 = t237 * t193;
t89 = pkin(2) * t147 - pkin(5) * t146 + t296;
t278 = pkin(3) * t147 - qJD(5) * t345 + t246 * t89 - t289;
t274 = -mrSges(6,1) * t231 - mrSges(6,2) * t233;
t272 = t362 + t443;
t270 = Ifges(3,5) * t248 - Ifges(3,6) * t243;
t267 = -t100 * t241 + t101 * t246;
t264 = t302 + t426;
t173 = -t234 * t329 + t326;
t171 = t234 * t331 + t324;
t262 = -t146 * t246 + t193 * t309;
t261 = t243 * (Ifges(3,1) * t248 - t362);
t110 = pkin(3) * t191 + t344;
t257 = qJD(5) * t110 + t241 * t89 + t288;
t254 = -qJD(4) * t268 + t451;
t253 = t370 + (mrSges(4,1) + t463) * t232;
t228 = Ifges(3,4) * t311;
t226 = -pkin(2) - t380;
t213 = pkin(5) * t334;
t211 = pkin(5) * t335;
t203 = -t225 - t380;
t185 = t266 * pkin(5);
t184 = t192 * pkin(5);
t179 = Ifges(3,1) * t312 + Ifges(3,5) * qJD(2) + t228;
t178 = Ifges(3,6) * qJD(2) + qJD(1) * t272;
t174 = t234 * t324 + t331;
t172 = -t234 * t326 + t329;
t157 = -mrSges(4,2) * t237 + mrSges(4,3) * t180;
t156 = t266 * t224;
t155 = t192 * t224;
t153 = t229 + t303;
t138 = -mrSges(4,1) * t180 + mrSges(4,2) * t181;
t135 = t266 * t193;
t134 = t192 * t193;
t79 = mrSges(6,1) * t170 - mrSges(6,3) * t99;
t78 = -mrSges(6,2) * t170 + mrSges(6,3) * t283;
t67 = t110 * t240 + t142 * t330;
t66 = t110 * t245 - t142 * t333;
t52 = t100 * t245 - t347;
t51 = -t100 * t240 - t346;
t38 = -t146 * t192 + t193 * t454;
t37 = -t145 * t193 - t146 * t266;
t30 = -mrSges(5,1) * t73 + mrSges(5,2) * t72;
t13 = -mrSges(6,2) * t109 + mrSges(6,3) * t20;
t12 = mrSges(6,1) * t109 - mrSges(6,3) * t19;
t11 = -t240 * t257 + t245 * t278;
t10 = t240 * t278 + t245 * t257;
t9 = -mrSges(6,1) * t20 + mrSges(6,2) * t19;
t1 = [t263 * (-t85 / 0.2e1 + t401) + (t179 * t386 + t270 * qJD(2) / 0.2e1) * qJD(2) + (Ifges(6,5) * t37 + Ifges(6,6) * t38) * t393 + (-t172 * mrSges(5,1) - t165 * mrSges(6,1) - t171 * mrSges(5,2) - t164 * mrSges(6,2) + (-t305 + t428) * t249 + (-m(6) * (t279 - t201) - m(5) * (t279 - t222) - t427 - t453) * t244) * g(1) + m(5) * (qJD(4) * t267 + t241 * t26 + t246 * t27) * t142 + (t272 / 0.2e1 + Ifges(3,2) * t386 + t304) * t195 + (-Ifges(6,5) * t135 - Ifges(6,6) * t134) * t400 + (-Ifges(6,4) * t135 - Ifges(6,2) * t134) * t415 + (-Ifges(6,1) * t135 - Ifges(6,4) * t134) * t416 + t60 * (mrSges(6,1) * t134 - mrSges(6,2) * t135) + (-t134 * t3 + t135 * t4 - t37 * t42 + t38 * t43) * mrSges(6,3) + (t302 + t301) * t191 / 0.2e1 + (-m(6) * (pkin(3) * t331 + t225 * t334 + t314) - t167 * mrSges(6,1) - t166 * mrSges(6,2) - m(5) * (pkin(2) * t334 + t314) - t174 * mrSges(5,1) - t173 * mrSges(5,2) + t461 * t336 + t427 * t249 + t428 * t244) * g(2) + t177 * (-Ifges(5,5) * t262 - Ifges(5,6) * t263) / 0.2e1 + (t100 * t262 - t101 * t263 - t26 * t339 - t27 * t338) * mrSges(5,3) + (Ifges(6,4) * t37 + Ifges(6,2) * t38) * t405 + t261 * t284 + (Ifges(6,6) * t405 - Ifges(4,4) * t389 + Ifges(6,3) * t393 + Ifges(5,5) * t395 + Ifges(6,5) * t403 + t456 + t444 / 0.2e1 + t445 / 0.2e1 - t455 + t42 * mrSges(6,1) - t43 * mrSges(6,2) - t129 / 0.2e1 - t447 / 0.2e1 - t446 / 0.2e1 - t281 - mrSges(4,1) * t297 + t457 / 0.2e1) * t147 + (t290 + Ifges(4,1) * t389 + t130 / 0.2e1 + t175 / 0.2e1 + t448 / 0.2e1 - t280 - mrSges(4,2) * t297) * t146 + (Ifges(3,1) * t196 - t284 * t443 + Ifges(3,4) * t195 / 0.2e1 + Ifges(3,5) * qJDD(2)) * t243 + t151 * (-Ifges(5,4) * t262 - Ifges(5,2) * t263) / 0.2e1 + pkin(3) * t9 * t339 + t45 * t344 + t46 * t345 + t37 * t410 + t38 * t412 + t338 * t414 - t135 * t419 - t134 * t420 + (Ifges(6,1) * t37 + Ifges(6,4) * t38) * t403 + t200 * (mrSges(5,1) * t263 - mrSges(5,2) * t262) + t132 * (-mrSges(6,1) * t38 + mrSges(6,2) * t37) + t10 * t78 + t11 * t79 + t66 * t12 + t67 * t13 + t425 * t89 + (-mrSges(4,1) * t383 - mrSges(4,3) * t187 - Ifges(4,4) * t114 + Ifges(5,5) * t409 + Ifges(6,5) * t416 - Ifges(4,2) * t115 - Ifges(4,6) * t236 + Ifges(5,6) * t408 + Ifges(6,6) * t415 + Ifges(5,3) * t399 + Ifges(6,3) * t400 + t424 + t426) * t191 + t138 * t296 + t116 * t288 + (-mrSges(4,2) * t383 - t186 * mrSges(4,3) + Ifges(4,1) * t114 + Ifges(4,4) * t115 + Ifges(4,5) * t236 + t168 * t275 + t269 * t399 + t271 * t408 + t273 * t409 + t285 * t86) * t193 + (t361 * t284 + t430 / 0.2e1) * t248 + (Ifges(3,4) * t196 + t430) * t386 + t196 * t361 / 0.2e1 + (-Ifges(5,1) * t262 - Ifges(5,4) * t263) * t395 - (-mrSges(4,1) * t115 + mrSges(4,2) * t114) * t235 - t117 * t289 - t282 * t307 - t178 * t310 / 0.2e1 + Ifges(2,3) * qJDD(1) - t23 * t339 / 0.2e1 + m(6) * (t10 * t43 + t11 * t42 + t3 * t67 + t4 * t66 + (t132 * t263 + t339 * t60) * pkin(3)); (-m(5) * (t235 + t313) - m(6) * (t235 + t317) + t422 + t452) * g(3) + t423 - (-Ifges(3,2) * t312 + t179 + t228) * t311 / 0.2e1 + (-t261 / 0.2e1 + t282) * qJD(1) ^ 2 + (m(4) * t381 + mrSges(3,1) * t243 + mrSges(4,1) * t232 + mrSges(3,2) * t248 + t370) * (g(1) * t249 + g(2) * t244) + t449 * t79 + (-g(1) * t213 - g(2) * t211 + t437 * t132 - t155 * t4 - t156 * t3 + t203 * t60 + t449 * t42 + t450 * t43) * m(6) + t450 * t78 + (mrSges(4,1) * t236 - mrSges(4,3) * t114) * t380 + (-mrSges(4,2) * t236 + mrSges(4,3) * t115) * t382 - t155 * t12 - t156 * t13 - t425 * t127 + (t246 * t116 - t241 * t117 + t157) * t293 + t181 * t281 + t437 * t53 + t438 * t294 + (m(5) * t254 - t116 * t309 - t117 * t308 + t433) * t224 + (-m(5) * t211 + t244 * t436 + t434) * g(2) + (-m(5) * t213 + t249 * t436 + t435) * g(1) + m(5) * (t168 * t226 + (t200 * t242 + t247 * t267) * t355) - t308 * t366 - t309 * t365 + m(4) * (t186 * t247 + t187 * t242) * pkin(1) + Ifges(3,6) * t195 + Ifges(3,5) * t196 + t203 * t9 + t226 * t30 + Ifges(3,3) * qJDD(2) - t138 * t298 - t270 * t307 / 0.2e1 + t178 * t312 / 0.2e1; t441 * t78 + ((-pkin(5) * t117 - t366) * t246 + (-pkin(5) * t116 - t365 + t401) * t241) * qJD(4) - t153 * t53 - t118 * t117 - t119 * t116 + t440 * t79 + (-t213 * t462 + t249 * t253 + t435) * g(1) + (-t211 * t462 + t244 * t253 + t434) * g(2) - t168 * t418 + (t371 + t422) * g(3) + (-t157 * t247 + (t367 - t438) * t242) * t356 + t433 * pkin(5) - t184 * t12 - t185 * t13 - t225 * t9 - pkin(2) * t30 + (-t317 * g(3) - t184 * t4 - t185 * t3 - t225 * t60 + t441 * t43 + t440 * t42 + (-t153 + t292) * t132) * m(6) + (pkin(5) * t254 - g(3) * t313 - t100 * t118 - t101 * t119 - t200 * t229) * m(5) + t423; -(mrSges(6,1) * t132 + Ifges(6,4) * t404 + Ifges(6,2) * t406 + Ifges(6,6) * t394 - t385 + t413) * t99 + (-mrSges(6,2) * t132 + mrSges(6,3) * t42 + Ifges(6,1) * t404 + Ifges(6,4) * t406 + Ifges(6,5) * t394 + t411) * t283 + t424 + (t364 - t116) * t100 + (mrSges(5,2) * t174 - t173 * t442 - t322) * g(1) + t301 + (t275 - t274 + t305) * t376 + (t363 + t117) * t101 + (t240 * t3 + t245 * t4 + (-t240 * t42 + t245 * t43) * qJD(5)) * t417 + (Ifges(5,5) * t151 - Ifges(5,6) * t152) * t392 + t85 * t395 + (Ifges(5,1) * t151 - t359) * t396 + (-Ifges(5,2) * t152 + t150 + t86) * t397 - t52 * t78 - t51 * t79 + t264 + (-mrSges(5,2) * t172 + t171 * t442 - t323) * g(2) - t200 * (mrSges(5,1) * t152 + mrSges(5,2) * t151) - t53 * t379 - m(6) * (t132 * t379 + t42 * t51 + t43 * t52) + ((-t240 * t79 + t245 * t78) * qJD(5) + t12 * t245 + t13 * t240) * pkin(3); -t132 * (mrSges(6,1) * t99 + mrSges(6,2) * t283) + (Ifges(6,1) * t283 - t384) * t404 + t40 * t403 + (Ifges(6,5) * t283 - Ifges(6,6) * t99) * t394 - t42 * t78 + t43 * t79 - g(1) * t322 - g(2) * t323 - t274 * t376 + (t283 * t42 + t43 * t99) * mrSges(6,3) + t264 + (-Ifges(6,2) * t99 + t41 + t93) * t406;];
tau = t1;
