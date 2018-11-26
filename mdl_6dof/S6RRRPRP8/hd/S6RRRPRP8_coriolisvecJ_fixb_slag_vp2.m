% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2018-11-23 17:46
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:46:25
% EndTime: 2018-11-23 17:46:44
% DurationCPUTime: 18.62s
% Computational Cost: add. (7814->680), mult. (18667->872), div. (0->0), fcn. (11864->6), ass. (0->311)
t262 = sin(qJ(3));
t263 = sin(qJ(2));
t338 = qJD(1) * t263;
t318 = t262 * t338;
t265 = cos(qJ(3));
t329 = t265 * qJD(2);
t209 = t318 - t329;
t316 = t265 * t338;
t210 = qJD(2) * t262 + t316;
t261 = sin(qJ(5));
t264 = cos(qJ(5));
t135 = t209 * t264 - t210 * t261;
t450 = Ifges(6,1) + Ifges(7,1);
t281 = t209 * t261 + t210 * t264;
t427 = Ifges(6,4) + Ifges(7,4);
t461 = t427 * t281;
t470 = t450 * t135 - t461;
t448 = Ifges(6,2) + Ifges(7,2);
t464 = t427 * t135;
t469 = t281 * t448 - t464;
t426 = Ifges(6,5) + Ifges(7,5);
t425 = Ifges(6,6) + Ifges(7,6);
t405 = pkin(8) - pkin(9);
t232 = t405 * t265;
t377 = pkin(7) * t262;
t319 = -pkin(3) - t377;
t278 = (-pkin(4) + t319) * t263;
t266 = cos(qJ(2));
t349 = t265 * t266;
t272 = -pkin(9) * t349 + t278;
t304 = pkin(2) * t263 - pkin(8) * t266;
t219 = t304 * qJD(1);
t352 = t219 * t265;
t468 = -qJD(1) * t272 + qJD(3) * t232 + t352;
t190 = t262 * t219;
t247 = qJ(4) * t338;
t334 = qJD(3) * t262;
t350 = t263 * t265;
t351 = t262 * t266;
t467 = t190 + t247 + (-pkin(7) * t350 + pkin(9) * t351) * qJD(1) + t405 * t334;
t226 = -pkin(2) * t266 - pkin(8) * t263 - pkin(1);
t194 = t226 * qJD(1);
t337 = qJD(1) * t266;
t253 = pkin(7) * t337;
t230 = qJD(2) * pkin(8) + t253;
t139 = t265 * t194 - t262 * t230;
t413 = qJD(4) - t139;
t328 = qJD(1) * qJD(2);
t309 = t263 * t328;
t333 = qJD(3) * t263;
t314 = t262 * t333;
t315 = t266 * t329;
t327 = qJD(2) * qJD(3);
t162 = t265 * t327 + (-t314 + t315) * qJD(1);
t332 = qJD(3) * t265;
t335 = qJD(2) * t266;
t274 = t262 * t335 + t263 * t332;
t163 = qJD(1) * t274 + t262 * t327;
t47 = qJD(5) * t135 + t162 * t264 + t163 * t261;
t48 = -qJD(5) * t281 - t162 * t261 + t163 * t264;
t466 = -t448 * t48 / 0.2e1 - t427 * t47 / 0.2e1 + t425 * t309 / 0.2e1;
t465 = qJD(2) / 0.2e1;
t457 = Ifges(4,1) + Ifges(5,1);
t449 = Ifges(4,5) + Ifges(5,4);
t279 = t261 * t262 + t264 * t265;
t412 = qJD(3) - qJD(5);
t142 = t412 * t279;
t275 = t279 * t266;
t175 = qJD(1) * t275;
t346 = t142 - t175;
t280 = t261 * t265 - t262 * t264;
t143 = t412 * t280;
t276 = t280 * t266;
t174 = qJD(1) * t276;
t463 = t143 - t174;
t462 = -t262 * qJD(4) - t253 + (t265 * t337 - t332) * qJ(4);
t460 = -pkin(9) * t210 + t413;
t244 = qJD(3) - t337;
t267 = -pkin(3) - pkin(4);
t73 = t244 * t267 + t460;
t140 = t262 * t194 + t265 * t230;
t102 = pkin(9) * t209 + t140;
t234 = t244 * qJ(4);
t86 = t102 + t234;
t21 = -t261 * t86 + t264 * t73;
t416 = qJ(6) * t281;
t17 = t21 - t416;
t237 = qJD(5) - t244;
t16 = pkin(5) * t237 + t17;
t22 = t261 * t73 + t264 * t86;
t439 = qJ(6) * t135;
t18 = t22 + t439;
t373 = Ifges(6,3) + Ifges(7,3);
t385 = -t237 / 0.2e1;
t43 = Ifges(7,6) * t48;
t44 = Ifges(6,6) * t48;
t45 = Ifges(7,5) * t47;
t46 = Ifges(6,5) * t47;
t229 = -qJD(2) * pkin(2) + pkin(7) * t338;
t119 = t209 * pkin(3) - t210 * qJ(4) + t229;
t92 = -pkin(4) * t209 - t119;
t51 = -pkin(5) * t135 + qJD(6) + t92;
t459 = (t135 * t426 - t281 * t425) * t385 + (t135 * t21 + t22 * t281) * mrSges(6,3) + (t135 * t16 + t18 * t281) * mrSges(7,3) - t373 * t309 - t51 * (mrSges(7,1) * t281 + mrSges(7,2) * t135) - t92 * (mrSges(6,1) * t281 + mrSges(6,2) * t135) + t43 + t44 + t45 + t46;
t313 = Ifges(3,5) * t465;
t456 = Ifges(4,6) - Ifges(5,6);
t455 = -t426 * t309 + t427 * t48 + t450 * t47;
t445 = t135 * t448 + t237 * t425 + t461;
t444 = t426 * t237 + t281 * t450 + t464;
t231 = t405 * t262;
t330 = qJD(5) * t264;
t331 = qJD(5) * t261;
t442 = t231 * t330 - t232 * t331 + t261 * t468 - t264 * t467;
t159 = t261 * t231 + t264 * t232;
t441 = -qJD(5) * t159 + t261 * t467 + t264 * t468;
t317 = t262 * t337;
t320 = t267 * t262;
t438 = qJD(3) * t320 - t267 * t317 - t462;
t400 = t135 / 0.2e1;
t451 = -qJD(2) / 0.2e1;
t447 = -qJ(6) * t463 - qJD(6) * t279 + t442;
t446 = pkin(5) * t338 - qJ(6) * t346 + qJD(6) * t280 + t441;
t443 = t449 * t309 + (-Ifges(4,4) + Ifges(5,5)) * t163 + t457 * t162;
t440 = pkin(5) * t463 + t438;
t224 = t264 * qJ(4) + t261 * t267;
t418 = -qJD(5) * t224 - t264 * t102 - t261 * t460;
t223 = -qJ(4) * t261 + t264 * t267;
t417 = qJD(5) * t223 - t261 * t102 + t264 * t460;
t437 = (-t317 + t334) * pkin(3) + t462;
t361 = Ifges(5,5) * t265;
t367 = Ifges(4,4) * t265;
t436 = t262 * t457 - t361 + t367;
t362 = Ifges(5,5) * t262;
t368 = Ifges(4,4) * t262;
t435 = t265 * t457 + t362 - t368;
t250 = Ifges(3,4) * t337;
t202 = Ifges(5,5) * t210;
t120 = t244 * Ifges(5,6) + t209 * Ifges(5,3) + t202;
t357 = t210 * Ifges(4,4);
t123 = -t209 * Ifges(4,2) + t244 * Ifges(4,6) + t357;
t282 = t139 * t265 + t140 * t262;
t111 = -pkin(3) * t244 + t413;
t113 = t234 + t140;
t283 = t111 * t265 - t113 * t262;
t289 = Ifges(5,3) * t262 + t361;
t293 = -Ifges(4,2) * t262 + t367;
t298 = mrSges(5,1) * t262 - mrSges(5,3) * t265;
t300 = mrSges(4,1) * t262 + mrSges(4,2) * t265;
t359 = Ifges(5,6) * t262;
t360 = Ifges(4,6) * t262;
t363 = Ifges(4,5) * t265;
t366 = Ifges(5,4) * t265;
t378 = t265 / 0.2e1;
t380 = t262 / 0.2e1;
t381 = -t262 / 0.2e1;
t388 = t210 / 0.2e1;
t390 = t209 / 0.2e1;
t391 = -t209 / 0.2e1;
t203 = Ifges(4,4) * t209;
t358 = t209 * Ifges(5,5);
t415 = t210 * t457 + t244 * t449 - t203 + t358;
t428 = t244 / 0.2e1;
t268 = t283 * mrSges(5,2) - t282 * mrSges(4,3) + t119 * t298 + t120 * t380 + t123 * t381 + t229 * t300 + t289 * t390 + t293 * t391 + t435 * t388 + (t359 + t366 - t360 + t363) * t428 + t415 * t378;
t434 = t268 + Ifges(3,1) * t338 / 0.2e1 + t250 / 0.2e1 + t313;
t401 = -t135 / 0.2e1;
t398 = -t281 / 0.2e1;
t312 = Ifges(3,6) * t451;
t420 = -t439 + t418;
t419 = -t416 + t417;
t184 = t280 * t263;
t245 = pkin(7) * t351;
t259 = t266 * pkin(3);
t376 = pkin(9) * t263;
t126 = pkin(4) * t266 + t245 + t259 + (-t226 - t376) * t265;
t246 = pkin(7) * t349;
t173 = t262 * t226 + t246;
t160 = -qJ(4) * t266 + t173;
t138 = t262 * t376 + t160;
t64 = t261 * t126 + t264 * t138;
t414 = t262 * t449 + t265 * t456;
t322 = Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1;
t323 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t324 = Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t325 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t326 = Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t369 = Ifges(3,4) * t263;
t410 = -t323 * t135 + t324 * t209 + t326 * t210 - (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t237 + t322 * t244 + t325 * t281 + t113 * mrSges(5,3) + t139 * mrSges(4,1) + t18 * mrSges(7,2) + t22 * mrSges(6,2) + Ifges(4,6) * t391 + Ifges(5,6) * t390 + t312 - (t266 * Ifges(3,2) + t369) * qJD(1) / 0.2e1 - t111 * mrSges(5,1) - t140 * mrSges(4,2) - t16 * mrSges(7,1) - t21 * mrSges(6,1) + (Ifges(4,3) + Ifges(5,2)) * t428 + t425 * t401 + t426 * t398 + t449 * t388 + t373 * t385;
t409 = t47 / 0.2e1;
t408 = t48 / 0.2e1;
t404 = pkin(1) * mrSges(3,1);
t403 = pkin(1) * mrSges(3,2);
t397 = t281 / 0.2e1;
t396 = t162 / 0.2e1;
t395 = -t163 / 0.2e1;
t394 = t163 / 0.2e1;
t389 = -t210 / 0.2e1;
t384 = t237 / 0.2e1;
t383 = -t244 / 0.2e1;
t379 = -t265 / 0.2e1;
t36 = mrSges(7,2) * t309 + mrSges(7,3) * t48;
t37 = mrSges(6,2) * t309 + mrSges(6,3) * t48;
t372 = t36 + t37;
t371 = mrSges(4,3) * t209;
t370 = mrSges(4,3) * t210;
t353 = qJD(2) * mrSges(3,2);
t256 = t262 * qJ(4);
t166 = -mrSges(4,2) * t244 - t371;
t169 = -mrSges(5,2) * t209 + mrSges(5,3) * t244;
t344 = -t166 - t169;
t167 = mrSges(4,1) * t244 - t370;
t168 = -mrSges(5,1) * t244 + mrSges(5,2) * t210;
t343 = -t167 + t168;
t222 = t304 * qJD(2);
t342 = t262 * t222 + t226 * t332;
t146 = t210 * pkin(3) + t209 * qJ(4);
t340 = qJ(4) * t315 + qJD(4) * t350;
t336 = qJD(2) * t263;
t225 = -t265 * pkin(3) - pkin(2) - t256;
t14 = -t48 * mrSges(7,1) + t47 * mrSges(7,2);
t63 = t264 * t126 - t138 * t261;
t158 = t264 * t231 - t232 * t261;
t172 = t226 * t265 - t245;
t201 = t265 * pkin(4) - t225;
t308 = pkin(7) * t309;
t307 = -pkin(7) + t320;
t306 = m(4) * t229 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t209 + mrSges(4,2) * t210 + mrSges(3,3) * t338;
t108 = -pkin(4) * t210 - t146;
t305 = t319 * t263;
t196 = qJD(1) * t222;
t303 = t194 * t334 - t196 * t265 + t230 * t332;
t302 = qJD(3) * t246 - t222 * t265 + t226 * t334;
t301 = mrSges(4,1) * t265 - mrSges(4,2) * t262;
t299 = mrSges(5,1) * t265 + mrSges(5,3) * t262;
t292 = Ifges(4,2) * t265 + t368;
t288 = -Ifges(5,3) * t265 + t362;
t75 = t194 * t332 + t262 * t196 - t230 * t334 - t265 * t308;
t59 = qJ(4) * t309 + t244 * qJD(4) + t75;
t277 = qJD(2) * t305;
t70 = qJD(1) * t277 + t303;
t285 = t262 * t70 + t265 * t59;
t76 = t262 * t308 - t303;
t284 = -t262 * t76 + t265 * t75;
t165 = -pkin(7) * t316 + t190;
t132 = -mrSges(5,1) * t309 + t162 * mrSges(5,2);
t31 = pkin(9) * t163 + t59;
t32 = -pkin(9) * t162 + t278 * t328 + t303;
t3 = t261 * t32 + t264 * t31 + t73 * t330 - t331 * t86;
t65 = pkin(9) * t314 + qJD(2) * t272 + t302;
t249 = qJ(4) * t336;
t67 = t249 + (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t350 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t262) * t266 + t342;
t7 = t126 * t330 - t138 * t331 + t261 * t65 + t264 * t67;
t243 = qJ(4) * t350;
t157 = t263 * t307 + t243;
t66 = t163 * pkin(3) - t162 * qJ(4) + qJD(2) * t253 - t210 * qJD(4);
t4 = -qJD(5) * t22 - t261 * t31 + t264 * t32;
t1 = -pkin(5) * t309 - qJ(6) * t47 - qJD(6) * t281 + t4;
t2 = qJ(6) * t48 + qJD(6) * t135 + t3;
t273 = t4 * mrSges(6,1) + t1 * mrSges(7,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t33 = -pkin(4) * t163 - t66;
t8 = -qJD(5) * t64 - t261 * t67 + t264 * t65;
t97 = (-t263 * t329 - t266 * t334) * pkin(7) + t342;
t271 = -t76 * mrSges(4,1) + t70 * mrSges(5,1) + t75 * mrSges(4,2) - t59 * mrSges(5,3) + t273;
t74 = (t265 * t267 - t256) * t333 + t307 * t335 + t340;
t242 = Ifges(5,2) * t309;
t241 = Ifges(4,3) * t309;
t228 = mrSges(3,3) * t337 - t353;
t218 = -pkin(5) + t223;
t185 = t279 * t263;
t178 = -t243 + (pkin(3) * t262 + pkin(7)) * t263;
t164 = pkin(7) * t318 + t352;
t161 = -t172 + t259;
t155 = Ifges(5,4) * t162;
t154 = Ifges(4,5) * t162;
t153 = Ifges(4,6) * t163;
t152 = Ifges(5,6) * t163;
t149 = pkin(5) * t279 + t201;
t147 = mrSges(5,1) * t209 - mrSges(5,3) * t210;
t145 = qJD(1) * t305 - t352;
t144 = t165 + t247;
t133 = -mrSges(4,2) * t309 - mrSges(4,3) * t163;
t131 = mrSges(4,1) * t309 - mrSges(4,3) * t162;
t130 = -mrSges(5,2) * t163 + mrSges(5,3) * t309;
t115 = -qJ(6) * t279 + t159;
t114 = qJ(6) * t280 + t158;
t107 = pkin(5) * t184 + t157;
t106 = mrSges(6,1) * t237 - mrSges(6,3) * t281;
t105 = mrSges(7,1) * t237 - mrSges(7,3) * t281;
t104 = -mrSges(6,2) * t237 + mrSges(6,3) * t135;
t103 = -mrSges(7,2) * t237 + mrSges(7,3) * t135;
t98 = t336 * t377 - t302;
t96 = pkin(3) * t274 + pkin(7) * t335 + qJ(4) * t314 - t340;
t91 = t277 + t302;
t90 = mrSges(4,1) * t163 + mrSges(4,2) * t162;
t89 = mrSges(5,1) * t163 - mrSges(5,3) * t162;
t85 = -qJD(4) * t266 + t249 + t97;
t80 = t162 * Ifges(4,4) - t163 * Ifges(4,2) + Ifges(4,6) * t309;
t79 = t162 * Ifges(5,5) + Ifges(5,6) * t309 + t163 * Ifges(5,3);
t78 = qJD(2) * t275 + t184 * t412;
t77 = -qJD(2) * t276 + t142 * t263;
t69 = -mrSges(6,1) * t135 + mrSges(6,2) * t281;
t68 = -mrSges(7,1) * t135 + mrSges(7,2) * t281;
t58 = -pkin(5) * t281 + t108;
t41 = -qJ(6) * t184 + t64;
t38 = pkin(5) * t266 - qJ(6) * t185 + t63;
t35 = -mrSges(6,1) * t309 - mrSges(6,3) * t47;
t34 = -mrSges(7,1) * t309 - mrSges(7,3) * t47;
t23 = -pkin(5) * t77 + t74;
t15 = -mrSges(6,1) * t48 + mrSges(6,2) * t47;
t13 = -pkin(5) * t48 + t33;
t6 = qJ(6) * t77 - qJD(6) * t184 + t7;
t5 = -pkin(5) * t336 - qJ(6) * t78 - qJD(6) * t185 + t8;
t9 = [(-t325 * t47 + ((0.3e1 / 0.2e1 * Ifges(3,4) * t266 - 0.2e1 * t403) * qJD(1) + t306 * pkin(7) + t313 + t434) * qJD(2) - t324 * t163 - t326 * t162 + t46 / 0.2e1 + t45 / 0.2e1 + t44 / 0.2e1 + t43 / 0.2e1 + t271 - t241 / 0.2e1 - t242 / 0.2e1 + t323 * t48 - t152 / 0.2e1 + t153 / 0.2e1 - t154 / 0.2e1 - t155 / 0.2e1) * t266 + (t289 * t394 + t293 * t395 + pkin(7) * t90 + t79 * t380 + t80 * t381 + t66 * t298 + (-t262 * t75 - t265 * t76) * mrSges(4,3) + (-t262 * t59 + t265 * t70) * mrSges(5,2) + (t123 * t379 + t119 * t299 + t288 * t391 + t292 * t390 + t229 * t301 + (t139 * t262 - t140 * t265) * mrSges(4,3) + (-t111 * t262 - t113 * t265) * mrSges(5,2) + t436 * t389 + t414 * t383 + t415 * t381) * qJD(3) + (((-0.3e1 / 0.2e1 * Ifges(3,4) + t363 / 0.2e1 - t360 / 0.2e1 + t366 / 0.2e1 + t359 / 0.2e1) * t263 - 0.2e1 * t404 + t325 * t185 + t323 * t184 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + (m(4) * pkin(7) + t300) * pkin(7) - t322 - t373) * t266) * qJD(1) - pkin(7) * t228 + t312 + t410) * qJD(2) + t435 * t396 + (qJD(3) * t120 + t443) * t378) * t263 + t184 * t466 + m(7) * (t1 * t38 + t107 * t13 + t16 * t5 + t18 * t6 + t2 * t41 + t23 * t51) + m(6) * (t157 * t33 + t21 * t8 + t22 * t7 + t3 * t64 + t4 * t63 + t74 * t92) + m(5) * (t111 * t91 + t113 * t85 + t119 * t96 + t160 * t59 + t161 * t70 + t178 * t66) + t6 * t103 + t7 * t104 + t5 * t105 + t8 * t106 + t107 * t14 + m(4) * (t139 * t98 + t140 * t97 + t172 * t76 + t173 * t75) + (t425 * t77 + t426 * t78) * t384 + (-t184 * t448 + t185 * t427) * t408 + (t427 * t78 + t448 * t77) * t400 + (t427 * t77 + t450 * t78) * t397 + (-t184 * t427 + t185 * t450) * t409 + t92 * (-mrSges(6,1) * t77 + mrSges(6,2) * t78) + t51 * (-mrSges(7,1) * t77 + mrSges(7,2) * t78) + t74 * t69 + t23 * t68 + t63 * t35 + t64 * t37 + t38 * t34 + t41 * t36 + t33 * (mrSges(6,1) * t184 + mrSges(6,2) * t185) + t13 * (mrSges(7,1) * t184 + mrSges(7,2) * t185) + (-t1 * t185 - t16 * t78 + t18 * t77 - t184 * t2) * mrSges(7,3) + (-t184 * t3 - t185 * t4 - t21 * t78 + t22 * t77) * mrSges(6,3) + t444 * t78 / 0.2e1 + t445 * t77 / 0.2e1 + t455 * t185 / 0.2e1 + t96 * t147 + t157 * t15 + t160 * t130 + t161 * t132 + t97 * t166 + t98 * t167 + t91 * t168 + t85 * t169 + t172 * t131 + t173 * t133 + t178 * t89; t440 * t68 + t441 * t106 + t442 * t104 + (t158 * t4 + t159 * t3 + t201 * t33 + t21 * t441 + t22 * t442 + t438 * t92) * m(6) + t443 * t380 + t444 * (t142 / 0.2e1 - t175 / 0.2e1) + t445 * (-t143 / 0.2e1 + t174 / 0.2e1) + t446 * t105 + (t1 * t114 + t115 * t2 + t13 * t149 + t16 * t446 + t18 * t447 + t440 * t51) * m(7) + t447 * t103 + (-t21 * t346 - t22 * t463 - t279 * t3 + t280 * t4) * mrSges(6,3) + (t1 * t280 - t16 * t346 - t18 * t463 - t2 * t279) * mrSges(7,3) + (mrSges(6,1) * t463 + mrSges(6,2) * t346) * t92 + (mrSges(7,1) * t463 + mrSges(7,2) * t346) * t51 + (t142 * t427 - t143 * t448) * t400 + (-t174 * t448 + t175 * t427) * t401 + (-t279 * t448 - t280 * t427) * t408 + (t142 * t450 - t143 * t427) * t397 + (-t174 * t427 + t175 * t450) * t398 + (-t279 * t427 - t280 * t450) * t409 + t114 * t34 + t115 * t36 + t284 * mrSges(4,3) + t285 * mrSges(5,2) + t279 * t466 + (m(5) * t285 + m(4) * t284 + (t130 + t133) * t265 + (-t131 + t132) * t262 + (-m(4) * t282 + m(5) * t283 + t262 * t344 + t265 * t343) * qJD(3)) * pkin(8) + (t142 * t426 - t143 * t425) * t384 + (-t174 * t425 + t175 * t426) * t385 + t268 * qJD(3) + t33 * (mrSges(6,1) * t279 - mrSges(6,2) * t280) + t13 * (mrSges(7,1) * t279 - mrSges(7,2) * t280) + t436 * t396 + (-t111 * t145 - t113 * t144 + t119 * t437 + t225 * t66) * m(5) + t437 * t147 + t438 * t69 - pkin(2) * t90 - m(4) * (t139 * t164 + t140 * t165) + ((t312 + (t228 + t353) * pkin(7) + (t369 / 0.2e1 + t404) * qJD(1) + (-t279 * t425 - t280 * t426) * t451 + t414 * t465 - t410) * t263 + (((-m(4) * pkin(2) - mrSges(3,1) - t301) * qJD(2) - t306) * pkin(7) + (t403 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t263) * qJD(1) - t250 / 0.2e1 + t313 - t434) * t266) * qJD(1) - t455 * t280 / 0.2e1 + t201 * t15 + t225 * t89 + t149 * t14 + t158 * t35 + t159 * t37 - t165 * t166 - t164 * t167 - t145 * t168 - t144 * t169 + t80 * t378 + t79 * t379 + t288 * t394 + t292 * t395 - t66 * t299; qJ(4) * t130 - pkin(3) * t132 + (t111 * t209 + t113 * t210) * mrSges(5,2) - t459 - t108 * t69 + (-pkin(3) * t70 + qJ(4) * t59 - t111 * t140 + t113 * t413 - t119 * t146) * m(5) + (-t209 * t457 + t120 + t202 - t357) * t389 - t271 + (-Ifges(4,2) * t210 - t203 + t415) * t390 + t241 + t242 - t58 * t68 + (-t449 * t209 - t210 * t456) * t383 + t417 * t104 + t418 * t106 + (-t108 * t92 + t21 * t418 + t22 * t417 + t223 * t4 + t224 * t3) * m(6) + t419 * t103 + t420 * t105 + (t1 * t218 + t16 * t420 + t18 * t419 + t2 * t224 - t51 * t58) * m(7) + t469 * t401 + (t445 - t470) * t398 - t119 * (mrSges(5,1) * t210 + mrSges(5,3) * t209) + t152 - t153 + t154 + t155 + t218 * t34 + t223 * t35 - t229 * (mrSges(4,1) * t210 - mrSges(4,2) * t209) - t146 * t147 + qJD(4) * t169 + t444 * t400 + t123 * t388 + (Ifges(5,3) * t210 - t358) * t391 + (-t343 + t370) * t140 + (t344 - t371) * t139 + t372 * t224; -t244 * t169 + (t147 - t68 - t69) * t210 + (t34 + t35 + t237 * (t103 + t104)) * t264 + (t372 - t237 * (t105 + t106)) * t261 + t132 + (t1 * t264 + t2 * t261 - t210 * t51 + t237 * (-t16 * t261 + t18 * t264)) * m(7) + (-t210 * t92 + t261 * t3 + t264 * t4 + t237 * (-t21 * t261 + t22 * t264)) * m(6) + (-t113 * t244 + t119 * t210 + t70) * m(5); t273 - t17 * t103 - t21 * t104 + t18 * t105 + t22 * t106 + (-t281 * t68 + t34) * pkin(5) + (-(-t16 + t17) * t18 + (-t281 * t51 + t1) * pkin(5)) * m(7) + t470 * t398 + t445 * t397 + (t444 - t469) * t401 + t459; -t135 * t103 + t281 * t105 + 0.2e1 * (t13 / 0.2e1 + t18 * t401 + t16 * t397) * m(7) + t14;];
tauc  = t9(:);
