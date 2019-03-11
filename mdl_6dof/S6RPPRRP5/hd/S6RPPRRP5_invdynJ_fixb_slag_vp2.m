% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:42
% EndTime: 2019-03-09 02:08:01
% DurationCPUTime: 12.10s
% Computational Cost: add. (3309->498), mult. (6124->625), div. (0->0), fcn. (3189->6), ass. (0->224)
t363 = Ifges(6,4) + Ifges(7,4);
t364 = Ifges(6,1) + Ifges(7,1);
t342 = Ifges(6,5) + Ifges(7,5);
t362 = Ifges(6,2) + Ifges(7,2);
t341 = Ifges(6,6) + Ifges(7,6);
t134 = -qJ(6) - pkin(8);
t372 = m(7) * t134 - mrSges(6,3) - mrSges(7,3);
t140 = cos(qJ(5));
t371 = t363 * t140;
t137 = sin(qJ(5));
t370 = t363 * t137;
t138 = sin(qJ(4));
t141 = cos(qJ(4));
t237 = qJD(1) * qJD(4);
t103 = qJDD(1) * t141 - t138 * t237;
t244 = qJD(4) * t140;
t249 = qJD(1) * t141;
t97 = -t137 * t249 + t244;
t42 = qJD(5) * t97 + qJDD(4) * t137 + t103 * t140;
t305 = t42 / 0.2e1;
t246 = qJD(4) * t137;
t98 = t140 * t249 + t246;
t43 = -qJD(5) * t98 + qJDD(4) * t140 - t103 * t137;
t304 = t43 / 0.2e1;
t104 = -qJDD(1) * t138 - t141 * t237;
t93 = qJDD(5) - t104;
t303 = t93 / 0.2e1;
t369 = -t97 / 0.2e1;
t368 = -t98 / 0.2e1;
t250 = qJD(1) * t138;
t119 = qJD(5) + t250;
t367 = -t119 / 0.2e1;
t122 = pkin(5) * t140 + pkin(4);
t366 = m(7) * t122;
t365 = qJD(4) / 0.2e1;
t361 = Ifges(6,3) + Ifges(7,3);
t183 = -mrSges(7,1) * t140 + mrSges(7,2) * t137;
t185 = -mrSges(6,1) * t140 + mrSges(6,2) * t137;
t360 = m(6) * pkin(4) - t183 - t185 + t366;
t325 = -t341 * t137 + t342 * t140;
t323 = -t362 * t137 + t371;
t321 = t364 * t140 - t370;
t359 = -m(6) * pkin(8) + t372;
t358 = t103 / 0.2e1;
t357 = t104 / 0.2e1;
t356 = t342 * t303 + t363 * t304 + t305 * t364;
t355 = t363 * t97;
t272 = Ifges(5,4) * t138;
t181 = t141 * Ifges(5,1) - t272;
t337 = t342 * t119 + t364 * t98 + t355;
t354 = Ifges(5,5) * qJD(4) + qJD(1) * t181 + t140 * t337;
t332 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t97 - mrSges(6,2) * t98 - mrSges(5,3) * t249;
t240 = qJD(5) * t140;
t241 = qJD(5) * t137;
t136 = pkin(1) + qJ(3);
t238 = qJD(1) * qJD(3);
t100 = qJDD(1) * t136 - qJDD(2) + t238;
t28 = -pkin(4) * t104 - pkin(8) * t103 + t100;
t123 = qJD(1) * qJ(2) + qJD(3);
t115 = -qJD(1) * pkin(7) + t123;
t243 = qJD(4) * t141;
t130 = qJD(1) * qJD(2);
t330 = qJDD(1) * qJ(2) + t130;
t113 = qJDD(3) + t330;
t99 = -qJDD(1) * pkin(7) + t113;
t54 = t115 * t243 + t138 * t99;
t46 = qJDD(4) * pkin(8) + t54;
t190 = pkin(4) * t138 - pkin(8) * t141;
t106 = t190 + t136;
t65 = qJD(1) * t106 - qJD(2);
t105 = t138 * t115;
t76 = qJD(4) * pkin(8) + t105;
t3 = t137 * t28 + t140 * t46 + t65 * t240 - t241 * t76;
t25 = t137 * t65 + t140 * t76;
t4 = -qJD(5) * t25 - t137 * t46 + t140 * t28;
t189 = -t137 * t4 + t140 * t3;
t24 = -t137 * t76 + t140 * t65;
t353 = -t24 * t240 - t25 * t241 + t189;
t135 = qJ(2) - pkin(7);
t317 = qJD(2) * t138 + t135 * t243;
t352 = t363 * t98;
t139 = sin(qJ(1));
t142 = cos(qJ(1));
t315 = g(1) * t142 + g(2) * t139;
t271 = Ifges(5,4) * t141;
t176 = -t138 * Ifges(5,2) + t271;
t351 = Ifges(5,6) * t365 + qJD(1) * t176 / 0.2e1 + t342 * t368 + t341 * t369 + t361 * t367;
t2 = qJ(6) * t43 + qJD(6) * t97 + t3;
t350 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t306 = m(7) * pkin(5);
t298 = -m(5) - m(4);
t349 = -m(6) - m(7);
t348 = t341 * t93 + t362 * t43 + t363 * t42;
t346 = -mrSges(7,1) - mrSges(6,1);
t345 = mrSges(3,2) - mrSges(4,3);
t344 = -mrSges(4,2) - mrSges(3,3);
t343 = mrSges(6,2) + mrSges(7,2);
t12 = -mrSges(6,1) * t43 + mrSges(6,2) * t42;
t340 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t103 - t12;
t338 = t341 * t119 + t362 * t97 + t352;
t202 = qJD(5) * t134;
t225 = t137 * t250;
t191 = pkin(4) * t141 + pkin(8) * t138;
t102 = t191 * qJD(1);
t253 = t140 * t141;
t50 = t137 * t102 + t115 * t253;
t336 = -qJ(6) * t225 + qJD(6) * t140 + t137 * t202 - t50;
t256 = t138 * t140;
t258 = t137 * t141;
t49 = t140 * t102 - t115 * t258;
t335 = -qJD(6) * t137 + t140 * t202 - (pkin(5) * t141 + qJ(6) * t256) * qJD(1) - t49;
t334 = mrSges(7,1) + t306;
t333 = -t105 + (t225 + t241) * pkin(5);
t108 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t250;
t288 = mrSges(7,3) * t98;
t61 = mrSges(7,1) * t119 - t288;
t290 = mrSges(6,3) * t98;
t62 = mrSges(6,1) * t119 - t290;
t274 = t61 + t62;
t289 = mrSges(7,3) * t97;
t59 = -mrSges(7,2) * t119 + t289;
t291 = mrSges(6,3) * t97;
t60 = -mrSges(6,2) * t119 + t291;
t275 = t59 + t60;
t154 = t137 * t274 - t140 * t275;
t331 = -t108 + t154;
t329 = -t138 * t325 + t141 * t361;
t328 = -t138 * t323 + t141 * t341;
t327 = -t138 * t321 + t141 * t342;
t326 = t137 * t342 + t140 * t341;
t324 = t140 * t362 + t370;
t322 = t137 * t364 + t371;
t319 = t341 * t43 + t342 * t42 + t361 * t93;
t245 = qJD(4) * t138;
t247 = qJD(2) * t141;
t318 = t135 * t245 - t247;
t53 = -t115 * t245 + t141 * t99;
t316 = -t138 * t54 - t141 * t53;
t314 = -m(5) + t349;
t313 = mrSges(6,1) + t334;
t310 = mrSges(2,2) + mrSges(5,3) + t344;
t1 = pkin(5) * t93 - qJ(6) * t42 - qJD(6) * t98 + t4;
t17 = -qJ(6) * t98 + t24;
t14 = pkin(5) * t119 + t17;
t18 = qJ(6) * t97 + t25;
t309 = -t1 * t137 - t14 * t240 + t140 * t2 - t18 * t241;
t186 = t138 * mrSges(5,1) + t141 * mrSges(5,2);
t307 = t138 * t366 + t372 * t141 + mrSges(2,1) + t186 - t345;
t143 = qJD(1) ^ 2;
t301 = t97 / 0.2e1;
t299 = t98 / 0.2e1;
t295 = t119 / 0.2e1;
t19 = mrSges(7,1) * t93 - mrSges(7,3) * t42;
t20 = mrSges(6,1) * t93 - mrSges(6,3) * t42;
t277 = -t19 - t20;
t21 = -mrSges(7,2) * t93 + mrSges(7,3) * t43;
t22 = -mrSges(6,2) * t93 + mrSges(6,3) * t43;
t276 = t21 + t22;
t273 = mrSges(7,2) * t140;
t58 = t137 * t106 + t135 * t256;
t262 = qJDD(1) * pkin(1);
t261 = t115 * t141;
t259 = t137 * t138;
t257 = t137 * t142;
t255 = t139 * t137;
t254 = t139 * t140;
t252 = t142 * t140;
t251 = t142 * pkin(1) + t139 * qJ(2);
t242 = qJD(5) * t135;
t239 = qJD(5) * t141;
t51 = -mrSges(7,1) * t97 + mrSges(7,2) * t98;
t234 = t51 - t332;
t232 = t298 + t349;
t226 = t142 * qJ(3) + t251;
t222 = t138 * t242;
t221 = t140 * t239;
t220 = t137 * t245;
t11 = -t43 * mrSges(7,1) + t42 * mrSges(7,2);
t201 = -t237 / 0.2e1;
t200 = (t138 ^ 2 + t141 ^ 2) * t115;
t94 = qJD(4) * t191 + qJD(3);
t198 = t106 * t240 + t137 * t94 + t140 * t317;
t188 = -t104 * mrSges(5,1) + t103 * mrSges(5,2);
t187 = mrSges(5,1) * t141 - mrSges(5,2) * t138;
t184 = mrSges(6,1) * t137 + mrSges(6,2) * t140;
t182 = mrSges(7,1) * t137 + t273;
t171 = -Ifges(5,5) * t138 - Ifges(5,6) * t141;
t166 = t137 * t25 + t140 * t24;
t77 = -qJD(4) * pkin(4) - t261;
t116 = qJD(1) * t136 - qJD(2);
t164 = qJD(3) * t116 + t100 * t136;
t45 = -qJDD(4) * pkin(4) - t53;
t80 = -t138 * t257 - t254;
t78 = t138 * t255 - t252;
t162 = t116 * t187;
t161 = t138 * (-Ifges(5,2) * t141 - t272);
t160 = t141 * (-Ifges(5,1) * t138 - t271);
t156 = t137 * t239 + t138 * t244;
t155 = t220 - t221;
t153 = -t137 * t275 - t140 * t274;
t146 = -qJD(6) * t141 + (qJ(6) * qJD(4) - t242) * t138;
t127 = t142 * qJ(2);
t124 = qJDD(2) - t262;
t112 = t134 * t140;
t110 = t134 * t137;
t101 = t186 * qJD(1);
t92 = (pkin(5) * t137 - t135) * t141;
t88 = t140 * t106;
t86 = t184 * t141;
t81 = t138 * t252 - t255;
t79 = -t138 * t254 - t257;
t71 = t140 * t94;
t67 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t104;
t57 = -t135 * t259 + t88;
t55 = -pkin(5) * t155 + t318;
t48 = -qJ(6) * t258 + t58;
t47 = -pkin(5) * t97 + qJD(6) + t77;
t44 = -qJ(6) * t253 + t88 + (-t135 * t137 + pkin(5)) * t138;
t16 = -t140 * t222 + t71 + (-qJD(5) * t106 - t317) * t137;
t15 = -t137 * t222 + t198;
t13 = -pkin(5) * t43 + qJDD(6) + t45;
t10 = -qJ(6) * t221 + t137 * t146 + t198;
t5 = pkin(5) * t243 + t71 + t146 * t140 + ((qJ(6) * t141 - t106) * qJD(5) - t317) * t137;
t6 = [(-t262 + t124) * mrSges(3,2) + 0.2e1 * t330 * mrSges(3,3) + (t100 + t238) * mrSges(4,3) + (-Ifges(5,4) * t103 / 0.2e1 - Ifges(5,2) * t104 / 0.2e1 + t361 * t303 + t319 / 0.2e1 + t1 * mrSges(7,1) - Ifges(5,6) * qJDD(4) + t341 * t304 + t342 * t305 + t350) * t138 + (t171 * t365 + t329 * t295 + t327 * t299 + t328 * t301 + t162) * qJD(4) + m(3) * (-pkin(1) * t124 + (t330 + t130) * qJ(2)) + (mrSges(4,3) * t136 + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + (Ifges(5,1) * t358 + Ifges(5,4) * t357 + Ifges(5,5) * qJDD(4) + t13 * t182 + t325 * t303 + t323 * t304 + t321 * t305) * t141 + t253 * t356 + t176 * t357 + t181 * t358 - t332 * t318 - t354 * t245 / 0.2e1 + (-m(5) * t316 + t340 * t141 + m(6) * (-t141 * t45 + t245 * t77) + t138 * t67) * t135 + m(5) * (qJD(2) * t200 + t164) + m(6) * (t15 * t25 + t16 * t24 - t77 * t247 + t3 * t58 + t4 * t57) + (-t295 * t326 - t299 * t322 - t301 * t324) * t239 + (t113 + t330) * mrSges(4,2) + t316 * mrSges(5,3) + t317 * t108 + t160 * t237 / 0.2e1 + qJD(3) * t101 + t100 * t186 + t45 * t86 + t92 * t11 + t55 * t51 + t57 * t20 + t58 * t22 + t10 * t59 + t15 * t60 + t5 * t61 + t16 * t62 + t44 * t19 + t48 * t21 - (t137 * t337 + t140 * t338) * t239 / 0.2e1 + (t155 * t25 + t156 * t24 - t253 * t4 - t258 * t3) * mrSges(6,3) + (-t1 * t253 + t14 * t156 + t155 * t18 - t2 * t258) * mrSges(7,3) + t77 * (-mrSges(6,1) * t155 - mrSges(6,2) * t156) + t47 * (-mrSges(7,1) * t155 - mrSges(7,2) * t156) - t348 * t258 / 0.2e1 + (t257 * t306 + t346 * t79 - t343 * t78 + t314 * (-pkin(7) * t142 + t127) + (-m(4) - m(3)) * t127 + t310 * t142 + (m(3) * pkin(1) + m(6) * t106 + (m(7) - t298) * t136 + t307) * t139) * g(1) + (t255 * t306 - m(3) * t251 - m(4) * t226 + t346 * t81 - t343 * t80 + t314 * (-pkin(7) * t139 + t226) + t310 * t139 + (-m(6) * t190 - t307) * t142) * g(2) + m(7) * (t1 * t44 + t10 * t18 + t13 * t92 + t14 * t5 + t2 * t48 + t47 * t55) + t161 * t201 + t338 * t220 / 0.2e1 + m(4) * (qJ(2) * t113 + qJD(2) * t123 + t164) + t136 * t188 + (t24 * mrSges(6,1) + t14 * mrSges(7,1) - t25 * mrSges(6,2) - t18 * mrSges(7,2) - t351) * t243; t277 * t140 - t276 * t137 + t345 * qJDD(1) + (-m(3) * qJ(2) + t344) * t143 + t154 * qJD(5) + m(6) * (-t137 * t3 - t140 * t4 + (t137 * t24 - t140 * t25) * qJD(5)) + m(7) * (-t1 * t140 - t137 * t2 + (t14 * t137 - t18 * t140) * qJD(5)) + m(3) * t124 - m(5) * t100 + (t234 * t141 + t331 * t138 - m(7) * (-t14 * t259 - t141 * t47 + t18 * t256) - m(6) * (-t141 * t77 - t24 * t259 + t25 * t256) - m(5) * t200) * qJD(1) - t188 + (-g(1) * t139 + g(2) * t142) * (m(3) - t232) + (-qJD(1) * t123 - t100) * m(4); -t143 * mrSges(4,3) + m(4) * t113 + qJDD(1) * mrSges(4,2) + (-t101 + t298 * t116 - m(6) * t166 - m(7) * (t137 * t18 + t14 * t140) + t153) * qJD(1) + (-t11 - t331 * qJD(4) + m(6) * (-t24 * t246 + t244 * t25 - t45) + m(7) * (-t14 * t246 + t18 * t244 - t13) + m(5) * t53 + t340) * t141 + (t67 + t276 * t140 + t277 * t137 + t234 * qJD(4) + t153 * qJD(5) + m(6) * (qJD(4) * t77 + t353) + m(7) * (qJD(4) * t47 + t309) + m(5) * t54) * t138 + t315 * t232; (t138 * t360 + t141 * t359 + t186) * g(3) + t315 * (t138 * t359 - t141 * t360 - t187) + t137 * t356 + (-t225 / 0.2e1 - t241 / 0.2e1) * t338 + t354 * t250 / 0.2e1 + (-t160 / 0.2e1 + t161 / 0.2e1) * t143 + (-t14 * (mrSges(7,1) * t141 + mrSges(7,3) * t256) - t24 * (mrSges(6,1) * t141 + mrSges(6,3) * t256) - t18 * (-mrSges(7,2) * t141 + mrSges(7,3) * t259) - t25 * (-mrSges(6,2) * t141 + mrSges(6,3) * t259) - t162) * qJD(1) + (-pkin(4) * t45 - t105 * t77 - t24 * t49 - t25 * t50) * m(6) + t332 * t105 + t333 * t51 + t322 * t305 + t324 * t304 + t326 * t303 + t309 * mrSges(7,3) - t108 * t261 - t122 * t11 + Ifges(5,6) * t104 + t110 * t19 - t112 * t21 + Ifges(5,5) * t103 + (t119 * t325 + t321 * t98 + t323 * t97) * qJD(5) / 0.2e1 - (t119 * t329 + t327 * t98 + t328 * t97) * qJD(1) / 0.2e1 + t53 * mrSges(5,1) - t54 * mrSges(5,2) - t50 * t60 - t49 * t62 - pkin(4) * t12 + t119 * (t182 * t47 + t184 * t77) + t353 * mrSges(6,3) + t348 * t140 / 0.2e1 + t171 * t201 + t335 * t61 + t336 * t59 + (t1 * t110 - t112 * t2 - t122 * t13 + t14 * t335 + t18 * t336 + t333 * t47) * m(7) + t337 * t240 / 0.2e1 + t13 * t183 + t45 * t185 + (t140 * t22 - t137 * t20 - t62 * t240 - t60 * t241 + m(6) * (-qJD(5) * t166 + t189)) * pkin(8) + t351 * t249 + Ifges(5,3) * qJDD(4); (t313 * t78 - t343 * t79) * g(2) + (-t313 * t80 + t343 * t81) * g(1) + t319 - t77 * (mrSges(6,1) * t98 + mrSges(6,2) * t97) - t47 * (mrSges(7,1) * t98 + mrSges(7,2) * t97) + (-t341 * t98 + t342 * t97) * t367 + (t364 * t97 - t352) * t368 + (-t362 * t98 + t337 + t355) * t369 + (t86 - (-t137 * t334 - t273) * t141) * g(3) + t334 * t1 + t338 * t299 - t17 * t59 + (t62 + t290) * t25 + (-t60 + t291) * t24 + t14 * t289 + (-m(7) * (-t14 + t17) + t61 + t288) * t18 + ((-m(7) * t47 - t51) * t98 + t19) * pkin(5) + t350; -t97 * t59 + t98 * t61 + (-g(3) * t138 + t14 * t98 + t141 * t315 - t18 * t97 + t13) * m(7) + t11;];
tau  = t6;
