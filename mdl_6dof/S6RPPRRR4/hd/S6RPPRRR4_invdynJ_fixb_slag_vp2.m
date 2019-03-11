% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:29
% EndTime: 2019-03-09 02:25:56
% DurationCPUTime: 16.49s
% Computational Cost: add. (7925->655), mult. (14524->881), div. (0->0), fcn. (8796->12), ass. (0->302)
t225 = sin(qJ(4));
t228 = cos(qJ(4));
t313 = qJD(1) * qJD(4);
t173 = -t228 * qJDD(1) + t225 * t313;
t454 = -t173 / 0.2e1;
t174 = -qJDD(1) * t225 - t228 * t313;
t453 = -t174 / 0.2e1;
t221 = sin(pkin(10));
t343 = cos(pkin(10));
t375 = sin(qJ(1));
t376 = cos(qJ(1));
t160 = -t221 * t375 - t343 * t376;
t452 = t160 * g(2);
t226 = cos(qJ(6));
t212 = t228 * qJD(1);
t196 = t212 + qJD(5);
t227 = cos(qJ(5));
t224 = sin(qJ(5));
t315 = t224 * qJD(4);
t326 = qJD(1) * t225;
t163 = t227 * t326 - t315;
t230 = -pkin(1) - pkin(2);
t194 = qJD(1) * t230 + qJD(2);
t277 = qJD(1) * t343;
t138 = qJ(2) * t277 + t221 * t194;
t128 = -qJD(1) * pkin(7) + t138;
t324 = qJD(3) * t225;
t106 = t128 * t228 + t324;
t97 = qJD(4) * pkin(8) + t106;
t327 = qJD(1) * t221;
t137 = -qJ(2) * t327 + t194 * t343;
t127 = qJD(1) * pkin(3) - t137;
t270 = t228 * pkin(4) + t225 * pkin(8);
t98 = qJD(1) * t270 + t127;
t44 = -t224 * t97 + t227 * t98;
t37 = pkin(9) * t163 + t44;
t34 = pkin(5) * t196 + t37;
t223 = sin(qJ(6));
t322 = qJD(4) * t227;
t162 = t224 * t326 + t322;
t45 = t224 * t98 + t227 * t97;
t38 = pkin(9) * t162 + t45;
t350 = t223 * t38;
t13 = t226 * t34 - t350;
t323 = qJD(4) * t225;
t191 = qJDD(1) * t230 + qJDD(2);
t314 = qJD(1) * qJD(2);
t195 = qJDD(1) * qJ(2) + t314;
t121 = t221 * t191 + t343 * t195;
t439 = -qJDD(1) * pkin(7) + qJD(3) * qJD(4) + t121;
t49 = t225 * qJDD(3) - t128 * t323 + t228 * t439;
t46 = qJDD(4) * pkin(8) + t49;
t120 = t191 * t343 - t221 * t195;
t113 = qJDD(1) * pkin(3) - t120;
t60 = -t173 * pkin(4) - t174 * pkin(8) + t113;
t12 = -qJD(5) * t45 - t224 * t46 + t227 * t60;
t159 = qJDD(5) - t173;
t92 = qJD(5) * t162 + qJDD(4) * t224 + t174 * t227;
t6 = pkin(5) * t159 - pkin(9) * t92 + t12;
t318 = qJD(5) * t227;
t320 = qJD(5) * t224;
t11 = t224 * t60 + t227 * t46 + t98 * t318 - t320 * t97;
t93 = qJD(5) * t163 + qJDD(4) * t227 - t174 * t224;
t9 = pkin(9) * t93 + t11;
t2 = qJD(6) * t13 + t223 * t6 + t226 * t9;
t347 = t226 * t38;
t14 = t223 * t34 + t347;
t3 = -qJD(6) * t14 - t223 * t9 + t226 * t6;
t451 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t450 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t229 = -pkin(9) - pkin(8);
t449 = -m(7) * t229 + mrSges(6,3) + mrSges(7,3);
t309 = qJD(5) + qJD(6);
t193 = t212 + t309;
t378 = t193 / 0.2e1;
t102 = t162 * t223 - t163 * t226;
t385 = t102 / 0.2e1;
t275 = t226 * t162 + t163 * t223;
t387 = t275 / 0.2e1;
t448 = Ifges(7,5) * t385 + Ifges(7,6) * t387 + Ifges(7,3) * t378;
t447 = t162 / 0.2e1;
t381 = -t163 / 0.2e1;
t446 = t196 / 0.2e1;
t368 = -qJD(1) / 0.2e1;
t295 = t224 * t212;
t296 = qJD(5) * t229;
t105 = qJD(3) * t228 - t225 * t128;
t269 = -pkin(4) * t225 + pkin(8) * t228;
t170 = t269 * qJD(1);
t65 = t227 * t105 + t224 * t170;
t445 = -pkin(9) * t295 + t224 * t296 - t65;
t331 = t227 * t228;
t254 = -pkin(5) * t225 + pkin(9) * t331;
t64 = -t105 * t224 + t227 * t170;
t444 = -qJD(1) * t254 + t227 * t296 - t64;
t148 = t221 * t331 - t224 * t343;
t279 = t228 * t343;
t292 = t225 * t315;
t419 = -qJD(5) * t148 + t221 * t292 - (t221 * t227 - t224 * t279) * qJD(1);
t334 = t224 * t228;
t147 = -t221 * t334 - t227 * t343;
t290 = t225 * t322;
t418 = qJD(5) * t147 - t221 * t290 - (t221 * t224 + t227 * t279) * qJD(1);
t219 = qJ(5) + qJ(6);
t213 = sin(t219);
t214 = cos(t219);
t267 = -mrSges(6,1) * t227 + mrSges(6,2) * t224;
t205 = pkin(5) * t227 + pkin(4);
t434 = m(7) * t205;
t443 = m(6) * pkin(4) + mrSges(7,1) * t214 - mrSges(7,2) * t213 - t267 + t434;
t442 = -m(6) * pkin(8) - t449;
t302 = mrSges(5,3) * t326;
t420 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t162 - mrSges(6,2) * t163 - t302;
t256 = t11 * t227 - t12 * t224;
t438 = -t44 * t318 - t45 * t320 + t256;
t437 = t213 * mrSges(7,1) + t214 * mrSges(7,2);
t436 = -t225 * t442 + t228 * t443;
t154 = qJDD(6) + t159;
t384 = t154 / 0.2e1;
t26 = -qJD(6) * t102 - t223 * t92 + t226 * t93;
t396 = t26 / 0.2e1;
t25 = qJD(6) * t275 + t223 * t93 + t226 * t92;
t397 = t25 / 0.2e1;
t399 = Ifges(7,4) * t397 + Ifges(7,2) * t396 + Ifges(7,6) * t384;
t361 = Ifges(5,4) * t225;
t262 = -t228 * Ifges(5,2) - t361;
t435 = t13 * mrSges(7,1) - Ifges(5,6) * qJD(4) / 0.2e1 + t262 * t368 + Ifges(6,5) * t381 + Ifges(6,6) * t447 + Ifges(6,3) * t446 - t14 * mrSges(7,2) + t448;
t398 = m(7) * pkin(5);
t433 = mrSges(3,1) + mrSges(2,1);
t432 = mrSges(4,2) - mrSges(5,3);
t431 = -mrSges(3,3) + mrSges(2,2);
t187 = t229 * t224;
t188 = t229 * t227;
t118 = t187 * t223 - t188 * t226;
t430 = -qJD(6) * t118 - t223 * t445 + t226 * t444;
t117 = t187 * t226 + t188 * t223;
t429 = qJD(6) * t117 + t223 * t444 + t226 * t445;
t80 = t147 * t226 - t148 * t223;
t428 = qJD(6) * t80 + t223 * t419 + t226 * t418;
t81 = t147 * t223 + t148 * t226;
t427 = -qJD(6) * t81 - t223 * t418 + t226 * t419;
t266 = t224 * mrSges(6,1) + t227 * mrSges(6,2);
t96 = -qJD(4) * pkin(4) - t105;
t425 = t266 * t96;
t186 = mrSges(5,1) * t228 - mrSges(5,2) * t225;
t424 = -mrSges(4,1) - t186;
t43 = -mrSges(6,1) * t93 + mrSges(6,2) * t92;
t423 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t174 + t43;
t422 = -t398 - mrSges(6,1);
t373 = pkin(5) * t224;
t421 = -t324 - (-qJD(1) * t373 + t128) * t228 + pkin(5) * t320;
t167 = t223 * t227 + t224 * t226;
t140 = t167 * t225;
t175 = -t221 * qJ(2) + t230 * t343;
t164 = pkin(3) - t175;
t129 = t164 + t270;
t176 = t343 * qJ(2) + t221 * t230;
t165 = -pkin(7) + t176;
t139 = t165 * t331;
t72 = t224 * t129 + t139;
t360 = Ifges(5,4) * t228;
t265 = -Ifges(5,1) * t225 - t360;
t158 = Ifges(6,4) * t162;
t85 = -t163 * Ifges(6,1) + t196 * Ifges(6,5) + t158;
t417 = Ifges(5,5) * qJD(4) + qJD(1) * t265 + t227 * t85;
t415 = t225 * (-Ifges(5,1) * t228 + t361) + t228 * (Ifges(5,2) * t225 - t360);
t321 = qJD(4) * t228;
t289 = t227 * t321;
t319 = qJD(5) * t225;
t243 = -t224 * t319 + t289;
t276 = qJD(2) * t343;
t271 = t225 * t276;
t413 = t165 * t321 + t271;
t50 = qJDD(3) * t228 - t128 * t321 - t225 * t439;
t411 = -t225 * t50 + t228 * t49;
t61 = mrSges(6,1) * t159 - mrSges(6,3) * t92;
t62 = -mrSges(6,2) * t159 + mrSges(6,3) * t93;
t410 = -t224 * t61 + t227 * t62;
t409 = -m(7) - m(6) - m(5);
t10 = -mrSges(7,1) * t26 + mrSges(7,2) * t25;
t408 = t10 + t423;
t407 = m(4) - t409;
t66 = -pkin(5) * t162 + t96;
t405 = -mrSges(7,1) * t66 + t14 * mrSges(7,3);
t404 = mrSges(7,2) * t66 - t13 * mrSges(7,3);
t402 = -m(6) * t96 - t420;
t142 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t173;
t301 = mrSges(5,3) * t212;
t185 = -qJD(4) * mrSges(5,2) - t301;
t47 = -qJDD(4) * pkin(4) - t50;
t401 = m(5) * ((-t105 * t228 - t106 * t225) * qJD(4) + t411) + m(6) * (t225 * t47 + t321 * t96) - t185 * t323 + t142 * t228;
t48 = -mrSges(7,1) * t275 + mrSges(7,2) * t102;
t400 = -m(7) * t66 + t402 - t48;
t231 = qJD(1) ^ 2;
t356 = Ifges(7,4) * t102;
t40 = Ifges(7,2) * t275 + Ifges(7,6) * t193 + t356;
t395 = -t40 / 0.2e1;
t394 = t40 / 0.2e1;
t94 = Ifges(7,4) * t275;
t41 = Ifges(7,1) * t102 + Ifges(7,5) * t193 + t94;
t393 = -t41 / 0.2e1;
t392 = t41 / 0.2e1;
t390 = t92 / 0.2e1;
t389 = t93 / 0.2e1;
t388 = -t275 / 0.2e1;
t386 = -t102 / 0.2e1;
t383 = t159 / 0.2e1;
t380 = t163 / 0.2e1;
t379 = -t193 / 0.2e1;
t374 = pkin(5) * t163;
t367 = qJD(5) / 0.2e1;
t161 = t221 * t376 - t343 * t375;
t337 = t214 * t228;
t338 = t213 * t228;
t366 = (-t160 * t214 + t161 * t338) * mrSges(7,1) + (t160 * t213 + t161 * t337) * mrSges(7,2);
t78 = t160 * t338 + t161 * t214;
t79 = -t160 * t337 + t161 * t213;
t365 = t78 * mrSges(7,1) - t79 * mrSges(7,2);
t364 = mrSges(4,1) * t221;
t363 = mrSges(6,3) * t162;
t362 = mrSges(6,3) * t163;
t359 = Ifges(6,4) * t163;
t358 = Ifges(6,4) * t224;
t357 = Ifges(6,4) * t227;
t341 = t161 * t224;
t335 = t224 * t225;
t333 = t225 * t227;
t332 = t226 * t227;
t330 = t437 * t225;
t329 = t376 * pkin(1) + t375 * qJ(2);
t328 = mrSges(4,1) * qJDD(1);
t325 = qJD(2) * t221;
t317 = qJDD(1) * mrSges(3,1);
t316 = qJDD(1) * mrSges(4,2);
t308 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t154;
t307 = Ifges(6,5) * t92 + Ifges(6,6) * t93 + Ifges(6,3) * t159;
t298 = t223 * t335;
t297 = t376 * pkin(2) + t329;
t291 = t228 * t315;
t281 = t318 / 0.2e1;
t278 = -t313 / 0.2e1;
t274 = -pkin(1) * t375 + t376 * qJ(2);
t272 = t228 * t276;
t268 = mrSges(5,1) * t225 + mrSges(5,2) * t228;
t264 = Ifges(6,1) * t227 - t358;
t263 = Ifges(6,1) * t224 + t357;
t261 = -Ifges(6,2) * t224 + t357;
t260 = Ifges(6,2) * t227 + t358;
t259 = -Ifges(5,5) * t228 + Ifges(5,6) * t225;
t258 = Ifges(6,5) * t227 - Ifges(6,6) * t224;
t257 = Ifges(6,5) * t224 + Ifges(6,6) * t227;
t123 = t227 * t129;
t59 = pkin(9) * t333 + t123 + (-t165 * t224 + pkin(5)) * t228;
t63 = pkin(9) * t335 + t72;
t27 = -t223 * t63 + t226 * t59;
t28 = t223 * t59 + t226 * t63;
t255 = t223 * t224 - t332;
t253 = t308 + t451;
t89 = t160 * t334 + t161 * t227;
t249 = t127 * t268;
t246 = -t137 * t221 + t138 * t343;
t242 = t225 * t318 + t291;
t241 = -pkin(2) * t375 + t274;
t144 = qJD(4) * t269 + t325;
t240 = t227 * t144 + t165 * t292 - t224 * t272;
t238 = -Ifges(6,5) * t225 - t228 * t264;
t237 = -Ifges(6,6) * t225 - t228 * t261;
t236 = -Ifges(6,3) * t225 - t228 * t258;
t234 = -t105 * t225 * t343 + t106 * t279 + t127 * t221;
t32 = t129 * t318 + t224 * t144 + t227 * t272 + (-t228 * t320 - t290) * t165;
t206 = -qJDD(1) * pkin(1) + qJDD(2);
t168 = t186 * qJD(1);
t157 = t266 * t225;
t141 = t225 * t332 - t298;
t131 = t255 * t212;
t130 = t167 * t212;
t126 = (t165 - t373) * t225;
t125 = mrSges(6,1) * t196 + t362;
t124 = -mrSges(6,2) * t196 + t363;
t112 = -mrSges(5,1) * t173 + mrSges(5,2) * t174;
t90 = -t160 * t331 + t341;
t84 = t162 * Ifges(6,2) + t196 * Ifges(6,6) - t359;
t73 = -pkin(5) * t242 + t413;
t71 = -t165 * t334 + t123;
t68 = mrSges(7,1) * t193 - mrSges(7,3) * t102;
t67 = -mrSges(7,2) * t193 + mrSges(7,3) * t275;
t53 = -qJD(6) * t298 + (t309 * t333 + t291) * t226 + t243 * t223;
t52 = t140 * t309 + t223 * t291 - t226 * t289;
t36 = t92 * Ifges(6,1) + t93 * Ifges(6,4) + t159 * Ifges(6,5);
t35 = t92 * Ifges(6,4) + t93 * Ifges(6,2) + t159 * Ifges(6,6);
t33 = -qJD(5) * t72 + t240;
t29 = -pkin(5) * t93 + t47;
t24 = pkin(9) * t242 + t32;
t21 = t254 * qJD(4) + (-t139 + (-pkin(9) * t225 - t129) * t224) * qJD(5) + t240;
t18 = -mrSges(7,2) * t154 + mrSges(7,3) * t26;
t17 = mrSges(7,1) * t154 - mrSges(7,3) * t25;
t16 = t226 * t37 - t350;
t15 = -t223 * t37 - t347;
t8 = t25 * Ifges(7,1) + t26 * Ifges(7,4) + t154 * Ifges(7,5);
t5 = -qJD(6) * t28 + t21 * t226 - t223 * t24;
t4 = qJD(6) * t27 + t21 * t223 + t226 * t24;
t1 = [(Ifges(7,1) * t52 + Ifges(7,4) * t53) * t385 + (-Ifges(7,1) * t141 + Ifges(7,4) * t140) * t397 + (-t44 * mrSges(6,1) + t45 * mrSges(6,2) + t106 * mrSges(5,3) - t435 - t448) * t323 + m(7) * (t126 * t29 + t13 * t5 + t14 * t4 + t2 * t28 + t27 * t3 + t66 * t73) + (t105 * t321 - t411) * mrSges(5,3) + (qJD(1) * t276 + t121) * mrSges(4,2) + m(6) * (t11 * t72 + t12 * t71 + t271 * t96 + t45 * t32 + t44 * t33) + m(5) * (qJD(2) * t234 + t113 * t164) + t401 * t165 + t185 * t272 + t168 * t325 + t176 * t316 + pkin(1) * t317 - t36 * t333 / 0.2e1 + m(3) * (-pkin(1) * t206 + (t195 + t314) * qJ(2)) + (qJD(5) * t85 + t35) * t335 / 0.2e1 + (t257 * t446 + t260 * t447 + t263 * t381) * t319 + (-t341 * t398 - m(3) * t329 - m(4) * t297 - t90 * mrSges(6,1) - t79 * mrSges(7,1) - t89 * mrSges(6,2) - t78 * mrSges(7,2) - t433 * t376 + t431 * t375 + t409 * (-t160 * pkin(3) + pkin(7) * t161 + t297) + t432 * t161 + (m(6) * t270 - t424) * t160) * g(2) + (-t13 * t52 + t14 * t53 + t140 * t2 + t141 * t3) * mrSges(7,3) + (t236 * t446 + t237 * t447 + t238 * t381 - t249 + t259 * qJD(4) / 0.2e1) * qJD(4) + (Ifges(7,5) * t52 + Ifges(7,6) * t53) * t378 + (-Ifges(7,5) * t141 + Ifges(7,6) * t140) * t384 + (Ifges(7,4) * t52 + Ifges(7,2) * t53) * t387 + (-Ifges(7,4) * t141 + Ifges(7,2) * t140) * t396 + (Ifges(4,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1) + (-m(3) * t274 - m(4) * t241 + t431 * t376 + t433 * t375 + t409 * (t161 * pkin(3) + t241) + (t424 - t436) * t161 + (-m(7) * t373 + pkin(7) * t409 - t266 + t432 - t437) * t160) * g(1) + t420 * t413 - t417 * t321 / 0.2e1 + t52 * t392 + t415 * t278 + (Ifges(5,4) * t453 + Ifges(5,2) * t454 + t434 * t452 + t308 / 0.2e1 + t307 / 0.2e1 + Ifges(6,3) * t383 + Ifges(7,3) * t384 + Ifges(6,6) * t389 + Ifges(6,5) * t390 - Ifges(5,6) * qJDD(4) + Ifges(7,6) * t396 + Ifges(7,5) * t397 + t450 + t451) * t228 + (Ifges(5,1) * t453 + Ifges(5,4) * t454 - Ifges(5,5) * qJDD(4) + t423 * t165 - t258 * t383 - t261 * t389 - t264 * t390 + t84 * t281 + t449 * t452) * t225 + m(4) * (qJD(2) * t246 + t120 * t175 + t121 * t176) + t96 * (-mrSges(6,1) * t242 - mrSges(6,2) * t243) + t173 * t262 / 0.2e1 + t174 * t265 / 0.2e1 + (t11 * t335 + t12 * t333 + t242 * t45 + t243 * t44) * mrSges(6,3) + t84 * t291 / 0.2e1 + t314 * t364 - t206 * mrSges(3,1) + 0.2e1 * t195 * mrSges(3,3) + t113 * t186 + t164 * t112 - t47 * t157 + t29 * (-mrSges(7,1) * t140 - mrSges(7,2) * t141) - t141 * t8 / 0.2e1 - t120 * mrSges(4,1) + t32 * t124 + t33 * t125 + t126 * t10 + t4 * t67 + t5 * t68 + t71 * t61 + t72 * t62 + t73 * t48 + t66 * (-mrSges(7,1) * t53 + mrSges(7,2) * t52) + t53 * t394 + t140 * t399 - t175 * t328 + t27 * t17 + t28 * t18; m(3) * t206 + t147 * t61 + t148 * t62 - t168 * t327 + t80 * t17 + t81 * t18 - t317 + t427 * t68 + t428 * t67 + t419 * t125 + t418 * t124 + (-m(4) * t246 - m(5) * t234) * qJD(1) + (t13 * t427 + t14 * t428 + t2 * t81 + t3 * t80) * m(7) + (t11 * t148 + t12 * t147 + t418 * t45 + t419 * t44) * m(6) + (m(4) * t120 - m(5) * t113 - t112 - t328) * t343 + (-m(3) * qJ(2) - mrSges(4,2) * t343 - mrSges(3,3) - t364) * t231 + (-g(1) * t375 + g(2) * t376) * (m(3) + t407) + (-t185 * t228 + t225 * t400) * t277 + (t408 * t225 + (t48 + t420) * t321 + m(4) * t121 + t316 + m(7) * (t225 * t29 + t321 * t66) + t401) * t221; m(4) * qJDD(3) - t140 * t17 + t141 * t18 - t52 * t67 - t53 * t68 + m(7) * (-t13 * t53 - t14 * t52 - t140 * t3 + t141 * t2) + t407 * g(3) + ((t124 * t227 - t125 * t224 + t185) * qJD(4) + m(5) * (qJD(4) * t106 + t50) - m(7) * t29 + m(6) * (-t315 * t44 + t322 * t45 - t47) - t408) * t228 + (t142 + (-t224 * t124 - t227 * t125) * qJD(5) + m(5) * t49 + m(6) * t438 + (-m(5) * t105 - t400) * qJD(4) + t410) * t225; (t186 + t436) * g(3) + (-(Ifges(7,4) * t385 + Ifges(7,2) * t387 + Ifges(7,6) * t378 + t394 + t405) * t309 + t29 * mrSges(7,2) + Ifges(7,4) * t396 + Ifges(7,1) * t397 + Ifges(7,5) * t384 - t3 * mrSges(7,3) + t8 / 0.2e1) * t167 + (t236 * t368 + t258 * t367) * t196 + t429 * t67 + (t117 * t3 + t118 * t2 + t13 * t430 + t14 * t429 - t205 * t29 + t421 * t66) * m(7) + t430 * t68 + (-(Ifges(7,1) * t385 + Ifges(7,4) * t387 + Ifges(7,5) * t378 + t392 + t404) * t309 + mrSges(7,1) * t29 - mrSges(7,3) * t2 - 0.2e1 * t399) * t255 + (Ifges(7,5) * t131 + Ifges(7,6) * t130) * t379 + t421 * t48 + (t264 * t381 + t425) * qJD(5) + (-t45 * (mrSges(6,2) * t225 + mrSges(6,3) * t334) - t44 * (-mrSges(6,1) * t225 + mrSges(6,3) * t331) + t249 + t238 * t380) * qJD(1) - (t320 + t295) * t84 / 0.2e1 + (-t301 - t185) * t105 + (g(1) * t160 + g(2) * t161) * (-t225 * t443 - t228 * t442 - t268) + (-t302 + t402) * t106 + (Ifges(7,4) * t131 + Ifges(7,2) * t130) * t388 + (t237 * t368 + t261 * t367) * t162 + t438 * mrSges(6,3) + (t13 * t131 - t130 * t14) * mrSges(7,3) + t257 * t383 + t260 * t389 + t263 * t390 + t131 * t393 + t85 * t281 + t259 * t278 + (-t125 * t318 - t124 * t320 + m(6) * ((-t224 * t45 - t227 * t44) * qJD(5) + t256) + t410) * pkin(8) + t415 * t231 / 0.2e1 + (Ifges(7,1) * t131 + Ifges(7,4) * t130) * t386 + t47 * t267 + (-Ifges(7,5) * t386 - Ifges(7,6) * t388 - Ifges(7,3) * t379 + t435) * t326 + (t425 + t417 / 0.2e1) * t212 + (-pkin(4) * t47 - t44 * t64 - t45 * t65) * m(6) + t227 * t35 / 0.2e1 + t224 * t36 / 0.2e1 - t205 * t10 + Ifges(5,5) * t174 + Ifges(5,6) * t173 + Ifges(5,3) * qJDD(4) - t66 * (-mrSges(7,1) * t130 + mrSges(7,2) * t131) + t117 * t17 + t118 * t18 - t65 * t124 - t64 * t125 + t130 * t395 - pkin(4) * t43 - t49 * mrSges(5,2) + t50 * mrSges(5,1); (-t362 + t125) * t45 + (-t366 - (t160 * t224 + t161 * t331) * mrSges(6,2) + t422 * (-t160 * t227 + t161 * t334)) * g(2) + (mrSges(6,2) * t90 + t422 * t89 - t365) * g(1) + (-t335 * t398 - t157 - t330) * g(3) - (Ifges(6,2) * t163 + t158 + t85) * t162 / 0.2e1 + (t363 - t124) * t44 + t307 - m(7) * (t13 * t15 + t14 * t16 - t374 * t66) + (Ifges(7,1) * t386 + Ifges(7,4) * t388 + Ifges(7,5) * t379 + t393 - t404) * t275 - (Ifges(7,4) * t386 + Ifges(7,2) * t388 + Ifges(7,6) * t379 + t395 - t405) * t102 + t450 + t253 + (Ifges(6,1) * t162 + t359) * t380 + t84 * t381 + t48 * t374 + ((-t223 * t68 + t226 * t67) * qJD(6) + t226 * t17 + t223 * t18) * pkin(5) - t196 * (Ifges(6,5) * t162 + Ifges(6,6) * t163) / 0.2e1 - t96 * (-mrSges(6,1) * t163 + mrSges(6,2) * t162) - t16 * t67 - t15 * t68 + (t2 * t223 + t226 * t3 + (-t13 * t223 + t14 * t226) * qJD(6)) * t398; -t66 * (mrSges(7,1) * t102 + mrSges(7,2) * t275) + (Ifges(7,1) * t275 - t356) * t386 + t40 * t385 + (Ifges(7,5) * t275 - Ifges(7,6) * t102) * t379 - t13 * t67 + t14 * t68 - g(1) * t365 - g(2) * t366 - g(3) * t330 + (t102 * t14 + t13 * t275) * mrSges(7,3) + t253 + (-Ifges(7,2) * t102 + t41 + t94) * t388;];
tau  = t1;
