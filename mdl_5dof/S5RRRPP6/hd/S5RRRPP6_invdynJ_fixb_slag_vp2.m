% Calculate vector of inverse dynamics joint torques for
% S5RRRPP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:59:59
% EndTime: 2019-12-31 21:00:33
% DurationCPUTime: 19.23s
% Computational Cost: add. (5256->597), mult. (11912->785), div. (0->0), fcn. (7843->10), ass. (0->280)
t445 = mrSges(5,1) + mrSges(6,1);
t444 = -mrSges(5,2) + mrSges(6,3);
t224 = sin(qJ(3));
t225 = sin(qJ(2));
t308 = qJD(1) * t225;
t287 = t224 * t308;
t227 = cos(qJ(3));
t305 = qJD(2) * t227;
t171 = -t287 + t305;
t285 = t227 * t308;
t172 = qJD(2) * t224 + t285;
t222 = sin(pkin(8));
t328 = cos(pkin(8));
t105 = -t328 * t171 + t172 * t222;
t368 = -t105 / 0.2e1;
t228 = cos(qJ(2));
t307 = qJD(1) * t228;
t199 = qJD(3) - t307;
t356 = -t199 / 0.2e1;
t355 = t199 / 0.2e1;
t244 = t222 * t171 + t172 * t328;
t365 = -t244 / 0.2e1;
t364 = t244 / 0.2e1;
t415 = Ifges(6,2) + Ifges(5,3);
t436 = Ifges(4,3) + t415;
t266 = pkin(2) * t225 - pkin(7) * t228;
t173 = t266 * qJD(1);
t119 = pkin(6) * t287 + t227 * t173;
t314 = t227 * t228;
t248 = pkin(3) * t225 - qJ(4) * t314;
t223 = -qJ(4) - pkin(7);
t272 = qJD(3) * t223;
t443 = -qJD(1) * t248 - qJD(4) * t224 + t227 * t272 - t119;
t151 = t224 * t173;
t300 = qJD(4) * t227;
t319 = t225 * t227;
t321 = t224 * t228;
t442 = t151 + (-pkin(6) * t319 - qJ(4) * t321) * qJD(1) - t224 * t272 - t300;
t212 = pkin(6) * t308;
t187 = -qJD(2) * pkin(2) + t212;
t118 = -pkin(3) * t171 + qJD(4) + t187;
t267 = pkin(2) * t228 + pkin(7) * t225;
t182 = -pkin(1) - t267;
t158 = t182 * qJD(1);
t213 = pkin(6) * t307;
t188 = qJD(2) * pkin(7) + t213;
t111 = t227 * t158 - t188 * t224;
t79 = -qJ(4) * t172 + t111;
t71 = pkin(3) * t199 + t79;
t112 = t158 * t224 + t188 * t227;
t80 = qJ(4) * t171 + t112;
t75 = t328 * t80;
t24 = t222 * t71 + t75;
t22 = qJ(5) * t199 + t24;
t31 = pkin(4) * t105 - qJ(5) * t244 + t118;
t441 = Ifges(5,4) * t364 + Ifges(6,5) * t365 + Ifges(5,6) * t355 + Ifges(6,6) * t356 + (Ifges(5,2) + Ifges(6,3)) * t368 - t118 * mrSges(5,1) - t31 * mrSges(6,1) + t22 * mrSges(6,2) + t24 * mrSges(5,3);
t299 = qJD(1) * qJD(2);
t176 = qJDD(1) * t225 + t228 * t299;
t100 = qJD(3) * t171 + qJDD(2) * t224 + t176 * t227;
t101 = -qJD(3) * t172 + qJDD(2) * t227 - t176 * t224;
t53 = t100 * t222 - t101 * t328;
t375 = -t53 / 0.2e1;
t55 = t100 * t328 + t222 * t101;
t373 = t55 / 0.2e1;
t423 = -m(6) - m(5);
t175 = qJDD(1) * t228 - t225 * t299;
t165 = qJDD(3) - t175;
t360 = t165 / 0.2e1;
t163 = t175 * pkin(6);
t440 = -mrSges(5,3) - mrSges(6,2);
t418 = Ifges(5,1) + Ifges(6,1);
t417 = Ifges(5,4) - Ifges(6,5);
t416 = Ifges(6,4) + Ifges(5,5);
t414 = Ifges(5,6) - Ifges(6,6);
t439 = -pkin(3) * t423 + mrSges(4,1);
t367 = t105 / 0.2e1;
t437 = Ifges(5,2) * t368 - Ifges(6,3) * t367 + t417 * t364 + t441;
t221 = qJ(3) + pkin(8);
t215 = sin(t221);
t216 = cos(t221);
t262 = -mrSges(4,1) * t227 + mrSges(4,2) * t224;
t435 = m(4) * pkin(2) + t444 * t215 + t445 * t216 - t262;
t410 = t442 * t222 + t443 * t328;
t409 = t443 * t222 - t442 * t328;
t304 = qJD(2) * t228;
t284 = t224 * t304;
t301 = qJD(3) * t227;
t239 = t225 * t301 + t284;
t226 = sin(qJ(1));
t229 = cos(qJ(1));
t313 = t228 * t229;
t149 = -t224 * t313 + t226 * t227;
t286 = t224 * t307;
t303 = qJD(3) * t224;
t402 = -t213 + (-t286 + t303) * pkin(3);
t432 = g(1) * t229 + g(2) * t226;
t335 = t222 * t80;
t23 = t328 * t71 - t335;
t21 = -t199 * pkin(4) + qJD(5) - t23;
t431 = -mrSges(6,2) * t21 + mrSges(5,3) * t23;
t429 = -Ifges(5,2) * t367 + Ifges(6,3) * t368 - t356 * t414 - t365 * t417 + t441;
t411 = -t105 * t417 + t199 * t416 + t244 * t418;
t428 = -t411 / 0.2e1 + t431;
t327 = qJDD(1) * pkin(1);
t113 = -pkin(2) * t175 - pkin(7) * t176 - t327;
t139 = qJDD(2) * pkin(7) + t163;
t35 = -t112 * qJD(3) + t227 * t113 - t139 * t224;
t17 = pkin(3) * t165 - qJ(4) * t100 - qJD(4) * t172 + t35;
t34 = t224 * t113 + t227 * t139 + t158 * t301 - t188 * t303;
t19 = qJ(4) * t101 + qJD(4) * t171 + t34;
t3 = t17 * t328 - t222 * t19;
t2 = -t165 * pkin(4) + qJDD(5) - t3;
t374 = t53 / 0.2e1;
t164 = t176 * pkin(6);
t140 = -qJDD(2) * pkin(2) + t164;
t70 = -pkin(3) * t101 + qJDD(4) + t140;
t5 = pkin(4) * t53 - qJ(5) * t55 - qJD(5) * t244 + t70;
t427 = mrSges(5,2) * t70 + mrSges(6,2) * t2 - mrSges(5,3) * t3 - mrSges(6,3) * t5 + Ifges(6,5) * t374 + 0.2e1 * t360 * t416 + 0.2e1 * t373 * t418 + (Ifges(5,4) + t417) * t375;
t391 = mrSges(5,2) * t118 - mrSges(6,3) * t31;
t425 = Ifges(5,4) * t368 + Ifges(6,5) * t367 + t355 * t416 + t364 * t418 + t391 - t428;
t424 = -m(4) - m(3);
t370 = t100 / 0.2e1;
t369 = t101 / 0.2e1;
t421 = t175 / 0.2e1;
t420 = t176 / 0.2e1;
t214 = pkin(6) * t304;
t419 = -mrSges(3,3) + mrSges(2,2);
t413 = -qJ(5) * t308 + t409;
t412 = pkin(4) * t308 - t410;
t271 = t328 * t224;
t167 = t222 * t227 + t271;
t124 = t167 * t307;
t270 = t328 * t227;
t324 = t222 * t224;
t243 = t270 - t324;
t236 = t228 * t243;
t125 = qJD(1) * t236;
t143 = t167 * qJD(3);
t144 = t243 * qJD(3);
t408 = -qJD(5) * t167 + t402 + (-t144 + t125) * qJ(5) + (-t124 + t143) * pkin(4);
t81 = -mrSges(5,2) * t199 - mrSges(5,3) * t105;
t84 = -mrSges(6,2) * t105 + mrSges(6,3) * t199;
t407 = -t81 - t84;
t82 = mrSges(5,1) * t199 - mrSges(5,3) * t244;
t83 = -mrSges(6,1) * t199 + mrSges(6,2) * t244;
t406 = t82 - t83;
t211 = Ifges(3,4) * t307;
t162 = Ifges(4,4) * t171;
t94 = t172 * Ifges(4,1) + t199 * Ifges(4,5) + t162;
t405 = Ifges(3,1) * t308 + Ifges(3,5) * qJD(2) + t227 * t94 + t211;
t404 = -qJD(2) * mrSges(3,1) - mrSges(4,1) * t171 + mrSges(4,2) * t172 + mrSges(3,3) * t308;
t203 = pkin(6) * t314;
t127 = t224 * t182 + t203;
t302 = qJD(3) * t225;
t403 = t224 * t302 - t227 * t304;
t401 = t163 * t228 + t164 * t225;
t400 = t440 * t225;
t399 = -t224 * t35 + t227 * t34;
t398 = t172 * Ifges(4,5) + t171 * Ifges(4,6) - t105 * t414 + t199 * t436 + t244 * t416;
t264 = mrSges(3,1) * t228 - mrSges(3,2) * t225;
t397 = t225 * mrSges(4,3) + mrSges(2,1) + t264;
t396 = Ifges(4,5) * t100 + Ifges(4,6) * t101 + t165 * t436 - t414 * t53 + t416 * t55;
t395 = m(6) * pkin(4) + t445;
t393 = m(6) * qJ(5) + t444;
t342 = Ifges(3,4) * t225;
t256 = t228 * Ifges(3,2) + t342;
t389 = t21 * mrSges(6,1) + t24 * mrSges(5,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t256 / 0.2e1 - t22 * mrSges(6,3) - t23 * mrSges(5,1);
t4 = t222 * t17 + t328 * t19;
t1 = qJ(5) * t165 + qJD(5) * t199 + t4;
t388 = -t35 * mrSges(4,1) - t3 * mrSges(5,1) + t2 * mrSges(6,1) + t34 * mrSges(4,2) + t4 * mrSges(5,2) - t1 * mrSges(6,3);
t386 = Ifges(5,4) * t367 + Ifges(6,5) * t368 + t356 * t416 + t365 * t418 - t391;
t383 = mrSges(5,1) * t70 + mrSges(6,1) * t5 - mrSges(6,2) * t1 - mrSges(5,3) * t4 + 0.2e1 * Ifges(6,3) * t374 - t55 * Ifges(5,4) / 0.2e1 - t165 * Ifges(5,6) / 0.2e1 + Ifges(6,6) * t360 + (-t417 + Ifges(6,5)) * t373 + (-t375 + t374) * Ifges(5,2);
t379 = m(5) * pkin(3);
t378 = Ifges(4,1) * t370 + Ifges(4,4) * t369 + Ifges(4,5) * t360;
t336 = t172 * Ifges(4,4);
t93 = t171 * Ifges(4,2) + t199 * Ifges(4,6) + t336;
t371 = -t93 / 0.2e1;
t357 = t172 / 0.2e1;
t351 = pkin(3) * t172;
t350 = pkin(3) * t222;
t347 = g(3) * t225;
t217 = t225 * pkin(6);
t174 = t266 * qJD(2);
t306 = qJD(2) * t225;
t292 = pkin(6) * t306;
t310 = t227 * t174 + t224 * t292;
t57 = -t225 * t300 + t248 * qJD(2) + (-t203 + (qJ(4) * t225 - t182) * t224) * qJD(3) + t310;
t311 = t224 * t174 + t182 * t301;
t62 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t319 + (-qJD(4) * t225 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t228) * t224 + t311;
t16 = t222 * t57 + t328 * t62;
t341 = Ifges(3,4) * t228;
t340 = Ifges(4,4) * t224;
t339 = Ifges(4,4) * t227;
t338 = t111 * mrSges(4,3);
t337 = t112 * mrSges(4,3);
t169 = t227 * t182;
t109 = -qJ(4) * t319 + t169 + (-pkin(6) * t224 - pkin(3)) * t228;
t322 = t224 * t225;
t115 = -qJ(4) * t322 + t127;
t66 = t222 * t109 + t328 * t115;
t323 = t223 * t225;
t320 = t224 * t229;
t318 = t225 * t229;
t317 = t226 * t224;
t315 = t226 * t228;
t208 = pkin(3) * t227 + pkin(2);
t189 = t228 * t208;
t312 = t229 * t215;
t200 = pkin(3) * t322;
t177 = t217 + t200;
t309 = t229 * pkin(1) + t226 * pkin(6);
t291 = m(4) * pkin(7) + mrSges(4,3);
t123 = t239 * pkin(3) + t214;
t288 = t328 * pkin(3);
t13 = t53 * mrSges(5,1) + t55 * mrSges(5,2);
t12 = t53 * mrSges(6,1) - t55 * mrSges(6,3);
t30 = -t165 * mrSges(6,1) + t55 * mrSges(6,2);
t269 = t299 / 0.2e1;
t263 = mrSges(3,1) * t225 + mrSges(3,2) * t228;
t261 = mrSges(4,1) * t224 + mrSges(4,2) * t227;
t258 = Ifges(4,1) * t227 - t340;
t257 = Ifges(4,1) * t224 + t339;
t255 = -Ifges(4,2) * t224 + t339;
t254 = Ifges(4,2) * t227 + t340;
t253 = Ifges(3,5) * t228 - Ifges(3,6) * t225;
t252 = Ifges(4,5) * t227 - Ifges(4,6) * t224;
t251 = Ifges(4,5) * t224 + Ifges(4,6) * t227;
t250 = pkin(4) * t216 + qJ(5) * t215;
t247 = pkin(1) * t263;
t147 = t224 * t315 + t227 * t229;
t15 = -t222 * t62 + t328 * t57;
t245 = t225 * (Ifges(3,1) * t228 - t342);
t65 = t109 * t328 - t222 * t115;
t235 = Ifges(4,5) * t225 + t228 * t258;
t234 = Ifges(4,6) * t225 + t228 * t255;
t233 = Ifges(4,3) * t225 + t228 * t252;
t219 = t229 * pkin(6);
t207 = -t288 - pkin(4);
t205 = qJ(5) + t350;
t186 = t223 * t227;
t184 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t307;
t159 = t261 * t225;
t150 = t227 * t313 + t317;
t148 = -t226 * t314 + t320;
t136 = t243 * t225;
t135 = t167 * t225;
t133 = t226 * t215 + t216 * t313;
t132 = -t226 * t216 + t228 * t312;
t131 = t216 * t315 - t312;
t130 = t215 * t315 + t216 * t229;
t126 = -pkin(6) * t321 + t169;
t122 = mrSges(4,1) * t199 - mrSges(4,3) * t172;
t121 = -mrSges(4,2) * t199 + mrSges(4,3) * t171;
t120 = -pkin(6) * t285 + t151;
t117 = -t186 * t328 + t223 * t324;
t116 = -t186 * t222 - t223 * t271;
t98 = -pkin(4) * t243 - qJ(5) * t167 - t208;
t86 = qJD(2) * t236 - t167 * t302;
t85 = t222 * t403 - t270 * t302 - t271 * t304;
t78 = -qJD(3) * t127 + t310;
t77 = (-t225 * t305 - t228 * t303) * pkin(6) + t311;
t74 = pkin(4) * t135 - qJ(5) * t136 + t177;
t73 = -mrSges(4,2) * t165 + mrSges(4,3) * t101;
t72 = mrSges(4,1) * t165 - mrSges(4,3) * t100;
t64 = mrSges(5,1) * t105 + mrSges(5,2) * t244;
t63 = mrSges(6,1) * t105 - mrSges(6,3) * t244;
t61 = t228 * pkin(4) - t65;
t60 = -qJ(5) * t228 + t66;
t58 = -mrSges(4,1) * t101 + mrSges(4,2) * t100;
t36 = pkin(4) * t244 + qJ(5) * t105 + t351;
t32 = t100 * Ifges(4,4) + t101 * Ifges(4,2) + t165 * Ifges(4,6);
t29 = mrSges(5,1) * t165 - mrSges(5,3) * t55;
t28 = -mrSges(5,2) * t165 - mrSges(5,3) * t53;
t27 = -mrSges(6,2) * t53 + mrSges(6,3) * t165;
t26 = t328 * t79 - t335;
t25 = t222 * t79 + t75;
t20 = -pkin(4) * t85 - qJ(5) * t86 - qJD(5) * t136 + t123;
t11 = -pkin(4) * t306 - t15;
t10 = qJ(5) * t306 - qJD(5) * t228 + t16;
t6 = [(t176 * t217 + t401) * mrSges(3,3) + (-mrSges(3,1) * t217 + Ifges(3,5) * t225) * qJDD(2) + (t163 * mrSges(3,3) - t396 / 0.2e1 - t436 * t360 - Ifges(4,6) * t369 - Ifges(4,5) * t370 - Ifges(6,6) * t374 - Ifges(5,6) * t375 + (-Ifges(3,2) * t225 + t341) * t269 + Ifges(3,4) * t420 + Ifges(3,2) * t421 - t416 * t373 + t388 + (-mrSges(3,2) * pkin(6) + Ifges(3,6)) * qJDD(2)) * t228 + t427 * t136 + t425 * t86 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t401) + t404 * t214 + t405 * t304 / 0.2e1 + (-t150 * mrSges(4,1) - t149 * mrSges(4,2) + t440 * t318 + t424 * t309 + t423 * (pkin(3) * t317 + t208 * t313 - t223 * t318 + t309) + t419 * t226 - t395 * t133 - t393 * t132 + (-m(4) * t267 - t397) * t229) * g(2) - t32 * t322 / 0.2e1 + m(6) * (t1 * t60 + t10 * t22 + t11 * t21 + t2 * t61 + t20 * t31 + t5 * t74) + m(5) * (t118 * t123 + t15 * t23 + t16 * t24 + t177 * t70 + t3 * t65 + t4 * t66) + (t111 * t403 - t112 * t239 - t319 * t35 - t322 * t34) * mrSges(4,3) + t187 * (mrSges(4,1) * t239 - mrSges(4,2) * t403) + t383 * t135 - t184 * t292 + (m(4) * t140 * pkin(6) + Ifges(3,1) * t176 + Ifges(3,4) * t421 + t255 * t369 + t258 * t370) * t225 + (-t148 * mrSges(4,1) - t147 * mrSges(4,2) + t423 * (pkin(3) * t320 + t226 * t323 + t219) + t419 * t229 + t424 * t219 + t395 * t131 + t393 * t130 + (m(3) * pkin(1) - m(4) * t182 + t423 * (-pkin(1) - t189) + t397 - t400) * t226) * g(1) + (qJD(2) * t233 - t251 * t302 + t414 * t85) * t355 + t341 * t420 + t256 * t421 + m(4) * (t111 * t78 + t112 * t77 + t126 * t35 + t127 * t34 + t187 * t214) + t437 * t85 + t171 * (qJD(2) * t234 - t254 * t302) / 0.2e1 - t247 * t299 - (t224 * t94 + t227 * t93) * t302 / 0.2e1 + (t398 / 0.2e1 + t415 * t355 + t111 * mrSges(4,1) - t112 * mrSges(4,2) + Ifges(5,6) * t368 + Ifges(6,6) * t367 + t364 * t416 - t389) * t306 + t245 * t269 + qJD(2) ^ 2 * t253 / 0.2e1 + t60 * t27 + t61 * t30 + t20 * t63 + t65 * t29 + t66 * t28 + t58 * t217 + (qJD(2) * t235 - t257 * t302) * t357 + t284 * t371 + t319 * t378 + t74 * t12 + t16 * t81 + t15 * t82 + t11 * t83 + t10 * t84 + (-t414 * t135 + t225 * t252) * t360 + t77 * t121 + t78 * t122 + t123 * t64 + t264 * t327 + t126 * t72 + t127 * t73 + t140 * t159 - pkin(1) * (-mrSges(3,1) * t175 + mrSges(3,2) * t176) + t177 * t13 + Ifges(2,3) * qJDD(1); t429 * t124 + (t386 + t428) * t125 + t427 * t167 + t425 * t144 + (m(4) * ((-t111 * t227 - t112 * t224) * qJD(3) + t399) - t224 * t72 + t227 * t73 - t122 * t301 - t121 * t303) * pkin(7) + t399 * mrSges(4,3) + t402 * t64 - t404 * t213 - (-Ifges(3,2) * t308 + t211 + t405) * t307 / 0.2e1 - t398 * t308 / 0.2e1 + t199 * t187 * t261 + (t27 + t28) * t117 + t408 * t63 + t409 * t81 + t410 * t82 + (-t116 * t3 + t117 * t4 + t118 * t402 - t208 * t70 + t23 * t410 + t24 * t409) * m(5) + t412 * t83 + (t1 * t117 + t116 * t2 + t21 * t412 + t22 * t413 + t31 * t408 + t5 * t98) * m(6) + t413 * t84 + t140 * t262 + (-t111 * (mrSges(4,1) * t225 - mrSges(4,3) * t314) - t112 * (-mrSges(4,2) * t225 - mrSges(4,3) * t321)) * qJD(1) + (Ifges(5,6) * t367 + Ifges(6,6) * t368 + t356 * t415 + t365 * t416 + t389) * t308 + (-t338 + t94 / 0.2e1) * t301 + (-t225 * t291 - t264 + t423 * (t189 - t323) + (-m(6) * t250 - t435) * t228 + t400) * g(3) + t432 * (t263 + (-t423 * t223 - t291 + t440) * t228 + (m(5) * t208 - m(6) * (-t208 - t250) + t435) * t225) + (t30 - t29) * t116 + (-t337 + t371) * t303 + (-pkin(2) * t140 - t111 * t119 - t112 * t120 - t187 * t213) * m(4) + t93 * t286 / 0.2e1 + (-t414 * t355 - t437) * t143 - t253 * t299 / 0.2e1 + (t171 * t255 + t172 * t258 + t199 * t252) * qJD(3) / 0.2e1 - (t171 * t234 + t172 * t235 + t199 * t233) * qJD(1) / 0.2e1 + (-t245 / 0.2e1 + t247) * qJD(1) ^ 2 + t184 * t212 - pkin(2) * t58 + t251 * t360 + t254 * t369 + t257 * t370 + t224 * t378 + t98 * t12 - t120 * t121 - t119 * t122 - (-t360 * t414 + t383) * t243 - t163 * mrSges(3,2) - t164 * mrSges(3,1) + Ifges(3,6) * t175 + Ifges(3,5) * t176 - t208 * t13 + t227 * t32 / 0.2e1 + Ifges(3,3) * qJDD(2); -(t386 + t431) * t105 + t429 * t244 - t388 + (m(5) * t23 + t406) * t25 + (-m(5) * t24 + t407) * t26 + (-t21 * t25 - t31 * t36 + t1 * t205 + t2 * t207 - g(3) * (-t200 + (-pkin(4) * t215 + qJ(5) * t216) * t225) + (-t26 + qJD(5)) * t22) * m(6) + (-m(5) * t118 - t64) * t351 + t396 - t172 * (Ifges(4,1) * t171 - t336) / 0.2e1 + (-t148 * mrSges(4,2) + t130 * t395 - t131 * t393 + t147 * t439) * g(2) + (t150 * mrSges(4,2) + t395 * t132 - t393 * t133 - t149 * t439) * g(1) + t411 * t367 - (-Ifges(4,2) * t172 + t162 + t94) * t171 / 0.2e1 + (-(-t215 * mrSges(6,1) + t216 * mrSges(6,3)) * t225 + t159) * g(3) - (-mrSges(5,1) * t215 - mrSges(5,2) * t216 - t224 * t379) * t347 - t36 * t63 + t172 * t337 + t171 * t338 + t28 * t350 + (Ifges(4,5) * t171 - Ifges(4,6) * t172) * t356 + t93 * t357 + (t222 * t4 + t3 * t328) * t379 + qJD(5) * t84 - t111 * t121 + t112 * t122 + t29 * t288 - t187 * (mrSges(4,1) * t172 + mrSges(4,2) * t171) + t205 * t27 + t207 * t30; -t407 * t105 + t406 * t244 + t12 + t13 + (-t228 * g(3) + t225 * t432) * t423 + (t105 * t22 - t21 * t244 + t5) * m(6) + (t105 * t24 + t23 * t244 + t70) * m(5); t244 * t63 - t199 * t84 + (-g(1) * t132 - g(2) * t130 - t22 * t199 - t215 * t347 + t244 * t31 + t2) * m(6) + t30;];
tau = t6;
