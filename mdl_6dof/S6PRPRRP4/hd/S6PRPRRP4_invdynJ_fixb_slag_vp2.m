% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:22
% EndTime: 2019-03-08 20:09:47
% DurationCPUTime: 16.80s
% Computational Cost: add. (6819->617), mult. (16246->823), div. (0->0), fcn. (12769->14), ass. (0->273)
t221 = sin(pkin(11));
t224 = cos(pkin(11));
t228 = sin(qJ(4));
t231 = cos(qJ(4));
t185 = t221 * t228 - t231 * t224;
t223 = sin(pkin(6));
t232 = cos(qJ(2));
t322 = t223 * t232;
t235 = t185 * t322;
t358 = pkin(8) + qJ(3);
t198 = t358 * t221;
t199 = t358 * t224;
t406 = -t231 * t198 - t199 * t228;
t415 = qJD(1) * t235 - t185 * qJD(3) + qJD(4) * t406;
t180 = t185 * qJD(4);
t186 = t221 * t231 + t224 * t228;
t181 = t186 * qJD(4);
t229 = sin(qJ(2));
t316 = qJD(1) * t223;
t290 = t229 * t316;
t435 = pkin(4) * t181 + pkin(9) * t180 - t290;
t426 = Ifges(6,1) + Ifges(7,1);
t434 = -Ifges(6,4) + Ifges(7,5);
t425 = Ifges(7,4) + Ifges(6,5);
t424 = Ifges(6,6) - Ifges(7,6);
t423 = Ifges(6,3) + Ifges(7,2);
t433 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t178 = t185 * qJD(2);
t431 = qJD(5) + t178;
t416 = mrSges(4,3) * (t221 ^ 2 + t224 ^ 2);
t120 = -qJD(2) * t180 + qJDD(2) * t186;
t227 = sin(qJ(5));
t230 = cos(qJ(5));
t179 = t186 * qJD(2);
t244 = t230 * qJD(4) - t179 * t227;
t432 = qJD(5) * t244;
t57 = qJDD(4) * t227 + t120 * t230 + t432;
t383 = t57 / 0.2e1;
t141 = qJD(4) * t227 + t179 * t230;
t58 = qJD(5) * t141 - t230 * qJDD(4) + t120 * t227;
t381 = t58 / 0.2e1;
t121 = -qJD(2) * t181 - qJDD(2) * t185;
t111 = qJDD(5) - t121;
t380 = t111 / 0.2e1;
t213 = pkin(3) * t224 + pkin(2);
t114 = pkin(4) * t185 - pkin(9) * t186 - t213;
t135 = -t198 * t228 + t199 * t231;
t407 = t227 * t114 + t230 * t135;
t418 = -qJD(5) * t407 - t227 * t415 + t435 * t230;
t308 = qJD(5) * t230;
t309 = qJD(5) * t227;
t417 = t114 * t308 - t135 * t309 + t435 * t227 + t230 * t415;
t430 = -m(6) - m(7);
t429 = t425 * t111 + t426 * t57 + t434 * t58;
t428 = qJ(6) * t181 + qJD(6) * t185 + t417;
t427 = mrSges(6,3) + mrSges(7,2);
t22 = mrSges(6,1) * t58 + mrSges(6,2) * t57;
t422 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t120 + t22;
t33 = -mrSges(7,2) * t58 + mrSges(7,3) * t111;
t36 = -mrSges(6,2) * t111 - mrSges(6,3) * t58;
t357 = t33 + t36;
t34 = mrSges(6,1) * t111 - mrSges(6,3) * t57;
t35 = -t111 * mrSges(7,1) + t57 * mrSges(7,2);
t356 = t35 - t34;
t421 = t141 * t425 + t244 * t424 + t423 * t431;
t139 = Ifges(6,4) * t244;
t342 = t244 * Ifges(7,5);
t420 = t141 * t426 + t425 * t431 + t139 - t342;
t81 = mrSges(7,2) * t244 + mrSges(7,3) * t431;
t351 = mrSges(6,3) * t244;
t82 = -mrSges(6,2) * t431 + t351;
t355 = t81 + t82;
t350 = mrSges(6,3) * t141;
t83 = mrSges(6,1) * t431 - t350;
t84 = -mrSges(7,1) * t431 + mrSges(7,2) * t141;
t354 = t83 - t84;
t419 = -pkin(5) * t181 - t418;
t352 = mrSges(5,3) * t179;
t414 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t244 - mrSges(6,2) * t141 - t352;
t254 = pkin(5) * t227 - qJ(6) * t230;
t192 = qJD(2) * qJ(3) + t290;
t225 = cos(pkin(6));
t315 = qJD(1) * t225;
t209 = t224 * t315;
t343 = pkin(8) * qJD(2);
t130 = t209 + (-t192 - t343) * t221;
t144 = t224 * t192 + t221 * t315;
t131 = t224 * t343 + t144;
t64 = t130 * t228 + t131 * t231;
t413 = -qJD(6) * t227 + t254 * t431 - t64;
t302 = qJDD(2) * t224;
t303 = qJDD(2) * t221;
t184 = -mrSges(4,1) * t302 + mrSges(4,2) * t303;
t42 = -t121 * mrSges(5,1) + t120 * mrSges(5,2);
t412 = t184 + t42;
t411 = Ifges(5,5) * qJD(4);
t410 = Ifges(5,6) * qJD(4);
t246 = t144 * t224 - (-t192 * t221 + t209) * t221;
t409 = t232 * t246;
t60 = qJD(4) * pkin(9) + t64;
t253 = -t232 * t316 + qJD(3);
t166 = -qJD(2) * t213 + t253;
t71 = pkin(4) * t178 - pkin(9) * t179 + t166;
t26 = -t227 * t60 + t230 * t71;
t408 = qJD(6) - t26;
t267 = -mrSges(4,1) * t224 + mrSges(4,2) * t221;
t405 = -mrSges(5,1) * t178 - mrSges(5,2) * t179 - t267 * qJD(2);
t262 = mrSges(7,1) * t227 - mrSges(7,3) * t230;
t264 = mrSges(6,1) * t227 + mrSges(6,2) * t230;
t63 = t130 * t231 - t228 * t131;
t59 = -qJD(4) * pkin(4) - t63;
t30 = -pkin(5) * t244 - qJ(6) * t141 + t59;
t404 = t30 * t262 + t59 * t264;
t403 = -t227 * t424 + t230 * t425;
t345 = Ifges(7,5) * t227;
t347 = Ifges(6,4) * t227;
t402 = t230 * t426 + t345 - t347;
t240 = -t180 * t227 + t186 * t308;
t401 = t111 * t423 - t424 * t58 + t425 * t57;
t333 = t178 * t230;
t400 = t308 + t333;
t334 = t178 * t227;
t399 = -t309 - t334;
t283 = qJD(2) * t316;
t203 = t232 * t283;
t306 = qJDD(1) * t223;
t169 = t229 * t306 + t203;
t142 = t169 + t433;
t305 = qJDD(1) * t225;
t207 = t224 * t305;
t107 = -t142 * t221 + t207;
t108 = t224 * t142 + t221 * t305;
t398 = -t107 * t221 + t108 * t224;
t311 = qJD(4) * t231;
t312 = qJD(4) * t228;
t91 = t207 + (-pkin(8) * qJDD(2) - t142) * t221;
t92 = pkin(8) * t302 + t108;
t15 = t130 * t311 - t131 * t312 + t228 * t91 + t231 * t92;
t13 = qJDD(4) * pkin(9) + t15;
t202 = t229 * t283;
t168 = t232 * t306 - t202;
t243 = qJDD(3) - t168;
t132 = -qJDD(2) * t213 + t243;
t39 = -pkin(4) * t121 - pkin(9) * t120 + t132;
t3 = t230 * t13 + t227 * t39 + t71 * t308 - t309 * t60;
t27 = t227 * t71 + t230 * t60;
t4 = -qJD(5) * t27 - t13 * t227 + t230 * t39;
t397 = -t227 * t4 + t230 * t3;
t1 = qJ(6) * t111 + qJD(6) * t431 + t3;
t2 = -pkin(5) * t111 + qJDD(6) - t4;
t396 = t1 * t230 + t2 * t227;
t392 = Ifges(7,5) * t383 + Ifges(7,6) * t380 - t57 * Ifges(6,4) / 0.2e1 - t111 * Ifges(6,6) / 0.2e1 + (Ifges(7,3) + Ifges(6,2)) * t381;
t284 = m(4) * qJ(3) + mrSges(4,3);
t391 = -mrSges(5,3) + mrSges(3,2) - t284;
t255 = pkin(5) * t230 + qJ(6) * t227;
t263 = -t230 * mrSges(7,1) - t227 * mrSges(7,3);
t265 = mrSges(6,1) * t230 - mrSges(6,2) * t227;
t390 = -m(7) * t255 - mrSges(5,1) + t263 - t265;
t237 = m(4) * pkin(2) - t267;
t220 = pkin(11) + qJ(4);
t214 = sin(t220);
t215 = cos(t220);
t266 = t215 * mrSges(5,1) - t214 * mrSges(5,2);
t389 = t266 + mrSges(3,1) + t237;
t272 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t387 = -m(6) * t59 + t414;
t269 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t341 = cos(pkin(10));
t274 = t341 * t229;
t222 = sin(pkin(10));
t324 = t222 * t232;
t175 = t225 * t274 + t324;
t275 = t223 * t341;
t123 = t175 * t215 - t214 * t275;
t273 = t341 * t232;
t325 = t222 * t229;
t177 = -t225 * t325 + t273;
t326 = t222 * t223;
t125 = t177 * t215 + t214 * t326;
t323 = t223 * t229;
t158 = t214 * t225 + t215 * t323;
t386 = -g(1) * t125 - g(2) * t123 - g(3) * t158;
t385 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t382 = -t58 / 0.2e1;
t379 = t244 / 0.2e1;
t378 = -t244 / 0.2e1;
t377 = -t141 / 0.2e1;
t376 = t141 / 0.2e1;
t375 = -t431 / 0.2e1;
t373 = t178 / 0.2e1;
t372 = -t179 / 0.2e1;
t371 = t179 / 0.2e1;
t368 = t227 / 0.2e1;
t367 = pkin(4) * t215;
t364 = g(3) * t223;
t353 = mrSges(5,3) * t178;
t349 = Ifges(5,4) * t179;
t348 = Ifges(6,4) * t141;
t346 = Ifges(6,4) * t230;
t344 = Ifges(7,5) * t230;
t112 = pkin(4) * t179 + pkin(9) * t178;
t32 = t227 * t112 + t230 * t63;
t174 = -t225 * t273 + t325;
t336 = t174 * t214;
t176 = t225 * t324 + t274;
t335 = t176 * t214;
t328 = t215 * t227;
t327 = t215 * t230;
t321 = t230 * t180;
t320 = t230 * t232;
t319 = -t174 * t213 + t175 * t358;
t318 = -t176 * t213 + t177 * t358;
t314 = qJD(2) * t229;
t73 = -mrSges(7,1) * t244 - mrSges(7,3) * t141;
t299 = -t73 + t414;
t296 = m(4) + m(5) - t430;
t295 = t214 * t322;
t294 = t223 * t320;
t205 = t227 * t322;
t138 = Ifges(7,5) * t141;
t43 = Ifges(7,6) * t431 - Ifges(7,3) * t244 + t138;
t293 = t43 * t368;
t289 = t223 * t314;
t280 = -t309 / 0.2e1;
t279 = t308 / 0.2e1;
t259 = -Ifges(6,2) * t227 + t346;
t256 = Ifges(7,3) * t227 + t344;
t31 = t112 * t230 - t227 * t63;
t54 = t114 * t230 - t135 * t227;
t172 = -t221 * t323 + t224 * t225;
t173 = t221 * t225 + t224 * t323;
t245 = t231 * t172 - t173 * t228;
t96 = t172 * t228 + t173 * t231;
t16 = -t130 * t312 - t131 * t311 - t228 * t92 + t231 * t91;
t75 = t227 * t96 + t294;
t239 = t186 * t309 + t321;
t236 = t186 * t322;
t14 = -qJDD(4) * pkin(4) - t16;
t86 = qJD(3) * t186 + qJD(4) * t135;
t233 = qJD(2) ^ 2;
t194 = -pkin(4) - t255;
t189 = t213 * t322;
t188 = -qJD(2) * pkin(2) + t253;
t167 = Ifges(5,4) * t178;
t157 = -t214 * t323 + t215 * t225;
t155 = -qJD(4) * mrSges(5,2) - t353;
t153 = -qJDD(2) * pkin(2) + t243;
t126 = t158 * t227 + t294;
t124 = -t177 * t214 + t215 * t326;
t122 = -t175 * t214 - t215 * t275;
t101 = t179 * Ifges(5,1) - t167 + t411;
t100 = -t178 * Ifges(5,2) + t349 + t410;
t99 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t121;
t76 = t230 * t96 - t205;
t72 = pkin(5) * t141 - qJ(6) * t244;
t69 = t125 * t227 - t176 * t230;
t67 = t123 * t227 - t174 * t230;
t65 = t186 * t254 - t406;
t62 = qJD(2) * t236 + qJD(4) * t96;
t61 = -qJD(2) * t235 + qJD(4) * t245;
t46 = Ifges(6,2) * t244 + Ifges(6,6) * t431 + t348;
t41 = -pkin(5) * t185 - t54;
t40 = qJ(6) * t185 + t407;
t29 = -pkin(5) * t179 - t31;
t28 = qJ(6) * t179 + t32;
t25 = -qJD(5) * t205 + t227 * t61 - t230 * t289 + t308 * t96;
t24 = -qJD(5) * t75 + t227 * t289 + t230 * t61;
t23 = -t254 * t180 + (qJD(5) * t255 - qJD(6) * t230) * t186 + t86;
t21 = mrSges(7,1) * t58 - mrSges(7,3) * t57;
t20 = qJ(6) * t431 + t27;
t19 = -pkin(5) * t431 + t408;
t5 = pkin(5) * t58 - qJ(6) * t57 - qJD(6) * t141 + t14;
t6 = [t61 * t155 + t96 * t99 + t357 * t76 + t356 * t75 - t354 * t25 + t355 * t24 + (-t172 * t221 + t173 * t224) * qJDD(2) * mrSges(4,3) - (t21 + t422) * t245 - t299 * t62 + (-m(2) - m(3) - t296) * g(3) + m(5) * (t15 * t96 + t16 * t245 + t61 * t64 - t62 * t63) + m(4) * (t107 * t172 + t108 * t173) + m(7) * (t1 * t76 + t19 * t25 + t2 * t75 + t20 * t24 - t245 * t5 + t30 * t62) + m(6) * (-t14 * t245 + t24 * t27 - t25 * t26 + t3 * t76 - t4 * t75 + t59 * t62) + ((mrSges(3,1) * qJDD(2) - t412) * t232 + (-qJDD(2) * mrSges(3,2) - qJD(2) * t405) * t229 + m(5) * (-t132 * t232 + t166 * t314) + m(4) * (qJD(2) * t409 - t153 * t232 + t188 * t314) + m(3) * (t168 * t232 + t169 * t229) + (-t229 * mrSges(3,1) + (-mrSges(3,2) + t416) * t232) * t233) * t223 + (m(3) * t225 ^ 2 + m(2)) * qJDD(1); (t229 * t364 - t169 + t203) * mrSges(3,2) + (-m(5) * t318 + t430 * (-pkin(9) * t335 - t176 * t367 + t318) - t272 * (-t176 * t327 + t177 * t227) + t269 * (-t176 * t328 - t177 * t230) + t427 * t335 + t391 * t177 + t389 * t176) * g(1) + (-m(5) * t319 + t430 * (-pkin(9) * t336 - t174 * t367 + t319) - t272 * (-t174 * t327 + t175 * t227) + t269 * (-t174 * t328 - t175 * t230) + t427 * t336 + t391 * t175 + t389 * t174) * g(2) + t428 * t81 + (t1 * t40 + t19 * t419 + t2 * t41 + t20 * t428 + t23 * t30 + t5 * t65) * m(7) + ((mrSges(7,2) * t2 - mrSges(6,3) * t4 + t429 / 0.2e1) * t230 + t420 * t280 + (-t1 * mrSges(7,2) - t3 * mrSges(6,3) + t392) * t227 + t132 * mrSges(5,2) - t16 * mrSges(5,3) + Ifges(5,1) * t120 + Ifges(5,4) * t121 + Ifges(5,5) * qJDD(4) + t14 * t264 + t256 * t381 + t259 * t382 + t262 * t5 + t279 * t43 + t403 * t380 + t402 * t383) * t186 + (-t203 + t433) * t416 - t422 * t406 + (-t14 * t406 + t26 * t418 + t27 * t417 + t3 * t407 + t4 * t54 + t59 * t86) * m(6) + (-t132 * t213 + t135 * t15 + t16 * t406 - t166 * t290 + t415 * t64 - t63 * t86) * m(5) + (t181 * t423 - t239 * t425 - t240 * t424) * t431 / 0.2e1 + (t425 * t181 - t426 * t239 + t240 * t434) * t376 + t27 * (-mrSges(6,2) * t181 - mrSges(6,3) * t240) + t20 * (-mrSges(7,2) * t240 + mrSges(7,3) * t181) + t26 * (mrSges(6,1) * t181 + mrSges(6,3) * t239) + t19 * (-mrSges(7,1) * t181 - mrSges(7,2) * t239) + t30 * (mrSges(7,1) * t240 + mrSges(7,3) * t239) + t59 * (mrSges(6,1) * t240 - mrSges(6,2) * t239) + (-Ifges(7,5) * t239 + Ifges(7,6) * t181 + Ifges(7,3) * t240) * t378 + (-Ifges(6,4) * t239 - Ifges(6,2) * t240 + Ifges(6,6) * t181) * t379 + (-Ifges(5,1) * t180 - Ifges(5,4) * t181) * t371 + (t180 * t63 - t181 * t64) * mrSges(5,3) + qJD(4) * (-Ifges(5,5) * t180 - Ifges(5,6) * t181) / 0.2e1 - t178 * (-Ifges(5,4) * t180 - Ifges(5,2) * t181) / 0.2e1 + t166 * (mrSges(5,1) * t181 - mrSges(5,2) * t180) + t417 * t82 + t418 * t83 + t419 * t84 - t420 * t321 / 0.2e1 + t421 * t181 / 0.2e1 + (m(5) * t63 - m(7) * t30 + t387 - t73) * qJD(1) * t236 + (-(t188 * t229 + t409) * t316 - pkin(2) * t153 + t246 * qJD(3) + t398 * qJ(3)) * m(4) + t407 * t36 + (-m(5) * t189 - t427 * t295 + t430 * (pkin(9) * t295 + t322 * t367 + t323 * t358 + t189) + t269 * (t205 * t215 - t230 * t323) + (-t266 * t232 - (m(5) * t358 + mrSges(5,3)) * t229 - t272 * (t215 * t320 + t227 * t229)) * t223) * g(3) + (-t232 * t364 + t168 + t202) * mrSges(3,1) - t240 * t46 / 0.2e1 + t405 * t290 + t398 * mrSges(4,3) + t153 * t267 + (Ifges(4,4) * t221 + Ifges(4,2) * t224) * t302 + (Ifges(4,1) * t221 + Ifges(4,4) * t224) * t303 + (t132 * mrSges(5,1) - t15 * mrSges(5,3) - Ifges(5,4) * t120 - Ifges(5,2) * t121 - Ifges(5,6) * qJDD(4) + Ifges(6,6) * t382 + Ifges(7,6) * t381 + t423 * t380 + t425 * t383 + t385 + t401 / 0.2e1) * t185 - t414 * t86 + t415 * t155 + Ifges(3,3) * qJDD(2) - (t229 * t284 + t232 * t237) * t364 - t213 * t42 - pkin(2) * t184 - t181 * t100 / 0.2e1 - t180 * t101 / 0.2e1 + t135 * t99 + t23 * t73 + t65 * t21 + t54 * t34 + t41 * t35 + t40 * t33 - t180 * t293; t178 * t155 - t233 * t416 + t299 * t179 + (t355 * t431 - t356) * t230 + (-t354 * t431 + t357) * t227 + (t1 * t227 - t179 * t30 - t2 * t230 + t431 * (t19 * t227 + t20 * t230)) * m(7) + (-t179 * t59 + t3 * t227 + t4 * t230 + t431 * (-t26 * t227 + t27 * t230)) * m(6) + (t178 * t64 + t179 * t63 + t132) * m(5) + (-qJD(2) * t246 + t153) * m(4) + (-g(1) * t176 - g(2) * t174 + g(3) * t322) * t296 + t412; -t354 * pkin(9) * t308 - t355 * pkin(9) * t309 + t356 * pkin(9) * t227 + t357 * pkin(9) * t230 + t429 * t368 + (mrSges(5,2) * t158 + t430 * (t157 * pkin(4) + pkin(9) * t158) + t390 * t157) * g(3) + (mrSges(5,2) * t125 + t430 * (t124 * pkin(4) + pkin(9) * t125) + t390 * t124) * g(1) + (mrSges(5,2) * t123 + t430 * (t122 * pkin(4) + pkin(9) * t123) + t390 * t122) * g(2) + (Ifges(6,2) * t382 - Ifges(7,3) * t381 + t380 * t424 - t392) * t230 + (-Ifges(5,2) * t373 + Ifges(6,6) * t378 + Ifges(7,6) * t379 + t27 * mrSges(6,2) - t20 * mrSges(7,3) - t26 * mrSges(6,1) + t19 * mrSges(7,1) - t166 * mrSges(5,1) + t410 / 0.2e1 + t425 * t377 + t423 * t375) * t179 + (t380 * t425 + t383 * t426) * t227 + (t352 + t387) * t64 + (-t353 - t155) * t63 - (-t259 / 0.2e1 + t256 / 0.2e1) * t432 + t100 * t371 - (t46 / 0.2e1 - t43 / 0.2e1) * t334 + (t141 * t402 + t403 * t431) * qJD(5) / 0.2e1 + t345 * t381 + t347 * t382 + (-t349 + t421) * t372 + (-t344 + t346) * t383 + (-t167 + t101) * t373 + (t19 * t400 + t20 * t399 + t386 + t396) * mrSges(7,2) + (-t26 * t400 + t27 * t399 + t386 + t397) * mrSges(6,3) + (t293 + t404) * qJD(5) + (-t26 * t31 - t27 * t32 - pkin(4) * t14 + ((-t27 * t227 - t26 * t230) * qJD(5) + t397) * pkin(9)) * m(6) - t14 * t265 + t5 * t263 + (t333 / 0.2e1 + t279) * t420 - (Ifges(5,1) * t372 + t259 * t378 + t256 * t379 - t166 * mrSges(5,2) - t411 / 0.2e1 + t402 * t377 + t403 * t375 - t404) * t178 + t413 * t73 + (-t19 * t29 - t20 * t28 + t194 * t5 + ((t19 * t230 - t20 * t227) * qJD(5) + t396) * pkin(9) + t413 * t30) * m(7) + Ifges(5,3) * qJDD(4) + t194 * t21 + Ifges(5,5) * t120 + Ifges(5,6) * t121 - t28 * t81 - t32 * t82 - t31 * t83 - t29 * t84 - pkin(4) * t22 + t16 * mrSges(5,1) - t15 * mrSges(5,2) + t46 * t280; (t141 * t20 - t19 * t244) * mrSges(7,2) - t59 * (mrSges(6,1) * t141 + mrSges(6,2) * t244) - t30 * (mrSges(7,1) * t141 - mrSges(7,3) * t244) + (-t141 * t424 + t244 * t425) * t375 + (t244 * t426 + t138 - t348 + t43) * t377 + (Ifges(7,3) * t141 + t342) * t379 + (t351 - t355) * t26 + (-Ifges(6,2) * t141 + t139 + t420) * t378 + t46 * t376 + t385 + (t350 + t354) * t27 + qJD(6) * t81 - t72 * t73 + qJ(6) * t33 - pkin(5) * t35 + (t269 * (t123 * t230 + t174 * t227) + t272 * t67) * g(2) + (t269 * (t125 * t230 + t176 * t227) + t272 * t69) * g(1) + (t269 * (t158 * t230 - t205) + t272 * t126) * g(3) + (-pkin(5) * t2 + qJ(6) * t1 - t19 * t27 + t20 * t408 - t30 * t72) * m(7) + t401; t141 * t73 - t431 * t81 + (-g(1) * t69 - g(2) * t67 - g(3) * t126 + t30 * t141 - t20 * t431 + t2) * m(7) + t35;];
tau  = t6;
