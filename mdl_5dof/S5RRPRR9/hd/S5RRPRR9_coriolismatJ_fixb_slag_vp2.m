% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:07
% EndTime: 2019-12-31 20:20:19
% DurationCPUTime: 5.76s
% Computational Cost: add. (15736->431), mult. (31000->597), div. (0->0), fcn. (34407->8), ass. (0->241)
t363 = sin(pkin(9));
t364 = cos(pkin(9));
t409 = sin(qJ(2));
t410 = cos(qJ(2));
t245 = -t363 * t410 - t364 * t409;
t417 = -t245 / 0.2e1;
t387 = pkin(4) * qJD(4);
t324 = t409 * pkin(2);
t460 = m(4) * t324;
t274 = cos(qJ(5));
t272 = sin(qJ(5));
t243 = t363 * t409 - t364 * t410;
t266 = -pkin(2) * t410 - pkin(1);
t195 = t243 * pkin(3) + t245 * pkin(7) + t266;
t273 = sin(qJ(4));
t275 = cos(qJ(4));
t323 = t409 * pkin(6);
t254 = -qJ(3) * t409 - t323;
t325 = t410 * pkin(6);
t255 = qJ(3) * t410 + t325;
t443 = t363 * t254 + t364 * t255;
t113 = t273 * t195 + t275 * t443;
t349 = t245 * t273;
t88 = pkin(8) * t349 + t113;
t374 = t272 * t88;
t112 = t275 * t195 - t273 * t443;
t348 = t245 * t275;
t87 = pkin(8) * t348 + t112;
t69 = t243 * pkin(4) + t87;
t46 = t274 * t69 - t374;
t53 = t274 * t87 - t374;
t459 = t46 - t53;
t309 = t363 * pkin(2);
t262 = t309 + pkin(7);
t398 = pkin(8) + t262;
t234 = t398 * t273;
t235 = t398 * t275;
t194 = -t234 * t272 + t235 * t274;
t303 = -t274 * t234 - t235 * t272;
t295 = t272 * t275 + t274 * t273;
t237 = Ifges(6,6) * t295;
t441 = -t272 * t273 + t274 * t275;
t238 = Ifges(6,5) * t441;
t328 = t238 - t237;
t27 = -t194 * mrSges(6,1) - t303 * mrSges(6,2) + t328;
t458 = t27 * qJD(5);
t176 = t441 * t243;
t423 = -t176 / 0.2e1;
t174 = t295 * t243;
t425 = t174 / 0.2e1;
t292 = Ifges(6,5) * t423 + Ifges(6,6) * t425;
t196 = -t245 * pkin(3) + t243 * pkin(7) + t324;
t210 = -t364 * t254 + t255 * t363;
t116 = t275 * t196 + t210 * t273;
t350 = t243 * t275;
t70 = -t245 * pkin(4) + pkin(8) * t350 + t116;
t117 = t273 * t196 - t210 * t275;
t351 = t243 * t273;
t89 = pkin(8) * t351 + t117;
t48 = -t272 * t89 + t274 * t70;
t49 = t272 * t70 + t274 * t89;
t440 = t49 * mrSges(6,2) / 0.2e1 - t48 * mrSges(6,1) / 0.2e1 - t292;
t438 = Ifges(6,3) * t417 - t440;
t457 = -t245 * mrSges(4,1) - t243 * mrSges(4,2);
t173 = t441 * t245;
t386 = t173 * mrSges(6,3);
t129 = mrSges(6,1) * t243 + t386;
t449 = -t129 / 0.2e1;
t154 = -pkin(4) * t349 + t210;
t456 = t154 / 0.2e1;
t416 = t441 / 0.2e1;
t370 = t274 * t88;
t47 = t272 * t69 + t370;
t52 = -t272 * t87 - t370;
t455 = t47 + t52;
t377 = t441 * mrSges(6,3);
t236 = t295 * mrSges(6,1);
t305 = mrSges(6,2) * t441 + t236;
t451 = qJD(5) * t305;
t391 = Ifges(6,4) * t295;
t450 = t194 * mrSges(6,3) - t391;
t267 = Ifges(5,5) * t275;
t389 = Ifges(5,6) * t273;
t293 = -t267 / 0.2e1 + t389 / 0.2e1;
t446 = Ifges(4,4) + t293;
t268 = Ifges(5,4) * t275;
t257 = t273 * Ifges(5,1) + t268;
t270 = t273 ^ 2;
t271 = t275 ^ 2;
t444 = t270 + t271;
t393 = Ifges(5,4) * t273;
t256 = t275 * Ifges(5,2) + t393;
t412 = -t275 / 0.2e1;
t413 = t273 / 0.2e1;
t442 = t256 * t413 + t257 * t412;
t297 = Ifges(5,2) * t273 - t268;
t175 = t295 * t245;
t383 = t175 * mrSges(6,3);
t127 = -mrSges(6,2) * t243 + t383;
t439 = t127 * t416 + t295 * t449;
t382 = t176 * mrSges(6,2);
t384 = t174 * mrSges(6,1);
t330 = t384 / 0.2e1 + t382 / 0.2e1;
t437 = 2 * m(6);
t436 = m(5) / 0.2e1;
t435 = -m(6) / 0.2e1;
t434 = m(6) / 0.2e1;
t433 = pkin(4) / 0.2e1;
t432 = m(4) * pkin(2);
t431 = -mrSges(5,1) / 0.2e1;
t430 = mrSges(5,2) / 0.2e1;
t428 = t127 / 0.2e1;
t427 = -t173 / 0.2e1;
t426 = t173 / 0.2e1;
t424 = t175 / 0.2e1;
t422 = -t303 / 0.2e1;
t421 = t303 / 0.2e1;
t420 = -t194 / 0.2e1;
t419 = -t243 / 0.2e1;
t418 = t243 / 0.4e1;
t415 = t295 / 0.2e1;
t414 = -t273 / 0.2e1;
t411 = t275 / 0.2e1;
t407 = t275 * pkin(4);
t406 = t46 * mrSges(6,2);
t405 = t47 * mrSges(6,1);
t402 = t52 * mrSges(6,1);
t401 = t53 * mrSges(6,2);
t399 = Ifges(6,1) - Ifges(6,2);
t395 = mrSges(4,3) * t243;
t394 = mrSges(4,3) * t245;
t169 = Ifges(6,4) * t175;
t392 = Ifges(6,4) * t441;
t385 = t173 * Ifges(6,4);
t376 = t295 * mrSges(6,3);
t375 = t272 * t49;
t373 = t273 * mrSges(5,1);
t371 = t274 * t48;
t369 = t275 * mrSges(5,2);
t126 = mrSges(6,2) * t245 + t174 * mrSges(6,3);
t128 = -mrSges(6,1) * t245 + t176 * mrSges(6,3);
t140 = -Ifges(5,6) * t245 + t243 * t297;
t298 = Ifges(5,1) * t275 - t393;
t142 = -Ifges(5,5) * t245 - t243 * t298;
t153 = -pkin(4) * t351 + t443;
t299 = t369 + t373;
t189 = t299 * t243;
t190 = t299 * t245;
t199 = mrSges(5,2) * t245 + mrSges(5,3) * t351;
t322 = mrSges(5,3) * t349;
t200 = -mrSges(5,2) * t243 + t322;
t201 = -t245 * mrSges(5,1) + mrSges(5,3) * t350;
t202 = t243 * mrSges(5,1) + mrSges(5,3) * t348;
t143 = Ifges(5,5) * t243 - t245 * t298;
t335 = t275 * t143;
t141 = Ifges(5,6) * t243 + t245 * t297;
t340 = t273 * t141;
t73 = -Ifges(6,4) * t176 + Ifges(6,2) * t174 - Ifges(6,6) * t245;
t74 = t175 * Ifges(6,2) + t243 * Ifges(6,6) - t385;
t75 = -Ifges(6,1) * t176 + Ifges(6,4) * t174 - Ifges(6,5) * t245;
t76 = -t173 * Ifges(6,1) + t243 * Ifges(6,5) + t169;
t92 = -t382 - t384;
t93 = -mrSges(6,1) * t175 - mrSges(6,2) * t173;
t3 = t75 * t427 + t76 * t423 + t73 * t424 + t74 * t425 + (-t409 ^ 2 + t410 ^ 2) * Ifges(3,4) + (-Ifges(3,2) + Ifges(3,1)) * t410 * t409 + m(6) * (t153 * t154 + t46 * t48 + t47 * t49) + (mrSges(4,1) * t324 + t340 / 0.2e1 - t335 / 0.2e1 + t292 + t446 * t243) * t243 + (Ifges(6,5) * t426 - Ifges(6,6) * t175 / 0.2e1 + t140 * t413 + t142 * t412 - mrSges(4,2) * t324 - t446 * t245 + (-Ifges(5,3) - Ifges(6,3) + Ifges(4,1) - Ifges(4,2)) * t243) * t245 - t210 * t189 + t113 * t199 + t117 * t200 + t112 * t201 + t116 * t202 + t153 * t93 + t154 * t92 + m(5) * (t112 * t116 + t113 * t117 + t210 * t443) - pkin(1) * (mrSges(3,1) * t409 + mrSges(3,2) * t410) - t443 * t190 + t47 * t126 + t49 * t127 + t46 * t128 + t48 * t129 + (t457 + t460) * t266;
t368 = t3 * qJD(1);
t253 = -mrSges(5,1) * t275 + t273 * mrSges(5,2);
t329 = Ifges(6,5) * t175 + Ifges(6,6) * t173;
t91 = -mrSges(6,1) * t173 + mrSges(6,2) * t175;
t289 = t47 * t386 + t154 * t91 + t243 * t329 / 0.2e1;
t296 = Ifges(5,5) * t273 + Ifges(5,6) * t275;
t352 = t210 * t245;
t361 = t113 * t275;
t94 = t173 * Ifges(6,2) + t169;
t95 = t175 * Ifges(6,1) + t385;
t4 = pkin(4) * t93 * t348 - m(6) * (t46 * t52 + t47 * t53) - t53 * t127 - t52 * t129 - t253 * t352 + t113 * t202 + (-t74 / 0.2e1 + t95 / 0.2e1) * t173 + (-t94 / 0.2e1 - t76 / 0.2e1 + t46 * mrSges(6,3)) * t175 + (t143 * t414 + t141 * t412 + m(6) * t154 * t407 - mrSges(5,3) * t361 + t296 * t419 + (t256 * t414 + t257 * t411) * t245) * t245 - t289 + (-t200 + t322) * t112;
t367 = t4 * qJD(1);
t7 = -t129 * t47 + t74 * t426 + t95 * t427 + (t127 - t383) * t46 + t289 + (t76 + t94) * t424;
t366 = t7 * qJD(1);
t327 = m(6) * t433;
t282 = (t369 / 0.2e1 + t373 / 0.2e1) * t243 + (t174 * t274 - t176 * t272) * t327 + t330;
t334 = t275 * t200;
t339 = t273 * t202;
t283 = (-t459 * t295 + t455 * t441) * t435 + t339 / 0.2e1 - t334 / 0.2e1;
t354 = t175 * t441;
t355 = t173 * t295;
t8 = (-t355 / 0.2e1 + t354 / 0.2e1) * mrSges(6,3) + t282 + t283 - t439;
t365 = t8 * qJD(1);
t204 = -mrSges(6,1) * t441 + mrSges(6,2) * t295;
t310 = t364 * pkin(2);
t263 = -t310 - pkin(3);
t252 = t263 - t407;
t278 = (-t243 * t262 * t444 - t245 * t263) * t436 + (t174 * t303 - t176 * t194 - t245 * t252) * t434 + (-t243 * t363 + t245 * t364) * t432 / 0.2e1 - t376 * t425 + t377 * t423 + (t204 + t253) * t417 + t444 * mrSges(5,3) * t419;
t281 = (t116 * t275 + t273 * t117) * t436 + (t295 * t49 + t441 * t48) * t434 + t128 * t416 + t126 * t415 + t199 * t413 + t201 * t411 + t460 / 0.2e1;
t12 = t278 - t281 - t457;
t362 = qJD(1) * t12;
t13 = -t176 * t127 + t174 * t129 + (t190 - t93 + t394) * t245 + (-t334 + t339 + t395) * t243 + m(6) * (-t154 * t245 + t174 * t46 - t176 * t47) + m(5) * (-t352 + (t112 * t273 - t361) * t243) + m(4) * (-t243 * t443 - t352);
t358 = t13 * qJD(1);
t357 = t154 * t273;
t302 = t439 - (t354 - t355) * mrSges(6,3) / 0.2e1;
t16 = t302 - t330;
t356 = t16 * qJD(1);
t345 = t262 * t273;
t344 = t262 * t275;
t343 = t272 * t126;
t342 = t272 * t173;
t338 = t274 * t128;
t337 = t274 * t175;
t332 = t194 * t376 - t252 * t305;
t319 = -Ifges(6,1) / 0.4e1 + Ifges(6,2) / 0.4e1;
t318 = Ifges(5,2) / 0.4e1 - Ifges(5,1) / 0.4e1;
t317 = -t46 / 0.2e1 + t53 / 0.2e1;
t316 = t47 / 0.2e1 + t52 / 0.2e1;
t311 = t93 * t413;
t304 = t267 - t389;
t300 = t328 * t418;
t291 = t116 * t431 + t117 * t430;
t290 = t175 * t422 + t194 * t426;
t280 = -(t74 / 0.4e1 - t95 / 0.4e1 - t385 / 0.2e1 + t319 * t175) * t295 + (mrSges(6,2) * t456 + t94 / 0.4e1 + t76 / 0.4e1 + t169 / 0.2e1 + t319 * t173) * t441 + t236 * t456 + t252 * t91 / 0.2e1;
t276 = (t200 * t414 + t202 * t412) * t262 + (-t295 * t316 + t317 * t441 + t290) * mrSges(6,3) + t194 * t449 + t127 * t421 + t210 * t299 / 0.2e1 - t340 / 0.4e1 + t335 / 0.4e1 + t280;
t279 = (t270 / 0.2e1 + t271 / 0.2e1) * t262 * mrSges(5,3) + (t263 * t430 + t257 / 0.4e1 + t268 / 0.4e1 - t318 * t273) * t273 + (0.3e1 / 0.4e1 * t393 + t263 * t431 + t256 / 0.4e1 + t318 * t275 + (-t204 / 0.2e1 + t252 * t435) * pkin(4)) * t275;
t286 = -t459 * t194 + t455 * t303;
t1 = (-0.3e1 / 0.4e1 * t389 + 0.3e1 / 0.4e1 * t267 + t238 / 0.4e1 - t237 / 0.4e1) * t243 + t276 + t286 * t434 + (-t343 / 0.2e1 + t311 - t338 / 0.2e1 + (t357 / 0.4e1 - t375 / 0.4e1 - t371 / 0.4e1) * t437) * pkin(4) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + t279) * t245 + t291 + t440;
t22 = -t263 * t299 + t298 * t414 - t297 * t412 - t450 * t295 + (-t295 * t399 - t392) * t441 + t332 + t442 + (-m(6) * t252 - t204) * pkin(4) * t273;
t288 = -t1 * qJD(1) + t22 * qJD(2);
t23 = -Ifges(6,4) * t441 ^ 2 - (t399 * t441 + t450) * t295 + t332;
t277 = mrSges(6,3) * t290 + t129 * t420 + t303 * t428 + t280 + t300;
t6 = t277 - t438;
t287 = -t6 * qJD(1) + t23 * qJD(2);
t284 = (t274 * t428 + t272 * t449 + (t342 / 0.2e1 - t337 / 0.2e1) * mrSges(6,3)) * pkin(4);
t11 = -mrSges(6,1) * t316 + mrSges(6,2) * t317 + t284;
t251 = (mrSges(6,1) * t272 + mrSges(6,2) * t274) * pkin(4);
t29 = (t422 + t421) * mrSges(6,2) + (t420 + t194 / 0.2e1) * mrSges(6,1);
t285 = -t11 * qJD(1) - t29 * qJD(2) + t251 * qJD(4);
t246 = t251 * qJD(5);
t15 = t302 + t330;
t14 = t278 + t281;
t10 = -t406 / 0.2e1 - t405 / 0.2e1 - t401 / 0.2e1 + t402 / 0.2e1 + t284 + t329;
t9 = t282 - t283 + t302;
t5 = t277 + t438;
t2 = t276 + (pkin(4) * t357 + t286) * t434 + t279 * t245 + t304 * t418 + Ifges(5,3) * t417 + t300 + (t371 + t375) * t327 + pkin(4) * t311 - t291 + (t343 + t338) * t433 + t293 * t243 + t438;
t17 = [qJD(2) * t3 + qJD(3) * t13 - qJD(4) * t4 + qJD(5) * t7, t14 * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t368 + ((-Ifges(4,5) + t442) * t243 + t75 * t415 + t73 * t416 + t140 * t411 + t142 * t413 + t309 * t394 + t310 * t395 - t48 * t376 + t49 * t377 + t199 * t344 - t263 * t189 + t252 * t92 + Ifges(4,6) * t245 + t153 * t204 + t194 * t126 - mrSges(3,1) * t325 - t201 * t345 + mrSges(3,2) * t323 - Ifges(3,6) * t409 + Ifges(3,5) * t410 + (m(5) * t262 + mrSges(5,3)) * (-t116 * t273 + t117 * t275) + (Ifges(6,2) * t441 + t391) * t425 + (Ifges(6,5) * t295 + Ifges(6,6) * t441 + t296) * t417 + (m(5) * t263 - t364 * t432 - mrSges(4,1) + t253) * t443 + (-t363 * t432 + mrSges(4,2)) * t210 + (Ifges(6,1) * t295 + t392) * t423 + m(6) * (t153 * t252 + t194 * t49 + t303 * t48) + t303 * t128) * qJD(2), t358 + t14 * qJD(2) + (t174 * t441 - t176 * t295) * m(6) * qJD(3) + t9 * qJD(4) + t15 * qJD(5), -t367 + t2 * qJD(2) + t9 * qJD(3) + (-t113 * mrSges(5,1) - t112 * mrSges(5,2) + Ifges(5,5) * t349 + Ifges(5,6) * t348 + t329 - t401 + t402) * qJD(4) + t10 * qJD(5) + (m(6) * (t272 * t53 + t274 * t52) + (-t337 + t342) * mrSges(6,3)) * t387, t366 + t5 * qJD(2) + t15 * qJD(3) + t10 * qJD(4) + (t329 - t405 - t406) * qJD(5); qJD(3) * t12 + qJD(4) * t1 + qJD(5) * t6 - t368, -qJD(4) * t22 - qJD(5) * t23, t362, (-mrSges(5,1) * t344 + mrSges(5,2) * t345 + t27 + t304) * qJD(4) + t458 + (m(6) * (-t194 * t274 + t272 * t303) + (-t272 * t295 - t274 * t441) * mrSges(6,3)) * t387 - t288, t27 * qJD(4) - t287 + t458; -qJD(2) * t12 - qJD(4) * t8 + qJD(5) * t16 - t358, -t362, 0, -t365 + (-t299 - t305) * qJD(4) - t451 + (t272 * t441 - t274 * t295) * t437 * t387 / 0.2e1, -qJD(4) * t305 + t356 - t451; -qJD(2) * t1 + qJD(3) * t8 + qJD(5) * t11 + t367, t29 * qJD(5) + t288, t365, -t246, -t246 - t285; -qJD(2) * t6 - qJD(3) * t16 - qJD(4) * t11 - t366, -t29 * qJD(4) + t287, -t356, t285, 0;];
Cq = t17;
