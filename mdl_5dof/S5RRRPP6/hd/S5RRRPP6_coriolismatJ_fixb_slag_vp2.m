% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:00:02
% EndTime: 2019-12-31 21:00:18
% DurationCPUTime: 7.02s
% Computational Cost: add. (8209->473), mult. (18321->638), div. (0->0), fcn. (17595->6), ass. (0->257)
t274 = sin(qJ(3));
t417 = t274 / 0.2e1;
t370 = sin(pkin(8));
t339 = t370 * t274;
t371 = cos(pkin(8));
t411 = cos(qJ(3));
t227 = -t371 * t411 + t339;
t262 = -pkin(3) * t411 - pkin(2);
t341 = t371 * t274;
t284 = t370 * t411 + t341;
t122 = t227 * pkin(4) - qJ(5) * t284 + t262;
t144 = mrSges(6,1) * t227 - mrSges(6,3) * t284;
t479 = m(6) * t122 + t144;
t145 = mrSges(5,1) * t227 + mrSges(5,2) * t284;
t478 = m(5) * t262 + t145;
t402 = Ifges(4,4) * t274;
t295 = -Ifges(4,2) * t411 - t402;
t269 = Ifges(4,4) * t411;
t245 = Ifges(4,1) * t274 + t269;
t358 = t411 / 0.2e1;
t317 = t245 * t358;
t477 = t295 * t417 + t317;
t474 = Ifges(5,1) + Ifges(6,1);
t473 = Ifges(6,4) + Ifges(5,5);
t472 = Ifges(5,6) - Ifges(6,6);
t275 = sin(qJ(2));
t352 = t275 * t411;
t320 = -t352 / 0.2e1;
t468 = mrSges(5,3) + mrSges(6,2);
t205 = t284 * t275;
t194 = Ifges(5,4) * t205;
t204 = t227 * t275;
t276 = cos(qJ(2));
t397 = Ifges(6,5) * t205;
t465 = -t474 * t204 - t473 * t276 - t194 + t397;
t206 = t284 * t276;
t207 = t276 * t227;
t464 = t473 * t275 - t474 * t207 + (-Ifges(5,4) + Ifges(6,5)) * t206;
t224 = Ifges(5,4) * t227;
t396 = Ifges(6,5) * t227;
t463 = t474 * t284 - t224 + t396;
t268 = t276 * mrSges(6,3);
t386 = t205 * mrSges(6,2);
t161 = -t268 - t386;
t385 = t205 * mrSges(5,3);
t162 = mrSges(5,2) * t276 - t385;
t462 = -t161 - t162;
t387 = t204 * mrSges(5,3);
t164 = -mrSges(5,1) * t276 + t387;
t388 = t204 * mrSges(6,2);
t165 = mrSges(6,1) * t276 - t388;
t461 = t164 - t165;
t460 = -Ifges(4,2) * t274 + t269;
t459 = -t473 * t227 - t472 * t284;
t293 = t411 * mrSges(4,1) - t274 * mrSges(4,2);
t457 = t472 * t204 - t473 * t205;
t456 = t164 / 0.2e1 - t165 / 0.2e1;
t455 = -t162 / 0.2e1 - t161 / 0.2e1;
t223 = Ifges(6,5) * t284;
t146 = Ifges(6,3) * t227 + t223;
t400 = Ifges(5,4) * t284;
t454 = -t474 * t227 + t146 + t223 - t400;
t453 = -Ifges(5,2) * t284 - t224 + t463;
t452 = Ifges(5,2) * t204 - t194 + t465;
t191 = Ifges(6,5) * t204;
t101 = -t276 * Ifges(6,6) + Ifges(6,3) * t205 - t191;
t401 = Ifges(5,4) * t204;
t451 = -t474 * t205 + t101 - t191 + t401;
t348 = t370 * pkin(3);
t256 = t348 + qJ(5);
t234 = m(6) * t256 + mrSges(6,3);
t267 = t411 * qJ(4);
t272 = t411 * pkin(7);
t362 = t272 + t267;
t405 = -qJ(4) - pkin(7);
t155 = -t405 * t341 + t362 * t370;
t450 = t405 * t339 + t362 * t371;
t349 = t371 * pkin(3);
t261 = -t349 - pkin(4);
t442 = m(5) * pkin(3);
t448 = m(6) * t261 - t371 * t442 - mrSges(5,1) - mrSges(6,1);
t447 = t370 * t442 - mrSges(5,2) + t234;
t445 = -m(5) / 0.2e1;
t444 = -m(6) / 0.2e1;
t443 = m(6) / 0.2e1;
t441 = -mrSges(5,2) / 0.2e1;
t440 = mrSges(6,3) / 0.2e1;
t408 = t275 * pkin(2);
t248 = -pkin(7) * t276 + t408;
t367 = t274 * t275;
t195 = pkin(6) * t367 + t411 * t248;
t141 = t275 * pkin(3) - t267 * t276 + t195;
t196 = -pkin(6) * t352 + t274 * t248;
t366 = t274 * t276;
t153 = -qJ(4) * t366 + t196;
t55 = t141 * t371 - t153 * t370;
t49 = -t275 * pkin(4) - t55;
t439 = -t49 / 0.2e1;
t407 = t275 * pkin(7);
t243 = -pkin(2) * t276 - pkin(1) - t407;
t351 = t276 * t411;
t324 = pkin(6) * t351;
t152 = t324 + (-qJ(4) * t275 + t243) * t274;
t139 = t371 * t152;
t231 = t411 * t243;
t301 = -t267 * t275 + t231;
t360 = pkin(6) * t366;
t151 = t301 - t360;
t57 = t151 * t370 + t139;
t438 = m(6) * t57;
t433 = -t204 / 0.2e1;
t432 = t204 / 0.2e1;
t430 = -t205 / 0.2e1;
t429 = t205 / 0.2e1;
t427 = -t206 / 0.2e1;
t426 = t206 / 0.2e1;
t425 = -t207 / 0.2e1;
t424 = -t227 / 0.2e1;
t423 = t227 / 0.2e1;
t421 = t284 / 0.2e1;
t420 = -t284 / 0.2e1;
t416 = -t275 / 0.2e1;
t415 = t275 / 0.2e1;
t414 = -t276 / 0.2e1;
t413 = t276 / 0.2e1;
t410 = pkin(3) * t274;
t409 = pkin(7) * t274;
t271 = t275 * pkin(6);
t273 = t276 * pkin(6);
t403 = Ifges(3,4) * t275;
t399 = Ifges(6,4) * t207;
t398 = Ifges(5,5) * t207;
t394 = Ifges(6,2) * t275;
t393 = Ifges(5,6) * t206;
t392 = Ifges(6,6) * t206;
t391 = Ifges(4,3) * t275;
t390 = Ifges(5,3) * t275;
t136 = (-pkin(6) * t274 - pkin(3)) * t276 + t301;
t54 = t370 * t136 + t139;
t46 = -qJ(5) * t276 + t54;
t340 = t370 * t152;
t53 = t136 * t371 - t340;
t47 = t276 * pkin(4) - t53;
t9 = t462 * t205 + t461 * t204 + m(5) * (t204 * t53 - t205 * t54) + m(6) * (-t204 * t47 - t205 * t46);
t389 = qJD(1) * t9;
t384 = t206 * mrSges(5,1);
t383 = t206 * mrSges(6,1);
t382 = t207 * mrSges(5,2);
t381 = t207 * mrSges(6,2);
t380 = t207 * mrSges(6,3);
t379 = t227 * mrSges(6,2);
t378 = t227 * mrSges(5,3);
t377 = t284 * mrSges(6,2);
t376 = t284 * mrSges(5,3);
t374 = t275 * mrSges(6,1);
t102 = -Ifges(6,5) * t207 + Ifges(6,6) * t275 + Ifges(6,3) * t206;
t103 = -Ifges(5,2) * t205 - t276 * Ifges(5,6) - t401;
t104 = -Ifges(5,4) * t207 - Ifges(5,2) * t206 + Ifges(5,6) * t275;
t110 = mrSges(6,1) * t205 + mrSges(6,3) * t204;
t111 = mrSges(5,1) * t205 - mrSges(5,2) * t204;
t112 = t380 + t383;
t113 = -t382 + t384;
t160 = -mrSges(6,2) * t206 + mrSges(6,3) * t275;
t163 = -mrSges(5,2) * t275 - mrSges(5,3) * t206;
t166 = mrSges(5,1) * t275 + mrSges(5,3) * t207;
t167 = -t374 - t381;
t185 = t231 - t360;
t186 = t274 * t243 + t324;
t199 = -Ifges(4,6) * t276 + t275 * t460;
t200 = Ifges(4,6) * t275 + t276 * t460;
t299 = Ifges(4,1) * t411 - t402;
t289 = t299 * t275;
t201 = -Ifges(4,5) * t276 + t289;
t202 = Ifges(4,5) * t275 + t276 * t299;
t292 = t274 * mrSges(4,1) + mrSges(4,2) * t411;
t225 = t292 * t276;
t357 = mrSges(4,3) * t367;
t237 = mrSges(4,2) * t276 - t357;
t238 = -t275 * mrSges(4,2) - mrSges(4,3) * t366;
t321 = mrSges(4,3) * t352;
t239 = -t276 * mrSges(4,1) - t321;
t240 = t275 * mrSges(4,1) - mrSges(4,3) * t351;
t241 = pkin(3) * t367 + t271;
t242 = pkin(3) * t366 + t273;
t287 = t292 * t271;
t296 = Ifges(4,5) * t411 - Ifges(4,6) * t274;
t290 = t276 * t296;
t318 = t351 / 0.2e1;
t319 = t352 / 0.2e1;
t343 = -t366 / 0.2e1;
t345 = -t367 / 0.2e1;
t365 = t275 * t276;
t56 = t370 * t141 + t371 * t153;
t48 = qJ(5) * t275 + t56;
t78 = pkin(4) * t205 + qJ(5) * t204 + t241;
t79 = pkin(4) * t206 + qJ(5) * t207 + t242;
t3 = (t275 * t296 - t403 - t472 * t205 - t473 * t204 + (Ifges(3,1) - Ifges(6,2) - Ifges(4,3) - Ifges(5,3)) * t276) * t415 + t102 * t429 + t104 * t430 + (Ifges(3,2) * t276 + t403) * t416 + t101 * t426 + t103 * t427 + t225 * t271 + t276 * t287 + t464 * t433 + t465 * t425 + m(4) * (pkin(6) ^ 2 * t365 + t185 * t195 + t186 * t196) + t200 * t345 + t199 * t343 + t201 * t318 + t202 * t319 + (t392 + t394 - t399 + t390 - t393 - t398 + t290 + t391) * t414 + (0.2e1 * Ifges(3,4) * t276 + (Ifges(3,1) - Ifges(3,2)) * t275) * t413 - pkin(1) * (t275 * mrSges(3,1) + mrSges(3,2) * t276) + t195 * t239 + t185 * t240 + t241 * t113 + t242 * t111 + t196 * t237 + t186 * t238 + t47 * t167 + t46 * t160 + t48 * t161 + t56 * t162 + t54 * t163 + t55 * t164 + t49 * t165 + t53 * t166 + t79 * t110 + t78 * t112 + m(5) * (t241 * t242 + t53 * t55 + t54 * t56) + m(6) * (t46 * t48 + t47 * t49 + t78 * t79);
t373 = t3 * qJD(1);
t294 = -Ifges(4,5) * t274 - Ifges(4,6) * t411;
t311 = -Ifges(6,3) * t204 - t397;
t325 = pkin(3) * t352;
t336 = -t204 * mrSges(5,1) - t205 * mrSges(5,2);
t338 = -t204 * mrSges(6,1) + t205 * mrSges(6,3);
t58 = t151 * t371 - t340;
t95 = -t204 * pkin(4) + t205 * qJ(5) + t325;
t4 = t311 * t429 + t103 * t432 + t53 * t385 + t54 * t387 + t46 * t388 - t47 * t386 + t111 * t325 - t294 * t365 / 0.2e1 + t201 * t345 + t199 * t320 + t95 * t110 + (m(6) * t95 + t338) * t78 + (m(5) * t325 + t336) * t241 + (m(5) * t54 + m(6) * t46 - t462) * t58 + (-m(5) * t53 + m(6) * t47 - t461) * t57 + t457 * t414 + (-t239 - t321) * t186 + (t357 + t237) * t185 + (pkin(6) * t293 - t477) * t275 ^ 2 + t451 * t433 + t452 * t430;
t372 = t4 * qJD(1);
t17 = t276 * t161 - m(6) * (t78 * t204 - t276 * t46) - t204 * t110;
t369 = qJD(1) * t17;
t361 = t442 / 0.2e1;
t356 = mrSges(6,2) / 0.2e1 + mrSges(5,3) / 0.2e1;
t342 = t227 * t413;
t337 = mrSges(6,1) * t284 + t227 * mrSges(6,3);
t335 = mrSges(5,1) * t284 - t227 * mrSges(5,2);
t328 = -t155 * t204 - t205 * t450;
t323 = mrSges(5,3) * t349;
t322 = mrSges(5,3) * t348;
t310 = Ifges(6,3) * t284 - t396;
t128 = pkin(4) * t284 + qJ(5) * t227 + t410;
t147 = -Ifges(5,2) * t227 + t400;
t308 = t155 * t57 + t450 * t58;
t277 = t468 * (t155 * t429 + t57 * t420 + t58 * t423 + t433 * t450) + t459 * t276 / 0.4e1 + (pkin(3) * t145 + t295) * t320 - t287 / 0.2e1 + (t122 * t95 + t128 * t78 + t308) * t444 + ((t241 * t274 + t262 * t352) * pkin(3) + t308) * t445 - (t444 * t46 + t445 * t54 + t455) * t155 + (t460 + 0.2e1 * t245) * t367 / 0.4e1 + (t452 / 0.4e1 - t311 / 0.4e1) * t227 + (t453 / 0.4e1 - t310 / 0.4e1) * t205 + (t454 / 0.4e1 - t147 / 0.4e1) * t204 - t111 * t410 / 0.2e1 + t293 * t408 / 0.2e1 + t237 * t409 / 0.2e1 + t46 * t377 / 0.2e1 - t53 * t378 / 0.2e1 + t47 * t379 / 0.2e1 + t54 * t376 / 0.2e1 + t239 * t272 / 0.2e1 - t78 * t337 / 0.2e1 - t122 * t338 / 0.2e1 - t241 * t335 / 0.2e1 - t262 * t336 / 0.2e1 + t290 / 0.4e1 - (t289 + t201) * t411 / 0.4e1 + t274 * t199 / 0.4e1 - t95 * t144 / 0.2e1 - t128 * t110 / 0.2e1 + (t274 ^ 2 + t411 ^ 2) * mrSges(4,3) * t407 / 0.2e1 + (t444 * t47 - t445 * t53 + t456) * t450 + (-t451 / 0.4e1 + t103 / 0.4e1) * t284;
t278 = (t256 * t48 + t261 * t49) * t443 - t399 / 0.2e1 - t398 / 0.2e1 + t394 / 0.2e1 - t393 / 0.2e1 + t392 / 0.2e1 + t391 / 0.2e1 + t390 / 0.2e1 + t195 * mrSges(4,1) / 0.2e1 - t196 * mrSges(4,2) / 0.2e1 + t256 * t160 / 0.2e1 + t261 * t167 / 0.2e1 + t48 * t440 + mrSges(6,1) * t439 + t55 * mrSges(5,1) / 0.2e1 + t56 * t441 + (t370 * t56 + t371 * t55) * t361 + Ifges(4,5) * t318 + Ifges(4,6) * t343 + t166 * t349 / 0.2e1 + t163 * t348 / 0.2e1;
t1 = t278 + t277;
t5 = -pkin(2) * t292 + t122 * t337 + t479 * t128 + t147 * t420 + t262 * t335 + t299 * t417 + t310 * t423 + t460 * t358 + t478 * t410 + t454 * t421 + t453 * t424 + t477;
t309 = -t1 * qJD(1) + t5 * qJD(2);
t15 = (t227 ^ 2 + t284 ^ 2) * t468 + (m(5) + m(6)) * (t155 * t284 - t227 * t450);
t279 = (t205 * t356 + t455) * t227 - (t204 * t356 + t456) * t284 + m(5) * (-t227 * t54 - t284 * t53 + t328) / 0.2e1 + (-t227 * t46 + t284 * t47 + t328) * t443;
t303 = t242 * t445 + t444 * t79;
t7 = -(t440 + t441) * t207 + (-mrSges(6,1) / 0.2e1 - mrSges(5,1) / 0.2e1) * t206 + t279 + t303;
t307 = qJD(1) * t7 + qJD(2) * t15;
t282 = (t122 * t204 - t276 * t450 - t284 * t78) * t443 + t144 * t432 + t110 * t420;
t300 = m(6) * t439 + t374 / 0.2e1;
t12 = (t342 + t207 / 0.2e1) * mrSges(6,2) + t282 + t300;
t27 = t479 * t284;
t306 = qJD(1) * t12 - qJD(2) * t27;
t281 = (-t204 * t261 - t205 * t256) * t443 + (t204 * t371 - t205 * t370) * t361;
t285 = t319 * t442 + t443 * t95;
t16 = -t281 + t285 + t336 + t338;
t280 = (-t227 * t256 + t261 * t284) * t443 + (-t227 * t370 - t284 * t371) * t361;
t291 = t128 * t443 + t274 * t361;
t20 = -t280 + t291 + t335 + t337;
t305 = qJD(1) * t16 + qJD(2) * t20;
t125 = m(6) * t284;
t98 = m(6) * t204;
t304 = qJD(1) * t98 - qJD(2) * t125;
t283 = -t268 + ((-qJ(5) - t256) * t276 + t54) * t443;
t19 = -t438 / 0.2e1 + t283;
t288 = qJD(1) * t19 + qJD(3) * t234;
t28 = 0.2e1 * t450 * t443 - t379;
t26 = t280 + t291;
t21 = t281 + t285;
t18 = -t386 + t438 / 0.2e1 + t283;
t11 = mrSges(6,2) * t342 - t381 / 0.2e1 + t282 - t300;
t6 = t380 / 0.2e1 + t383 / 0.2e1 - t382 / 0.2e1 + t384 / 0.2e1 + t279 - t303;
t2 = t278 - t277;
t8 = [qJD(2) * t3 + qJD(3) * t4 + qJD(4) * t9 - qJD(5) * t17, t373 + ((-t227 * t472 + t284 * t473) * t415 + t294 * t416 + t202 * t417 + t102 * t423 + t104 * t424 + t146 * t426 + t147 * t427 + mrSges(3,2) * t271 + t463 * t425 + t464 * t421 - t55 * t376 + t238 * t272 + t200 * t358 + t478 * t242 - t240 * t409 - t48 * t379 + t49 * t377 - t56 * t378 - Ifges(3,6) * t275 + t262 * t113 - pkin(2) * t225 + t79 * t144 + (m(6) * t79 + t112) * t122 + (-m(5) * t55 + m(6) * t49 - t166 + t167) * t155 + (-m(4) * pkin(2) - mrSges(3,1) - t293) * t273 + (t317 + Ifges(3,5)) * t276 + (m(4) * pkin(7) + mrSges(4,3)) * (-t195 * t274 + t411 * t196) + (m(5) * t56 + m(6) * t48 + t160 + t163) * t450 - t295 * t343) * qJD(2) + t2 * qJD(3) + t6 * qJD(4) + t11 * qJD(5), t372 + t2 * qJD(2) + (-t186 * mrSges(4,1) - t185 * mrSges(4,2) - Ifges(4,5) * t367 - Ifges(4,6) * t352 + t204 * t322 + t205 * t323 + t256 * t388 - t261 * t386 + t447 * t58 + t448 * t57 + t457) * qJD(3) + t21 * qJD(4) + t18 * qJD(5), qJD(2) * t6 + qJD(3) * t21 + t389, qJD(2) * t11 + qJD(3) * t18 - t369; -qJD(3) * t1 + qJD(4) * t7 + qJD(5) * t12 - t373, qJD(3) * t5 + qJD(4) * t15 - qJD(5) * t27, (-pkin(7) * t293 - t155 * t447 + t227 * t323 - t256 * t377 - t261 * t379 - t284 * t322 + t448 * t450 + t296 + t459) * qJD(3) + t26 * qJD(4) + t28 * qJD(5) + t309, qJD(3) * t26 + t307, qJD(3) * t28 + t306; qJD(2) * t1 - qJD(4) * t16 + qJD(5) * t19 - t372, -qJD(4) * t20 - t309, t234 * qJD(5), -t305, t288; -qJD(2) * t7 + qJD(3) * t16 + qJD(5) * t98 - t389, qJD(3) * t20 - qJD(5) * t125 - t307, t305, 0, t304; -qJD(2) * t12 - qJD(3) * t19 - qJD(4) * t98 + t369, qJD(4) * t125 - t306, -t288, -t304, 0;];
Cq = t8;
