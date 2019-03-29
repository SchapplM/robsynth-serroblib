% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:25:53
% EndTime: 2019-03-29 15:26:02
% DurationCPUTime: 4.75s
% Computational Cost: add. (6658->358), mult. (18360->505), div. (0->0), fcn. (17975->8), ass. (0->246)
t227 = cos(qJ(3));
t466 = t227 ^ 2;
t224 = sin(qJ(3));
t412 = sin(qJ(4));
t413 = cos(qJ(4));
t255 = t224 * t413 + t227 * t412;
t398 = Ifges(6,3) * t255;
t464 = t398 / 0.2e1;
t189 = t412 * t224 - t413 * t227;
t141 = mrSges(5,1) * t189 + mrSges(5,2) * t255;
t409 = pkin(2) * t224;
t431 = t141 * t409 + Ifges(4,4) * t466 + (-Ifges(4,4) * t224 + (-Ifges(4,2) + Ifges(4,1)) * t227) * t224;
t225 = sin(qJ(2));
t411 = pkin(1) * t225;
t165 = t189 * t411;
t462 = t165 / 0.2e1;
t460 = Ifges(5,1) - Ifges(5,2);
t223 = sin(qJ(5));
t219 = t223 ^ 2;
t226 = cos(qJ(5));
t221 = t226 ^ 2;
t443 = t219 + t221;
t438 = mrSges(6,3) * t443;
t391 = t226 * mrSges(6,1);
t394 = t223 * mrSges(6,2);
t191 = -t391 + t394;
t459 = mrSges(5,1) - t191;
t215 = Ifges(6,5) * t226;
t399 = Ifges(6,6) * t223;
t274 = t399 - t215;
t458 = t189 * t274;
t408 = t227 * pkin(2);
t228 = cos(qJ(2));
t410 = pkin(1) * t228;
t210 = -t408 - t410;
t112 = t165 * t223 + t210 * t226;
t113 = -t165 * t226 + t210 * t223;
t248 = pkin(1) * t255;
t166 = t225 * t248;
t367 = t166 * t223;
t135 = t226 * t409 + t367;
t366 = t166 * t226;
t137 = t223 * t409 - t366;
t426 = m(6) * pkin(2);
t344 = t426 / 0.2e1;
t229 = pkin(2) ^ 2;
t352 = t224 * t229;
t359 = t255 * t223;
t130 = -mrSges(6,2) * t189 - mrSges(6,3) * t359;
t374 = t137 * t130;
t358 = t255 * t226;
t132 = mrSges(6,1) * t189 - mrSges(6,3) * t358;
t375 = t135 * t132;
t376 = t132 * t226;
t441 = t130 * t223 + t376;
t445 = m(5) * t210;
t256 = (t224 * mrSges(4,1) + t227 * mrSges(4,2)) * t410;
t446 = -t256 / 0.2e1;
t457 = t375 / 0.2e1 + t374 / 0.2e1 + ((t112 * t224 - t135 * t227) * t226 + (t113 * t224 - t137 * t227) * t223) * t344 + t446 - m(5) * t227 * t352 / 0.2e1 + (t441 + t445) * t409 / 0.2e1 + t431;
t454 = -t413 / 0.2e1;
t401 = Ifges(6,4) * t223;
t195 = Ifges(6,2) * t226 + t401;
t121 = t195 * t255;
t200 = Ifges(6,1) * t226 - t401;
t76 = Ifges(6,5) * t189 + t200 * t255;
t453 = t76 - t121;
t217 = Ifges(6,4) * t226;
t199 = Ifges(6,1) * t223 + t217;
t275 = Ifges(6,2) * t223 - t217;
t452 = t199 - t275;
t361 = t189 * t223;
t129 = -mrSges(6,2) * t255 + mrSges(6,3) * t361;
t321 = t226 * t412;
t107 = pkin(2) * t129 * t321;
t325 = t255 * t412;
t285 = mrSges(5,3) * t325;
t173 = pkin(2) * t285;
t120 = t255 * (Ifges(6,5) * t223 + Ifges(6,6) * t226);
t360 = t189 * t226;
t301 = -t360 / 0.2e1;
t302 = t361 / 0.2e1;
t414 = t226 / 0.2e1;
t415 = t223 / 0.2e1;
t73 = Ifges(6,6) * t255 + t189 * t275;
t75 = Ifges(6,5) * t255 - t189 * t200;
t257 = t195 * t302 + t199 * t301 + t73 * t414 + t75 * t415 + t120 / 0.2e1 - Ifges(5,6) * t255 - Ifges(5,5) * t189;
t131 = mrSges(6,1) * t255 + mrSges(6,3) * t360;
t315 = t412 * t131;
t390 = t226 * mrSges(6,2);
t395 = t223 * mrSges(6,1);
t193 = t390 + t395;
t118 = t193 * t189;
t320 = t413 * t118;
t451 = Ifges(4,5) * t227 - Ifges(4,6) * t224 + t107 - t173 + t257 + (mrSges(5,3) * t413 * t189 - t223 * t315 + t320) * pkin(2);
t450 = -mrSges(4,1) * t227 + mrSges(4,2) * t224;
t343 = t223 * t408;
t370 = t166 * t118;
t449 = -t129 * t343 / 0.2e1 - t370 / 0.2e1;
t122 = t199 * t255;
t389 = t226 * t76;
t74 = Ifges(6,6) * t189 - t255 * t275;
t393 = t223 * t74;
t448 = -t226 * t121 / 0.4e1 - t223 * t122 / 0.4e1 - t458 / 0.4e1 + t389 / 0.4e1 - t393 / 0.4e1 - t452 * t359 / 0.4e1;
t119 = t193 * t255;
t140 = mrSges(5,1) * t255 - mrSges(5,2) * t189;
t357 = t210 * t140;
t383 = t113 * t129;
t385 = t112 * t131;
t447 = -t357 / 0.2e1 + t119 * t462 - t383 / 0.2e1 - t385 / 0.2e1;
t444 = t459 * t165;
t442 = t389 / 0.2e1 - t393 / 0.2e1;
t377 = t131 * t226;
t440 = t140 / 0.2e1 + t377 / 0.2e1;
t439 = -mrSges(5,2) * t413 - t412 * t459;
t350 = t226 * t130;
t291 = t350 / 0.2e1;
t342 = t412 / 0.2e1;
t436 = t413 * t291 + t285 / 0.2e1 + t119 * t342 + t320 / 0.2e1;
t117 = t191 * t255;
t340 = pkin(2) * t412;
t283 = -t340 / 0.2e1;
t270 = t223 * t283;
t341 = pkin(2) * t413;
t284 = -t341 / 0.2e1;
t427 = pkin(2) / 0.2e1;
t434 = -t325 * t427 * t438 - t117 * t284 + t130 * t270 + t283 * t376;
t168 = t189 * t410;
t136 = t168 * t223 + t226 * t411;
t138 = -t168 * t226 + t223 * t411;
t167 = t228 * t248;
t433 = t167 * t119 + t138 * t130 + t136 * t132 + (t167 * t255 + t168 * t189) * mrSges(5,3);
t396 = t255 * mrSges(5,3);
t123 = t165 * t396;
t115 = t123 / 0.2e1;
t432 = t115 + t447;
t416 = -t223 / 0.2e1;
t238 = (-mrSges(5,1) / 0.2e1 + t191 / 0.2e1) * t167 + (t136 * t416 + t138 * t414) * mrSges(6,3) + t168 * mrSges(5,2) / 0.2e1;
t180 = Ifges(5,4) * t189;
t429 = (t464 + t458 / 0.2e1 + t180) * t189 + ((Ifges(5,2) / 0.2e1 - t460 / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(6,3) / 0.2e1) * t189 + (-Ifges(5,4) - t274 / 0.2e1) * t255) * t255;
t296 = t358 / 0.2e1;
t298 = -t359 / 0.2e1;
t249 = t75 * t296 + t73 * t298 + t76 * t301 + t74 * t302 + t429;
t335 = -t408 / 0.2e1;
t287 = t226 * t335;
t428 = -t123 / 0.2e1 + t115 + t140 * t335 + t249 + t131 * t287 - t447 + t449;
t424 = -mrSges(6,1) / 0.2e1;
t423 = mrSges(6,1) / 0.2e1;
t422 = -mrSges(6,2) / 0.2e1;
t421 = mrSges(6,2) / 0.2e1;
t419 = -t166 / 0.2e1;
t417 = t193 / 0.2e1;
t405 = mrSges(6,3) * t255;
t397 = t166 * mrSges(5,2);
t235 = t249 - t123 + (-t119 + t396) * t165 + t383 + t385 + t357;
t99 = t166 * t165;
t3 = t235 - t256 + t375 + t374 - t370 + m(6) * (t112 * t135 + t113 * t137 - t99) + t409 * t445 + t431;
t387 = t3 * qJD(1);
t381 = t113 * t226;
t384 = t112 * t223;
t266 = -t381 + t384;
t355 = t223 * t132;
t6 = t235 + (-t118 - t350 + t355) * t166 + m(6) * (t166 * t266 - t99);
t386 = t6 * qJD(1);
t382 = t113 * t132;
t380 = t129 * t223;
t294 = -t358 / 0.2e1;
t363 = t189 * t120;
t277 = -t122 * t296 + t74 * t294 - t363 / 0.2e1 + t453 * t298;
t371 = t166 * t117;
t13 = t112 * t130 + t266 * t405 + t277 - t371 - t382;
t379 = t13 * qJD(1);
t345 = t224 ^ 2 + t466;
t276 = mrSges(4,3) * t345 - mrSges(3,2);
t334 = -mrSges(3,1) + t141 + t450;
t368 = t166 * t167;
t16 = t334 * t411 + m(6) * (t112 * t136 + t113 * t138 + t368) + m(5) * (t165 * t168 + t210 * t411 + t368) + (t276 + m(4) * (-0.1e1 + t345) * t411) * t410 + t433;
t373 = t16 * qJD(1);
t339 = t135 * t223 * mrSges(6,3);
t338 = t137 * t226 * mrSges(6,3);
t331 = mrSges(5,3) * t462;
t326 = Ifges(6,5) * t301 + Ifges(6,6) * t302 + t464;
t324 = t223 * t413;
t323 = t223 * t412;
t322 = t226 * t413;
t319 = t413 * t165;
t318 = t413 * t167;
t317 = t413 * t193;
t313 = t412 * t219;
t312 = t412 * t221;
t304 = t166 * t417;
t303 = t367 / 0.2e1;
t300 = t359 / 0.2e1;
t295 = t358 / 0.4e1;
t293 = -t358 / 0.4e1;
t290 = t107 / 0.2e1 - t173 / 0.2e1;
t56 = t195 * t416 + t200 * t415 + t414 * t452;
t279 = t132 * t454;
t232 = m(5) * (-t168 * t412 - t318) * t427 + (-t136 * t323 + t138 * t321 - t318) * t344 + t446 + t238;
t1 = t442 * t189 - t255 * t331 + t75 * t294 + t73 * t300 + t440 * t408 + t232 - t429 + t432 - t449 - t457;
t260 = t140 + t377 + t380;
t8 = t249 + t441 * t409 + ((-m(6) * t443 - m(5)) * t352 - t260 * pkin(2)) * t227 + t431;
t273 = -t1 * qJD(1) + t8 * qJD(2);
t263 = t215 / 0.2e1 - t399 / 0.2e1;
t236 = t75 * t414 + t73 * t416 - (Ifges(5,4) - t263) * t255 + (Ifges(6,3) - t460) * t189;
t254 = t263 * t189;
t237 = (-t180 + t254 + t442) * t189;
t14 = -t236 * t255 + t260 * t408 + t237;
t4 = (t291 - t355 / 0.2e1 + t118 / 0.2e1) * t166 + (t380 / 0.2e1 + t440) * t408 + t237 - (t331 + t236) * t255 + t238 + t432;
t269 = -t4 * qJD(1) - t14 * qJD(2);
t106 = t132 * t343;
t24 = t350 * t408 - t106 + t74 * t296 - t122 * t294 + t363 / 0.2e1 + t453 * t300;
t234 = (t287 + t112 / 0.2e1) * t130 - (t381 / 0.2e1 - t384 / 0.2e1) * t405 + t106 / 0.2e1 - t382 / 0.2e1 - t371 / 0.2e1 + t277;
t261 = t136 * t423 + t138 * t422;
t9 = t234 - t261;
t268 = t9 * qJD(1) - t24 * qJD(2);
t262 = t135 * t424 + t137 * t421;
t247 = t195 * t293 + t200 * t295 + t448;
t23 = -t398 / 0.2e1 + t254 + t247;
t12 = t23 + t304 + t262 + t434;
t17 = t200 * t293 + t195 * t295 + (t117 * t454 + (t132 * t342 + t224 * t423) * t226 + (t130 * t342 + t224 * t422) * t223 - (-t312 / 0.2e1 - t313 / 0.2e1) * t405) * pkin(2) + t326 - t448;
t243 = (-t199 / 0.2e1 + t275 / 0.2e1) * t226 + (-t200 / 0.2e1 + t195 / 0.2e1) * t223;
t38 = pkin(2) * t317 + t243;
t253 = t12 * qJD(1) - t17 * qJD(2) - t38 * qJD(3);
t231 = (-t112 * t324 + t113 * t322 + t319 + (-t313 - t312 + t412) * t166) * t344 + t131 * t270 + t290 + t419 * t438 + (t223 * t279 + t436) * pkin(2);
t250 = t339 / 0.2e1 - t338 / 0.2e1;
t15 = t231 + t250;
t25 = ((t279 - t315 / 0.2e1) * t223 + t436) * pkin(2) + t290;
t33 = t341 * t438 + m(6) * (-0.1e1 + t443) * t229 * t412 * t413 + t439 * pkin(2);
t252 = t15 * qJD(1) + t25 * qJD(2) + t33 * qJD(3);
t246 = t257 + t444;
t245 = t322 * t422 + t324 * t424;
t20 = (t417 - t390 / 0.2e1 - t395 / 0.2e1) * t166 + t23;
t31 = (t317 / 0.2e1 + t245) * pkin(2) + t243;
t244 = -qJD(1) * t20 - qJD(2) * t23 + qJD(3) * t31 - qJD(4) * t56;
t22 = t247 + t326;
t239 = t22 + t434;
t32 = pkin(2) * t245 + t193 * t284 + t56;
t21 = t25 + t257;
t19 = mrSges(6,1) * t303 + t366 * t421 + t22 + t304;
t18 = t239 + (-t394 / 0.2e1 + t391 / 0.2e1) * t409;
t11 = t239 + t304 - t262;
t10 = t234 + t261;
t7 = t231 + t246 - t250 + t397;
t5 = t132 * t303 + t350 * t419 + t238 + t428;
t2 = t232 + t428 + t457;
t26 = [qJD(2) * t16 + qJD(3) * t3 + qJD(4) * t6 + qJD(5) * t13, t2 * qJD(3) + t5 * qJD(4) + t10 * qJD(5) + t373 + (m(6) * (-t136 * t226 - t138 * t223) * t408 + t433 + (t276 * t228 + (-m(5) * t408 + t334) * t225) * pkin(1)) * qJD(2), t2 * qJD(2) + t7 * qJD(4) + t11 * qJD(5) + t387 + ((-t135 * t323 + t137 * t321 + t319) * t426 + t338 - t339 + t397 + m(5) * (pkin(2) * t319 - t166 * t340) + t444 + t450 * t411 + t451) * qJD(3), t5 * qJD(2) + t7 * qJD(3) + t19 * qJD(5) + t386 + (t246 + (mrSges(5,2) - t438) * t166) * qJD(4), t379 + t10 * qJD(2) + t11 * qJD(3) + t19 * qJD(4) + (-mrSges(6,1) * t113 - mrSges(6,2) * t112 - t120) * qJD(5); -qJD(3) * t1 - qJD(4) * t4 + qJD(5) * t9 - t373, qJD(3) * t8 - qJD(4) * t14 - qJD(5) * t24, qJD(3) * t451 + t21 * qJD(4) + t18 * qJD(5) + t273, t21 * qJD(3) + qJD(4) * t257 + t22 * qJD(5) + t269, t18 * qJD(3) + t22 * qJD(4) + (t193 * t408 - t120) * qJD(5) + t268; qJD(2) * t1 + qJD(4) * t15 + qJD(5) * t12 - t387, qJD(4) * t25 - qJD(5) * t17 - t273, qJD(4) * t33 - qJD(5) * t38, t32 * qJD(5) + (t413 * t438 + t439) * qJD(4) * pkin(2) + t252, t32 * qJD(4) + ((-mrSges(6,1) * t321 + mrSges(6,2) * t323) * pkin(2) - t274) * qJD(5) + t253; qJD(2) * t4 - qJD(3) * t15 + qJD(5) * t20 - t386, -qJD(3) * t25 + qJD(5) * t23 - t269, -qJD(5) * t31 - t252, t56 * qJD(5), -qJD(5) * t274 - t244; -qJD(2) * t9 - qJD(3) * t12 - qJD(4) * t20 - t379, qJD(3) * t17 - qJD(4) * t23 - t268, qJD(4) * t31 - t253, t244, 0;];
Cq  = t26;
