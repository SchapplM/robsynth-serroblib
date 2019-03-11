% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:24
% EndTime: 2019-03-09 08:21:34
% DurationCPUTime: 4.78s
% Computational Cost: add. (7404->536), mult. (15569->719), div. (0->0), fcn. (14280->6), ass. (0->258)
t447 = Ifges(4,1) + Ifges(6,3);
t446 = Ifges(6,1) + Ifges(5,3);
t445 = Ifges(5,4) + Ifges(6,6);
t444 = Ifges(4,6) + Ifges(6,4);
t443 = pkin(4) + qJ(3);
t288 = cos(pkin(9));
t287 = sin(pkin(9));
t325 = qJ(4) * t287 + pkin(2);
t236 = -pkin(3) * t288 - t325;
t239 = mrSges(5,2) * t288 - mrSges(5,3) * t287;
t442 = m(5) * t236 + t239;
t292 = cos(qJ(2));
t369 = t292 * mrSges(7,2);
t289 = sin(qJ(6));
t291 = cos(qJ(6));
t221 = -t289 * t287 + t288 * t291;
t290 = sin(qJ(2));
t187 = t221 * t290;
t381 = t187 * mrSges(7,3);
t142 = t369 + t381;
t440 = -t381 / 0.2e1 + t142 / 0.2e1;
t439 = -m(6) - m(5);
t438 = -t290 / 0.2e1;
t397 = t290 / 0.2e1;
t437 = -t292 / 0.2e1;
t436 = t292 / 0.2e1;
t319 = t288 * mrSges(6,1) - t287 * mrSges(6,3);
t204 = t319 * t290;
t222 = t287 * t291 + t288 * t289;
t188 = t222 * t290;
t97 = -mrSges(7,1) * t187 + mrSges(7,2) * t188;
t434 = -t204 + t97;
t392 = pkin(3) + qJ(5);
t433 = t288 * t392;
t215 = Ifges(7,4) * t222;
t127 = -t221 * Ifges(7,1) + t215;
t432 = Ifges(7,2) * t221 + t127 + t215;
t358 = t287 * t290;
t225 = -t292 * mrSges(6,1) + mrSges(6,2) * t358;
t234 = mrSges(5,1) * t358 + t292 * mrSges(5,3);
t431 = t234 - t225;
t125 = -mrSges(7,1) * t222 - mrSges(7,2) * t221;
t391 = -pkin(5) - qJ(4);
t165 = t287 * t391 - pkin(2) - t433;
t241 = mrSges(6,1) * t287 + mrSges(6,3) * t288;
t430 = m(7) * t165 + t125 - t241;
t355 = t288 * t290;
t429 = pkin(4) * t355 + t292 * qJ(5);
t346 = t291 * t187;
t351 = t289 * t188;
t428 = -t351 / 0.2e1 - t346 / 0.2e1;
t359 = t188 * t291;
t360 = t187 * t289;
t427 = -t359 / 0.2e1 + t360 / 0.2e1;
t383 = Ifges(5,6) * t288;
t385 = Ifges(6,5) * t288;
t389 = Ifges(4,4) * t288;
t426 = -(Ifges(4,6) / 0.2e1 + Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1) * t290 + (-t287 * Ifges(4,2) + t389) * t437 + Ifges(5,5) * t397 + t444 * t438 + (t446 * t287 - t383 + t385) * t436;
t421 = m(4) / 0.2e1;
t420 = -m(5) / 0.2e1;
t419 = m(5) / 0.2e1;
t418 = -m(6) / 0.2e1;
t417 = m(6) / 0.2e1;
t416 = -m(7) / 0.2e1;
t415 = m(7) / 0.2e1;
t414 = m(4) * pkin(7);
t238 = t443 * t287;
t216 = pkin(8) * t287 + t238;
t240 = t443 * t288;
t217 = pkin(8) * t288 + t240;
t113 = t216 * t291 + t217 * t289;
t413 = -t113 / 0.2e1;
t411 = t187 / 0.2e1;
t409 = t188 / 0.2e1;
t189 = t221 * t292;
t408 = t189 / 0.2e1;
t190 = t222 * t292;
t407 = t190 / 0.2e1;
t406 = -t221 / 0.2e1;
t405 = t221 / 0.2e1;
t403 = t222 / 0.2e1;
t401 = -t287 / 0.2e1;
t399 = -t288 / 0.2e1;
t398 = t288 / 0.2e1;
t396 = -t291 / 0.2e1;
t394 = m(6) * t288;
t284 = t290 * pkin(7);
t390 = Ifges(4,4) * t287;
t388 = Ifges(7,4) * t188;
t387 = Ifges(7,4) * t221;
t386 = Ifges(6,5) * t287;
t384 = Ifges(5,6) * t287;
t237 = -pkin(2) * t292 - t290 * qJ(3) - pkin(1);
t357 = t287 * t292;
t259 = pkin(7) * t357;
t152 = t237 * t288 - t259;
t285 = t292 * pkin(3);
t140 = -t152 + t285;
t109 = t140 + t429;
t336 = -pkin(7) * t287 - pkin(3);
t354 = t288 * t292;
t306 = pkin(4) * t354 + (-qJ(5) + t336) * t290;
t367 = qJ(3) * t292;
t242 = t290 * pkin(2) - t367;
t356 = t288 * t242;
t111 = t306 - t356;
t341 = pkin(3) * t358 + t284;
t366 = qJ(5) * t287;
t117 = (t288 * t391 + t366) * t290 + t341;
t255 = qJ(4) * t354;
t322 = t287 * t392 + pkin(7);
t118 = -t255 + (-pkin(5) * t288 + t322) * t292;
t153 = pkin(7) * t354 + t287 * t237;
t139 = qJ(4) * t292 - t153;
t120 = -pkin(4) * t358 - t139;
t223 = t287 * t242;
t164 = -pkin(7) * t355 + t223;
t280 = t290 * qJ(4);
t141 = -t164 - t280;
t122 = -pkin(4) * t357 - t141;
t137 = (qJ(4) * t288 - t366) * t290 - t341;
t138 = t292 * t322 - t255;
t370 = t292 * mrSges(7,1);
t380 = t188 * mrSges(7,3);
t143 = -t370 - t380;
t144 = -mrSges(7,2) * t290 + mrSges(7,3) * t189;
t145 = mrSges(7,1) * t290 - mrSges(7,3) * t190;
t146 = t290 * t336 - t356;
t163 = pkin(7) * t358 + t356;
t181 = -t280 * t288 + t341;
t182 = -t255 + (pkin(3) * t287 + pkin(7)) * t292;
t205 = t319 * t292;
t318 = -t287 * mrSges(5,2) - t288 * mrSges(5,3);
t206 = t318 * t290;
t207 = t318 * t292;
t208 = (t287 * mrSges(4,1) + t288 * mrSges(4,2)) * t292;
t224 = t292 * mrSges(4,2) - mrSges(4,3) * t358;
t226 = -t290 * mrSges(4,2) - mrSges(4,3) * t357;
t256 = mrSges(6,2) * t357;
t375 = t290 * mrSges(6,1);
t227 = t256 + t375;
t228 = -t292 * mrSges(4,1) - mrSges(4,3) * t355;
t229 = t290 * mrSges(4,1) - mrSges(4,3) * t354;
t373 = t290 * mrSges(6,3);
t230 = -mrSges(6,2) * t354 - t373;
t231 = mrSges(5,1) * t357 - t290 * mrSges(5,3);
t257 = mrSges(5,1) * t354;
t374 = t290 * mrSges(5,2);
t232 = t257 + t374;
t233 = -mrSges(6,2) * t355 + t292 * mrSges(6,3);
t235 = mrSges(5,1) * t355 - t292 * mrSges(5,2);
t301 = Ifges(4,5) * t397 + (-t288 * Ifges(5,2) + t384) * t437 + (Ifges(4,5) / 0.2e1 - Ifges(6,6) / 0.2e1 - Ifges(5,4) / 0.2e1) * t290 + t445 * t438 + (t447 * t288 + t386 - t390) * t436;
t309 = -Ifges(7,5) * t190 / 0.2e1 - Ifges(7,6) * t189 / 0.2e1;
t89 = t259 + t285 + (pkin(8) * t290 - t237) * t288 + t429;
t338 = t287 * (-pkin(4) - pkin(8));
t90 = t290 * t338 + t292 * t391 + t153;
t41 = -t289 * t89 + t291 * t90;
t42 = t289 * t90 + t291 * t89;
t91 = (pkin(8) * t292 - t242) * t288 + t306;
t92 = t223 + t280 + (-pkin(7) * t288 + pkin(5)) * t290 + t292 * t338;
t43 = -t289 * t91 + t291 * t92;
t44 = t289 * t92 + t291 * t91;
t85 = Ifges(7,2) * t187 - t292 * Ifges(7,6) + t388;
t86 = Ifges(7,4) * t190 + Ifges(7,2) * t189 + t290 * Ifges(7,6);
t174 = Ifges(7,4) * t187;
t87 = Ifges(7,1) * t188 - t292 * Ifges(7,5) + t174;
t88 = Ifges(7,1) * t190 + Ifges(7,4) * t189 + t290 * Ifges(7,5);
t378 = t190 * mrSges(7,2);
t379 = t189 * mrSges(7,1);
t98 = t378 - t379;
t1 = (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t290 + Ifges(7,5) * t409 + Ifges(7,6) * t411 + pkin(7) * t208 + t287 * t426 + t301 * t288) * t290 + (Ifges(3,4) * t292 - pkin(1) * mrSges(3,2) + (-Ifges(4,5) + t445) * t354 + (-Ifges(5,5) + t444) * t357 + (-Ifges(3,2) + Ifges(3,1) + m(4) * pkin(7) ^ 2 - Ifges(7,3) - Ifges(6,2) - Ifges(5,1) - Ifges(4,3) + (pkin(7) * mrSges(4,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t288) * t288 + (pkin(7) * mrSges(4,1) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t287 + (-Ifges(4,4) + Ifges(6,5) - Ifges(5,6)) * t288) * t287) * t290 + t309) * t292 + m(4) * (t152 * t163 + t153 * t164) + m(5) * (t139 * t141 + t140 * t146 + t181 * t182) + m(6) * (t109 * t111 + t120 * t122 - t137 * t138) + m(7) * (t117 * t118 + t41 * t43 + t42 * t44) + t87 * t407 + t85 * t408 + t88 * t409 + t86 * t411 + t117 * t98 + t118 * t97 + t44 * t142 + t43 * t143 + t42 * t144 + t41 * t145 - t138 * t204 + t137 * t205 + t182 * t206 + t181 * t207 + t164 * t224 + t122 * t225 + t153 * t226 + t120 * t227 + t163 * t228 + t152 * t229 + t109 * t230 + t139 * t231 + t140 * t232 + t111 * t233 + t141 * t234 + t146 * t235;
t382 = t1 * qJD(1);
t377 = t289 * mrSges(7,1);
t376 = t289 * mrSges(7,2);
t372 = t291 * mrSges(7,1);
t371 = t291 * mrSges(7,2);
t100 = Ifges(7,1) * t187 - t388;
t343 = Ifges(7,5) * t187 - Ifges(7,6) * t188;
t96 = t188 * mrSges(7,1) + mrSges(7,2) * t187;
t99 = -Ifges(7,2) * t188 + t174;
t4 = t117 * t96 + t41 * t142 - t42 * t143 + t343 * t437 + (-t42 * mrSges(7,3) + t100 / 0.2e1 - t85 / 0.2e1) * t188 + (-t41 * mrSges(7,3) + t87 / 0.2e1 + t99 / 0.2e1) * t187;
t368 = t4 * qJD(1);
t11 = t187 * t142 - t188 * t143 + m(7) * (t187 * t42 - t188 * t41) + ((-t228 + t233 + t235) * t288 + (-t224 + t431) * t287 + m(6) * (t109 * t288 - t120 * t287) + m(5) * (t139 * t287 + t140 * t288) + m(4) * (-t152 * t288 - t153 * t287)) * t290;
t365 = qJD(1) * t11;
t347 = t291 * t143;
t353 = t289 * t142;
t12 = (m(5) * t181 - m(6) * t137 + m(7) * t117 + t206 + t434) * t355 + (t347 + t353 - m(7) * (-t289 * t42 - t291 * t41) + m(6) * t120 - m(5) * t139 - t431) * t292;
t364 = qJD(1) * t12;
t348 = t291 * t142;
t352 = t289 * t143;
t15 = t434 * t358 + (t233 + t348 - t352) * t292 + m(7) * (t117 * t358 + (-t289 * t41 + t291 * t42) * t292) + m(6) * (t109 * t292 - t137 * t358);
t363 = qJD(1) * t15;
t303 = t369 / 0.2e1 + t440;
t304 = -t370 / 0.2e1 + t380 / 0.2e1 + t143 / 0.2e1;
t17 = t289 * t303 + t291 * t304;
t362 = t17 * qJD(1);
t18 = t289 * t304 - t291 * t303;
t361 = t18 * qJD(1);
t350 = t289 * t221;
t349 = t289 * t222;
t345 = t291 * t221;
t344 = t291 * t222;
t342 = Ifges(7,5) * t222 + Ifges(7,6) * t221;
t106 = (-m(6) + m(7) * (-t289 ^ 2 - t291 ^ 2)) * t292;
t340 = t106 * qJD(1);
t339 = mrSges(5,2) / 0.2e1 - mrSges(6,3) / 0.2e1;
t337 = -t97 / 0.2e1 + t204 / 0.2e1;
t335 = t358 / 0.2e1;
t333 = -t355 / 0.2e1;
t329 = t125 / 0.2e1 - t241 / 0.2e1;
t321 = t329 * t287;
t317 = Ifges(7,1) * t222 + t387;
t124 = mrSges(7,1) * t221 - t222 * mrSges(7,2);
t126 = Ifges(7,2) * t222 - t387;
t14 = t165 * t124 + t126 * t406 + t317 * t405 - t432 * t222 / 0.2e1;
t112 = -t216 * t289 + t217 * t291;
t294 = t143 * t413 - t117 * t124 / 0.2e1 + t165 * t96 / 0.2e1 - t292 * t342 / 0.4e1 + t432 * t187 / 0.4e1 + (t87 + t99) * t222 / 0.4e1 + (-t100 / 0.4e1 + t85 / 0.4e1) * t221 + t440 * t112 + (t413 * mrSges(7,3) - t126 / 0.4e1 + t317 / 0.4e1) * t188;
t299 = Ifges(7,3) * t397 + t43 * mrSges(7,1) / 0.2e1 - t44 * mrSges(7,2) / 0.2e1 - t309;
t2 = t294 - t299;
t316 = t2 * qJD(1) - t14 * qJD(2);
t22 = -m(7) * (t112 * t221 + t113 * t222) - m(6) * (t238 * t287 + t240 * t288) + (-t221 ^ 2 - t222 ^ 2) * mrSges(7,3) + (-mrSges(6,2) - 0.4e1 * (-m(5) / 0.4e1 - m(4) / 0.4e1) * qJ(3) + mrSges(4,3) + mrSges(5,1)) * (-t287 ^ 2 - t288 ^ 2);
t293 = (t187 * t403 + t188 * t406) * mrSges(7,3) + (-t234 / 0.2e1 + t225 / 0.2e1 + t224 / 0.2e1) * t288 + (t235 / 0.2e1 + t233 / 0.2e1 - t228 / 0.2e1) * t287 + (-t152 * t287 + t153 * t288) * t421 + (-t139 * t288 + t140 * t287) * t419 + ((t238 * t290 + t120) * t288 + (-t240 * t290 + t109) * t287) * t417 + (-t112 * t188 + t113 * t187 + t221 * t41 + t222 * t42) * t415 + t143 * t405 + t142 * t403;
t300 = t182 * t420 + t138 * t418 + t118 * t416 + t379 / 0.2e1 - t378 / 0.2e1;
t6 = t293 + (-t414 / 0.2e1 + (mrSges(6,1) / 0.2e1 - mrSges(4,2) / 0.2e1 + mrSges(5,3) / 0.2e1) * t288 + (-mrSges(4,1) / 0.2e1 + t339) * t287) * t292 + t300;
t315 = qJD(1) * t6 - qJD(2) * t22;
t202 = t325 + t433;
t32 = (-m(6) * t202 + t430 + t442) * t287;
t295 = (-t206 / 0.2e1 + t337) * t287 + (-t287 * t181 + (-t236 * t290 - t367) * t288) * t419 + (t287 * t137 + t202 * t355 - t240 * t292) * t417 + (-t165 * t355 - t287 * t117 + (-t112 * t291 - t113 * t289) * t292) * t415;
t297 = t146 * t420 + t111 * t418 + (-t289 * t43 + t291 * t44) * t416 + t289 * t145 / 0.2e1 + t144 * t396;
t307 = -t349 / 0.2e1 - t345 / 0.2e1;
t302 = t307 * mrSges(7,3);
t308 = (-t239 / 0.2e1 - t329) * t288;
t7 = -t257 + (t288 * mrSges(6,2) + t302) * t292 + (t308 - t339) * t290 + t295 + t297;
t314 = -qJD(1) * t7 + qJD(2) * t32;
t45 = -t202 * t394 + t288 * t430;
t296 = t337 * t288 + (t344 / 0.2e1 - t350 / 0.2e1) * t292 * mrSges(7,3) + (t288 * t137 - t202 * t358 + t238 * t292) * t417 + (t165 * t358 - t288 * t117 + (-t112 * t289 + t113 * t291) * t292) * t415;
t298 = t122 * t418 + (t289 * t44 + t291 * t43) * t416 - t289 * t144 / 0.2e1 + t145 * t396;
t9 = -t256 + (-mrSges(6,1) / 0.2e1 + t321) * t290 + t296 + t298;
t313 = -qJD(1) * t9 + qJD(2) * t45;
t51 = (t333 + t428) * m(7) + t439 * t355;
t305 = (t344 - t350) * t415;
t69 = t305 + (t415 - t439) * t287;
t312 = qJD(1) * t51 - qJD(2) * t69;
t54 = m(6) * t358 + (t335 - t427) * m(7);
t78 = t394 + (t398 - t307) * m(7);
t311 = qJD(1) * t54 - qJD(2) * t78;
t310 = qJD(1) * t96 - qJD(2) * t124;
t94 = m(7) * t399 + (t345 + t349) * t415;
t93 = m(7) * t401 + t305;
t64 = m(7) * t333 + (t346 + t351) * t415;
t63 = m(7) * t335 + (-t359 + t360) * t415;
t20 = t348 / 0.2e1 - t352 / 0.2e1 + (-t371 / 0.2e1 - t377 / 0.2e1) * t292 + t428 * mrSges(7,3);
t19 = -t353 / 0.2e1 - t347 / 0.2e1 + (t376 / 0.2e1 - t372 / 0.2e1) * t292 + t427 * mrSges(7,3);
t10 = t375 / 0.2e1 + t290 * t321 + t296 - t298;
t8 = -t373 / 0.2e1 + t374 / 0.2e1 + t292 * t302 + t290 * t308 + t295 - t297;
t5 = t293 + t414 * t436 + mrSges(4,2) * t354 / 0.2e1 - mrSges(5,2) * t357 / 0.2e1 - t300 + (mrSges(4,1) + mrSges(6,3)) * t357 / 0.2e1 - (mrSges(6,1) + mrSges(5,3)) * t354 / 0.2e1;
t3 = t294 + t299;
t13 = [qJD(2) * t1 + qJD(3) * t11 - qJD(4) * t12 + qJD(5) * t15 + qJD(6) * t4, t5 * qJD(3) + t8 * qJD(4) + t10 * qJD(5) + t3 * qJD(6) + t382 + (t442 * t182 + ((Ifges(4,2) * t288 + t390) * t401 + (-Ifges(5,2) * t287 - t383) * t399 + Ifges(3,5) + (-m(4) * pkin(2) - mrSges(4,1) * t288 + mrSges(4,2) * t287 - mrSges(3,1)) * pkin(7) + (-t446 * t288 - t384 + t386) * t287 / 0.2e1 + (t287 * t447 - t385 + t389) * t398) * t292 + (-t163 * mrSges(4,3) + t146 * mrSges(5,1) - t111 * mrSges(6,2) + (-t229 + t232) * qJ(3) + t301) * t287 + (t164 * mrSges(4,3) - t141 * mrSges(5,1) - t122 * mrSges(6,2) + (t226 - t231) * qJ(3) - t426) * t288 + t86 * t403 + t88 * t406 + t127 * t407 + t126 * t408 + 0.2e1 * (t112 * t43 + t113 * t44 + t118 * t165) * t415 + 0.2e1 * (t111 * t238 + t122 * t240 - t138 * t202) * t417 + mrSges(3,2) * t284 + (-Ifges(7,5) * t221 + Ifges(7,6) * t222) * t397 + 0.2e1 * ((-t163 * t287 + t164 * t288) * t421 + (-t141 * t288 + t146 * t287) * t419) * qJ(3) + (t221 * t43 + t222 * t44) * mrSges(7,3) + t118 * t125 + t113 * t144 + t112 * t145 + t165 * t98 + t202 * t205 - pkin(2) * t208 + t236 * t207 + t238 * t230 + t240 * t227 - t138 * t241 - Ifges(3,6) * t290) * qJD(2), qJD(2) * t5 + qJD(4) * t64 + qJD(5) * t63 + t365, qJD(2) * t8 + qJD(3) * t64 + qJD(6) * t19 - t364, qJD(2) * t10 + qJD(3) * t63 + qJD(6) * t20 + t363, t368 + t3 * qJD(2) + t19 * qJD(4) + t20 * qJD(5) + (-mrSges(7,1) * t42 - mrSges(7,2) * t41 + t343) * qJD(6); qJD(3) * t6 + qJD(4) * t7 + qJD(5) * t9 + qJD(6) * t2 - t382, -qJD(3) * t22 - qJD(4) * t32 - qJD(5) * t45 - qJD(6) * t14, qJD(4) * t93 + qJD(5) * t94 + t315, qJD(3) * t93 - t314, qJD(3) * t94 - t313 (-mrSges(7,1) * t113 - mrSges(7,2) * t112 + t342) * qJD(6) + t316; -qJD(2) * t6 + qJD(4) * t51 + qJD(5) * t54 + qJD(6) * t96 - t365, -qJD(4) * t69 - qJD(5) * t78 - qJD(6) * t124 - t315, 0, t312, t311, t310; -qJD(2) * t7 - qJD(3) * t51 - qJD(5) * t106 - qJD(6) * t17 + t364, qJD(3) * t69 + t314, -t312, 0, -t340, -t362 + (-t372 + t376) * qJD(6); -qJD(2) * t9 - qJD(3) * t54 + qJD(4) * t106 - qJD(6) * t18 - t363, qJD(3) * t78 + t313, -t311, t340, 0, -t361 + (-t371 - t377) * qJD(6); -qJD(2) * t2 - qJD(3) * t96 + qJD(4) * t17 + qJD(5) * t18 - t368, qJD(3) * t124 - t316, -t310, t362, t361, 0;];
Cq  = t13;
