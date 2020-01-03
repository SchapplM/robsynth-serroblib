% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:47
% EndTime: 2019-12-31 19:09:58
% DurationCPUTime: 5.04s
% Computational Cost: add. (14862->410), mult. (29995->571), div. (0->0), fcn. (33428->8), ass. (0->235)
t262 = sin(pkin(9));
t263 = cos(pkin(9));
t266 = sin(qJ(3));
t397 = cos(qJ(3));
t239 = t397 * t262 + t266 * t263;
t405 = t239 / 0.2e1;
t372 = pkin(4) * qJD(4);
t237 = t262 * t266 - t397 * t263;
t264 = sin(qJ(5));
t265 = sin(qJ(4));
t267 = cos(qJ(5));
t268 = cos(qJ(4));
t291 = t264 * t268 + t267 * t265;
t168 = t291 * t239;
t369 = t168 * mrSges(6,3);
t126 = -mrSges(6,2) * t237 - t369;
t417 = t126 / 0.2e1;
t403 = t291 / 0.2e1;
t305 = -pkin(2) * t263 - pkin(1);
t183 = pkin(3) * t237 - pkin(7) * t239 + t305;
t384 = pkin(6) + qJ(2);
t245 = t384 * t263;
t301 = t384 * t262;
t201 = t397 * t245 - t266 * t301;
t103 = t265 * t183 + t201 * t268;
t338 = t239 * t265;
t83 = -pkin(8) * t338 + t103;
t360 = t264 * t83;
t102 = t268 * t183 - t265 * t201;
t337 = t239 * t268;
t82 = -pkin(8) * t337 + t102;
t69 = t237 * pkin(4) + t82;
t40 = t267 * t69 - t360;
t49 = t267 * t82 - t360;
t446 = t40 - t49;
t418 = -pkin(8) - pkin(7);
t250 = t418 * t265;
t251 = t418 * t268;
t209 = t250 * t264 - t251 * t267;
t298 = t267 * t250 + t251 * t264;
t230 = Ifges(6,6) * t291;
t290 = t264 * t265 - t267 * t268;
t231 = Ifges(6,5) * t290;
t320 = -t231 - t230;
t28 = -t209 * mrSges(6,1) - t298 * mrSges(6,2) + t320;
t445 = t28 * qJD(5);
t169 = t290 * t237;
t412 = t169 / 0.2e1;
t167 = t291 * t237;
t414 = t167 / 0.2e1;
t288 = Ifges(6,5) * t412 + Ifges(6,6) * t414;
t393 = pkin(7) * t237;
t396 = pkin(3) * t239;
t194 = t393 + t396;
t200 = t245 * t266 + t397 * t301;
t115 = t268 * t194 + t265 * t200;
t339 = t237 * t268;
t72 = t239 * pkin(4) + pkin(8) * t339 + t115;
t116 = t265 * t194 - t268 * t200;
t340 = t237 * t265;
t89 = pkin(8) * t340 + t116;
t44 = -t264 * t89 + t267 * t72;
t45 = t264 * t72 + t267 * t89;
t428 = t45 * mrSges(6,2) / 0.2e1 - t44 * mrSges(6,1) / 0.2e1 - t288;
t426 = Ifges(6,3) * t405 - t428;
t147 = pkin(4) * t338 + t200;
t444 = t147 / 0.2e1;
t356 = t267 * t83;
t41 = t264 * t69 + t356;
t48 = -t264 * t82 - t356;
t443 = t41 + t48;
t229 = t291 * mrSges(6,1);
t300 = -mrSges(6,2) * t290 + t229;
t440 = qJD(5) * t300;
t166 = t290 * t239;
t371 = t166 * mrSges(6,3);
t128 = mrSges(6,1) * t237 + t371;
t439 = t128 * t403 + t290 * t417;
t376 = Ifges(6,4) * t291;
t438 = t209 * mrSges(6,3) - t376;
t437 = -t128 / 0.2e1;
t256 = Ifges(5,5) * t268;
t374 = Ifges(5,6) * t265;
t289 = -t256 / 0.2e1 + t374 / 0.2e1;
t434 = Ifges(4,4) + t289;
t257 = Ifges(5,4) * t268;
t249 = Ifges(5,1) * t265 + t257;
t399 = -t268 / 0.2e1;
t433 = t249 * t399;
t293 = -Ifges(5,2) * t265 + t257;
t431 = -t115 * t265 + t116 * t268;
t430 = -mrSges(5,1) * t268 + t265 * mrSges(5,2);
t342 = t168 * t290;
t343 = t166 * t291;
t429 = (-t343 / 0.2e1 + t342 / 0.2e1) * mrSges(6,3) + t439;
t427 = -t439 - (t342 - t343) * mrSges(6,3) / 0.2e1;
t425 = 2 * m(6);
t424 = m(5) / 0.2e1;
t423 = m(6) / 0.2e1;
t422 = pkin(4) / 0.2e1;
t421 = -mrSges(5,1) / 0.2e1;
t420 = mrSges(5,2) / 0.2e1;
t416 = -t166 / 0.2e1;
t415 = t166 / 0.2e1;
t413 = -t168 / 0.2e1;
t195 = mrSges(6,1) * t290 + mrSges(6,2) * t291;
t411 = t195 / 0.2e1;
t410 = -t298 / 0.2e1;
t409 = t298 / 0.2e1;
t408 = -t209 / 0.2e1;
t407 = t237 / 0.2e1;
t406 = t237 / 0.4e1;
t404 = -t290 / 0.2e1;
t394 = pkin(4) * t268;
t255 = -pkin(3) - t394;
t402 = t255 / 0.2e1;
t401 = -t265 / 0.2e1;
t400 = t265 / 0.2e1;
t398 = t268 / 0.2e1;
t392 = t40 * mrSges(6,2);
t391 = t41 * mrSges(6,1);
t388 = t48 * mrSges(6,1);
t387 = t49 * mrSges(6,2);
t385 = Ifges(6,1) - Ifges(6,2);
t379 = Ifges(5,4) * t265;
t378 = Ifges(6,4) * t166;
t161 = Ifges(6,4) * t168;
t377 = Ifges(6,4) * t290;
t370 = t167 * mrSges(6,1);
t368 = t169 * mrSges(6,2);
t362 = t291 * mrSges(6,3);
t361 = t264 * t45;
t359 = t265 * mrSges(5,1);
t357 = t267 * t44;
t355 = t268 * mrSges(5,2);
t125 = -mrSges(6,2) * t239 + t167 * mrSges(6,3);
t127 = mrSges(6,1) * t239 - t169 * mrSges(6,3);
t140 = Ifges(5,6) * t239 - t237 * t293;
t294 = Ifges(5,1) * t268 - t379;
t142 = Ifges(5,5) * t239 - t237 * t294;
t148 = -pkin(4) * t340 + t201;
t247 = t355 + t359;
t181 = t247 * t237;
t182 = t247 * t239;
t189 = -mrSges(5,2) * t239 + mrSges(5,3) * t340;
t316 = mrSges(5,3) * t338;
t190 = -t237 * mrSges(5,2) - t316;
t191 = t239 * mrSges(5,1) + mrSges(5,3) * t339;
t192 = t237 * mrSges(5,1) - mrSges(5,3) * t337;
t227 = t237 * mrSges(4,2);
t143 = Ifges(5,5) * t237 + t239 * t294;
t326 = t268 * t143;
t141 = Ifges(5,6) * t237 + t239 * t293;
t331 = t265 * t141;
t73 = Ifges(6,4) * t169 + Ifges(6,2) * t167 + Ifges(6,6) * t239;
t74 = -Ifges(6,2) * t168 + Ifges(6,6) * t237 - t378;
t75 = Ifges(6,1) * t169 + Ifges(6,4) * t167 + Ifges(6,5) * t239;
t76 = -Ifges(6,1) * t166 + Ifges(6,5) * t237 - t161;
t91 = t368 - t370;
t92 = mrSges(6,1) * t168 - mrSges(6,2) * t166;
t3 = -t305 * t227 + t115 * t192 - t200 * t181 + t201 * t182 + t103 * t189 + t116 * t190 + t102 * t191 + t41 * t125 + t45 * t126 + t40 * t127 + t44 * t128 + t147 * t91 + t148 * t92 + t75 * t416 + t74 * t414 + t73 * t413 + t76 * t412 + m(6) * (t147 * t148 + t40 * t44 + t41 * t45) + m(5) * (t102 * t115 + t103 * t116 + t200 * t201) + (t305 * mrSges(4,1) + Ifges(6,5) * t416 + Ifges(6,6) * t413 + t140 * t401 + t142 * t398 - t239 * t434) * t239 + (t331 / 0.2e1 - t326 / 0.2e1 + (Ifges(5,3) + Ifges(6,3) - Ifges(4,1) + Ifges(4,2)) * t239 + t288 + t434 * t237) * t237;
t354 = t3 * qJD(1);
t248 = Ifges(5,2) * t268 + t379;
t323 = -Ifges(6,5) * t168 + Ifges(6,6) * t166;
t90 = -mrSges(6,1) * t166 - mrSges(6,2) * t168;
t283 = t147 * t90 + t323 * t407 + t41 * t371;
t292 = Ifges(5,5) * t265 + Ifges(5,6) * t268;
t341 = t200 * t239;
t349 = t103 * t268;
t93 = Ifges(6,2) * t166 - t161;
t94 = -Ifges(6,1) * t168 + t378;
t4 = -t49 * t126 - t48 * t128 - m(6) * (t40 * t48 + t41 * t49) - pkin(4) * t92 * t337 + t430 * t341 + t103 * t192 + (t94 / 0.2e1 - t74 / 0.2e1) * t166 - (t40 * mrSges(6,3) - t76 / 0.2e1 - t93 / 0.2e1) * t168 + (t143 * t400 + t141 * t398 - m(6) * t147 * t394 + t292 * t407 + mrSges(5,3) * t349 + (t248 * t401 - t433) * t239) * t239 - t283 + (-t190 - t316) * t102;
t353 = t4 * qJD(1);
t7 = -t128 * t41 + t74 * t415 + t94 * t416 + (t126 + t369) * t40 + t283 + (t76 + t93) * t413;
t352 = t7 * qJD(1);
t286 = t370 / 0.2e1 - t368 / 0.2e1;
t318 = m(6) * t422;
t275 = (t355 / 0.2e1 + t359 / 0.2e1) * t237 + (t167 * t267 + t169 * t264) * t318 + t286;
t325 = t268 * t190;
t330 = t265 * t192;
t276 = -m(6) * (-t443 * t290 - t446 * t291) / 0.2e1 + t330 / 0.2e1 - t325 / 0.2e1;
t8 = t275 + t276 + t429;
t351 = t8 * qJD(1);
t260 = t265 ^ 2;
t261 = t268 ^ 2;
t297 = mrSges(5,3) * (-t261 / 0.2e1 - t260 / 0.2e1);
t273 = (t169 * t404 - t291 * t414) * mrSges(6,3) + t237 * t297 + (-t396 + (-t260 - t261) * t393) * t424 + (t167 * t298 + t209 * t169 + t255 * t239) * t423;
t274 = (t268 * t115 + t265 * t116) * t424 + (-t290 * t44 + t291 * t45) * t423 + t127 * t404 + t125 * t403 + t189 * t400 + t191 * t398;
t302 = t411 + t430 / 0.2e1;
t13 = t227 + (-mrSges(4,1) + t302) * t239 + t273 - t274;
t350 = qJD(1) * t13;
t12 = t169 * t126 + t167 * t128 + (mrSges(4,3) * t239 + t182 + t92) * t239 + (mrSges(4,3) * t237 - t325 + t330) * t237 + m(6) * (t147 * t239 + t167 * t40 + t169 * t41) + m(5) * (t341 + (t102 * t265 - t349) * t237) + m(4) * (-t201 * t237 + t341) + (m(3) * qJ(2) + mrSges(3,3)) * (t262 ^ 2 + t263 ^ 2);
t346 = t12 * qJD(1);
t345 = t147 * t265;
t15 = t286 + t429;
t344 = t15 * qJD(1);
t334 = t264 * t125;
t333 = t264 * t166;
t329 = t267 * t127;
t328 = t267 * t168;
t322 = t209 * t362 - t255 * t300;
t313 = -Ifges(6,1) / 0.4e1 + Ifges(6,2) / 0.4e1;
t312 = -Ifges(5,2) / 0.4e1 + Ifges(5,1) / 0.4e1;
t311 = -t40 / 0.2e1 + t49 / 0.2e1;
t310 = t41 / 0.2e1 + t48 / 0.2e1;
t306 = t92 * t400;
t299 = t256 - t374;
t296 = t320 * t406;
t287 = t115 * t421 + t116 * t420;
t285 = -t168 * t410 + t209 * t415;
t284 = t248 * t400 + t433;
t272 = -(-t378 / 0.2e1 - t94 / 0.4e1 + t74 / 0.4e1 - t313 * t168) * t291 - (-t161 / 0.2e1 + mrSges(6,2) * t444 + t76 / 0.4e1 + t93 / 0.4e1 + t313 * t166) * t290 + t229 * t444 + t90 * t402;
t269 = (t190 * t401 + t192 * t399) * pkin(7) + (-t290 * t311 - t291 * t310 + t285) * mrSges(6,3) + t200 * t247 / 0.2e1 + t209 * t437 + t126 * t409 - t331 / 0.4e1 + t326 / 0.4e1 + t272;
t271 = pkin(7) * t297 + (pkin(3) * t420 - t249 / 0.4e1 - t257 / 0.4e1 - t312 * t265) * t265 + (pkin(3) * t421 - 0.3e1 / 0.4e1 * t379 - t248 / 0.4e1 + t312 * t268 + (m(6) * t402 + t411) * pkin(4)) * t268;
t279 = -t446 * t209 + t443 * t298;
t2 = t269 + t279 * t423 + (-t334 / 0.2e1 + t306 - t329 / 0.2e1 + (t345 / 0.4e1 - t361 / 0.4e1 - t357 / 0.4e1) * t425) * pkin(4) + (-Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + t271) * t239 + (-0.3e1 / 0.4e1 * t374 + 0.3e1 / 0.4e1 * t256 - t231 / 0.4e1 - t230 / 0.4e1) * t237 + t287 + t428;
t22 = t293 * t399 + pkin(3) * t247 + t294 * t401 - t438 * t291 - (-t291 * t385 + t377) * t290 + t284 + t322 + (-m(6) * t255 - t195) * pkin(4) * t265;
t281 = -t2 * qJD(1) + t22 * qJD(3);
t25 = -Ifges(6,4) * t290 ^ 2 - (-t290 * t385 + t438) * t291 + t322;
t270 = t285 * mrSges(6,3) + t128 * t408 + t298 * t417 + t272 + t296;
t6 = t270 - t426;
t280 = -t6 * qJD(1) + t25 * qJD(3);
t277 = (t267 * t417 + t264 * t437 + (t333 / 0.2e1 + t328 / 0.2e1) * mrSges(6,3)) * pkin(4);
t11 = -mrSges(6,1) * t310 + mrSges(6,2) * t311 + t277;
t244 = (mrSges(6,1) * t264 + mrSges(6,2) * t267) * pkin(4);
t35 = (t410 + t409) * mrSges(6,2) + (t408 + t209 / 0.2e1) * mrSges(6,1);
t278 = -t11 * qJD(1) - t35 * qJD(3) + t244 * qJD(4);
t240 = t244 * qJD(5);
t16 = t286 + t427;
t14 = t239 * t302 + t273 + t274;
t10 = -t392 / 0.2e1 - t391 / 0.2e1 - t387 / 0.2e1 + t388 / 0.2e1 + t277 + t323;
t9 = t275 - t276 + t427;
t5 = t270 + t426;
t1 = t269 + t271 * t239 + (t357 + t361) * t318 + t299 * t406 + Ifges(5,3) * t405 + t296 + (pkin(4) * t345 + t279) * t423 + pkin(4) * t306 - t287 + (t334 + t329) * t422 + t289 * t237 + t426;
t17 = [qJD(2) * t12 + qJD(3) * t3 - qJD(4) * t4 + qJD(5) * t7, t346 + (-t167 * t290 + t169 * t291) * m(6) * qJD(2) + t14 * qJD(3) + t9 * qJD(4) + t16 * qJD(5), t14 * qJD(2) + t1 * qJD(4) + t5 * qJD(5) + t354 + ((-Ifges(4,5) + t284) * t237 + (-Ifges(6,2) * t290 + t376) * t414 + 0.2e1 * (t148 * t255 + t209 * t45 + t298 * t44) * t423 + 0.2e1 * (-pkin(3) * t201 + pkin(7) * t431) * t424 + t75 * t403 + t73 * t404 + (Ifges(6,1) * t291 - t377) * t412 + t140 * t398 + t142 * t400 + t255 * t91 + t201 * t430 - Ifges(4,6) * t239 + t209 * t125 - t201 * mrSges(4,1) + t298 * t127 + t148 * t195 + t200 * mrSges(4,2) + pkin(3) * t181 - t265 * pkin(7) * t191 + t268 * pkin(7) * t189 - t44 * t362 - t45 * t290 * mrSges(6,3) + (Ifges(6,5) * t291 - Ifges(6,6) * t290 + t292) * t405 + t431 * mrSges(5,3)) * qJD(3), -t353 + t9 * qJD(2) + t1 * qJD(3) + (-t103 * mrSges(5,1) - t102 * mrSges(5,2) - Ifges(5,5) * t338 - Ifges(5,6) * t337 + t323 - t387 + t388) * qJD(4) + t10 * qJD(5) + (m(6) * (t264 * t49 + t267 * t48) + (t328 + t333) * mrSges(6,3)) * t372, t352 + t16 * qJD(2) + t5 * qJD(3) + t10 * qJD(4) + (t323 - t391 - t392) * qJD(5); -qJD(3) * t13 - qJD(4) * t8 - qJD(5) * t15 - t346, 0, -t350, -t351 + (-t247 - t300) * qJD(4) - t440 + (-t264 * t290 - t267 * t291) * t425 * t372 / 0.2e1, -qJD(4) * t300 - t344 - t440; qJD(2) * t13 + qJD(4) * t2 + qJD(5) * t6 - t354, t350, -qJD(4) * t22 - qJD(5) * t25, (pkin(7) * t430 + t28 + t299) * qJD(4) + t445 + (m(6) * (-t209 * t267 + t264 * t298) + (-t264 * t291 + t267 * t290) * mrSges(6,3)) * t372 - t281, t28 * qJD(4) - t280 + t445; qJD(2) * t8 - qJD(3) * t2 + qJD(5) * t11 + t353, t351, t35 * qJD(5) + t281, -t240, -t240 - t278; qJD(2) * t15 - qJD(3) * t6 - qJD(4) * t11 - t352, t344, -t35 * qJD(4) + t280, t278, 0;];
Cq = t17;
