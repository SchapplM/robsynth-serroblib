% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR7
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:15
% EndTime: 2019-12-31 19:03:25
% DurationCPUTime: 4.72s
% Computational Cost: add. (9037->381), mult. (18955->538), div. (0->0), fcn. (17973->8), ass. (0->214)
t253 = sin(qJ(3));
t397 = t253 / 0.2e1;
t251 = sin(qJ(5));
t252 = sin(qJ(4));
t254 = cos(qJ(5));
t255 = cos(qJ(4));
t331 = t254 * t255;
t218 = -t251 * t252 + t331;
t219 = -t251 * t255 - t254 * t252;
t135 = -mrSges(6,1) * t219 + mrSges(6,2) * t218;
t378 = Ifges(6,4) * t219;
t138 = Ifges(6,2) * t218 - t378;
t244 = -pkin(4) * t255 - pkin(3);
t296 = Ifges(6,1) * t218 + t378;
t401 = -t219 / 0.2e1;
t210 = Ifges(6,4) * t218;
t139 = -Ifges(6,1) * t219 + t210;
t442 = Ifges(6,2) * t219 + t139 + t210;
t22 = -t244 * t135 + t138 * t401 + t219 * t296 / 0.2e1 - t442 * t218 / 0.2e1;
t256 = cos(qJ(3));
t309 = -cos(pkin(9)) * pkin(1) - pkin(2);
t389 = t253 * pkin(7);
t207 = -pkin(3) * t256 + t309 - t389;
t196 = t255 * t207;
t241 = sin(pkin(9)) * pkin(1) + pkin(6);
t335 = t252 * t256;
t129 = -t241 * t335 + t196;
t334 = t253 * t255;
t321 = pkin(8) * t334;
t117 = t129 - t321;
t329 = t255 * t256;
t130 = t252 * t207 + t241 * t329;
t336 = t252 * t253;
t118 = -pkin(8) * t336 + t130;
t339 = t251 * t118;
t287 = t254 * t117 - t339;
t107 = -t321 + t196 + (-t241 * t252 - pkin(4)) * t256;
t48 = t107 * t254 - t339;
t445 = t287 - t48;
t193 = t219 * t253;
t175 = Ifges(6,4) * t193;
t192 = t251 * t336 - t253 * t331;
t110 = Ifges(6,2) * t192 + t175;
t99 = -Ifges(6,1) * t192 - t256 * Ifges(6,5) + t175;
t444 = t99 + t110;
t137 = Ifges(6,5) * t218 + Ifges(6,6) * t219;
t413 = -pkin(8) - pkin(7);
t230 = t413 * t252;
t231 = t413 * t255;
t151 = t230 * t251 - t231 * t254;
t304 = t254 * t230 + t231 * t251;
t26 = -t151 * mrSges(6,1) - t304 * mrSges(6,2) + t137;
t443 = t26 * qJD(5);
t194 = t219 * t256;
t195 = t218 * t256;
t281 = -Ifges(6,5) * t195 / 0.2e1 - Ifges(6,6) * t194 / 0.2e1;
t390 = t253 * pkin(3);
t391 = pkin(7) * t256;
t232 = t390 - t391;
t152 = t255 * t232 + t241 * t336;
t121 = t253 * pkin(4) - pkin(8) * t329 + t152;
t153 = t252 * t232 - t241 * t334;
t128 = -pkin(8) * t335 + t153;
t57 = t121 * t254 - t128 * t251;
t58 = t121 * t251 + t128 * t254;
t422 = Ifges(6,3) * t397 - t58 * mrSges(6,2) / 0.2e1 + t57 * mrSges(6,1) / 0.2e1 - t281;
t108 = t192 * mrSges(6,1) - t193 * mrSges(6,2);
t441 = qJD(5) * t108;
t415 = m(6) * pkin(4);
t438 = -t256 / 0.4e1;
t357 = t218 * mrSges(6,3);
t355 = t219 * mrSges(6,3);
t333 = t254 * t118;
t286 = -t251 * t117 - t333;
t49 = t107 * t251 + t333;
t436 = t286 + t49;
t155 = -t253 * mrSges(6,2) + mrSges(6,3) * t194;
t157 = t253 * mrSges(6,1) - mrSges(6,3) * t195;
t373 = Ifges(5,6) * t252;
t377 = Ifges(5,5) * t255;
t282 = t377 / 0.2e1 - t373 / 0.2e1;
t364 = t153 * mrSges(5,2);
t365 = t152 * mrSges(5,1);
t396 = t254 / 0.2e1;
t417 = m(6) / 0.2e1;
t434 = Ifges(5,3) * t397 + t365 / 0.2e1 - t364 / 0.2e1 + t422 + (t157 * t396 + (t251 * t58 + t254 * t57) * t417 + t251 * t155 / 0.2e1) * pkin(4) + t282 * t256;
t360 = t193 * mrSges(6,3);
t154 = mrSges(6,2) * t256 + t360;
t362 = t192 * mrSges(6,3);
t156 = -mrSges(6,1) * t256 + t362;
t409 = -t156 / 0.2e1;
t433 = t151 * t409 + t304 * t154 / 0.2e1;
t431 = 0.2e1 * t253;
t248 = t252 ^ 2;
t249 = t255 ^ 2;
t324 = t248 + t249;
t423 = t324 * mrSges(5,3);
t429 = t253 * t256;
t352 = t255 * mrSges(5,2);
t353 = t252 * mrSges(5,1);
t280 = -t352 / 0.2e1 - t353 / 0.2e1;
t428 = t280 * t256;
t424 = -mrSges(5,1) * t255 + mrSges(5,2) * t252;
t277 = t424 * t253;
t316 = mrSges(5,3) * t336;
t283 = mrSges(5,2) * t256 - t316;
t284 = -mrSges(5,1) * t256 - mrSges(5,3) * t334;
t425 = t252 * t283 + t255 * t284;
t285 = -t152 * t252 + t153 * t255;
t358 = t195 * mrSges(6,2);
t359 = t194 * mrSges(6,1);
t326 = t359 / 0.2e1 - t358 / 0.2e1;
t379 = Ifges(6,4) * t192;
t111 = Ifges(6,1) * t193 + t379;
t392 = pkin(4) * t252;
t306 = t241 + t392;
t200 = t306 * t253;
t97 = Ifges(6,2) * t193 - t256 * Ifges(6,6) - t379;
t420 = t151 * t362 / 0.2e1 - t304 * t360 / 0.2e1 + t137 * t438 - t244 * t108 / 0.2e1 + t200 * t135 / 0.2e1 + t442 * t193 / 0.4e1 + t444 * t218 / 0.4e1 + (-t111 / 0.4e1 + t97 / 0.4e1) * t219 + (-t296 / 0.4e1 + t138 / 0.4e1) * t192;
t419 = t253 ^ 2;
t418 = m(5) / 0.2e1;
t408 = -t192 / 0.2e1;
t407 = t192 / 0.2e1;
t406 = t193 / 0.2e1;
t405 = t194 / 0.2e1;
t404 = t195 / 0.2e1;
t300 = t352 + t353;
t206 = t256 * t300;
t403 = -t206 / 0.2e1;
t400 = -t252 / 0.2e1;
t398 = t252 / 0.2e1;
t395 = t255 / 0.2e1;
t393 = t256 / 0.2e1;
t388 = t48 * mrSges(6,2);
t387 = t49 * mrSges(6,1);
t384 = Ifges(5,1) - Ifges(5,2);
t381 = Ifges(5,4) * t252;
t380 = Ifges(5,4) * t255;
t376 = Ifges(5,5) * t256;
t374 = Ifges(5,2) * t252;
t372 = Ifges(5,6) * t256;
t363 = t192 * mrSges(6,2);
t361 = t193 * mrSges(6,1);
t351 = t424 - mrSges(4,1);
t13 = t108 * t393 + t154 * t406 + t156 * t407 - (t192 ^ 2 + t193 ^ 2) * mrSges(6,3) / 0.2e1;
t322 = pkin(4) * t334;
t10 = -m(6) * (-t445 * t192 + t436 * t193 - t256 * t322) / 0.2e1 - t277 * t393 - t13 + t425 * t397 + t419 * t423 / 0.2e1;
t350 = t10 * qJD(1);
t349 = t13 * qJD(1);
t341 = t241 * t253;
t338 = t251 * t192;
t224 = t253 * mrSges(5,1) - mrSges(5,3) * t329;
t337 = t252 * t224;
t332 = t254 * t193;
t223 = -t253 * mrSges(5,2) - mrSges(5,3) * t335;
t330 = t255 * t223;
t328 = t256 * t135;
t325 = Ifges(6,5) * t193 + Ifges(6,6) * t192;
t323 = qJD(3) * t256;
t312 = -t355 / 0.2e1;
t311 = t355 / 0.2e1;
t310 = -t328 / 0.2e1 + (t311 + t312) * t192;
t302 = t49 * t312;
t299 = -t361 - t363;
t298 = Ifges(5,1) * t255 - t381;
t297 = Ifges(5,1) * t252 + t380;
t295 = -t374 + t380;
t294 = Ifges(5,2) * t255 + t381;
t293 = -t373 + t377;
t100 = Ifges(6,1) * t195 + Ifges(6,4) * t194 + Ifges(6,5) * t253;
t109 = t358 - t359;
t190 = Ifges(5,6) * t253 + t295 * t256;
t191 = Ifges(5,5) * t253 + t298 * t256;
t201 = t306 * t256;
t98 = Ifges(6,4) * t195 + Ifges(6,2) * t194 + Ifges(6,6) * t253;
t4 = t129 * t224 + t130 * t223 + t99 * t404 + t200 * t109 + t201 * t299 + t100 * t408 + t98 * t406 + t97 * t405 + t57 * t156 + t48 * t157 + t58 * t154 + t49 * t155 + m(5) * (t129 * t152 + t130 * t153) + m(6) * (t200 * t201 + t48 * t57 + t49 * t58) + (t191 * t395 + t190 * t400 + t241 * t206 + Ifges(6,5) * t408 + Ifges(6,6) * t406 + t309 * mrSges(4,1) + (-Ifges(4,4) + t282) * t253 + (-t152 * t255 - t153 * t252) * mrSges(5,3)) * t253 + (-t365 + t364 + t309 * mrSges(4,2) + (Ifges(5,1) * t249 / 0.2e1 + Ifges(4,1) - Ifges(4,2) - Ifges(6,3) - Ifges(5,3) + (m(5) * t241 + t352) * t241 + (t241 * mrSges(5,1) - t380 + t374 / 0.2e1) * t252) * t253 + t281 + (Ifges(4,4) - t293) * t256) * t256;
t9 = t156 * t405 + t157 * t406 + t154 * t404 + t155 * t408 + (t330 / 0.2e1 - t337 / 0.2e1 - t361 / 0.2e1 - t363 / 0.2e1 - t280 * t253) * t253 + (-t109 / 0.2e1 + t403 - t428) * t256 + (-t58 * t192 + t57 * t193 + t48 * t194 + t49 * t195 + t200 * t253 - t201 * t256) * t417 + ((-t129 * t252 + t130 * t255 - t241 * t256) * t256 + (t285 + t341) * t253) * t418;
t292 = t4 * qJD(1) + t9 * qJD(2);
t274 = t49 * t362 - t200 * t108 - t256 * t325 / 0.2e1;
t279 = pkin(4) * t299;
t5 = t419 * t241 * t424 - t287 * t154 - t286 * t156 - m(6) * (t48 * t286 + t49 * t287) - t279 * t334 + t130 * t284 + (t111 / 0.2e1 - t97 / 0.2e1) * t192 + (t48 * mrSges(6,3) - t99 / 0.2e1 - t110 / 0.2e1) * t193 + ((-Ifges(5,4) * t336 - t376) * t252 + (t130 * mrSges(5,3) + Ifges(5,4) * t334 - t200 * t415 + t384 * t336 - t372) * t255) * t253 - t274 + (-t283 - t316) * t129;
t290 = -t5 * qJD(1) - t10 * qJD(2);
t8 = -t156 * t49 + t97 * t407 + t111 * t408 + (t154 - t360) * t48 + t274 + t444 * t406;
t289 = t8 * qJD(1) + t13 * qJD(2);
t35 = m(5) * (-0.1e1 + t324) * t429 + m(6) * (-t192 * t195 + t193 * t194 - t429);
t288 = t9 * qJD(1) + t35 * qJD(2);
t273 = t287 * mrSges(6,2);
t272 = t286 * mrSges(6,1);
t136 = -mrSges(6,1) * t218 - mrSges(6,2) * t219;
t14 = pkin(3) * t353 + Ifges(5,4) * t248 + (pkin(3) * mrSges(5,2) - t384 * t252 - t380) * t255 + (-m(6) * t244 - t136) * t392 + t22;
t257 = ((t200 * t252 + t244 * t334) * pkin(4) + t436 * t304 + t445 * t151) * t417 + t293 * t438 + pkin(3) * t277 / 0.2e1 + t279 * t398 - t297 * t336 / 0.2e1 + t300 * t341 / 0.2e1 - t294 * t334 / 0.2e1 + t286 * t311 + t136 * t322 / 0.2e1 - (t295 * t431 - t372) * t252 / 0.4e1 + (t298 * t431 - t376) * t255 / 0.4e1 - t425 * pkin(7) / 0.2e1 + (-t48 / 0.2e1 + t287 / 0.2e1) * t357 - t389 * t423 / 0.2e1 + t433;
t2 = -t257 + t302 - t420 + t434;
t263 = t335 * t415;
t266 = (t194 * t254 + t195 * t251) * t415 / 0.2e1 + t326;
t20 = t328 / 0.2e1 + t263 / 0.2e1 + t266;
t268 = -t2 * qJD(1) - t20 * qJD(2) - t14 * qJD(3);
t25 = t310 - t326;
t264 = t49 * t311 + t420;
t258 = t264 + t302 + t433;
t6 = t258 - t422;
t267 = t6 * qJD(1) + t25 * qJD(2) - t22 * qJD(3);
t259 = (t154 * t396 + t251 * t409 + (t338 / 0.2e1 - t332 / 0.2e1) * mrSges(6,3)) * pkin(4) - t388 / 0.2e1 - t387 / 0.2e1;
t261 = t273 / 0.2e1 - t272 / 0.2e1;
t12 = t259 + t261;
t226 = (mrSges(6,1) * t251 + mrSges(6,2) * t254) * pkin(4);
t265 = -t12 * qJD(1) + t226 * qJD(4);
t217 = t226 * qJD(5);
t24 = t310 + t326;
t21 = -t263 / 0.2e1 + t403 + t428 + t266 + t310;
t11 = t259 - t261 + t325;
t7 = t258 + t422;
t3 = qJD(3) * t9 - qJD(4) * t10 + qJD(5) * t13;
t1 = t257 + t264 + t434;
t15 = [qJD(3) * t4 - qJD(4) * t5 + qJD(5) * t8, t3, t1 * qJD(4) + t7 * qJD(5) + (t297 * t395 + t294 * t400 + Ifges(4,5) + (-m(5) * pkin(3) + t351) * t241) * t323 + t292 + (mrSges(4,2) * t341 + m(6) * (t151 * t58 + t201 * t244 + t304 * t57) + t57 * t355 + t58 * t357 + t190 * t395 - Ifges(4,6) * t253 + t191 * t398 + t244 * t109 + t218 * t98 / 0.2e1 + t100 * t401 - pkin(3) * t206 + t139 * t404 + t201 * t136 + t138 * t405 + t304 * t157 + t151 * t155 + (m(5) * t285 + t330 - t337) * pkin(7) + (Ifges(5,5) * t252 - Ifges(6,5) * t219 + Ifges(5,6) * t255 + Ifges(6,6) * t218) * t397 + t285 * mrSges(5,3)) * qJD(3), t1 * qJD(3) + (-Ifges(5,5) * t336 - Ifges(5,6) * t334 - t273 + t272 + (-t251 ^ 2 - t254 ^ 2) * t118 * t415 - t129 * mrSges(5,2) - t130 * mrSges(5,1) + t325 + (t338 - t332) * mrSges(6,3) * pkin(4)) * qJD(4) + t11 * qJD(5) + t290, t7 * qJD(3) + t11 * qJD(4) + (t325 - t387 - t388) * qJD(5) + t289; t3, t35 * qJD(3), t21 * qJD(4) + t24 * qJD(5) + (-mrSges(4,2) + t423) * t323 + t288 + ((t194 * t219 + t195 * t218) * mrSges(6,3) + (t136 + t351) * t253 + 0.2e1 * (t151 * t195 + t194 * t304 + t244 * t253) * t417 + 0.2e1 * (t324 * t391 - t390) * t418) * qJD(3), -t350 + t21 * qJD(3) + (t277 + t108 + (t192 * t254 + t193 * t251) * t415) * qJD(4) + t441, qJD(3) * t24 + qJD(4) * t108 + t349 + t441; -qJD(4) * t2 + qJD(5) * t6 - t292, -qJD(4) * t20 + qJD(5) * t25 - t288, -qJD(4) * t14 - qJD(5) * t22, t443 + t268 + (t293 + (m(6) * (-t151 * t254 + t251 * t304) + (-t254 * t218 + t251 * t219) * mrSges(6,3)) * pkin(4) + t424 * pkin(7) + t26) * qJD(4), t26 * qJD(4) + t267 + t443; qJD(3) * t2 + qJD(5) * t12 - t290, t20 * qJD(3) + t350, -t268, -t217, -t217 - t265; -qJD(3) * t6 - qJD(4) * t12 - t289, -t25 * qJD(3) - t349, -t267, t265, 0;];
Cq = t15;
