% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR13_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:35
% EndTime: 2019-12-31 19:14:46
% DurationCPUTime: 5.08s
% Computational Cost: add. (8771->445), mult. (18187->620), div. (0->0), fcn. (17203->6), ass. (0->254)
t255 = sin(qJ(3));
t254 = sin(qJ(4));
t256 = cos(qJ(4));
t402 = sin(qJ(5));
t403 = cos(qJ(5));
t448 = t402 * t254 - t403 * t256;
t186 = t448 * t255;
t257 = cos(qJ(3));
t145 = mrSges(6,1) * t257 - mrSges(6,3) * t186;
t233 = pkin(3) * t257 + pkin(7) * t255;
t215 = t256 * t233;
t258 = -pkin(1) - pkin(6);
t340 = t257 * t258;
t156 = -t254 * t340 + t215;
t157 = t254 * t233 + t256 * t340;
t328 = pkin(4) * t403;
t292 = t328 / 0.2e1;
t387 = Ifges(5,5) * t255;
t321 = -t387 / 0.2e1;
t327 = pkin(4) * t402;
t437 = m(6) * pkin(4);
t332 = t437 / 0.2e1;
t349 = t254 * t255;
t404 = t257 / 0.2e1;
t276 = t403 * t254 + t402 * t256;
t184 = t276 * t255;
t143 = -mrSges(6,2) * t257 + mrSges(6,3) * t184;
t423 = t143 / 0.2e1;
t436 = -mrSges(5,2) / 0.2e1;
t416 = t186 / 0.2e1;
t419 = t184 / 0.2e1;
t279 = Ifges(6,5) * t416 + Ifges(6,6) * t419;
t298 = -t254 * t258 + pkin(4);
t347 = t255 * t256;
t119 = pkin(8) * t347 + t298 * t257 + t215;
t138 = pkin(8) * t349 + t157;
t57 = t403 * t119 - t402 * t138;
t58 = t402 * t119 + t403 * t138;
t441 = Ifges(6,3) * t404 + t57 * mrSges(6,1) / 0.2e1 - t58 * mrSges(6,2) / 0.2e1 + t279;
t468 = Ifges(5,3) * t404 + t156 * mrSges(5,1) / 0.2e1 + t157 * t436 + (t402 * t58 + t403 * t57) * t332 + t256 * t321 + Ifges(5,6) * t349 / 0.2e1 + t145 * t292 + t327 * t423 + t441;
t371 = t256 * mrSges(5,1);
t372 = t254 * mrSges(5,2);
t374 = t276 * mrSges(6,2);
t376 = t448 * mrSges(6,1);
t452 = -t376 / 0.2e1 - t374 / 0.2e1;
t467 = t372 / 0.2e1 - t371 / 0.2e1 - (t276 * t402 - t403 * t448) * t332 - t452;
t398 = pkin(7) * t257;
t401 = pkin(3) * t255;
t223 = qJ(2) - t398 + t401;
t207 = t256 * t223;
t342 = t256 * t257;
t331 = pkin(8) * t342;
t116 = t298 * t255 + t207 - t331;
t346 = t255 * t258;
t154 = t254 * t223 + t256 * t346;
t348 = t254 * t257;
t125 = -pkin(8) * t348 + t154;
t312 = t402 * t125;
t53 = t403 * t116 - t312;
t153 = -t254 * t346 + t207;
t124 = t153 - t331;
t60 = t403 * t124 - t312;
t466 = t53 - t60;
t252 = t254 ^ 2;
t253 = t256 ^ 2;
t334 = t252 + t253;
t465 = mrSges(5,3) * t334;
t432 = -pkin(8) - pkin(7);
t232 = t432 * t254;
t234 = t432 * t256;
t139 = -t402 * t232 + t403 * t234;
t284 = t403 * t232 + t402 * t234;
t335 = -Ifges(6,5) * t448 - Ifges(6,6) * t276;
t28 = t139 * mrSges(6,1) - t284 * mrSges(6,2) + t335;
t464 = t28 * qJD(5);
t183 = t448 * t257;
t185 = t276 * t257;
t105 = -mrSges(6,1) * t183 - mrSges(6,2) * t185;
t144 = -mrSges(6,2) * t255 - t185 * mrSges(6,3);
t146 = mrSges(6,1) * t255 + t183 * mrSges(6,3);
t357 = t186 * t183;
t358 = t184 * t185;
t405 = -t257 / 0.2e1;
t420 = -t184 / 0.2e1;
t263 = t105 * t405 + t144 * t420 + t146 * t416 + (-t357 / 0.2e1 - t358 / 0.2e1) * mrSges(6,3);
t296 = t186 * mrSges(6,1) + t184 * mrSges(6,2);
t463 = qJD(5) * t296;
t126 = mrSges(6,1) * t276 - mrSges(6,2) * t448;
t400 = pkin(4) * t254;
t297 = -t258 + t400;
t209 = t297 * t257;
t407 = t255 / 0.4e1;
t399 = pkin(4) * t256;
t245 = -pkin(3) - t399;
t410 = t245 / 0.2e1;
t424 = t139 / 0.2e1;
t461 = t284 / 0.2e1;
t462 = t146 * t424 + t144 * t461 + t209 * t126 / 0.2e1 + t105 * t410 + t335 * t407;
t315 = t403 * t125;
t54 = t402 * t116 + t315;
t59 = -t402 * t124 - t315;
t460 = t59 + t54;
t446 = t372 - t371;
t197 = t446 * t257;
t438 = m(6) / 0.2e1;
t458 = (-t257 ^ 2 * t399 - t460 * t184 + t186 * t466) * t438 - t197 * t405 + t263;
t251 = Ifges(5,4) * t256;
t228 = Ifges(5,1) * t254 + t251;
t201 = t228 * t257;
t250 = Ifges(5,5) * t256;
t385 = Ifges(5,6) * t254;
t295 = t250 - t385;
t388 = Ifges(6,4) * t276;
t129 = -Ifges(6,2) * t448 + t388;
t130 = -Ifges(6,1) * t448 - t388;
t300 = t129 / 0.4e1 - t130 / 0.4e1;
t217 = -mrSges(5,2) * t255 - mrSges(5,3) * t348;
t352 = t254 * t217;
t304 = -t352 / 0.2e1;
t384 = Ifges(5,6) * t255;
t451 = -Ifges(5,2) * t254 + t251;
t180 = t257 * t451 + t384;
t353 = t254 * t180;
t377 = t185 * mrSges(6,1);
t378 = t183 * mrSges(6,2);
t107 = t377 - t378;
t430 = t107 / 0.2e1;
t456 = pkin(3) * t197 / 0.2e1 - t353 / 0.4e1 - t254 * t201 / 0.4e1 + t295 * t407 + t400 * t430 + pkin(7) * t304 + (t139 * t466 + t209 * t400 + t460 * t284) * t438 + t300 * t183 + t462;
t435 = -mrSges(6,3) / 0.2e1;
t373 = t276 * mrSges(6,3);
t454 = t255 * t257;
t165 = Ifges(6,4) * t185;
t100 = -Ifges(6,1) * t183 + Ifges(6,5) * t255 - t165;
t108 = Ifges(6,2) * t183 - t165;
t389 = Ifges(6,4) * t183;
t109 = -Ifges(6,1) * t185 + t389;
t98 = -Ifges(6,2) * t185 + t255 * Ifges(6,6) - t389;
t288 = (-t109 / 0.4e1 + t98 / 0.4e1) * t276;
t450 = -(t100 / 0.4e1 + t108 / 0.4e1) * t448 - t288;
t361 = t139 * t183;
t362 = t284 * t185;
t449 = -t361 / 0.2e1 + t362 / 0.2e1;
t280 = -t156 * t254 + t157 * t256;
t390 = Ifges(5,4) * t254;
t226 = Ifges(5,2) * t256 + t390;
t229 = Ifges(5,1) * t256 - t390;
t445 = t229 / 0.4e1 - t226 / 0.4e1;
t219 = mrSges(5,1) * t255 - mrSges(5,3) * t342;
t343 = t256 * t219;
t302 = -t343 / 0.2e1;
t444 = t405 * t465 + t302;
t337 = -t377 / 0.2e1 + t378 / 0.2e1;
t439 = m(5) / 0.2e1;
t434 = t60 / 0.2e1;
t127 = t374 + t376;
t428 = t127 / 0.2e1;
t425 = -t284 / 0.2e1;
t421 = -t183 / 0.2e1;
t418 = -t185 / 0.2e1;
t225 = t254 * mrSges(5,1) + t256 * mrSges(5,2);
t199 = t257 * t225;
t415 = t199 / 0.2e1;
t200 = t257 * t226;
t414 = -t200 / 0.4e1;
t409 = -t254 / 0.2e1;
t408 = t254 / 0.2e1;
t406 = t256 / 0.2e1;
t397 = t53 * mrSges(6,2);
t396 = t54 * mrSges(6,1);
t393 = t59 * mrSges(6,1);
t392 = t60 * mrSges(6,2);
t375 = t448 * mrSges(6,3);
t106 = -mrSges(6,1) * t184 + mrSges(6,2) * t186;
t179 = Ifges(5,6) * t257 - t255 * t451;
t181 = Ifges(5,5) * t257 - t229 * t255;
t198 = t225 * t255;
t208 = t297 * t255;
t216 = -mrSges(5,2) * t257 + mrSges(5,3) * t349;
t218 = mrSges(5,1) * t257 + mrSges(5,3) * t347;
t277 = t250 / 0.2e1 - t385 / 0.2e1 - Ifges(4,4);
t182 = t257 * t229 + t387;
t345 = t256 * t182;
t97 = Ifges(6,4) * t186 + Ifges(6,2) * t184 + Ifges(6,6) * t257;
t99 = Ifges(6,1) * t186 + Ifges(6,4) * t184 + Ifges(6,5) * t257;
t3 = t153 * t218 + t156 * t219 + t154 * t216 + t157 * t217 - t208 * t107 + t209 * t106 + t100 * t416 + t54 * t143 + t58 * t144 + t53 * t145 + t57 * t146 + t99 * t421 + t98 * t419 + t97 * t418 + m(5) * (t153 * t156 + t154 * t157) + m(6) * (-t208 * t209 + t53 * t57 + t54 * t58) + (qJ(2) * mrSges(4,1) + Ifges(6,5) * t421 + Ifges(6,6) * t418 + t179 * t409 + t181 * t406 + t258 * t198 + t277 * t257) * t257 + (t258 * t199 - t345 / 0.2e1 + t353 / 0.2e1 - qJ(2) * mrSges(4,2) - t277 * t255 + (-m(5) * t258 ^ 2 - Ifges(4,1) + Ifges(4,2) + Ifges(5,3) + Ifges(6,3)) * t257 + t279) * t255;
t370 = t3 * qJD(1);
t369 = t53 * t185;
t368 = t53 * t448;
t336 = -Ifges(6,5) * t185 + Ifges(6,6) * t183;
t265 = t209 * t105 - (t100 / 0.2e1 + t108 / 0.2e1) * t185 + (-t109 / 0.2e1 + t98 / 0.2e1 + t54 * mrSges(6,3)) * t183 + mrSges(6,3) * t369 + t255 * t336 / 0.2e1;
t6 = t59 * t146 + m(6) * (t53 * t59 + t54 * t60) + t60 * t144 + t153 * t217 - t154 * t219 + (t258 * t197 + (t321 + t153 * mrSges(5,3) - t182 / 0.2e1 + t200 / 0.2e1) * t254 + (-t384 / 0.2e1 - t154 * mrSges(5,3) - t201 / 0.2e1 - t180 / 0.2e1 + (m(6) * t209 + t107) * pkin(4)) * t256) * t257 + t265;
t367 = t6 * qJD(1);
t7 = t53 * t144 - t54 * t146 + t265;
t366 = t7 * qJD(1);
t365 = t446 - mrSges(4,1);
t11 = (t304 + t444) * t255 + t458 + t467;
t364 = t11 * qJD(1);
t13 = t263 - t452;
t363 = t13 * qJD(1);
t356 = t186 * t276;
t25 = t255 * mrSges(4,1) + t257 * mrSges(4,2) + t276 * t144 - t448 * t146 + t352 + t343 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(6) * (t276 * t54 - t368) + m(5) * (t153 * t256 + t154 * t254);
t354 = t25 * qJD(1);
t351 = t254 * t218;
t350 = t254 * t226;
t344 = t256 * t216;
t341 = t257 * t126;
t333 = qJD(3) * t255;
t324 = t59 / 0.2e1 + t54 / 0.2e1;
t323 = t434 - t53 / 0.2e1;
t317 = -t373 / 0.2e1;
t316 = -t341 / 0.2e1 - t186 * t317 + t356 * t435;
t314 = t403 * t185;
t311 = t402 * t146;
t310 = t402 * t183;
t301 = -t258 * t225 / 0.2e1;
t204 = Ifges(6,4) * t448;
t128 = -Ifges(6,2) * t276 - t204;
t131 = Ifges(6,1) * t276 - t204;
t299 = t131 / 0.4e1 + t128 / 0.4e1;
t294 = mrSges(6,3) * t328;
t293 = mrSges(6,3) * t327;
t289 = mrSges(5,3) * (-t253 / 0.2e1 - t252 / 0.2e1);
t286 = t185 * t294;
t285 = t183 * t293;
t35 = m(5) * (-0.1e1 + t334) * t454 + (t357 + t358 - t454) * m(6);
t259 = (-t106 / 0.2e1 + t217 * t406 + t219 * t409 + t198 / 0.2e1) * t257 + (t430 + t344 / 0.2e1 - t351 / 0.2e1 + t415) * t255 + ((-t153 * t254 + t154 * t256) * t257 + (t280 - 0.2e1 * t340) * t255) * t439 + (-t54 * t183 - t184 * t57 - t186 * t58 + t257 * t208 + t209 * t255 - t369) * t438 + t144 * t421 + t145 * t420 + t146 * t418 - t186 * t423;
t278 = m(6) * (-t139 * t276 - t284 * t448);
t9 = -t278 / 0.2e1 + t259;
t281 = t9 * qJD(1) + t35 * qJD(2);
t20 = t245 * t126 - (-t130 / 0.2e1 + t129 / 0.2e1) * t276 - (t131 / 0.2e1 + t128 / 0.2e1) * t448;
t14 = -pkin(3) * t225 + t229 * t408 - t350 / 0.2e1 + (t228 / 0.2e1 + t451 / 0.2e1) * t256 + t20 + (m(6) * t245 + t127) * t400;
t266 = (-t314 - t310) * t332 - mrSges(5,1) * t348 / 0.2e1 + t342 * t436 + t337;
t267 = t332 * t348 + t415;
t18 = t341 / 0.2e1 + t266 + t267;
t2 = -(t108 + t100) * t448 / 0.4e1 - (t131 + t128) * t185 / 0.4e1 - (t228 + t451) * t348 / 0.4e1 + ((t245 * t438 + t428) * pkin(4) + t445) * t342 - t375 * t434 - t368 * t435 + t256 * t414 + t460 * t317 + t444 * pkin(7) + t449 * mrSges(6,3) + t257 * t301 - t288 + t345 / 0.4e1 + t456 - t468;
t275 = t2 * qJD(1) - t18 * qJD(2) + t14 * qJD(3);
t27 = t316 - t337;
t260 = (t139 * t435 + t300) * t183 - (mrSges(6,3) * t425 + t299) * t185 + t450 + t462;
t5 = t260 - t441;
t274 = t5 * qJD(1) + t27 * qJD(2) + t20 * qJD(3);
t16 = -t323 * mrSges(6,2) + t324 * mrSges(6,1) + (-t403 * t144 / 0.2e1 + t311 / 0.2e1 + (-t310 / 0.2e1 - t314 / 0.2e1) * mrSges(6,3)) * pkin(4);
t220 = (mrSges(6,1) * t402 + mrSges(6,2) * t403) * pkin(4);
t29 = (t425 + t461) * mrSges(6,2) + (t424 - t139 / 0.2e1) * mrSges(6,1);
t269 = t16 * qJD(1) - t29 * qJD(3) + t220 * qJD(4);
t210 = t220 * qJD(5);
t26 = t316 + t337;
t19 = t266 - t267 + t316;
t15 = -t396 / 0.2e1 + t144 * t292 + t285 / 0.2e1 - pkin(4) * t311 / 0.2e1 + t286 / 0.2e1 - t397 / 0.2e1 + t393 / 0.2e1 - t392 / 0.2e1 + t336;
t12 = t263 + t452;
t10 = (t257 * t289 + t302 + t304) * t255 + t458 - t467;
t8 = t278 / 0.2e1 + t259;
t4 = t260 + t441;
t1 = (-t276 * t324 - t323 * t448 + t449) * mrSges(6,3) + (-pkin(7) * t219 / 0.2e1 + t414 + t182 / 0.4e1) * t256 - t299 * t185 + (t301 + (-t228 / 0.4e1 - t451 / 0.4e1) * t254 + pkin(7) * t289 + ((m(6) * t410 + t428) * pkin(4) + t445) * t256) * t257 + t450 + t456 + t468;
t17 = [qJD(2) * t25 + qJD(3) * t3 + qJD(4) * t6 + qJD(5) * t7, t354 + m(6) * (t184 * t448 - t356) * qJD(2) + t8 * qJD(3) + t10 * qJD(4) + t12 * qJD(5), t370 + t8 * qJD(2) + t1 * qJD(4) + t4 * qJD(5) + (-t256 * t228 / 0.2e1 + t350 / 0.2e1 - Ifges(4,5) + (-m(5) * pkin(3) + t365) * t258) * t333 + (-mrSges(4,2) * t340 + m(6) * (-t139 * t58 - t208 * t245 + t284 * t57) - t58 * t375 - t57 * t373 + t179 * t406 - Ifges(4,6) * t257 + t181 * t408 + t245 * t106 - t448 * t97 / 0.2e1 + t276 * t99 / 0.2e1 - t208 * t127 + pkin(3) * t198 + t131 * t416 - t139 * t143 + t284 * t145 + t129 * t419 + (m(5) * t280 + t344 - t351) * pkin(7) + (Ifges(5,5) * t254 + Ifges(6,5) * t276 + Ifges(5,6) * t256 - Ifges(6,6) * t448) * t404 + t280 * mrSges(5,3)) * qJD(3), t367 + t10 * qJD(2) + t1 * qJD(3) + (-Ifges(5,5) * t348 - Ifges(5,6) * t342 + t393 + (t402 * t60 + t403 * t59) * t437 + t285 + t286 - t392 - t153 * mrSges(5,2) - t154 * mrSges(5,1) + t336) * qJD(4) + t15 * qJD(5), t366 + t12 * qJD(2) + t4 * qJD(3) + t15 * qJD(4) + (t336 - t396 - t397) * qJD(5); qJD(3) * t9 + qJD(4) * t11 + qJD(5) * t13 - t354, t35 * qJD(3), t19 * qJD(4) + t26 * qJD(5) + (t127 + t365) * t333 + t281 + ((t183 * t448 + t185 * t276) * mrSges(6,3) + (-mrSges(4,2) + t465) * t257 + 0.2e1 * (t245 * t255 + t361 - t362) * t438 + 0.2e1 * (t334 * t398 - t401) * t439) * qJD(3), t364 + t19 * qJD(3) + (m(6) * (-t184 * t327 + t186 * t328) + t446 * t255 + t296) * qJD(4) + t463, t26 * qJD(3) + qJD(4) * t296 + t363 + t463; -qJD(2) * t9 + qJD(4) * t2 + qJD(5) * t5 - t370, -qJD(4) * t18 + qJD(5) * t27 - t281, qJD(4) * t14 + qJD(5) * t20, (-t276 * t293 + t448 * t294 + (t139 * t403 + t284 * t402) * t437 + t295 + t446 * pkin(7) + t28) * qJD(4) + t464 + t275, t28 * qJD(4) + t274 + t464; -qJD(2) * t11 - qJD(3) * t2 - qJD(5) * t16 - t367, t18 * qJD(3) - t364, qJD(5) * t29 - t275, -t210, -t210 - t269; -qJD(2) * t13 - qJD(3) * t5 + qJD(4) * t16 - t366, -t27 * qJD(3) - t363, -qJD(4) * t29 - t274, t269, 0;];
Cq = t17;
