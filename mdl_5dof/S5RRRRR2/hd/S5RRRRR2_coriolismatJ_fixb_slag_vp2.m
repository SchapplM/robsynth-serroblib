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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:52:59
% EndTime: 2019-12-05 18:53:11
% DurationCPUTime: 4.99s
% Computational Cost: add. (6658->360), mult. (18360->519), div. (0->0), fcn. (17975->8), ass. (0->247)
t458 = Ifges(5,1) - Ifges(5,2);
t465 = Ifges(6,3) - t458;
t225 = sin(qJ(3));
t228 = cos(qJ(3));
t407 = sin(qJ(4));
t408 = cos(qJ(4));
t190 = t407 * t225 - t408 * t228;
t259 = t408 * t225 + t407 * t228;
t143 = mrSges(5,1) * t190 + mrSges(5,2) * t259;
t396 = Ifges(4,4) * t225;
t404 = pkin(2) * t225;
t462 = t225 / 0.2e1;
t464 = t143 * t404 - t225 * (Ifges(4,2) * t228 + t396) / 0.2e1 + (Ifges(4,1) * t228 - t396) * t462 + (0.2e1 * Ifges(4,4) * t228 + (Ifges(4,1) - Ifges(4,2)) * t225) * t228 / 0.2e1;
t463 = t190 / 0.2e1;
t226 = sin(qJ(2));
t406 = pkin(1) * t226;
t167 = t190 * t406;
t403 = pkin(2) * t228;
t229 = cos(qJ(2));
t405 = pkin(1) * t229;
t211 = -t403 - t405;
t224 = sin(qJ(5));
t227 = cos(qJ(5));
t115 = t167 * t224 + t211 * t227;
t116 = -t167 * t227 + t211 * t224;
t252 = pkin(1) * t259;
t168 = t226 * t252;
t367 = t168 * t224;
t138 = t227 * t404 + t367;
t366 = t168 * t227;
t140 = t224 * t404 - t366;
t298 = m(5) * t211 * t404;
t426 = m(6) * pkin(2);
t347 = t426 / 0.2e1;
t230 = pkin(2) ^ 2;
t353 = t225 * t230;
t359 = t259 * t224;
t133 = -mrSges(6,2) * t190 - mrSges(6,3) * t359;
t371 = t140 * t133;
t358 = t259 * t227;
t135 = mrSges(6,1) * t190 - mrSges(6,3) * t358;
t372 = t138 * t135;
t373 = t135 * t227;
t436 = (t133 * t224 + t373) * t404;
t260 = (mrSges(4,1) * t225 + mrSges(4,2) * t228) * t405;
t443 = -t260 / 0.2e1;
t460 = t372 / 0.2e1 + t371 / 0.2e1 + ((t115 * t225 - t138 * t228) * t227 + (t116 * t225 - t140 * t228) * t224) * t347 + t443 + t298 / 0.2e1 - m(5) * t228 * t353 / 0.2e1 + t436 / 0.2e1 + t464;
t459 = -t190 / 0.2e1;
t220 = t224 ^ 2;
t222 = t227 ^ 2;
t441 = t220 + t222;
t438 = t441 * mrSges(6,3);
t384 = t227 * mrSges(6,1);
t386 = t224 * mrSges(6,2);
t192 = -t384 + t386;
t457 = mrSges(5,1) - t192;
t216 = Ifges(6,5) * t227;
t391 = Ifges(6,6) * t224;
t280 = t391 - t216;
t456 = t280 * t190;
t303 = t358 / 0.2e1;
t453 = -t359 / 0.2e1;
t452 = -t408 / 0.2e1;
t218 = Ifges(6,4) * t227;
t200 = Ifges(6,1) * t224 + t218;
t392 = Ifges(6,2) * t224;
t281 = t392 - t218;
t451 = -t281 + t200;
t361 = t190 * t224;
t132 = -mrSges(6,2) * t259 + mrSges(6,3) * t361;
t324 = t227 * t407;
t110 = pkin(2) * t132 * t324;
t328 = t259 * t407;
t292 = mrSges(5,3) * t328;
t175 = pkin(2) * t292;
t123 = t259 * (Ifges(6,5) * t224 + Ifges(6,6) * t227);
t393 = Ifges(6,4) * t224;
t196 = Ifges(6,2) * t227 + t393;
t248 = Ifges(6,6) * t259 + t281 * t190;
t201 = Ifges(6,1) * t227 - t393;
t249 = Ifges(6,5) * t259 - t201 * t190;
t360 = t190 * t227;
t308 = -t360 / 0.2e1;
t309 = t361 / 0.2e1;
t409 = t227 / 0.2e1;
t410 = t224 / 0.2e1;
t262 = t196 * t309 + t200 * t308 + t248 * t409 + t249 * t410 + t123 / 0.2e1 - Ifges(5,6) * t259 - Ifges(5,5) * t190;
t134 = mrSges(6,1) * t259 + mrSges(6,3) * t360;
t318 = t407 * t134;
t383 = t227 * mrSges(6,2);
t387 = t224 * mrSges(6,1);
t194 = t383 + t387;
t121 = t194 * t190;
t323 = t408 * t121;
t450 = Ifges(4,5) * t228 - Ifges(4,6) * t225 + t110 - t175 + t262 + (mrSges(5,3) * t408 * t190 - t224 * t318 + t323) * pkin(2);
t170 = t190 * t405;
t139 = t170 * t224 + t227 * t406;
t141 = -t170 * t227 + t224 * t406;
t169 = t229 * t252;
t254 = (-mrSges(5,1) / 0.2e1 + t192 / 0.2e1) * t169 + t170 * mrSges(5,2) / 0.2e1;
t411 = -t224 / 0.2e1;
t449 = t254 + (t139 * t411 + t141 * t409) * mrSges(6,3);
t448 = -mrSges(4,1) * t228 + mrSges(4,2) * t225;
t124 = t196 * t259;
t125 = t200 * t259;
t81 = Ifges(6,5) * t190 + t201 * t259;
t382 = t227 * t81;
t80 = Ifges(6,6) * t190 - t259 * t281;
t385 = t224 * t80;
t446 = -t227 * t124 / 0.4e1 - t224 * t125 / 0.4e1 - t456 / 0.4e1 + t382 / 0.4e1 - t385 / 0.4e1 - t451 * t359 / 0.4e1;
t445 = t123 * t459 + (t81 - t124) * t453;
t282 = mrSges(5,1) * t259 - mrSges(5,2) * t190;
t108 = t211 * t282;
t40 = t115 * t134;
t41 = t116 * t132;
t61 = t168 * t121;
t122 = t194 * t259;
t62 = t167 * t122;
t444 = -t108 / 0.2e1 - t40 / 0.2e1 - t41 / 0.2e1 + t61 / 0.2e1 + t62 / 0.2e1;
t442 = t457 * t167;
t440 = t382 / 0.2e1 - t385 / 0.2e1;
t439 = -mrSges(5,2) * t408 - t457 * t407;
t325 = t227 * t408;
t346 = t407 / 0.2e1;
t435 = t133 * t325 / 0.2e1 + t292 / 0.2e1 + t323 / 0.2e1 + t122 * t346;
t120 = t192 * t259;
t344 = pkin(2) * t407;
t290 = -t344 / 0.2e1;
t278 = t224 * t290;
t345 = pkin(2) * t408;
t291 = -t345 / 0.2e1;
t427 = pkin(2) / 0.2e1;
t433 = -t328 * t427 * t438 - t120 * t291 + t133 * t278 + t290 * t373;
t432 = t169 * t122 + t141 * t133 + t139 * t135 + (t169 * t259 + t170 * t190) * mrSges(5,3);
t389 = t167 * mrSges(5,3);
t126 = t259 * t389;
t117 = -t126 / 0.2e1;
t431 = t117 + t444;
t390 = Ifges(6,3) * t259;
t428 = t249 * t303 + t248 * t453 + (t390 + t456) * t463 - 0.2e1 * Ifges(5,4) * t190 * t459 + (t458 * t459 + t465 * t463 + (-Ifges(5,4) - t280 / 0.2e1) * t259) * t259;
t253 = t81 * t308 + t80 * t309 + t428;
t261 = t282 * t403;
t339 = -t403 / 0.2e1;
t296 = t227 * t339;
t375 = t132 * t224;
t429 = t134 * t296 + t117 + t126 / 0.2e1 - t261 / 0.2e1 + t253 + t339 * t375 - t444;
t424 = -mrSges(6,1) / 0.2e1;
t423 = -mrSges(6,2) / 0.2e1;
t422 = mrSges(6,2) / 0.2e1;
t416 = -t139 / 0.2e1;
t415 = t141 / 0.2e1;
t414 = -t168 / 0.2e1;
t412 = t194 / 0.2e1;
t398 = mrSges(6,3) * t259;
t388 = t168 * mrSges(5,2);
t239 = t253 + t464;
t271 = -t126 + t62 - t108 - t40 - t41 + t61;
t3 = m(6) * (t115 * t138 + t116 * t140 - t167 * t168) + t239 - t260 - t126 + t298 + t372 + t371 - t271;
t381 = t3 * qJD(1);
t236 = (-t222 * Ifges(6,1) / 0.2e1 + (t218 - t392 / 0.2e1) * t224 + t465) * t190 - (Ifges(5,4) + t280) * t259;
t269 = t216 / 0.2e1 - t391 / 0.2e1;
t238 = (-Ifges(5,4) + t269) * t190 + t440;
t377 = t116 * t227;
t379 = t115 * t224;
t274 = t377 - t379;
t356 = t224 * t135;
t60 = t133 * t366;
t6 = t60 + (-t356 - m(6) * (-t167 - t274)) * t168 + t238 * t190 - (t236 - t389) * t259 + t271;
t380 = t6 * qJD(1);
t378 = t116 * t135;
t301 = -t358 / 0.2e1;
t284 = -t125 * t303 + t80 * t301 + t445;
t369 = t168 * t120;
t13 = t115 * t133 - t274 * t398 + t284 - t369 - t378;
t376 = t13 * qJD(1);
t348 = t225 ^ 2 + t228 ^ 2;
t283 = t348 * mrSges(4,3) - mrSges(3,2);
t337 = -mrSges(3,1) + t143 + t448;
t368 = t168 * t169;
t16 = t337 * t406 + m(6) * (t115 * t139 + t116 * t141 + t368) + m(5) * (t167 * t170 + t211 * t406 + t368) + (t283 + m(4) * (-0.1e1 + t348) * t406) * t405 + t432;
t370 = t16 * qJD(1);
t343 = t138 * t224 * mrSges(6,3);
t342 = t140 * t227 * mrSges(6,3);
t338 = t403 / 0.2e1;
t334 = -t389 / 0.2e1;
t329 = Ifges(6,5) * t308 + Ifges(6,6) * t309 + t390 / 0.2e1;
t327 = t224 * t408;
t326 = t224 * t407;
t322 = t408 * t167;
t321 = t408 * t169;
t320 = t408 * t194;
t316 = t407 * t220;
t315 = t407 * t222;
t311 = t168 * t412;
t310 = t367 / 0.2e1;
t302 = t358 / 0.4e1;
t300 = -t358 / 0.4e1;
t299 = t110 / 0.2e1 - t175 / 0.2e1;
t295 = t132 * t338;
t294 = t134 * t338;
t63 = t196 * t411 + t201 * t410 + t451 * t409;
t286 = t135 * t452;
t233 = m(5) * (-t407 * t170 - t321) * t427 + (-t139 * t326 + t141 * t324 - t321) * t347 + t443 + t449;
t1 = t227 * t294 + t224 * t295 + t233 + t431 - t428 + t261 / 0.2e1 + t440 * t190 - t259 * t334 - t460;
t273 = t134 * t227 + t375;
t8 = t239 + t436 + ((-m(6) * t441 - m(5)) * t353 + (-t273 - t282) * pkin(2)) * t228;
t279 = -t1 * qJD(1) + t8 * qJD(2);
t14 = t273 * t403 + (-mrSges(5,2) * t403 + t238) * t190 - (-mrSges(5,1) * t403 + t236) * t259;
t4 = t60 / 0.2e1 + (mrSges(6,3) * t415 + t294) * t227 + (mrSges(6,3) * t416 + t135 * t414 + t295) * t224 + (mrSges(5,2) * t339 + t238) * t190 - (mrSges(5,1) * t339 + t236 + t334) * t259 + t254 + t431;
t277 = -t4 * qJD(1) - t14 * qJD(2);
t109 = t356 * t403;
t235 = (t296 + t115 / 0.2e1) * t133 - (-t379 / 0.2e1 + t377 / 0.2e1) * t398 + t109 / 0.2e1 - t378 / 0.2e1 - t369 / 0.2e1 + t284;
t267 = mrSges(6,1) * t416 + mrSges(6,2) * t415;
t10 = t235 + t267;
t24 = t227 * t133 * t403 - t125 * t301 + t80 * t303 - t109 - t445;
t275 = t10 * qJD(1) - t24 * qJD(2);
t268 = t138 * t424 + t140 * t422;
t251 = t196 * t300 + t201 * t302 + t446;
t23 = -t390 / 0.2e1 + t269 * t190 + t251;
t11 = t23 + t311 + t268 + t433;
t17 = t196 * t302 + t201 * t300 + (t120 * t452 + (mrSges(6,1) * t462 + t135 * t346) * t227 + (t133 * t346 + t225 * t423) * t224 - (-t315 / 0.2e1 - t316 / 0.2e1) * t398) * pkin(2) + t329 - t446;
t245 = (-t200 / 0.2e1 + t281 / 0.2e1) * t227 + (-t201 / 0.2e1 + t196 / 0.2e1) * t224;
t42 = pkin(2) * t320 + t245;
t258 = t11 * qJD(1) - t17 * qJD(2) - t42 * qJD(3);
t232 = (-t115 * t327 + t116 * t325 + t322 + (-t316 - t315 + t407) * t168) * t347 + t134 * t278 + t299 + t414 * t438 + (t224 * t286 + t435) * pkin(2);
t255 = t343 / 0.2e1 - t342 / 0.2e1;
t15 = t232 + t255;
t25 = ((t286 - t318 / 0.2e1) * t224 + t435) * pkin(2) + t299;
t35 = t345 * t438 + m(6) * (-0.1e1 + t441) * t230 * t407 * t408 + t439 * pkin(2);
t257 = t15 * qJD(1) + t25 * qJD(2) + t35 * qJD(3);
t250 = t262 + t442;
t247 = t325 * t423 + t327 * t424;
t19 = (t412 - t387 / 0.2e1 - t383 / 0.2e1) * t168 + t23;
t32 = (t320 / 0.2e1 + t247) * pkin(2) + t245;
t246 = -qJD(1) * t19 - qJD(2) * t23 + qJD(3) * t32 - qJD(4) * t63;
t22 = t251 + t329;
t237 = t22 + t433;
t33 = t247 * pkin(2) + t194 * t291 + t63;
t21 = t25 + t262;
t20 = mrSges(6,1) * t310 + t366 * t422 + t22 + t311;
t18 = t237 + (-t386 / 0.2e1 + t384 / 0.2e1) * t404;
t12 = t237 + t311 - t268;
t9 = t235 - t267;
t7 = t232 + t250 - t255 + t388;
t5 = t429 - t60 / 0.2e1 + t135 * t310 + t449;
t2 = t233 + t429 + t460;
t26 = [qJD(2) * t16 + qJD(3) * t3 - qJD(4) * t6 + qJD(5) * t13, t2 * qJD(3) + t5 * qJD(4) + t9 * qJD(5) + t370 + (m(6) * (-t139 * t227 - t141 * t224) * t403 + t432 + (t283 * t229 + (-m(5) * t403 + t337) * t226) * pkin(1)) * qJD(2), t2 * qJD(2) + t7 * qJD(4) + t12 * qJD(5) + t381 + (t342 - t343 + (-t138 * t326 + t140 * t324 + t322) * t426 + m(5) * (pkin(2) * t322 - t168 * t344) + t388 + t442 + t448 * t406 + t450) * qJD(3), t5 * qJD(2) + t7 * qJD(3) + t20 * qJD(5) - t380 + (t250 + (mrSges(5,2) - t438) * t168) * qJD(4), t376 + t9 * qJD(2) + t12 * qJD(3) + t20 * qJD(4) + (-mrSges(6,1) * t116 - mrSges(6,2) * t115 - t123) * qJD(5); -qJD(3) * t1 - qJD(4) * t4 + qJD(5) * t10 - t370, qJD(3) * t8 - qJD(4) * t14 - qJD(5) * t24, t450 * qJD(3) + t21 * qJD(4) + t18 * qJD(5) + t279, t21 * qJD(3) + qJD(4) * t262 + t22 * qJD(5) + t277, t18 * qJD(3) + t22 * qJD(4) + (t194 * t403 - t123) * qJD(5) + t275; qJD(2) * t1 + qJD(4) * t15 + qJD(5) * t11 - t381, qJD(4) * t25 - qJD(5) * t17 - t279, qJD(4) * t35 - qJD(5) * t42, t33 * qJD(5) + (t408 * t438 + t439) * qJD(4) * pkin(2) + t257, t33 * qJD(4) + ((-mrSges(6,1) * t324 + mrSges(6,2) * t326) * pkin(2) - t280) * qJD(5) + t258; qJD(2) * t4 - qJD(3) * t15 + qJD(5) * t19 + t380, -qJD(3) * t25 + qJD(5) * t23 - t277, -qJD(5) * t32 - t257, t63 * qJD(5), -qJD(5) * t280 - t246; -qJD(2) * t10 - qJD(3) * t11 - qJD(4) * t19 - t376, qJD(3) * t17 - qJD(4) * t23 - t275, qJD(4) * t32 - t258, t246, 0;];
Cq = t26;
