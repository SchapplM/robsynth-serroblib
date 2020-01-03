% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR11_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR11_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR11_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:32:29
% EndTime: 2019-12-31 21:32:47
% DurationCPUTime: 8.10s
% Computational Cost: add. (9432->528), mult. (19618->713), div. (0->0), fcn. (17760->6), ass. (0->265)
t337 = sin(qJ(3));
t511 = pkin(7) - pkin(8);
t299 = t511 * t337;
t340 = cos(qJ(3));
t300 = t511 * t340;
t336 = sin(qJ(5));
t339 = cos(qJ(5));
t163 = t299 * t339 - t300 * t336;
t366 = t299 * t336 + t300 * t339;
t413 = t337 * t339;
t269 = -t336 * t340 + t413;
t365 = t336 * t337 + t339 * t340;
t405 = -Ifges(6,5) * t365 - Ifges(6,6) * t269;
t531 = -t366 * mrSges(6,1) - t163 * mrSges(6,2) + t405;
t538 = t531 * qJD(5);
t341 = cos(qJ(2));
t224 = t269 * t341;
t338 = sin(qJ(2));
t170 = mrSges(6,2) * t338 + mrSges(6,3) * t224;
t461 = pkin(7) * t341;
t301 = pkin(2) * t338 - t461;
t411 = t338 * t340;
t198 = -pkin(6) * t411 + t337 * t301;
t172 = t338 * qJ(4) + t198;
t393 = -pkin(6) * t337 - pkin(3);
t417 = t301 * t340;
t173 = t338 * t393 - t417;
t414 = t337 * t338;
t197 = pkin(6) * t414 + t417;
t408 = t340 * t341;
t311 = mrSges(5,2) * t408;
t433 = t338 * mrSges(5,1);
t276 = t311 - t433;
t412 = t337 * t341;
t277 = -mrSges(5,2) * t412 + mrSges(5,3) * t338;
t486 = pkin(3) + pkin(4);
t279 = -t336 * qJ(4) - t339 * t486;
t280 = t339 * qJ(4) - t336 * t486;
t225 = t365 * t341;
t474 = t225 / 0.2e1;
t475 = t224 / 0.2e1;
t505 = Ifges(6,5) * t474 + Ifges(6,6) * t475;
t115 = (-pkin(8) * t341 - t301) * t340 + (-pkin(4) + t393) * t338;
t122 = pkin(8) * t412 + t172;
t52 = t115 * t339 - t122 * t336;
t524 = -t338 / 0.2e1;
t53 = t115 * t336 + t122 * t339;
t380 = Ifges(6,3) * t524 - t53 * mrSges(6,2) / 0.2e1 + t52 * mrSges(6,1) / 0.2e1 + t505;
t471 = -t280 / 0.2e1;
t171 = -mrSges(6,1) * t338 - mrSges(6,3) * t225;
t479 = -t171 / 0.2e1;
t489 = mrSges(5,1) / 0.2e1;
t491 = -m(6) / 0.2e1;
t493 = -m(5) / 0.2e1;
t537 = (-pkin(3) * t173 + qJ(4) * t172) * t493 + (t279 * t52 + t280 * t53) * t491 + pkin(3) * t276 / 0.2e1 - qJ(4) * t277 / 0.2e1 - t172 * mrSges(5,3) / 0.2e1 + t173 * t489 - t197 * mrSges(4,1) / 0.2e1 + t198 * mrSges(4,2) / 0.2e1 + t279 * t479 + t170 * t471 + t380;
t127 = mrSges(6,1) * t269 - mrSges(6,2) * t365;
t262 = Ifges(6,4) * t269;
t132 = -Ifges(6,2) * t365 + t262;
t261 = Ifges(6,4) * t365;
t133 = Ifges(6,2) * t269 + t261;
t135 = Ifges(6,1) * t269 - t261;
t136 = Ifges(6,1) * t365 + t262;
t428 = qJ(4) * t337;
t495 = -t340 * t486 - t428;
t257 = pkin(2) - t495;
t536 = -t257 * t127 + (t136 / 0.2e1 + t132 / 0.2e1) * t269 + (t135 / 0.2e1 - t133 / 0.2e1) * t365;
t377 = t340 * mrSges(5,1) + t337 * mrSges(5,3);
t535 = t377 / 0.2e1;
t510 = Ifges(4,5) + Ifges(5,4);
t128 = mrSges(6,1) * t365 + mrSges(6,2) * t269;
t534 = m(6) * t257 + t128;
t376 = t336 * mrSges(6,1) + t339 * mrSges(6,2);
t222 = t336 * t411 - t338 * t413;
t201 = Ifges(6,4) * t222;
t223 = t365 * t338;
t105 = Ifges(6,2) * t223 + t201;
t98 = Ifges(6,1) * t223 + t341 * Ifges(6,5) - t201;
t532 = t98 / 0.4e1 - t105 / 0.4e1;
t529 = t135 / 0.4e1 - t133 / 0.4e1;
t527 = t136 / 0.4e1 + t132 / 0.4e1;
t102 = mrSges(6,1) * t223 - mrSges(6,2) * t222;
t427 = qJ(4) * t340;
t265 = -t337 * t486 + t427;
t358 = -pkin(6) + t265;
t159 = t358 * t338;
t199 = Ifges(6,6) * t223;
t200 = Ifges(6,5) * t222;
t373 = t200 + t199;
t466 = t341 / 0.2e1;
t526 = t159 * t102 - t373 * t466;
t168 = -mrSges(6,2) * t341 - t222 * mrSges(6,3);
t480 = t168 / 0.2e1;
t525 = t341 * t405 / 0.4e1 + t257 * t102 / 0.2e1 + t159 * t127 / 0.2e1 + t163 * t480;
t468 = t338 / 0.2e1;
t523 = -t341 / 0.2e1;
t282 = -pkin(2) * t341 - pkin(7) * t338 - pkin(1);
t418 = t282 * t337;
t196 = pkin(6) * t408 + t418;
t470 = -t336 / 0.2e1;
t477 = t222 / 0.2e1;
t354 = (t223 * t470 + t339 * t477) * mrSges(6,3);
t403 = pkin(6) * t412 - t340 * t282;
t146 = pkin(8) * t411 - t403;
t314 = pkin(8) * t414;
t147 = t196 + t314;
t59 = -t146 * t336 + t147 * t339;
t60 = t146 * t339 + t147 * t336;
t518 = t196 * t493 + (t336 * t60 + t339 * t59) * t491 + t354;
t327 = Ifges(5,5) * t337;
t504 = Ifges(5,1) * t340 + t327;
t212 = -t341 * Ifges(5,4) + t338 * t504;
t450 = Ifges(4,4) * t337;
t294 = Ifges(4,1) * t340 - t450;
t214 = -t341 * Ifges(4,5) + t294 * t338;
t397 = Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1;
t382 = t397 * t341;
t517 = -t382 + t212 / 0.2e1 + t214 / 0.2e1;
t404 = t376 * t466;
t410 = t339 * t168;
t169 = mrSges(6,1) * t341 - t223 * mrSges(6,3);
t416 = t336 * t169;
t407 = -t416 / 0.2e1 + t410 / 0.2e1;
t516 = t404 - t407;
t515 = qJD(5) * t376;
t367 = t163 * t336 - t339 * t366;
t490 = m(6) / 0.2e1;
t513 = 0.2e1 * t490;
t512 = -Ifges(6,5) / 0.2e1;
t509 = Ifges(5,2) + Ifges(4,3);
t437 = t365 * mrSges(6,3);
t202 = Ifges(6,4) * t223;
t107 = Ifges(6,1) * t222 + t202;
t96 = -Ifges(6,2) * t222 + t341 * Ifges(6,6) + t202;
t507 = t96 + t107;
t503 = Ifges(5,6) * t337 + t510 * t340;
t330 = Ifges(4,4) * t340;
t502 = -Ifges(4,2) * t337 + t330;
t293 = Ifges(4,1) * t337 + t330;
t501 = -t197 * t337 + t198 * t340;
t500 = t172 * t340 + t173 * t337;
t499 = Ifges(5,3) * t340 - t327;
t498 = -pkin(3) * t340 - t428;
t281 = -pkin(2) + t498;
t497 = m(5) * t281 - t377;
t396 = Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t362 = -t200 / 0.2e1 - t199 / 0.2e1;
t449 = Ifges(5,5) * t340;
t288 = Ifges(5,3) * t337 + t449;
t496 = Ifges(4,6) * t524 + Ifges(5,6) * t468 + t288 * t466 + t502 * t523;
t492 = m(5) / 0.2e1;
t103 = mrSges(6,1) * t222 + mrSges(6,2) * t223;
t485 = t103 / 0.2e1;
t484 = t128 / 0.2e1;
t482 = t163 / 0.2e1;
t481 = -t366 / 0.2e1;
t478 = -t222 / 0.2e1;
t476 = t223 / 0.2e1;
t462 = pkin(6) * t340;
t333 = t341 * pkin(3);
t114 = pkin(4) * t341 - t146 + t333;
t166 = t418 + (-qJ(4) + t462) * t341;
t119 = t314 + t166;
t50 = t114 * t339 - t119 * t336;
t460 = t50 * mrSges(6,2);
t51 = t114 * t336 + t119 * t339;
t459 = t51 * mrSges(6,1);
t456 = t59 * mrSges(6,1);
t455 = t60 * mrSges(6,2);
t436 = t269 * mrSges(6,3);
t104 = -mrSges(6,1) * t224 + mrSges(6,2) * t225;
t160 = t358 * t341;
t167 = t333 + t403;
t213 = Ifges(5,4) * t338 + t341 * t504;
t215 = Ifges(4,5) * t338 + t294 * t341;
t284 = pkin(3) * t337 - t427;
t363 = pkin(6) + t284;
t216 = t363 * t338;
t217 = t363 * t341;
t285 = mrSges(5,1) * t337 - mrSges(5,3) * t340;
t250 = t285 * t338;
t251 = t285 * t341;
t286 = mrSges(4,1) * t337 + mrSges(4,2) * t340;
t252 = t341 * t286;
t271 = mrSges(4,2) * t341 - mrSges(4,3) * t414;
t272 = -mrSges(4,2) * t338 - mrSges(4,3) * t412;
t273 = -mrSges(4,1) * t341 - mrSges(4,3) * t411;
t274 = mrSges(5,1) * t341 + mrSges(5,2) * t411;
t275 = mrSges(4,1) * t338 - mrSges(4,3) * t408;
t325 = t341 * mrSges(5,3);
t400 = mrSges(5,2) * t414;
t278 = -t325 - t400;
t313 = Ifges(5,5) * t411;
t208 = -Ifges(5,6) * t341 + Ifges(5,3) * t414 + t313;
t431 = t341 * Ifges(4,6);
t210 = t338 * t502 - t431;
t388 = t208 / 0.2e1 - t210 / 0.2e1;
t97 = Ifges(6,4) * t225 + Ifges(6,2) * t224 - Ifges(6,6) * t338;
t99 = Ifges(6,1) * t225 + Ifges(6,4) * t224 - Ifges(6,5) * t338;
t3 = (Ifges(3,4) * t341 - pkin(1) * mrSges(3,2) + t517 * t340 + (-t341 * t396 + t388) * t337 + t505) * t341 + m(4) * (t196 * t198 - t197 * t403) - t403 * t275 + (pkin(6) * t252 - pkin(1) * mrSges(3,1) + t223 * t512 + Ifges(6,6) * t477 - Ifges(3,4) * t338 + (t215 / 0.2e1 + t213 / 0.2e1 + t397 * t338) * t340 + (t396 * t338 + t496) * t337 + (Ifges(3,1) - Ifges(3,2) - Ifges(6,3) + (m(4) * pkin(6) + t286) * pkin(6) - t509) * t341) * t338 + t99 * t476 + t97 * t478 + t98 * t474 + t96 * t475 + t198 * t271 + t196 * t272 + t197 * t273 + t173 * t274 + t167 * t276 + t166 * t277 + t172 * t278 + t217 * t250 + t216 * t251 + t50 * t171 + t53 * t168 + t52 * t169 + t51 * t170 + t159 * t104 + t160 * t103 + m(5) * (t166 * t172 + t167 * t173 + t216 * t217) + m(6) * (t159 * t160 + t50 * t52 + t51 * t53);
t435 = t3 * qJD(1);
t184 = t495 * t338;
t247 = t498 * t338;
t248 = t377 * t338;
t378 = t340 * mrSges(4,1) - t337 * mrSges(4,2);
t249 = t378 * t338;
t253 = t499 * t338;
t289 = Ifges(4,2) * t340 + t450;
t254 = t338 * t289;
t255 = -Ifges(5,1) * t414 + t313;
t256 = t338 * t293;
t312 = Ifges(5,6) * t411;
t4 = t312 * t523 + t216 * t248 - t247 * t250 + t98 * t477 + t105 * t478 + t184 * t103 + t60 * t168 + t59 * t169 + (-t273 + t274) * t196 - (t271 + t278) * t403 + (-t222 * t50 + t223 * t51) * mrSges(6,3) + m(5) * (-t166 * t403 + t167 * t196 - t216 * t247) + m(6) * (t159 * t184 + t50 * t59 + t51 * t60) + (pkin(6) * t249 + (t255 / 0.2e1 - t256 / 0.2e1 + t431 / 0.2e1 - t166 * mrSges(5,2) - t196 * mrSges(4,3) + t388) * t340 + (t253 / 0.2e1 + t254 / 0.2e1 - t167 * mrSges(5,2) - t403 * mrSges(4,3) - t517) * t337) * t338 + t507 * t476 - t526;
t430 = t4 * qJD(1);
t7 = t50 * t168 - t51 * t169 + (-t51 * mrSges(6,3) - t107 / 0.2e1 - t96 / 0.2e1) * t223 + (t50 * mrSges(6,3) - t98 / 0.2e1 + t105 / 0.2e1) * t222 + t526;
t429 = t7 * qJD(1);
t18 = (t103 - t250) * t411 + (-t278 - t410 + t416) * t341 + m(6) * (t159 * t411 + (t336 * t50 - t339 * t51) * t341) + m(5) * (-t166 * t341 - t216 * t411);
t426 = qJD(1) * t18;
t49 = t280 * mrSges(6,1) + mrSges(6,2) * t279;
t425 = qJD(5) * t49;
t420 = t279 * t222;
t419 = t280 * t223;
t415 = t336 * t269;
t409 = t339 * t365;
t395 = t50 / 0.2e1 + t60 / 0.2e1;
t394 = t51 / 0.2e1 - t59 / 0.2e1;
t385 = -t499 / 0.2e1 - t289 / 0.2e1;
t291 = Ifges(5,1) * t337 - t449;
t384 = t291 / 0.2e1 + t293 / 0.2e1;
t381 = (t484 + t535) * t340;
t10 = -t284 * t377 - pkin(2) * t286 + (m(5) * t284 + t285) * t281 + (-t288 / 0.2e1 + t502 / 0.2e1 + t384) * t340 + (t294 / 0.2e1 + t504 / 0.2e1 + t385) * t337 + t534 * t265 + t536;
t343 = t529 * t222 + t527 * t223 + (t163 * t478 + t269 * t394) * mrSges(6,3) + (t216 * t284 - t247 * t281) * t492 + (t159 * t265 + t184 * t257 + (-t51 + t59) * t163) * t490 - pkin(2) * t249 / 0.2e1 + t184 * t484 + t216 * t285 / 0.2e1 + t247 * t535 + t265 * t485 + t281 * t248 / 0.2e1 + t284 * t250 / 0.2e1 + t507 * t269 / 0.4e1 - t503 * t341 / 0.4e1 - (t395 * mrSges(6,3) - t532) * t365 + (t476 * mrSges(6,3) + (t50 + t60) * t490 + t169 / 0.2e1) * t366 - t525;
t347 = pkin(6) * t286 / 0.2e1 + (t294 / 0.4e1 + t504 / 0.4e1 - t289 / 0.4e1 - t499 / 0.4e1) * t340 + (-t291 / 0.4e1 - t502 / 0.4e1 + t288 / 0.4e1 - t293 / 0.4e1) * t337 + (mrSges(5,2) + mrSges(4,3)) * pkin(7) * (-t340 ^ 2 / 0.2e1 - t337 ^ 2 / 0.2e1);
t349 = (t196 / 0.2e1 - t166 / 0.2e1) * mrSges(5,2) + ((-t166 + t196) * t492 - t278 / 0.2e1 - t271 / 0.2e1) * pkin(7) + t208 / 0.4e1 - t210 / 0.4e1 + t255 / 0.4e1 - t256 / 0.4e1;
t350 = -t254 / 0.4e1 - t253 / 0.4e1 + t214 / 0.4e1 + t212 / 0.4e1 + (-t403 / 0.2e1 + t167 / 0.2e1) * mrSges(5,2) + (t274 / 0.2e1 - t273 / 0.2e1 + (t167 - t403) * t492) * pkin(7);
t2 = (-Ifges(5,2) / 0.2e1 - Ifges(4,3) / 0.2e1 + t347) * t338 + t343 + ((0.3e1 / 0.4e1 * Ifges(4,6) - Ifges(5,6) / 0.2e1) * t341 + t349) * t337 + (-t382 + t350) * t340 + t537;
t370 = t2 * qJD(1) + t10 * qJD(2);
t344 = (-t107 / 0.4e1 - t96 / 0.4e1) * t269 - t532 * t365 + (mrSges(6,3) * t481 - t527) * t223 + (mrSges(6,3) * t482 - t529) * t222 + t169 * t481 + t525;
t5 = t344 - t380;
t369 = t5 * qJD(1) - qJD(2) * t536;
t345 = (t485 - t250 / 0.2e1) * t337 + (t409 / 0.2e1 - t415 / 0.2e1) * t341 * mrSges(6,3) + (-t216 * t337 + (-t281 * t338 - t461) * t340) * t492 + (t159 * t337 + t257 * t411 + t341 * t367) * t490;
t352 = t173 * t493 + (t336 * t53 + t339 * t52) * t491 + t170 * t470 + t339 * t479;
t11 = -t311 + (t489 + t381) * t338 + t345 + t352;
t40 = (-t497 + t534) * t337;
t368 = qJD(1) * t11 + qJD(2) * t40;
t22 = -t354 + t516;
t364 = t22 * qJD(1) - qJD(3) * t376;
t359 = m(6) * t367;
t357 = t367 * t491;
t20 = (t482 - t163 / 0.2e1) * mrSges(6,2);
t351 = (t420 / 0.2e1 - t419 / 0.2e1) * mrSges(6,3) + t279 * t480 + t169 * t471 - t362;
t8 = mrSges(6,1) * t394 + mrSges(6,2) * t395 + t351 + t362;
t356 = t8 * qJD(1) + t20 * qJD(2) + t49 * qJD(3);
t348 = -t325 + (t418 + (-0.2e1 * qJ(4) + t462) * t341) * t492 + ((-t280 * t341 + t51) * t339 + (t279 * t341 - t50) * t336) * t490 - t516;
t14 = t348 + t518;
t29 = t357 + t359 / 0.2e1;
t90 = mrSges(5,3) + m(5) * qJ(4) + m(6) * (-t279 * t336 + t280 * t339) + t376;
t355 = t14 * qJD(1) + t29 * qJD(2) + t90 * qJD(3);
t25 = t357 - t359 / 0.2e1 + (m(5) * pkin(7) + mrSges(5,2)) * t340 + (-t409 + t415) * mrSges(6,3);
t23 = t354 + t404 + t407;
t13 = t348 - t400 - t518;
t12 = -t433 / 0.2e1 + t338 * t381 + t345 - t352;
t9 = t460 / 0.2e1 + t459 / 0.2e1 - t455 / 0.2e1 + t456 / 0.2e1 + t351 - t362;
t6 = t344 + t380;
t1 = t347 * t338 + t343 + t350 * t340 + (t431 / 0.4e1 + t349) * t337 + t509 * t468 + t396 * t412 + t510 * t408 / 0.2e1 - t537;
t15 = [qJD(2) * t3 + qJD(3) * t4 + qJD(4) * t18 + qJD(5) * t7, t1 * qJD(3) + t12 * qJD(4) + t6 * qJD(5) + t435 + ((m(4) * t501 + m(5) * t500) * pkin(7) + (Ifges(3,5) + (-m(4) * pkin(2) - mrSges(3,1) - t378) * pkin(6)) * t341 - t52 * t436 + (t160 * t257 + t163 * t52 + t366 * t53) * t513 + t135 * t474 + t132 * t475 - t53 * t437 + t281 * t251 - t365 * t97 / 0.2e1 + t269 * t99 / 0.2e1 - pkin(2) * t252 + t257 * t104 + t163 * t171 + t366 * t170 + t160 * t128 + (t215 + t213) * t337 / 0.2e1 + ((-t275 + t276) * pkin(7) + t385 * t341 + t510 * t468) * t337 + ((t272 + t277) * pkin(7) + t384 * t341 + (Ifges(4,6) - Ifges(5,6)) * t468 - t496) * t340 + (pkin(6) * mrSges(3,2) - Ifges(3,6) + t269 * t512 + Ifges(6,6) * t365 / 0.2e1) * t338 + t497 * t217 + t501 * mrSges(4,3) + t500 * mrSges(5,2)) * qJD(2), t1 * qJD(2) + t13 * qJD(4) + t9 * qJD(5) + t430 + (-t196 * mrSges(4,1) - t196 * mrSges(5,1) + t403 * mrSges(4,2) - t403 * mrSges(5,3) + t312 - t373 + t455 - t456 + (t279 * t59 + t280 * t60) * t513 + 0.2e1 * (-pkin(3) * t196 - qJ(4) * t403) * t492 + ((-mrSges(5,2) * qJ(4) - Ifges(4,6)) * t340 + (mrSges(5,2) * pkin(3) - t510) * t337) * t338 + (t419 - t420) * mrSges(6,3)) * qJD(3), qJD(2) * t12 + qJD(3) * t13 + qJD(5) * t23 + t426, t429 + t6 * qJD(2) + t9 * qJD(3) + t23 * qJD(4) + (-t373 - t459 - t460) * qJD(5); qJD(3) * t2 + qJD(4) * t11 + qJD(5) * t5 - t435, qJD(3) * t10 + qJD(4) * t40 - qJD(5) * t536, t25 * qJD(4) - t538 + t370 + (t280 * t436 - t279 * t437 + m(6) * (-t163 * t280 + t279 * t366) - Ifges(4,6) * t337 + (m(5) * t498 - t377 - t378) * pkin(7) + t498 * mrSges(5,2) + t503 + t531) * qJD(3), qJD(3) * t25 + t368, -qJD(3) * t531 + t369 + t538; -qJD(2) * t2 + qJD(4) * t14 + qJD(5) * t8 - t430, qJD(4) * t29 + qJD(5) * t20 - t370, qJD(4) * t90 + t425, t355, t356 - t425; -qJD(2) * t11 - qJD(3) * t14 - qJD(5) * t22 - t426, -qJD(3) * t29 - t368, -t355 + t515, 0, -t364 - t515; -qJD(2) * t5 - qJD(3) * t8 + qJD(4) * t22 - t429, -qJD(3) * t20 - t369, -qJD(4) * t376 - t356, t364, 0;];
Cq = t15;
