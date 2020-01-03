% Calculate vector of inverse dynamics joint torques for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:10
% EndTime: 2019-12-31 22:09:07
% DurationCPUTime: 30.99s
% Computational Cost: add. (9453->749), mult. (23456->1017), div. (0->0), fcn. (18076->10), ass. (0->331)
t525 = Ifges(5,4) + Ifges(6,4);
t526 = Ifges(5,1) + Ifges(6,1);
t524 = -Ifges(6,5) - Ifges(5,5);
t523 = Ifges(5,2) + Ifges(6,2);
t522 = Ifges(5,6) + Ifges(6,6);
t270 = cos(pkin(5));
t385 = qJD(1) * t270;
t256 = qJD(2) + t385;
t273 = sin(qJ(3));
t277 = cos(qJ(3));
t274 = sin(qJ(2));
t269 = sin(pkin(5));
t386 = qJD(1) * t269;
t359 = t274 * t386;
t185 = t256 * t277 - t273 * t359;
t278 = cos(qJ(2));
t375 = qJD(1) * qJD(2);
t213 = (qJDD(1) * t274 + t278 * t375) * t269;
t373 = qJDD(1) * t270;
t255 = qJDD(2) + t373;
t108 = qJD(3) * t185 + t213 * t277 + t255 * t273;
t468 = t108 / 0.2e1;
t552 = Ifges(4,4) * t468;
t264 = t270 * t274 * pkin(1);
t323 = pkin(3) * t273 - pkin(9) * t277;
t399 = t269 * t278;
t551 = -(t264 + (pkin(7) + t323) * t399) * qJD(1) + t323 * qJD(3);
t186 = t256 * t273 + t277 * t359;
t109 = -qJD(3) * t186 - t213 * t273 + t255 * t277;
t104 = qJDD(4) - t109;
t358 = t278 * t386;
t240 = qJD(3) - t358;
t272 = sin(qJ(4));
t276 = cos(qJ(4));
t148 = t186 * t276 + t240 * t272;
t212 = (-qJDD(1) * t278 + t274 * t375) * t269;
t203 = qJDD(3) + t212;
t370 = pkin(1) * t385;
t374 = qJDD(1) * t269;
t539 = pkin(7) * t374 + qJD(2) * t370;
t540 = -pkin(7) * t269 * t375 + pkin(1) * t373;
t149 = t274 * t540 + t278 * t539;
t127 = pkin(8) * t255 + t149;
t135 = -pkin(1) * t374 + pkin(2) * t212 - pkin(8) * t213;
t388 = pkin(7) * t399 + t264;
t209 = t388 * qJD(1);
t172 = t256 * pkin(8) + t209;
t178 = (-pkin(2) * t278 - pkin(8) * t274 - pkin(1)) * t386;
t381 = qJD(3) * t277;
t382 = qJD(3) * t273;
t34 = t277 * t127 + t273 * t135 - t172 * t382 + t178 * t381;
t26 = pkin(9) * t203 + t34;
t150 = -t274 * t539 + t278 * t540;
t128 = -t255 * pkin(2) - t150;
t32 = -t109 * pkin(3) - t108 * pkin(9) + t128;
t206 = -pkin(7) * t359 + t278 * t370;
t171 = -t256 * pkin(2) - t206;
t85 = -t185 * pkin(3) - t186 * pkin(9) + t171;
t106 = t172 * t277 + t178 * t273;
t88 = pkin(9) * t240 + t106;
t37 = t272 * t85 + t276 * t88;
t4 = -qJD(4) * t37 - t26 * t272 + t276 * t32;
t147 = -t186 * t272 + t240 * t276;
t49 = qJD(4) * t147 + t108 * t276 + t203 * t272;
t1 = pkin(4) * t104 - qJ(5) * t49 - qJD(5) * t148 + t4;
t452 = t203 / 0.2e1;
t467 = t109 / 0.2e1;
t469 = t104 / 0.2e1;
t50 = -qJD(4) * t148 - t108 * t272 + t203 * t276;
t473 = t50 / 0.2e1;
t474 = t49 / 0.2e1;
t521 = -Ifges(6,3) - Ifges(5,3);
t378 = qJD(4) * t276;
t380 = qJD(4) * t272;
t3 = t276 * t26 + t272 * t32 + t85 * t378 - t380 * t88;
t2 = qJ(5) * t50 + qJD(5) * t147 + t3;
t532 = t4 * mrSges(5,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2);
t550 = t1 * mrSges(6,1) - 0.2e1 * Ifges(4,2) * t467 - 0.2e1 * Ifges(4,6) * t452 - t521 * t469 + t522 * t473 - t524 * t474 + t532 - t552;
t549 = t525 * t147;
t392 = t277 * t278;
t187 = (-t272 * t392 + t274 * t276) * t269;
t176 = qJD(1) * t187;
t355 = t272 * t381;
t503 = t273 * t378 + t176 + t355;
t408 = t185 * t272;
t548 = t380 - t408;
t547 = t525 * t148;
t180 = qJD(4) - t185;
t517 = t147 * t522 - t148 * t524 - t180 * t521;
t546 = t517 / 0.2e1;
t520 = -t104 * t524 + t49 * t526 + t50 * t525;
t545 = t520 / 0.2e1;
t529 = t104 * t522 + t49 * t525 + t50 * t523;
t544 = t529 / 0.2e1;
t516 = t147 * t523 + t180 * t522 + t547;
t515 = t148 * t526 - t524 * t180 + t549;
t293 = (pkin(2) * t274 - pkin(8) * t278) * t269;
t207 = qJD(1) * t293;
t139 = t277 * t206 + t273 * t207;
t115 = pkin(9) * t359 + t139;
t368 = pkin(8) * t382;
t543 = t551 * t276 + (t115 + t368) * t272;
t324 = pkin(3) * t277 + pkin(9) * t273;
t236 = -pkin(2) - t324;
t542 = -t276 * t115 + t236 * t378 + t272 * t551;
t267 = pkin(4) * t276 + pkin(3);
t271 = -qJ(5) - pkin(9);
t319 = mrSges(4,1) * t277 - mrSges(4,2) * t273;
t541 = -m(6) * (t267 * t277 - t271 * t273) - t273 * mrSges(6,3) - t319;
t138 = -t273 * t206 + t207 * t277;
t114 = -pkin(3) * t359 - t138;
t367 = pkin(8) * t381;
t538 = t367 - t114;
t537 = t525 * t276;
t536 = t525 * t272;
t447 = cos(qJ(1));
t361 = t447 * t274;
t275 = sin(qJ(1));
t394 = t275 * t278;
t223 = t270 * t361 + t394;
t362 = t269 * t447;
t161 = t223 * t277 - t273 * t362;
t360 = t447 * t278;
t395 = t274 * t275;
t222 = -t270 * t360 + t395;
t535 = t161 * t272 - t222 * t276;
t534 = -t161 * t276 - t222 * t272;
t533 = Ifges(4,1) * t468 + Ifges(4,5) * t452;
t475 = Ifges(4,4) * t467 + t533;
t476 = m(6) * pkin(4);
t530 = -t104 * t521 - t49 * t524 + t50 * t522;
t436 = -mrSges(5,1) - mrSges(6,1);
t528 = mrSges(5,2) + mrSges(6,2);
t527 = -mrSges(4,3) + mrSges(3,2);
t188 = (t272 * t274 + t276 * t392) * t269;
t177 = qJD(1) * t188;
t393 = t276 * t277;
t265 = pkin(8) * t393;
t332 = t273 * t358;
t377 = qJD(5) * t276;
t519 = -pkin(4) * t332 + t177 * qJ(5) - t273 * t377 + (pkin(4) * t273 - qJ(5) * t393) * qJD(3) + (-t265 + (qJ(5) * t273 - t236) * t272) * qJD(4) + t543;
t396 = t273 * t276;
t518 = -qJ(5) * t176 + (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t396 + (-qJD(5) * t273 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t277) * t272 + t542;
t514 = t240 * Ifges(4,5);
t513 = t240 * Ifges(4,6);
t442 = pkin(4) * t272;
t512 = -m(6) * t442 - mrSges(4,3);
t337 = qJD(4) * t271;
t407 = t185 * t276;
t105 = -t273 * t172 + t178 * t277;
t132 = pkin(3) * t186 - pkin(9) * t185;
t54 = -t105 * t272 + t276 * t132;
t511 = -pkin(4) * t186 + qJ(5) * t407 - qJD(5) * t272 + t276 * t337 - t54;
t510 = pkin(4) * t548 - t106;
t55 = t276 * t105 + t272 * t132;
t509 = qJ(5) * t408 + t272 * t337 + t377 - t55;
t508 = t476 + mrSges(6,1);
t197 = t272 * t236 + t265;
t507 = -qJD(4) * t197 + t543;
t506 = (-t276 * t382 - t277 * t380) * pkin(8) + t542;
t505 = pkin(4) * t503 + t538;
t402 = t269 * t274;
t257 = pkin(7) * t402;
t445 = pkin(1) * t278;
t198 = t257 + (-pkin(2) - t445) * t270;
t220 = -t270 * t277 + t273 * t402;
t400 = t269 * t277;
t221 = t270 * t273 + t274 * t400;
t111 = t220 * pkin(3) - t221 * pkin(9) + t198;
t199 = pkin(8) * t270 + t388;
t389 = pkin(2) * t399 + pkin(8) * t402;
t446 = pkin(1) * t269;
t200 = -t389 - t446;
t134 = t277 * t199 + t273 * t200;
t113 = -pkin(9) * t399 + t134;
t53 = t272 * t111 + t276 * t113;
t336 = mrSges(3,3) * t359;
t504 = -mrSges(3,1) * t256 - mrSges(4,1) * t185 + mrSges(4,2) * t186 + t336;
t379 = qJD(4) * t273;
t502 = t272 * t379 - t276 * t381 + t177;
t229 = t270 * t445 - t257;
t314 = -mrSges(6,1) * t276 + mrSges(6,2) * t272;
t317 = -mrSges(5,1) * t276 + mrSges(5,2) * t272;
t501 = m(5) * pkin(3) + m(6) * t267 - t314 - t317;
t313 = mrSges(6,1) * t272 + mrSges(6,2) * t276;
t316 = mrSges(5,1) * t272 + mrSges(5,2) * t276;
t87 = -pkin(3) * t240 - t105;
t56 = -pkin(4) * t147 + qJD(5) + t87;
t500 = t56 * t313 + t87 * t316;
t499 = -t272 * t524 + t276 * t522;
t498 = -t272 * t522 - t276 * t524;
t497 = t276 * t523 + t536;
t496 = -t272 * t523 + t537;
t495 = t272 * t526 + t537;
t494 = t276 * t526 - t536;
t493 = t332 - t382;
t492 = -m(5) * pkin(9) + m(6) * t271 - mrSges(5,3) - mrSges(6,3);
t491 = -t378 + t407;
t35 = -t273 * t127 + t135 * t277 - t172 * t381 - t178 * t382;
t489 = -t273 * t35 + t277 * t34;
t488 = -t272 * t4 + t276 * t3;
t376 = m(5) + m(4) + m(6);
t430 = Ifges(3,4) * t274;
t487 = -t274 * (Ifges(3,1) * t278 - t430) / 0.2e1 + pkin(1) * (mrSges(3,1) * t274 + mrSges(3,2) * t278);
t486 = -mrSges(5,1) - t508;
t485 = mrSges(4,1) + t501;
t284 = mrSges(4,2) + t492;
t482 = m(5) * t324 + t273 * mrSges(5,3);
t481 = pkin(8) * t376 - mrSges(3,2) - t512;
t480 = t482 + mrSges(3,1) - t541;
t479 = t35 * mrSges(4,1) - t34 * mrSges(4,2) + Ifges(4,5) * t108 + Ifges(4,6) * t109 + Ifges(4,3) * t203;
t478 = t150 * mrSges(3,1) - t149 * mrSges(3,2) + Ifges(3,5) * t213 - Ifges(3,6) * t212 + Ifges(3,3) * t255;
t429 = Ifges(4,4) * t186;
t98 = t185 * Ifges(4,2) + t429 + t513;
t470 = -t98 / 0.2e1;
t466 = -t147 / 0.2e1;
t465 = t147 / 0.2e1;
t464 = -t148 / 0.2e1;
t463 = t148 / 0.2e1;
t457 = -t180 / 0.2e1;
t456 = t180 / 0.2e1;
t455 = -t185 / 0.2e1;
t454 = -t186 / 0.2e1;
t453 = t186 / 0.2e1;
t158 = -t221 * t272 - t276 * t399;
t443 = pkin(4) * t158;
t434 = mrSges(5,3) * t147;
t433 = mrSges(5,3) * t148;
t432 = mrSges(6,3) * t147;
t431 = mrSges(6,3) * t148;
t428 = Ifges(4,4) * t273;
t427 = Ifges(4,4) * t277;
t422 = Ifges(3,6) * t256;
t421 = t105 * mrSges(4,3);
t420 = t106 * mrSges(4,3);
t417 = t185 * Ifges(4,6);
t416 = t186 * Ifges(4,5);
t415 = t240 * Ifges(4,3);
t414 = t256 * Ifges(3,5);
t27 = -pkin(3) * t203 - t35;
t413 = t27 * t273;
t404 = t223 * t272;
t225 = -t270 * t395 + t360;
t403 = t225 * t272;
t401 = t269 * t275;
t398 = t272 * t273;
t397 = t272 * t277;
t387 = t447 * pkin(1) + pkin(7) * t401;
t383 = qJD(2) * t269;
t363 = t225 * pkin(2) + t387;
t357 = t274 * t383;
t356 = t278 * t383;
t16 = -t50 * mrSges(6,1) + t49 * mrSges(6,2);
t347 = t381 / 0.2e1;
t342 = -t379 / 0.2e1;
t340 = -pkin(1) * t275 + pkin(7) * t362;
t36 = -t272 * t88 + t276 * t85;
t52 = t276 * t111 - t113 * t272;
t133 = -t273 * t199 + t200 * t277;
t335 = mrSges(3,3) * t358;
t326 = -t223 * pkin(2) + t340;
t112 = pkin(3) * t399 - t133;
t320 = mrSges(4,1) * t220 + mrSges(4,2) * t221;
t292 = -t221 * t276 + t272 * t399;
t318 = mrSges(5,1) * t158 + mrSges(5,2) * t292;
t315 = -t158 * mrSges(6,1) - mrSges(6,2) * t292;
t312 = Ifges(4,1) * t277 - t428;
t307 = -Ifges(4,2) * t273 + t427;
t302 = Ifges(4,5) * t277 - Ifges(4,6) * t273;
t165 = t225 * t277 + t273 * t401;
t224 = t270 * t394 + t361;
t120 = -t165 * t272 + t224 * t276;
t208 = qJD(2) * t293;
t210 = t229 * qJD(2);
t69 = -t199 * t381 - t200 * t382 + t208 * t277 - t273 * t210;
t19 = -qJ(5) * t148 + t36;
t68 = -t199 * t382 + t200 * t381 + t273 * t208 + t277 * t210;
t64 = pkin(9) * t357 + t68;
t156 = qJD(3) * t221 + t273 * t356;
t157 = -qJD(3) * t220 + t277 * t356;
t211 = t388 * qJD(2);
t78 = t156 * pkin(3) - t157 * pkin(9) + t211;
t14 = t111 * t378 - t113 * t380 + t272 * t78 + t276 * t64;
t160 = t223 * t273 + t277 * t362;
t65 = -pkin(3) * t357 - t69;
t15 = -qJD(4) * t53 - t272 * t64 + t276 * t78;
t252 = Ifges(3,4) * t358;
t242 = t271 * t276;
t241 = t271 * t272;
t235 = (pkin(8) + t442) * t273;
t233 = t276 * t236;
t226 = (-mrSges(3,1) * t278 + mrSges(3,2) * t274) * t269;
t205 = -t256 * mrSges(3,2) + t335;
t196 = -pkin(8) * t397 + t233;
t179 = Ifges(4,4) * t185;
t170 = Ifges(3,1) * t359 + t252 + t414;
t169 = t422 + (t278 * Ifges(3,2) + t430) * t386;
t168 = -qJ(5) * t398 + t197;
t164 = t225 * t273 - t275 * t400;
t153 = -qJ(5) * t396 + t233 + (-pkin(8) * t272 - pkin(4)) * t277;
t152 = mrSges(4,1) * t240 - mrSges(4,3) * t186;
t151 = -mrSges(4,2) * t240 + mrSges(4,3) * t185;
t121 = t165 * t276 + t224 * t272;
t99 = t186 * Ifges(4,1) + t179 + t514;
t97 = t415 + t416 + t417;
t92 = mrSges(5,1) * t180 - t433;
t91 = mrSges(6,1) * t180 - t431;
t90 = -mrSges(5,2) * t180 + t434;
t89 = -mrSges(6,2) * t180 + t432;
t84 = qJD(4) * t158 + t157 * t276 + t272 * t357;
t83 = qJD(4) * t292 - t157 * t272 + t276 * t357;
t77 = -mrSges(4,2) * t203 + mrSges(4,3) * t109;
t76 = mrSges(4,1) * t203 - mrSges(4,3) * t108;
t75 = -mrSges(5,1) * t147 + mrSges(5,2) * t148;
t74 = -mrSges(6,1) * t147 + mrSges(6,2) * t148;
t70 = t112 - t443;
t51 = -mrSges(4,1) * t109 + mrSges(4,2) * t108;
t38 = qJ(5) * t158 + t53;
t31 = pkin(4) * t220 + qJ(5) * t292 + t52;
t28 = -pkin(4) * t83 + t65;
t24 = -mrSges(5,2) * t104 + mrSges(5,3) * t50;
t23 = -mrSges(6,2) * t104 + mrSges(6,3) * t50;
t22 = mrSges(5,1) * t104 - mrSges(5,3) * t49;
t21 = mrSges(6,1) * t104 - mrSges(6,3) * t49;
t20 = qJ(5) * t147 + t37;
t18 = pkin(4) * t180 + t19;
t17 = -mrSges(5,1) * t50 + mrSges(5,2) * t49;
t13 = -pkin(4) * t50 + qJDD(5) + t27;
t6 = qJ(5) * t83 + qJD(5) * t158 + t14;
t5 = pkin(4) * t156 - qJ(5) * t84 + qJD(5) * t292 + t15;
t7 = [(-m(3) * t340 + t223 * mrSges(3,1) + t275 * mrSges(2,1) + t447 * mrSges(2,2) - m(5) * (-pkin(3) * t161 + t326) - m(4) * t326 + t161 * mrSges(4,1) - m(6) * (-t161 * t267 + t326) + t436 * t534 - t528 * t535 + t481 * t222 - t284 * t160) * g(1) + (-m(3) * t387 - t225 * mrSges(3,1) - t447 * mrSges(2,1) + t275 * mrSges(2,2) - m(5) * (pkin(3) * t165 + t363) - m(4) * t363 - t165 * mrSges(4,1) - m(6) * (t165 * t267 + t363) + t436 * t121 - t528 * t120 - t481 * t224 + t284 * t164) * g(2) + (-t552 + t530 / 0.2e1 - t34 * mrSges(4,3) + t550) * t220 + t388 * (-mrSges(3,2) * t255 - mrSges(3,3) * t212) + m(3) * (t149 * t388 + t150 * t229 - t206 * t211 + t209 * t210) + t156 * t546 + (-t156 * t521 + t522 * t83 - t524 * t84) * t456 + (t156 * t522 + t523 * t83 + t525 * t84) * t465 + (-t156 * t524 + t525 * t83 + t526 * t84) * t463 + t478 * t270 + (-t105 * t157 - t106 * t156) * mrSges(4,3) + (-t35 * mrSges(4,3) + 0.2e1 * t475) * t221 + ((-g(1) * t447 - g(2) * t275) * mrSges(3,3) + (-mrSges(3,1) * t212 - mrSges(3,2) * t213 + (m(3) * t446 - t226) * qJDD(1)) * pkin(1) + (-mrSges(3,3) * t150 + Ifges(3,1) * t213 - Ifges(3,4) * t212 + Ifges(3,5) * t255) * t274 + (t149 * mrSges(3,3) + Ifges(3,4) * t213 - Ifges(3,2) * t212 + Ifges(3,6) * t255 - t479) * t278 + ((t414 / 0.2e1 + t170 / 0.2e1 - t206 * mrSges(3,3)) * t278 + (-t422 / 0.2e1 + t97 / 0.2e1 - t169 / 0.2e1 - t209 * mrSges(3,3) - t106 * mrSges(4,2) + t105 * mrSges(4,1) + t417 / 0.2e1 + t416 / 0.2e1 + t415 / 0.2e1) * t274 + (t278 * (Ifges(3,4) * t278 - Ifges(3,2) * t274) / 0.2e1 - t487) * t386) * qJD(2)) * t269 + t504 * t211 + t185 * (Ifges(4,4) * t157 - Ifges(4,2) * t156) / 0.2e1 + t128 * t320 + t13 * t315 - t27 * t318 + m(6) * (t1 * t31 + t13 * t70 + t18 * t5 + t2 * t38 + t20 * t6 + t28 * t56) + m(5) * (t112 * t27 + t14 * t37 + t15 * t36 + t3 * t53 + t4 * t52 + t65 * t87) + m(4) * (t105 * t69 + t106 * t68 + t128 * t198 + t133 * t35 + t134 * t34 + t171 * t211) + t515 * t84 / 0.2e1 + t516 * t83 / 0.2e1 + t156 * t470 + (Ifges(4,1) * t157 - Ifges(4,4) * t156) * t453 + t240 * (Ifges(4,5) * t157 - Ifges(4,6) * t156) / 0.2e1 + Ifges(2,3) * qJDD(1) + (t1 * mrSges(6,3) + t4 * mrSges(5,3) + t524 * t469 - t525 * t473 - t526 * t474 - t520 / 0.2e1) * t292 + t31 * t21 + t38 * t23 + t52 * t22 + t53 * t24 + t70 * t16 + t28 * t74 + t65 * t75 + t56 * (-mrSges(6,1) * t83 + mrSges(6,2) * t84) + t87 * (-mrSges(5,1) * t83 + mrSges(5,2) * t84) + t6 * t89 + t14 * t90 + t5 * t91 + t15 * t92 + (t3 * mrSges(5,3) + t2 * mrSges(6,3) + t522 * t469 + t523 * t473 + t525 * t474 + t544) * t158 + t112 * t17 + t133 * t76 + t134 * t77 + t68 * t151 + t69 * t152 + t20 * (-mrSges(6,2) * t156 + mrSges(6,3) * t83) + t37 * (-mrSges(5,2) * t156 + mrSges(5,3) * t83) + t18 * (mrSges(6,1) * t156 - mrSges(6,3) * t84) + t36 * (mrSges(5,1) * t156 - mrSges(5,3) * t84) + t157 * t99 / 0.2e1 + t171 * (mrSges(4,1) * t156 + mrSges(4,2) * t157) + t198 * t51 + t210 * t205 + t229 * (mrSges(3,1) * t255 - mrSges(3,3) * t213); -((-Ifges(3,2) * t359 + t273 * t517 + t277 * t99 + t170 + t252) * t278 + t186 * (Ifges(4,5) * t274 + t278 * t312) + t185 * (Ifges(4,6) * t274 + t278 * t307) + t256 * (Ifges(3,5) * t278 - Ifges(3,6) * t274) + t274 * t97 + t240 * (Ifges(4,3) * t274 + t278 * t302)) * t386 / 0.2e1 + (t335 - t205) * t206 + t505 * t74 + t506 * t90 + (-t114 * t87 + t196 * t4 + t197 * t3 + (t381 * t87 + t413) * pkin(8) + t506 * t37 + t507 * t36) * m(5) + t507 * t92 + (-t367 - t138) * t152 + t487 * qJD(1) ^ 2 * t269 ^ 2 + t240 * t171 * (mrSges(4,1) * t273 + mrSges(4,2) * t277) + (t185 * t307 + t186 * t312 + t240 * t302) * qJD(3) / 0.2e1 + t396 * t545 + (-t495 * t379 + (-t273 * t524 + t277 * t494) * qJD(3)) * t463 + (t176 * t522 - t177 * t524 - t332 * t521) * t457 + (t176 * t523 + t177 * t525 + t332 * t522) * t466 + (t176 * t525 + t177 * t526 - t332 * t524) * t464 + t478 + (-t499 * t379 + (-t273 * t521 + t277 * t498) * qJD(3)) * t456 + (-t497 * t379 + (t273 * t522 + t277 * t496) * qJD(3)) * t465 + (-t3 * t398 - t4 * t396) * mrSges(5,3) + (-t1 * t396 - t2 * t398) * mrSges(6,3) + t169 * t359 / 0.2e1 + t316 * t413 + (-t76 + t17) * pkin(8) * t273 + (-m(4) * t171 + t336 - t504) * t209 + (-t105 * t138 - t106 * t139 - pkin(2) * t128 + ((-t105 * t277 - t106 * t273) * qJD(3) + t489) * pkin(8)) * m(4) + t489 * mrSges(4,3) + t99 * t347 + t427 * t468 - t128 * t319 + (-t105 * (mrSges(4,1) * t274 - mrSges(4,3) * t392) - t106 * (-mrSges(4,3) * t273 * t278 - mrSges(4,2) * t274)) * t386 + (pkin(8) * t77 - t550) * t277 + t518 * t89 + t519 * t91 + (t1 * t153 + t13 * t235 + t168 * t2 + t18 * t519 + t20 * t518 + t505 * t56) * m(6) + (-t368 - t139) * t151 - t381 * t421 + t98 * t332 / 0.2e1 + (t13 * t313 + t498 * t469 + t496 * t473 + t494 * t474 + t475 + t533) * t273 - t529 * t398 / 0.2e1 - t530 * t277 / 0.2e1 + (-t404 * t476 - t376 * (-t222 * pkin(2) + pkin(8) * t223) + t527 * t223 + t436 * (-t222 * t393 + t404) - t528 * (t222 * t397 + t223 * t276) + t480 * t222) * g(2) + (-t403 * t476 - t376 * (-t224 * pkin(2) + pkin(8) * t225) + t527 * t225 + t436 * (-t224 * t393 + t403) - t528 * (t224 * t397 + t225 * t276) + t480 * t224) * g(1) + t538 * t75 + (-t482 * t399 + t226 + (t512 * t274 + t278 * t541) * t269 + t436 * t188 - t528 * t187 - t376 * t389) * g(3) - pkin(2) * t51 + t428 * t467 + t515 * (t272 * t342 + t276 * t347 - t177 / 0.2e1) + t516 * (t276 * t342 - t355 / 0.2e1 - t176 / 0.2e1) + (t470 - t420 + t546) * t382 + t153 * t21 + t168 * t23 + t196 * t22 + t197 * t24 + t235 * t16 + (-mrSges(6,1) * t493 + mrSges(6,3) * t502) * t18 + (-mrSges(5,1) * t493 + mrSges(5,3) * t502) * t36 + (mrSges(5,1) * t503 - mrSges(5,2) * t502) * t87 + (mrSges(6,1) * t503 - mrSges(6,2) * t502) * t56 + (mrSges(6,2) * t493 - mrSges(6,3) * t503) * t20 + (mrSges(5,2) * t493 - mrSges(5,3) * t503) * t37; (t36 * t491 - t37 * t548 + t488) * mrSges(5,3) + (-t1 * t272 + t18 * t491 + t2 * t276 - t20 * t548) * mrSges(6,3) + t495 * t474 + t276 * t544 + t272 * t545 + (-t75 + t152) * t106 + (-pkin(3) * t27 - t106 * t87 - t36 * t54 - t37 * t55) * m(5) + (t420 - t18 * mrSges(6,1) - t36 * mrSges(5,1) + t20 * mrSges(6,2) + t37 * mrSges(5,2) - Ifges(4,2) * t455 - t171 * mrSges(4,1) + t513 / 0.2e1 + t522 * t466 - t524 * t464 - t521 * t457) * t186 + (-t380 / 0.2e1 + t408 / 0.2e1) * t516 + (t378 / 0.2e1 - t407 / 0.2e1) * t515 + t479 + t509 * t89 + t510 * t74 + t511 * t91 + (t1 * t241 - t13 * t267 + t18 * t511 - t2 * t242 + t20 * t509 + t510 * t56) * m(6) + (t421 + Ifges(4,1) * t454 - t171 * mrSges(4,2) - t514 / 0.2e1 + t496 * t466 + t494 * t464 + t498 * t457 - t500) * t185 + (t164 * t485 + t165 * t284) * g(1) + (t160 * t485 + t161 * t284) * g(2) + (-t90 * t380 - t92 * t378 + m(5) * ((-t272 * t37 - t276 * t36) * qJD(4) + t488) - t272 * t22 + t276 * t24) * pkin(9) + t13 * t314 + t27 * t317 + (t99 + t179) * t455 + (-t429 + t517) * t454 + t98 * t453 - pkin(3) * t17 + (t147 * t496 + t148 * t494 + t180 * t498) * qJD(4) / 0.2e1 - t55 * t90 - t54 * t92 - t105 * t151 + t241 * t21 - t242 * t23 - t267 * t16 + t497 * t473 + t499 * t469 + t500 * qJD(4) + (t220 * t501 + t221 * t492 + t320) * g(3); (-t486 * t535 - t528 * t534) * g(2) + t516 * t463 + (t433 + t92) * t37 + (t434 - t90) * t36 + t532 + (-m(6) * (-t18 + t19) + t431 + t91) * t20 + t18 * t432 + t508 * t1 + (-m(6) * t443 + t315 - t318) * g(3) - t19 * t89 - t56 * (mrSges(6,1) * t148 + mrSges(6,2) * t147) - t87 * (mrSges(5,1) * t148 + mrSges(5,2) * t147) + (t147 * t526 - t547) * t464 + (-t147 * t524 - t148 * t522) * t457 + (t120 * t486 + t121 * t528) * g(1) + (-t148 * t523 + t515 + t549) * t466 + t530 + ((-m(6) * t56 - t74) * t148 + t21) * pkin(4); -t147 * t89 + t148 * t91 + (-g(1) * t164 - g(2) * t160 - g(3) * t220 - t20 * t147 + t18 * t148 + t13) * m(6) + t16;];
tau = t7;
