% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:06:15
% EndTime: 2018-11-23 18:06:36
% DurationCPUTime: 21.76s
% Computational Cost: add. (15156->734), mult. (38129->982), div. (0->0), fcn. (26759->8), ass. (0->323)
t313 = sin(qJ(2));
t316 = cos(qJ(2));
t339 = pkin(2) * t313 - pkin(8) * t316;
t270 = t339 * qJD(1);
t315 = cos(qJ(3));
t312 = sin(qJ(3));
t372 = qJD(1) * t313;
t353 = t312 * t372;
t219 = pkin(7) * t353 + t315 * t270;
t378 = t315 * t316;
t328 = pkin(3) * t313 - pkin(9) * t378;
t442 = -pkin(9) - pkin(8);
t355 = qJD(3) * t442;
t512 = -qJD(1) * t328 + t315 * t355 - t219;
t249 = t312 * t270;
t379 = t313 * t315;
t380 = t312 * t316;
t511 = t249 + (-pkin(7) * t379 - pkin(9) * t380) * qJD(1) - t312 * t355;
t311 = sin(qJ(4));
t314 = cos(qJ(4));
t329 = t311 * t312 - t314 * t315;
t458 = qJD(3) + qJD(4);
t208 = t458 * t329;
t323 = t329 * t316;
t229 = qJD(1) * t323;
t510 = t208 - t229;
t266 = t311 * t315 + t312 * t314;
t209 = t458 * t266;
t324 = t266 * t316;
t228 = qJD(1) * t324;
t507 = t209 - t228;
t282 = t442 * t312;
t283 = t442 * t315;
t364 = qJD(4) * t314;
t365 = qJD(4) * t311;
t490 = t282 * t364 + t283 * t365 + t311 * t512 - t511 * t314;
t215 = t311 * t282 - t314 * t283;
t489 = -qJD(4) * t215 + t511 * t311 + t314 * t512;
t509 = t507 * qJ(5) + qJD(5) * t329 - t490;
t508 = -pkin(4) * t372 + t510 * qJ(5) - qJD(5) * t266 + t489;
t363 = qJD(2) * qJD(3);
t367 = qJD(3) * t312;
t368 = qJD(2) * t316;
t216 = t315 * t363 + (-t313 * t367 + t315 * t368) * qJD(1);
t366 = qJD(3) * t315;
t487 = t312 * t368 + t313 * t366;
t217 = -qJD(1) * t487 - t312 * t363;
t369 = qJD(2) * t315;
t263 = -t353 + t369;
t352 = t315 * t372;
t264 = qJD(2) * t312 + t352;
t343 = t314 * t263 - t264 * t311;
t111 = qJD(4) * t343 + t216 * t314 + t217 * t311;
t201 = t263 * t311 + t264 * t314;
t112 = -qJD(4) * t201 - t216 * t311 + t217 * t314;
t310 = sin(pkin(10));
t386 = cos(pkin(10));
t63 = t111 * t310 - t112 * t386;
t447 = t63 / 0.2e1;
t64 = t111 * t386 + t310 * t112;
t446 = t64 / 0.2e1;
t477 = Ifges(6,1) + Ifges(7,1);
t476 = -Ifges(6,4) + Ifges(7,5);
t475 = Ifges(7,4) + Ifges(6,5);
t140 = -t208 * t310 + t209 * t386;
t161 = t228 * t386 - t229 * t310;
t377 = t140 - t161;
t141 = -t208 * t386 - t310 * t209;
t162 = -t310 * t228 - t229 * t386;
t376 = t141 - t162;
t371 = qJD(1) * t316;
t306 = pkin(7) * t371;
t410 = pkin(3) * t312;
t258 = t371 * t410 + t306;
t506 = pkin(3) * t367 - t258;
t132 = t201 * t310 - t343 * t386;
t484 = t201 * t386 + t310 * t343;
t505 = pkin(5) * t484 + qJ(6) * t132;
t191 = Ifges(5,4) * t343;
t297 = qJD(3) - t371;
t285 = qJD(4) + t297;
t127 = Ifges(5,1) * t201 + Ifges(5,5) * t285 + t191;
t280 = -qJD(2) * pkin(2) + pkin(7) * t372;
t223 = -pkin(3) * t263 + t280;
t399 = Ifges(5,4) * t201;
t416 = -t285 / 0.2e1;
t425 = -t201 / 0.2e1;
t427 = -t343 / 0.2e1;
t431 = -t484 / 0.2e1;
t434 = t132 / 0.2e1;
t435 = -t132 / 0.2e1;
t370 = qJD(2) * t313;
t345 = qJD(1) * t370;
t275 = -pkin(2) * t316 - t313 * pkin(8) - pkin(1);
t255 = t275 * qJD(1);
t281 = qJD(2) * pkin(8) + t306;
t205 = t255 * t312 + t281 * t315;
t273 = t339 * qJD(2);
t256 = qJD(1) * t273;
t342 = pkin(7) * t345;
t148 = -qJD(3) * t205 + t315 * t256 + t312 * t342;
t105 = pkin(3) * t345 - pkin(9) * t216 + t148;
t147 = t255 * t366 + t312 * t256 - t281 * t367 - t315 * t342;
t118 = pkin(9) * t217 + t147;
t204 = t315 * t255 - t281 * t312;
t168 = -pkin(9) * t264 + t204;
t158 = pkin(3) * t297 + t168;
t169 = pkin(9) * t263 + t205;
t165 = t314 * t169;
t94 = t158 * t311 + t165;
t30 = -qJD(4) * t94 + t314 * t105 - t118 * t311;
t13 = pkin(4) * t345 - qJ(5) * t111 - qJD(5) * t201 + t30;
t29 = t311 * t105 + t314 * t118 + t158 * t364 - t169 * t365;
t15 = qJ(5) * t112 + qJD(5) * t343 + t29;
t6 = t310 * t13 + t386 * t15;
t2 = qJ(6) * t345 + qJD(6) * t285 + t6;
t5 = t13 * t386 - t310 * t15;
t3 = -pkin(5) * t345 - t5;
t456 = t30 * mrSges(5,1) + t5 * mrSges(6,1) - t3 * mrSges(7,1) - t29 * mrSges(5,2) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t474 = -Ifges(6,6) + Ifges(7,6);
t500 = Ifges(5,3) + Ifges(6,3) + Ifges(7,2);
t457 = Ifges(5,5) * t111 + Ifges(5,6) * t112 + t345 * t500 + t474 * t63 + t475 * t64;
t153 = -pkin(4) * t343 + qJD(5) + t223;
t466 = qJ(5) * t343;
t83 = t94 + t466;
t388 = t310 * t83;
t491 = qJ(5) * t201;
t163 = t311 * t169;
t93 = t314 * t158 - t163;
t82 = t93 - t491;
t77 = pkin(4) * t285 + t82;
t25 = t386 * t77 - t388;
t23 = -t285 * pkin(5) + qJD(6) - t25;
t45 = pkin(5) * t132 - qJ(6) * t484 + t153;
t473 = t476 * t132 + t475 * t285 + t477 * t484;
t501 = t153 * mrSges(6,2) + t23 * mrSges(7,2) - t25 * mrSges(6,3) - t45 * mrSges(7,3) + t473 / 0.2e1;
t504 = (-Ifges(6,4) * t434 - Ifges(7,5) * t435 - t475 * t416 - t431 * t477 + t501) * t132 + t456 + t457 + (-Ifges(5,2) * t201 + t127 + t191) * t427 - t223 * (mrSges(5,1) * t201 + mrSges(5,2) * t343) + (Ifges(5,1) * t343 - t399) * t425;
t349 = Ifges(3,5) * qJD(2) / 0.2e1;
t503 = t477 * t446 + t476 * t447 + t475 * t345 / 0.2e1;
t494 = t310 * t509 + t508 * t386;
t493 = t508 * t310 - t386 * t509;
t488 = pkin(4) * t507 + t506;
t409 = pkin(4) * t201;
t496 = pkin(5) * t372 - t494;
t495 = -qJ(6) * t372 + t493;
t197 = t266 * t386 - t310 * t329;
t492 = pkin(5) * t377 - t376 * qJ(6) - qJD(6) * t197 + t488;
t101 = -t168 * t311 - t165;
t327 = t101 - t466;
t362 = t310 * t311 * pkin(3);
t102 = t314 * t168 - t163;
t86 = t102 - t491;
t467 = -qJD(4) * t362 - t310 * t327 + (pkin(3) * t364 - t86) * t386;
t304 = Ifges(3,4) * t371;
t391 = t264 * Ifges(4,4);
t179 = t263 * Ifges(4,2) + t297 * Ifges(4,6) + t391;
t259 = Ifges(4,4) * t263;
t180 = t264 * Ifges(4,1) + t297 * Ifges(4,5) + t259;
t330 = t204 * t315 + t205 * t312;
t400 = Ifges(4,4) * t315;
t334 = -Ifges(4,2) * t312 + t400;
t401 = Ifges(4,4) * t312;
t336 = Ifges(4,1) * t315 - t401;
t337 = mrSges(4,1) * t312 + mrSges(4,2) * t315;
t397 = Ifges(4,6) * t312;
t398 = Ifges(4,5) * t315;
t412 = t315 / 0.2e1;
t413 = -t312 / 0.2e1;
t417 = t264 / 0.2e1;
t317 = -t330 * mrSges(4,3) + t280 * t337 + t263 * t334 / 0.2e1 + t336 * t417 + t297 * (-t397 + t398) / 0.2e1 + t179 * t413 + t180 * t412;
t486 = t317 + Ifges(3,1) * t372 / 0.2e1 + t304 / 0.2e1 + t349;
t485 = Ifges(5,5) * t343 - Ifges(5,6) * t201;
t80 = t386 * t83;
t26 = t310 * t77 + t80;
t24 = qJ(6) * t285 + t26;
t70 = Ifges(7,5) * t484 + Ifges(7,6) * t285 + Ifges(7,3) * t132;
t73 = Ifges(6,4) * t484 - Ifges(6,2) * t132 + Ifges(6,6) * t285;
t483 = -t153 * mrSges(6,1) - t45 * mrSges(7,1) + t24 * mrSges(7,2) + t26 * mrSges(6,3) - t70 / 0.2e1 + t73 / 0.2e1;
t480 = -Ifges(6,6) / 0.2e1;
t479 = Ifges(7,6) / 0.2e1;
t126 = Ifges(5,2) * t343 + Ifges(5,6) * t285 + t399;
t478 = t126 / 0.2e1;
t415 = t285 / 0.2e1;
t430 = t484 / 0.2e1;
t348 = -Ifges(3,6) * qJD(2) / 0.2e1;
t469 = qJD(6) + t467;
t344 = t386 * t311;
t468 = -t310 * t86 + t327 * t386 + (t310 * t314 + t344) * qJD(4) * pkin(3);
t241 = t329 * t313;
t34 = t386 * t82 - t388;
t463 = qJD(6) - t34;
t122 = mrSges(6,1) * t285 - mrSges(6,3) * t484;
t123 = -mrSges(7,1) * t285 + mrSges(7,2) * t484;
t462 = t123 - t122;
t262 = t315 * t275;
t407 = pkin(7) * t312;
t203 = -pkin(9) * t379 + t262 + (-pkin(3) - t407) * t316;
t299 = pkin(7) * t378;
t227 = t312 * t275 + t299;
t381 = t312 * t313;
t211 = -pkin(9) * t381 + t227;
t143 = t311 * t203 + t314 * t211;
t461 = Ifges(4,5) * t216 + Ifges(4,6) * t217;
t460 = -t148 * mrSges(4,1) + t147 * mrSges(4,2);
t459 = pkin(1) * mrSges(3,2) * qJD(1);
t341 = Ifges(5,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t358 = t479 + t480;
t359 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t402 = Ifges(3,4) * t313;
t455 = -t358 * t132 - t341 * t285 - t359 * t484 - t204 * mrSges(4,1) - t24 * mrSges(7,3) - t25 * mrSges(6,1) - t93 * mrSges(5,1) - t201 * Ifges(5,5) - t343 * Ifges(5,6) - t297 * Ifges(4,3) - t264 * Ifges(4,5) - t263 * Ifges(4,6) - t348 + (Ifges(3,2) * t316 + t402) * qJD(1) / 0.2e1 - Ifges(6,6) * t435 - Ifges(7,6) * t434 + t205 * mrSges(4,2) + t23 * mrSges(7,1) + t26 * mrSges(6,2) + t94 * mrSges(5,2) - t475 * t430 - t500 * t415;
t451 = -Ifges(6,2) * t434 + Ifges(7,3) * t435 + t431 * t476 + t483;
t450 = Ifges(7,5) * t446 + Ifges(7,3) * t447 + t345 * t479;
t449 = -t64 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t447 + t345 * t480;
t448 = -t63 / 0.2e1;
t441 = pkin(1) * mrSges(3,1);
t437 = t111 / 0.2e1;
t436 = t112 / 0.2e1;
t426 = t343 / 0.2e1;
t424 = t201 / 0.2e1;
t423 = t216 / 0.2e1;
t422 = t217 / 0.2e1;
t240 = t266 * t313;
t421 = -t240 / 0.2e1;
t420 = -t241 / 0.2e1;
t419 = -t263 / 0.2e1;
t418 = -t264 / 0.2e1;
t414 = -t297 / 0.2e1;
t408 = pkin(4) * t310;
t149 = -qJD(2) * t323 - t209 * t313;
t373 = t315 * t273 + t370 * t407;
t139 = t328 * qJD(2) + (-t299 + (pkin(9) * t313 - t275) * t312) * qJD(3) + t373;
t166 = t312 * t273 + t275 * t366 + (-t313 * t369 - t316 * t367) * pkin(7);
t146 = -pkin(9) * t487 + t166;
t47 = -qJD(4) * t143 + t314 * t139 - t146 * t311;
t31 = pkin(4) * t370 - qJ(5) * t149 + qJD(5) * t241 + t47;
t150 = -qJD(2) * t324 + t241 * t458;
t46 = t311 * t139 + t314 * t146 + t203 * t364 - t211 * t365;
t35 = qJ(5) * t150 - qJD(5) * t240 + t46;
t10 = t310 * t31 + t386 * t35;
t404 = mrSges(5,3) * t343;
t403 = mrSges(5,3) * t201;
t383 = qJD(2) * mrSges(3,2);
t142 = t314 * t203 - t311 * t211;
t113 = -pkin(4) * t316 + t241 * qJ(5) + t142;
t119 = -qJ(5) * t240 + t143;
t69 = t310 * t113 + t386 * t119;
t302 = pkin(3) * t314 + pkin(4);
t243 = pkin(3) * t344 + t310 * t302;
t274 = pkin(3) * t381 + t313 * pkin(7);
t356 = Ifges(4,3) * t345 + t461;
t224 = pkin(3) * t487 + pkin(7) * t368;
t303 = -pkin(3) * t315 - pkin(2);
t354 = t386 * pkin(4);
t22 = t63 * mrSges(6,1) + t64 * mrSges(6,2);
t21 = t63 * mrSges(7,1) - t64 * mrSges(7,3);
t195 = -pkin(3) * t217 + qJD(2) * t306;
t214 = t314 * t282 + t283 * t311;
t340 = m(4) * t280 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t263 + mrSges(4,2) * t264 + mrSges(3,3) * t372;
t206 = pkin(4) * t240 + t274;
t160 = pkin(3) * t264 + t409;
t338 = mrSges(4,1) * t315 - mrSges(4,2) * t312;
t335 = Ifges(4,1) * t312 + t400;
t333 = Ifges(4,2) * t315 + t401;
t332 = Ifges(4,5) * t312 + Ifges(4,6) * t315;
t331 = t147 * t315 - t148 * t312;
t230 = pkin(4) * t329 + t303;
t124 = -pkin(4) * t150 + t224;
t52 = -mrSges(7,1) * t345 + t64 * mrSges(7,2);
t89 = -pkin(4) * t112 + t195;
t326 = -qJ(5) * t266 + t214;
t9 = t31 * t386 - t310 * t35;
t68 = t113 * t386 - t310 * t119;
t242 = t302 * t386 - t362;
t301 = -t354 - pkin(5);
t300 = qJ(6) + t408;
t277 = mrSges(3,3) * t371 - t383;
t234 = -pkin(5) - t242;
t233 = qJ(6) + t243;
t226 = -pkin(7) * t380 + t262;
t222 = mrSges(4,1) * t297 - mrSges(4,3) * t264;
t221 = -mrSges(4,2) * t297 + mrSges(4,3) * t263;
t220 = -pkin(7) * t352 + t249;
t196 = t266 * t310 + t329 * t386;
t193 = -mrSges(4,2) * t345 + mrSges(4,3) * t217;
t192 = mrSges(4,1) * t345 - mrSges(4,3) * t216;
t177 = -qJ(5) * t329 + t215;
t173 = -t310 * t240 - t241 * t386;
t172 = t240 * t386 - t241 * t310;
t171 = mrSges(5,1) * t285 - t403;
t170 = -mrSges(5,2) * t285 + t404;
t167 = -qJD(3) * t227 + t373;
t159 = -mrSges(4,1) * t217 + mrSges(4,2) * t216;
t152 = t216 * Ifges(4,1) + t217 * Ifges(4,4) + Ifges(4,5) * t345;
t151 = t216 * Ifges(4,4) + t217 * Ifges(4,2) + Ifges(4,6) * t345;
t138 = -mrSges(5,1) * t343 + mrSges(5,2) * t201;
t121 = -mrSges(6,2) * t285 - mrSges(6,3) * t132;
t120 = -mrSges(7,2) * t132 + mrSges(7,3) * t285;
t117 = t177 * t386 + t310 * t326;
t116 = t177 * t310 - t326 * t386;
t115 = pkin(5) * t196 - qJ(6) * t197 + t230;
t99 = -mrSges(5,2) * t345 + mrSges(5,3) * t112;
t98 = mrSges(5,1) * t345 - mrSges(5,3) * t111;
t90 = pkin(5) * t172 - qJ(6) * t173 + t206;
t88 = t149 * t386 + t310 * t150;
t87 = t149 * t310 - t150 * t386;
t79 = mrSges(6,1) * t132 + mrSges(6,2) * t484;
t78 = mrSges(7,1) * t132 - mrSges(7,3) * t484;
t67 = -mrSges(5,1) * t112 + mrSges(5,2) * t111;
t66 = t409 + t505;
t65 = t316 * pkin(5) - t68;
t62 = -qJ(6) * t316 + t69;
t54 = t111 * Ifges(5,1) + t112 * Ifges(5,4) + Ifges(5,5) * t345;
t53 = t111 * Ifges(5,4) + t112 * Ifges(5,2) + Ifges(5,6) * t345;
t51 = mrSges(6,1) * t345 - mrSges(6,3) * t64;
t50 = -mrSges(6,2) * t345 - mrSges(6,3) * t63;
t49 = -mrSges(7,2) * t63 + mrSges(7,3) * t345;
t48 = t160 + t505;
t33 = t310 * t82 + t80;
t16 = pkin(5) * t87 - qJ(6) * t88 - qJD(6) * t173 + t124;
t11 = pkin(5) * t63 - qJ(6) * t64 - qJD(6) * t484 + t89;
t8 = -pkin(5) * t370 - t9;
t7 = qJ(6) * t370 - qJD(6) * t316 + t10;
t1 = [-(t356 + t457 + t461) * t316 / 0.2e1 + (pkin(7) * t159 + t151 * t413 + t152 * t412 + t336 * t423 + t334 * t422 + (-t147 * t312 - t148 * t315) * mrSges(4,3) + (t180 * t413 - t315 * t179 / 0.2e1 + t332 * t414 + t280 * t338 + t333 * t419 + t335 * t418 + (t204 * t312 - t205 * t315) * mrSges(4,3)) * qJD(3) + ((-0.2e1 * t441 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t398 / 0.2e1 - t397 / 0.2e1) * t313 + Ifges(5,5) * t420 + Ifges(5,6) * t421 + t359 * t173 + t358 * t172) * qJD(1) + t348 - pkin(7) * t277 - t455) * qJD(2)) * t313 + t173 * t503 + (Ifges(6,4) * t435 + Ifges(7,5) * t434 + t475 * t415 + t477 * t430 + t501) * t88 + t150 * t478 + (Ifges(6,4) * t173 - Ifges(6,2) * t172) * t448 + (t340 * pkin(7) + t349 - 0.2e1 * t459 + t486) * t368 + (Ifges(7,5) * t173 + Ifges(7,3) * t172) * t447 + (-t172 * t2 + t173 * t3) * mrSges(7,2) + m(5) * (t142 * t30 + t143 * t29 + t195 * t274 + t223 * t224 + t46 * t94 + t47 * t93) + m(7) * (t11 * t90 + t16 * t45 + t2 * t62 + t23 * t8 + t24 * t7 + t3 * t65) + m(6) * (t10 * t26 + t124 * t153 + t206 * t89 + t25 * t9 + t5 * t68 + t6 * t69) + (-Ifges(6,2) * t435 + Ifges(7,3) * t434 + t474 * t415 + t476 * t430 - t483) * t87 + t274 * t67 + t226 * t192 + t227 * t193 + t166 * t221 + t167 * t222 + t223 * (-mrSges(5,1) * t150 + mrSges(5,2) * t149) + t224 * t138 + t206 * t22 + t11 * (mrSges(7,1) * t172 - mrSges(7,3) * t173) + t89 * (mrSges(6,1) * t172 + mrSges(6,2) * t173) + t46 * t170 + t47 * t171 + t142 * t98 + t143 * t99 + t149 * t127 / 0.2e1 + t8 * t123 + t124 * t79 + t7 * t120 + t10 * t121 + t9 * t122 + t90 * t21 + t16 * t78 + t68 * t51 + t69 * t50 + t65 * t52 + t62 * t49 + m(4) * (t147 * t227 + t148 * t226 + t205 * t166 + t204 * t167) + (-t149 * t93 + t150 * t94 - t240 * t29 + t241 * t30) * mrSges(5,3) + (-Ifges(5,4) * t241 - Ifges(5,2) * t240) * t436 + t195 * (mrSges(5,1) * t240 - mrSges(5,2) * t241) + (-Ifges(5,1) * t241 - Ifges(5,4) * t240) * t437 + (-Ifges(5,5) * t437 - Ifges(5,6) * t436 - Ifges(6,6) * t448 - Ifges(7,6) * t447 - t475 * t446 + (0.3e1 / 0.2e1 * Ifges(3,4) * t368 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1 + (m(4) * pkin(7) + t337) * pkin(7) - t341) * t370) * qJD(1) - t456 + t460) * t316 + (t172 * t476 + t173 * t477) * t446 + (Ifges(5,5) * t149 + Ifges(5,6) * t150) * t415 + t54 * t420 + t53 * t421 + (Ifges(5,1) * t149 + Ifges(5,4) * t150) * t424 + (Ifges(5,4) * t149 + Ifges(5,2) * t150) * t426 + t172 * t449 + t172 * t450 + (-t172 * t6 - t173 * t5) * mrSges(6,3); t197 * t503 + (mrSges(5,1) * t507 - mrSges(5,2) * t510) * t223 + (-t266 * t30 - t29 * t329 - t507 * t94 + t510 * t93) * mrSges(5,3) + (m(4) * (-qJD(3) * t330 + t331) - t192 * t312 + t193 * t315 + (-t221 * t312 - t222 * t315) * qJD(3)) * pkin(8) + (-Ifges(5,4) * t229 - Ifges(5,2) * t228) * t427 + (t229 / 0.2e1 - t208 / 0.2e1) * t127 + (-Ifges(5,1) * t208 - Ifges(5,4) * t209) * t424 + (-Ifges(5,4) * t208 - Ifges(5,2) * t209) * t426 + (t228 / 0.2e1 - t209 / 0.2e1) * t126 + (-Ifges(5,1) * t229 - Ifges(5,4) * t228) * t425 + (t73 - t70) * (t161 / 0.2e1 - t140 / 0.2e1) - t329 * t53 / 0.2e1 + t195 * (mrSges(5,1) * t329 + mrSges(5,2) * t266) + (Ifges(5,4) * t266 - Ifges(5,2) * t329) * t436 + (Ifges(5,1) * t266 - Ifges(5,4) * t329) * t437 + (t195 * t303 + t214 * t30 + t215 * t29 + t223 * t506 + t489 * t93 + t490 * t94) * m(5) + (t138 * t410 + t317) * qJD(3) + t492 * t78 + t473 * (-t162 / 0.2e1 + t141 / 0.2e1) + t493 * t121 + t494 * t122 + (-t116 * t5 + t117 * t6 + t153 * t488 + t230 * t89 + t25 * t494 + t26 * t493) * m(6) + t331 * mrSges(4,3) - m(4) * (t204 * t219 + t205 * t220) + ((t349 - t304 / 0.2e1 + t459 + ((-m(4) * pkin(2) - mrSges(3,1) - t338) * qJD(2) - t340) * pkin(7) - t486) * t316 + ((t277 + t383) * pkin(7) + (t441 + t402 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t316) * qJD(1) + t348 + t455) * t313 + (Ifges(5,5) * t266 - Ifges(5,6) * t329 + t196 * t474 + t197 * t475 + t332) * t370 / 0.2e1) * qJD(1) + (Ifges(6,4) * t162 + Ifges(7,5) * t141 - Ifges(6,2) * t161 + Ifges(7,3) * t140) * t434 + (t52 - t51) * t116 + (Ifges(6,4) * t141 + Ifges(7,5) * t162 - Ifges(6,2) * t140 + Ifges(7,3) * t161) * t435 + (t140 * t476 + t141 * t477) * t430 + (t196 * t476 + t197 * t477) * t446 + (t50 + t49) * t117 + (mrSges(6,1) * t377 + mrSges(6,2) * t376) * t153 + (mrSges(7,1) * t377 - mrSges(7,3) * t376) * t45 + (-t196 * t2 + t197 * t3 + t23 * t376 - t24 * t377) * mrSges(7,2) + (-t196 * t6 - t197 * t5 - t25 * t376 - t26 * t377) * mrSges(6,3) + t495 * t120 + (t11 * t115 + t116 * t3 + t117 * t2 + t23 * t496 + t24 * t495 + t45 * t492) * m(7) + t496 * t123 + (-Ifges(5,5) * t208 - Ifges(5,6) * t209 + t140 * t474 + t141 * t475) * t415 + (-Ifges(5,5) * t229 - Ifges(5,6) * t228 + t161 * t474 + t162 * t475) * t416 + (t161 * t476 + t162 * t477) * t431 + t312 * t152 / 0.2e1 + t488 * t79 + t489 * t171 + t490 * t170 + t303 * t67 + t266 * t54 / 0.2e1 - t258 * t138 + t230 * t22 - t220 * t221 - t219 * t222 + t214 * t98 + t215 * t99 + t11 * (mrSges(7,1) * t196 - mrSges(7,3) * t197) + t89 * (mrSges(6,1) * t196 + mrSges(6,2) * t197) - pkin(2) * t159 + t115 * t21 + t151 * t412 + t333 * t422 + t335 * t423 + (Ifges(7,5) * t197 + Ifges(7,3) * t196) * t447 + (Ifges(6,4) * t197 - Ifges(6,2) * t196) * t448 + t196 * t449 + t196 * t450; t451 * t484 + (t474 * t484 + t485) * t416 + t467 * t121 + t462 * t468 + (-t153 * t160 + t242 * t5 + t243 * t6 - t25 * t468 + t26 * t467) * m(6) + (t2 * t233 + t23 * t468 + t234 * t3 + t24 * t469 - t45 * t48) * m(7) + t469 * t120 + t504 + t201 * t478 + (t201 * t94 + t343 * t93) * mrSges(5,3) + t356 + (t204 * t263 + t205 * t264) * mrSges(4,3) - t460 + (-t264 * t138 + t311 * t99 + t314 * t98 + (t170 * t314 - t171 * t311) * qJD(4) + (0.2e1 * t223 * t418 + t29 * t311 + t30 * t314 + t364 * t94 - t365 * t93) * m(5)) * pkin(3) - m(5) * (t101 * t93 + t102 * t94) + (-Ifges(4,2) * t264 + t180 + t259) * t419 - t280 * (mrSges(4,1) * t264 + mrSges(4,2) * t263) + t242 * t51 + t243 * t50 + t233 * t49 + t234 * t52 - t204 * t221 + t205 * t222 - t102 * t170 - t101 * t171 - t160 * t79 - t48 * t78 + (Ifges(4,5) * t263 - Ifges(4,6) * t264) * t414 + t179 * t417 + (Ifges(4,1) * t263 - t391) * t418; -t462 * t33 + t463 * t120 + (t2 * t300 - t23 * t33 + t24 * t463 + t3 * t301 - t45 * t66) * m(7) + (t171 + t403) * t94 + t51 * t354 + (-t170 + t404) * t93 + t300 * t49 + t301 * t52 + ((t310 * t6 + t386 * t5) * pkin(4) - t153 * t409 + t25 * t33 - t26 * t34) * m(6) - t34 * t121 - t66 * t78 + t485 * t416 + (t416 * t474 + t451) * t484 + t50 * t408 + t126 * t424 - t79 * t409 + t504; -(-t120 - t121) * t132 - t462 * t484 + t21 + t22 + (t132 * t24 - t23 * t484 + t11) * m(7) + (t132 * t26 + t25 * t484 + t89) * m(6); -t285 * t120 + t484 * t78 + 0.2e1 * (t3 / 0.2e1 + t45 * t430 + t24 * t416) * m(7) + t52;];
tauc  = t1(:);
