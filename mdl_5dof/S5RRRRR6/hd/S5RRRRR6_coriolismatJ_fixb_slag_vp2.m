% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:44
% EndTime: 2022-01-20 12:08:00
% DurationCPUTime: 6.77s
% Computational Cost: add. (17888->312), mult. (34216->409), div. (0->0), fcn. (35872->8), ass. (0->209)
t269 = sin(qJ(3));
t430 = -pkin(8) - pkin(7);
t253 = t430 * t269;
t273 = cos(qJ(3));
t254 = t430 * t273;
t268 = sin(qJ(4));
t272 = cos(qJ(4));
t200 = t253 * t268 - t272 * t254;
t237 = -t268 * t269 + t272 * t273;
t417 = pkin(9) * t237;
t158 = t200 + t417;
t267 = sin(qJ(5));
t271 = cos(qJ(5));
t238 = -t268 * t273 - t272 * t269;
t233 = t238 * pkin(9);
t447 = t272 * t253 + t268 * t254;
t461 = t447 + t233;
t494 = -t158 * t267 + t271 * t461;
t536 = t494 / 0.2e1;
t270 = sin(qJ(2));
t259 = pkin(1) * t270 + pkin(7);
t408 = pkin(8) + t259;
t234 = t408 * t269;
t235 = t408 * t273;
t176 = -t234 * t268 + t272 * t235;
t136 = t176 + t417;
t448 = -t272 * t234 - t268 * t235;
t460 = t448 + t233;
t495 = -t136 * t267 + t271 * t460;
t535 = t495 / 0.2e1;
t534 = qJD(3) + qJD(4);
t186 = t237 * t267 - t238 * t271;
t312 = t271 * t237 + t238 * t267;
t502 = Ifges(6,5) * t312 - Ifges(6,6) * t186;
t505 = t495 * mrSges(6,2);
t74 = t136 * t271 + t267 * t460;
t514 = t74 * mrSges(6,1);
t520 = -t514 / 0.2e1 - t505 / 0.2e1;
t528 = t502 + 0.2e1 * t520;
t533 = t528 * qJD(5);
t506 = t494 * mrSges(6,2);
t103 = t158 * t271 + t267 * t461;
t513 = t103 * mrSges(6,1);
t521 = -t513 / 0.2e1 - t506 / 0.2e1;
t529 = t502 + 0.2e1 * t521;
t532 = t529 * qJD(5);
t399 = Ifges(6,4) * t186;
t425 = t186 / 0.2e1;
t426 = -t186 / 0.2e1;
t427 = t312 / 0.2e1;
t278 = (Ifges(6,2) * t312 + t399) * t426 + (Ifges(6,1) * t312 - t399) * t425 + (0.2e1 * Ifges(6,4) * t312 + (Ifges(6,1) - Ifges(6,2)) * t186) * t427;
t296 = Ifges(5,4) * t237 ^ 2 + (-Ifges(5,4) * t238 + (Ifges(5,2) - Ifges(5,1)) * t237) * t238 + t278;
t187 = -mrSges(5,1) * t238 + mrSges(5,2) * t237;
t262 = -t273 * pkin(3) - pkin(2);
t342 = t262 * t187;
t204 = -t237 * pkin(4) + t262;
t458 = mrSges(6,1) * t186 + mrSges(6,2) * t312;
t356 = t204 * t458;
t531 = t296 + t342 + t356;
t274 = cos(qJ(2));
t421 = pkin(1) * t274;
t244 = t262 - t421;
t345 = t244 * t187;
t201 = t204 - t421;
t357 = t201 * t458;
t530 = t296 + t345 + t357;
t482 = -t200 * mrSges(5,1) - t447 * mrSges(5,2);
t497 = -t506 - t513;
t527 = t482 + t497;
t483 = -t176 * mrSges(5,1) - t448 * mrSges(5,2);
t496 = -t505 - t514;
t526 = t483 + t496;
t207 = t238 * t421;
t208 = t237 * t421;
t138 = t207 * t271 - t208 * t267;
t139 = t207 * t267 + t208 * t271;
t337 = t138 * mrSges(6,1) / 0.2e1 - t139 * mrSges(6,2) / 0.2e1;
t292 = t207 * mrSges(5,1) / 0.2e1 - t208 * mrSges(5,2) / 0.2e1 + t337;
t385 = t312 * mrSges(6,3);
t429 = -t494 / 0.2e1;
t432 = -t495 / 0.2e1;
t462 = -t296 + (t535 + t536) * t385;
t472 = (t201 / 0.2e1 + t204 / 0.2e1) * t458;
t524 = (t429 + t432) * t385 - t472 + (-t244 / 0.2e1 - t262 / 0.2e1) * t187 + t292 + t462;
t507 = (t429 + t536) * mrSges(6,2);
t523 = qJD(5) * t507;
t508 = (t432 + t535) * mrSges(6,2);
t522 = qJD(5) * t508;
t326 = mrSges(6,3) * t426;
t328 = mrSges(6,3) * t427;
t386 = t186 * mrSges(6,3);
t378 = t238 * mrSges(5,3);
t485 = (t200 + t176) * t378;
t471 = t485 / 0.2e1;
t519 = t357 / 0.2e1 + t356 / 0.2e1 + t345 / 0.2e1 + t342 / 0.2e1 + t292 + t471 - t462 + (t494 + t495) * t328 + (-t326 - t386 / 0.2e1) * (t74 + t103);
t518 = qJD(1) * t508 + qJD(2) * t507;
t418 = pkin(4) * t238;
t419 = pkin(3) * t269;
t211 = -t418 + t419;
t515 = m(6) * t211;
t509 = -t103 * t271 + t267 * t494;
t300 = m(6) * (t267 * t495 - t271 * t74);
t260 = pkin(3) * t272 + pkin(4);
t339 = t267 * t268;
t220 = -pkin(3) * t339 + t260 * t271;
t225 = (t271 * t272 - t339) * pkin(3);
t498 = -t220 + t225;
t492 = t273 / 0.2e1;
t489 = m(6) * t418;
t338 = t268 * t271;
t221 = pkin(3) * t338 + t260 * t267;
t224 = (-t267 * t272 - t338) * pkin(3);
t484 = t221 + t224;
t436 = -m(6) / 0.2e1;
t481 = m(5) * t419;
t470 = t269 ^ 2 + t273 ^ 2;
t468 = -Ifges(4,2) * t273 / 0.2e1 - Ifges(4,4) * t269 + Ifges(4,1) * t492;
t437 = m(5) / 0.2e1;
t459 = (t498 * t103 + t484 * t494) * t436 - t521;
t454 = (t244 + t262) * t419;
t446 = t470 * t274;
t445 = -mrSges(4,1) * t273 + mrSges(4,2) * t269;
t222 = t224 * mrSges(6,1);
t380 = t225 * mrSges(6,2);
t442 = (mrSges(5,1) * t268 + mrSges(5,2) * t272) * pkin(3) - t222 + t380;
t435 = m(6) / 0.2e1;
t441 = (t138 * t220 + t139 * t221) * t435 + (t207 * t272 + t208 * t268) * pkin(3) * t437;
t264 = Ifges(4,4) * t273;
t248 = -Ifges(4,2) * t269 + t264;
t249 = Ifges(4,1) * t269 + t264;
t440 = (t248 + t249) * t492 + t468 * t269;
t108 = -mrSges(6,1) * t312 + mrSges(6,2) * t186;
t188 = -mrSges(5,1) * t237 - mrSges(5,2) * t238;
t371 = t211 * t108 + t188 * t419;
t439 = t371 + t440;
t438 = (-t138 * t186 + t139 * t312) * mrSges(6,3) + (t207 * t238 + t208 * t237) * mrSges(5,3);
t434 = -pkin(4) / 0.2e1;
t433 = m(6) * pkin(4);
t423 = m(5) * t270;
t422 = m(6) * t270;
t375 = t273 * mrSges(4,2);
t376 = t269 * mrSges(4,1);
t246 = t375 + t376;
t420 = pkin(2) * t246;
t398 = pkin(3) * qJD(3);
t397 = pkin(4) * qJD(4);
t377 = t267 * mrSges(6,1);
t261 = -pkin(2) - t421;
t343 = t261 * t246;
t5 = t201 * t515 + t244 * t481 + t343 + t439 + t530;
t374 = t5 * qJD(1);
t348 = t238 * t108;
t333 = pkin(4) * t348;
t6 = -t201 * t489 - t333 + t530;
t373 = t6 * qJD(1);
t212 = t221 * mrSges(6,1);
t61 = t220 * mrSges(6,2) + t212;
t370 = qJD(5) * t61;
t15 = t278 + t357;
t364 = t15 * qJD(1);
t283 = (mrSges(4,3) * t470 - mrSges(3,2)) * t274 + (-mrSges(3,1) + t108 + t188 + t445) * t270;
t21 = m(6) * (t138 * t495 + t139 * t74) + m(5) * (t176 * t208 + t207 * t448) + (t201 * t422 + t244 * t423 + m(4) * (t446 * t259 + t261 * t270) + t283) * pkin(1) + t438;
t353 = t21 * qJD(1);
t352 = t220 * t312;
t167 = t267 * pkin(4) * t386;
t327 = -t385 / 0.2e1;
t336 = -t167 / 0.2e1 + t271 * pkin(4) * t327;
t334 = t201 + t204;
t332 = t271 * t385;
t331 = t272 * t237 * mrSges(5,3);
t330 = -t421 / 0.2e1;
t314 = t334 * t238;
t311 = Ifges(5,5) * t237 + Ifges(5,6) * t238 + t502;
t310 = -t220 / 0.2e1 + t271 * t434;
t280 = t334 * t211 * t436 - t371;
t2 = (mrSges(4,2) * t330 - t249 / 0.2e1 - t248 / 0.2e1) * t273 + (pkin(2) / 0.2e1 - t261 / 0.2e1) * t246 + (mrSges(4,1) * t330 - t468) * t269 + t280 - t454 * m(5) / 0.2e1 + t441 + t524;
t7 = t204 * t515 + t262 * t481 - t420 + t439 + t531;
t307 = -t2 * qJD(1) + t7 * qJD(2);
t306 = -t167 + t311;
t10 = -t204 * t489 - t333 + t531;
t297 = m(6) * (t138 * t271 + t139 * t267);
t4 = 0.2e1 * (t297 / 0.4e1 + m(6) * t314 / 0.4e1) * pkin(4) + t333 + t524;
t305 = -t4 * qJD(1) + t10 * qJD(2);
t18 = t278 + t356;
t276 = t278 + t472;
t8 = t276 - t337;
t304 = t8 * qJD(1) + t18 * qJD(2);
t302 = -t212 / 0.2e1 + t377 * t434;
t125 = t221 * t386;
t301 = t220 * t327 - t125 / 0.2e1 + t224 * t326 + t225 * t328;
t279 = (t484 * t495 + t498 * t74) * t435 + t301 + t520;
t12 = t279 + t167 / 0.2e1 + (-t300 / 0.2e1 + t332 / 0.2e1) * pkin(4) - t520;
t285 = t509 * t433 / 0.2e1 + t336 + t521;
t14 = t125 / 0.2e1 + (t352 / 0.2e1 - t225 * t312 / 0.2e1 + t224 * t425) * mrSges(6,3) + t285 + t459;
t60 = -m(6) * (t220 * t224 + t221 * t225) + t442;
t294 = t12 * qJD(1) - t14 * qJD(2) - t60 * qJD(3);
t293 = -t61 * qJD(3) + t518;
t243 = (t271 * mrSges(6,2) + t377) * pkin(4);
t49 = -t222 / 0.2e1 + (t225 / 0.2e1 + t310) * mrSges(6,2) + t302;
t289 = -qJD(3) * t49 + qJD(4) * t243 - t518;
t284 = t268 * pkin(3) * t378 - mrSges(6,3) * t352 + Ifges(4,5) * t273 - Ifges(4,6) * t269 - t125 + t311;
t236 = t243 * qJD(5);
t50 = -t380 / 0.2e1 + t222 / 0.2e1 + t310 * mrSges(6,2) + t302;
t13 = t301 + t285 + t311 - t459 + t482;
t11 = pkin(4) * t300 / 0.2e1 + t279 + t311 + t336 + t520 + t483;
t9 = t276 + t337;
t3 = -t485 / 0.2e1 + (-t314 * t435 - t348 + t297 / 0.2e1) * pkin(4) + t519;
t1 = -t280 + t454 * t437 + (-t375 / 0.2e1 - t376 / 0.2e1) * t421 - t420 / 0.2e1 + t343 / 0.2e1 - t471 + t440 + t441 + t519;
t16 = [qJD(2) * t21 + qJD(3) * t5 + qJD(4) * t6 + qJD(5) * t15, t1 * qJD(3) + t3 * qJD(4) + t9 * qJD(5) + t353 + (0.2e1 * (t103 * t139 + t138 * t494) * t435 + 0.2e1 * (t200 * t208 + t207 * t447) * t437 + (t204 * t422 + t262 * t423 + m(4) * (-pkin(2) * t270 + t446 * pkin(7)) + t283) * pkin(1) + t438) * qJD(2), t374 + t1 * qJD(2) + (t284 + m(6) * (-t220 * t74 + t221 * t495) + t445 * t259 + t526) * qJD(3) + t11 * qJD(4) + t533 + (-t331 + m(5) * (-t176 * t272 + t268 * t448)) * t398, t373 + t3 * qJD(2) + t11 * qJD(3) + (t306 + t526) * qJD(4) + t533 + (-t332 + t300) * t397, t364 + t9 * qJD(2) + (t502 + t496) * qJD(5) + t534 * t528; -qJD(3) * t2 - qJD(4) * t4 + qJD(5) * t8 - t353, qJD(3) * t7 + qJD(4) * t10 + qJD(5) * t18, (t284 + m(6) * (-t103 * t220 + t221 * t494) + t445 * pkin(7) + t527) * qJD(3) + t13 * qJD(4) + t532 + (-t331 + m(5) * (-t200 * t272 + t268 * t447)) * t398 + t307, t13 * qJD(3) + (t306 + t527) * qJD(4) + t532 + (m(6) * t509 - t332) * t397 + t305, (t502 + t497) * qJD(5) + t304 + t534 * t529; qJD(2) * t2 + qJD(4) * t12 - t374 + t522, -qJD(4) * t14 - t307 + t523, -qJD(4) * t60 - t370, ((t224 * t271 + t225 * t267) * t433 - t442) * qJD(4) + t50 * qJD(5) + t294, t50 * qJD(4) + t293 - t370; qJD(2) * t4 - qJD(3) * t12 - t373 + t522, qJD(3) * t14 - t305 + t523, qJD(5) * t49 - t294, -t236, -t236 - t289; -qJD(2) * t8 - t508 * t534 - t364, -t507 * t534 - t304, -qJD(4) * t49 - t293, t289, 0;];
Cq = t16;
