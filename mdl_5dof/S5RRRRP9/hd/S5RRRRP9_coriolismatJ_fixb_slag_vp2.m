% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:52
% EndTime: 2019-12-31 22:04:11
% DurationCPUTime: 9.00s
% Computational Cost: add. (12954->583), mult. (27881->776), div. (0->0), fcn. (26623->6), ass. (0->300)
t631 = Ifges(6,4) + Ifges(5,5);
t630 = Ifges(5,6) - Ifges(6,6);
t368 = cos(qJ(3));
t363 = t368 * pkin(7);
t464 = t368 * pkin(8) + t363;
t365 = sin(qJ(4));
t366 = sin(qJ(3));
t478 = t365 * t366;
t554 = cos(qJ(4));
t579 = -pkin(8) - pkin(7);
t591 = t464 * t554 + t579 * t478;
t575 = t591 / 0.2e1;
t549 = pkin(3) * t368;
t350 = -pkin(2) - t549;
t444 = t554 * t368;
t393 = t444 - t478;
t445 = t554 * t366;
t394 = t365 * t368 + t445;
t404 = -pkin(4) * t393 - qJ(5) * t394;
t175 = t350 + t404;
t201 = -mrSges(6,1) * t393 - mrSges(6,3) * t394;
t611 = m(6) * t175 + t201;
t629 = Ifges(5,3) + Ifges(6,2);
t369 = cos(qJ(2));
t290 = t394 * t369;
t451 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t628 = t451 * t290;
t417 = t631 * t393 - t630 * t394;
t511 = t393 * Ifges(6,5);
t210 = Ifges(6,1) * t394 - t511;
t312 = Ifges(5,4) * t393;
t212 = Ifges(5,1) * t394 + t312;
t627 = t212 + t210;
t367 = sin(qJ(2));
t476 = t366 * t367;
t288 = t365 * t476 - t367 * t444;
t271 = Ifges(6,5) * t288;
t289 = t394 * t367;
t153 = -Ifges(6,6) * t369 + Ifges(6,3) * t289 - t271;
t535 = Ifges(5,4) * t288;
t626 = t153 - t271 + t535 + (-Ifges(5,1) - Ifges(6,1)) * t289;
t532 = Ifges(6,5) * t289;
t157 = -Ifges(6,1) * t288 - t369 * Ifges(6,4) + t532;
t274 = Ifges(5,4) * t289;
t159 = -Ifges(5,1) * t288 - t369 * Ifges(5,5) - t274;
t625 = Ifges(5,2) * t288 + t157 + t159 - t274;
t235 = t365 * t464 - t579 * t445;
t609 = t288 * t575 - t235 * t289 / 0.2e1;
t199 = mrSges(6,1) * t394 - mrSges(6,3) * t393;
t200 = mrSges(5,1) * t394 + mrSges(5,2) * t393;
t203 = Ifges(6,3) * t394 + t511;
t207 = -Ifges(5,2) * t394 + t312;
t309 = Ifges(6,5) * t394;
t209 = Ifges(6,1) * t393 + t309;
t509 = t394 * Ifges(5,4);
t211 = Ifges(5,1) * t393 - t509;
t204 = -Ifges(6,3) * t393 + t309;
t208 = Ifges(5,2) * t393 + t509;
t593 = t204 / 0.2e1 - t208 / 0.2e1;
t624 = -(t203 / 0.2e1 - t212 / 0.2e1 - t207 / 0.2e1 - t210 / 0.2e1) * t393 - (-t211 / 0.2e1 - t209 / 0.2e1 - t593) * t394 + t175 * t199 + t350 * t200;
t583 = m(6) / 0.2e1;
t623 = 0.2e1 * t583;
t622 = -mrSges(6,1) / 0.2e1;
t291 = t393 * t369;
t617 = t291 / 0.2e1;
t558 = t367 / 0.2e1;
t238 = -mrSges(6,2) * t290 + mrSges(6,3) * t367;
t241 = -mrSges(5,2) * t367 - mrSges(5,3) * t290;
t244 = mrSges(5,1) * t367 - mrSges(5,3) * t291;
t508 = t367 * mrSges(6,1);
t513 = t291 * mrSges(6,2);
t245 = -t508 + t513;
t548 = t367 * pkin(2);
t336 = -pkin(7) * t369 + t548;
t267 = pkin(6) * t476 + t368 * t336;
t474 = t367 * t368;
t268 = -pkin(6) * t474 + t366 * t336;
t551 = pkin(3) * t365;
t345 = qJ(5) + t551;
t458 = t554 * pkin(3);
t349 = -t458 - pkin(4);
t472 = t368 * t369;
t185 = t367 * pkin(3) - pkin(8) * t472 + t267;
t475 = t366 * t369;
t227 = -pkin(8) * t475 + t268;
t96 = t365 * t185 + t554 * t227;
t84 = qJ(5) * t367 + t96;
t95 = t185 * t554 - t365 * t227;
t85 = -t367 * pkin(4) - t95;
t382 = t84 * mrSges(6,3) / 0.2e1 + t85 * t622 + t95 * mrSges(5,1) / 0.2e1 - t96 * mrSges(5,2) / 0.2e1 + t631 * t617 + t629 * t558 + t628;
t420 = t458 / 0.2e1;
t454 = t551 / 0.2e1;
t561 = t345 / 0.2e1;
t581 = m(5) * pkin(3);
t610 = (t345 * t84 + t349 * t85) * t583 + Ifges(4,3) * t558 + t267 * mrSges(4,1) / 0.2e1 - t268 * mrSges(4,2) / 0.2e1 + t238 * t561 + t349 * t245 / 0.2e1 + (t365 * t96 + t554 * t95) * t581 / 0.2e1 + Ifges(4,5) * t472 / 0.2e1 - Ifges(4,6) * t475 / 0.2e1 + t241 * t454 + t244 * t420 + t382;
t356 = m(6) * qJ(5) + mrSges(6,3);
t608 = qJD(4) * t356;
t607 = t356 * qJD(5);
t515 = t289 * mrSges(6,2);
t259 = -t515 / 0.2e1;
t605 = 0.2e1 * t259;
t603 = mrSges(6,1) + mrSges(5,1);
t358 = Ifges(4,5) * t368;
t529 = Ifges(4,6) * t366;
t601 = Ifges(3,4) - t358 / 0.2e1 + t529 / 0.2e1;
t471 = t369 * qJ(5);
t328 = -pkin(2) * t369 - t367 * pkin(7) - pkin(1);
t460 = pkin(6) * t472;
t224 = t460 + (-pkin(8) * t367 + t328) * t366;
t446 = t554 * t224;
t314 = t368 * t328;
t416 = -pkin(8) * t474 + t314;
t184 = (-pkin(6) * t366 - pkin(3)) * t369 + t416;
t480 = t365 * t184;
t90 = t446 + t480;
t76 = t90 - t471;
t357 = t369 * mrSges(6,3);
t239 = -t357 - t515;
t514 = t289 * mrSges(5,3);
t240 = mrSges(5,2) * t369 - t514;
t600 = t239 + t240;
t359 = Ifges(4,4) * t368;
t599 = -Ifges(4,2) * t366 + t359;
t332 = Ifges(4,1) * t366 + t359;
t461 = pkin(6) * t475;
t223 = t416 - t461;
t479 = t365 * t224;
t103 = t223 * t554 - t479;
t102 = t223 * t365 + t446;
t564 = t394 / 0.2e1;
t443 = t102 * t564;
t566 = t393 / 0.2e1;
t597 = t103 * t566 + t443;
t418 = t630 * t288 - t631 * t289;
t596 = -t267 * t366 + t268 * t368;
t595 = -mrSges(4,1) * t368 + mrSges(4,2) * t366;
t594 = t240 / 0.2e1 + t239 / 0.2e1;
t361 = t367 * pkin(6);
t326 = pkin(3) * t476 + t361;
t503 = qJ(5) * t288;
t406 = pkin(4) * t289 + t503;
t131 = t406 + t326;
t170 = mrSges(6,1) * t289 + mrSges(6,3) * t288;
t592 = -m(6) * t131 - t170;
t590 = -m(6) * t76 - t600;
t516 = t288 * mrSges(5,3);
t242 = -mrSges(5,1) * t369 + t516;
t517 = t288 * mrSges(6,2);
t243 = mrSges(6,1) * t369 - t517;
t89 = t184 * t554 - t479;
t79 = t369 * pkin(4) - t89;
t589 = -m(6) * t79 + t242 - t243;
t405 = -pkin(4) * t288 + qJ(5) * t289;
t462 = pkin(3) * t474;
t144 = t405 + t462;
t171 = mrSges(5,1) * t289 - mrSges(5,2) * t288;
t198 = pkin(4) * t394 - qJ(5) * t393;
t550 = pkin(3) * t366;
t179 = t198 + t550;
t456 = mrSges(4,3) * t476;
t321 = mrSges(4,2) * t369 - t456;
t323 = -mrSges(4,1) * t369 - mrSges(4,3) * t474;
t512 = t393 * mrSges(6,2);
t304 = t512 / 0.2e1;
t504 = t90 * t394;
t447 = -t504 / 0.2e1;
t510 = t394 * mrSges(6,2);
t449 = -t510 / 0.2e1;
t155 = -Ifges(5,2) * t289 - Ifges(5,6) * t369 - t535;
t407 = -Ifges(6,3) * t288 - t532;
t412 = -mrSges(6,1) * t288 + mrSges(6,3) * t289;
t413 = -mrSges(5,1) * t288 - mrSges(5,2) * t289;
t587 = -t350 * t413 / 0.2e1 - t175 * t412 / 0.2e1 - t326 * t200 / 0.2e1 - t289 * t203 / 0.4e1 - t288 * t208 / 0.4e1 - t131 * t199 / 0.2e1 + t417 * t369 / 0.4e1 + (t211 + t209 + t204) * t288 / 0.4e1 + (t207 + t627) * t289 / 0.4e1 + (t155 / 0.4e1 - t626 / 0.4e1) * t394 + (t407 / 0.4e1 - t625 / 0.4e1) * t393;
t375 = t79 * t304 + t76 * t449 - t587 + t609 * mrSges(6,2) + (t447 + t609) * mrSges(5,3);
t427 = t102 * t235 + t103 * t591;
t428 = t358 - t529;
t536 = Ifges(4,4) * t366;
t333 = t368 * Ifges(4,1) - t536;
t284 = -t369 * Ifges(4,5) + t367 * t333;
t473 = t368 * t284;
t282 = -t369 * Ifges(4,6) + t367 * t599;
t477 = t366 * t282;
t557 = -t368 / 0.2e1;
t560 = -t366 / 0.2e1;
t577 = t170 / 0.2e1;
t588 = (t131 * t179 + t144 * t175 - t235 * t76 + t591 * t79 + t427) * t583 - t369 * t428 / 0.4e1 + t473 / 0.4e1 - t477 / 0.4e1 + t144 * t201 / 0.2e1 + t179 * t577 + t171 * t550 / 0.2e1 + t375 + (t321 * t560 + t323 * t557) * pkin(7) + t597 * mrSges(6,2) + (-t242 / 0.2e1 + t243 / 0.2e1) * t591 - t594 * t235;
t585 = m(5) / 0.2e1;
t584 = -m(6) / 0.2e1;
t582 = -pkin(4) / 0.2e1;
t580 = m(6) * pkin(3);
t578 = -t155 / 0.2e1;
t202 = -mrSges(5,1) * t393 + mrSges(5,2) * t394;
t576 = t202 / 0.2e1;
t574 = -t591 / 0.2e1;
t569 = t288 / 0.2e1;
t565 = -t393 / 0.2e1;
t330 = Ifges(4,2) * t368 + t536;
t563 = -t330 / 0.2e1;
t562 = t333 / 0.4e1;
t559 = t366 / 0.2e1;
t556 = t368 / 0.2e1;
t555 = t369 / 0.2e1;
t553 = m(6) * t102;
t552 = m(6) * t591;
t364 = t369 * pkin(6);
t545 = t89 * mrSges(5,2);
t544 = t89 * mrSges(6,3);
t543 = t90 * mrSges(5,1);
t542 = t90 * mrSges(6,1);
t525 = t102 * mrSges(5,1);
t524 = t102 * mrSges(6,1);
t523 = t103 * mrSges(5,2);
t522 = t103 * mrSges(6,3);
t521 = t591 * mrSges(5,1);
t520 = t591 * mrSges(6,1);
t519 = t235 * mrSges(5,2);
t518 = t235 * mrSges(6,3);
t327 = pkin(3) * t475 + t364;
t132 = pkin(4) * t290 - qJ(5) * t291 + t327;
t154 = Ifges(6,5) * t291 + Ifges(6,6) * t367 + Ifges(6,3) * t290;
t156 = Ifges(5,4) * t291 - Ifges(5,2) * t290 + Ifges(5,6) * t367;
t158 = Ifges(6,1) * t291 + Ifges(6,4) * t367 + Ifges(6,5) * t290;
t160 = Ifges(5,1) * t291 - Ifges(5,4) * t290 + Ifges(5,5) * t367;
t172 = mrSges(6,1) * t290 - mrSges(6,3) * t291;
t173 = mrSges(5,1) * t290 + mrSges(5,2) * t291;
t265 = t314 - t461;
t266 = t366 * t328 + t460;
t283 = Ifges(4,6) * t367 + t369 * t599;
t285 = Ifges(4,5) * t367 + t333 * t369;
t329 = mrSges(4,1) * t366 + mrSges(4,2) * t368;
t306 = t369 * t329;
t322 = -t367 * mrSges(4,2) - mrSges(4,3) * t475;
t324 = t367 * mrSges(4,1) - mrSges(4,3) * t472;
t452 = -Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1;
t5 = (-pkin(1) * mrSges(3,2) - t477 / 0.2e1 + t473 / 0.2e1 + t452 * t291 - t628 + t601 * t369) * t369 + (t157 / 0.2e1 + t159 / 0.2e1) * t291 + m(4) * (t265 * t267 + t266 * t268) + (t154 / 0.2e1 - t156 / 0.2e1) * t289 + (-t158 / 0.2e1 - t160 / 0.2e1) * t288 + m(6) * (t131 * t132 + t76 * t84 + t79 * t85) + m(5) * (t326 * t327 + t89 * t95 + t90 * t96) + (t153 / 0.2e1 + t578) * t290 + t326 * t173 + t327 * t171 + t268 * t321 + t266 * t322 + t267 * t323 + t265 * t324 + t96 * t240 + t90 * t241 + t95 * t242 + t85 * t243 + t89 * t244 + t79 * t245 + t76 * t238 + t84 * t239 + t132 * t170 + t131 * t172 + (-pkin(1) * mrSges(3,1) + pkin(6) * t306 + t283 * t560 + t285 * t556 - t601 * t367 + t451 * t289 + t452 * t288 + (Ifges(3,1) - Ifges(3,2) - Ifges(4,3) + (m(4) * pkin(6) + t329) * pkin(6) - t629) * t369) * t367;
t507 = t5 * qJD(1);
t372 = -t131 * t412 + t288 * t578 - t326 * t413 + t418 * t555 - t89 * t514 + t79 * t515 - t90 * t516 - t76 * t517 + t626 * t569 + (-t407 / 0.2e1 + t625 / 0.2e1) * t289;
t408 = Ifges(4,5) * t366 + Ifges(4,6) * t368;
t436 = t330 * t560;
t6 = -t171 * t462 + t266 * t323 + t372 + (-t321 - t456) * t265 + t592 * t144 + (-m(5) * t90 + t590) * t103 + (m(5) * t89 + t589) * t102 + (-t369 * t408 / 0.2e1 + t284 * t559 + t282 * t556 - m(5) * t326 * t549 + t266 * t368 * mrSges(4,3) + (pkin(6) * t595 - t332 * t557 + t436) * t367) * t367;
t506 = t6 * qJD(1);
t7 = t592 * t405 + t589 * t90 + t590 * t89 + t372;
t505 = t7 * qJD(1);
t25 = m(6) * (t131 * t288 - t369 * t76) - t369 * t239 + t288 * t170;
t502 = qJD(1) * t25;
t483 = t326 * t366;
t482 = t345 * t394;
t481 = t349 * t393;
t463 = mrSges(5,3) * t551;
t459 = -pkin(7) * mrSges(4,3) / 0.2e1;
t457 = t345 * t517;
t455 = -t551 / 0.2e1;
t453 = pkin(6) * t329 / 0.2e1;
t450 = -t512 / 0.2e1;
t448 = t89 * t565;
t433 = t393 * t555;
t426 = t288 * t463;
t423 = mrSges(5,3) * t458;
t419 = 0.2e1 * t304;
t415 = t289 * t423;
t10 = -pkin(2) * t329 + (t599 / 0.2e1 + t332 / 0.2e1) * t368 + (t563 + t333 / 0.2e1 + pkin(3) * t202) * t366 + m(5) * t350 * t550 + t624 + t611 * t179;
t383 = -t235 * t90 - t591 * t89 + t427;
t2 = t588 + t462 * t576 + t474 * t562 + t474 * t563 + ((t350 * t474 + t483) * pkin(3) + t383) * t585 + t595 * t548 / 0.2e1 + (t448 + t597) * mrSges(5,3) + (t453 + (t366 ^ 2 + t368 ^ 2) * t459) * t367 - t610 - (0.2e1 * t332 + t599) * t476 / 0.4e1;
t401 = t2 * qJD(1) + t10 * qJD(2);
t11 = t611 * t198 + t624;
t373 = (t198 * t131 + t175 * t405 + (t89 + t79) * t591) * t584 - t198 * t170 / 0.2e1 + t242 * t575 + t243 * t574 - t405 * t201 / 0.2e1 + ((-t76 + t90) * t584 + t594) * t235 + (t448 + t447) * mrSges(6,2);
t377 = (-pkin(4) * t85 + qJ(5) * t84) * t583 + t245 * t582 + qJ(5) * t238 / 0.2e1 + t382;
t3 = t587 + t377 + t373 + t76 * t510 / 0.2e1 + t79 * t450 - (mrSges(6,2) + mrSges(5,3)) * t609;
t400 = -t3 * qJD(1) + t11 * qJD(2);
t379 = (-t131 * t394 + t175 * t288 - t369 * t591) * t583 + t201 * t569 - t394 * t577;
t397 = t85 * t584 + t508 / 0.2e1;
t21 = (-t433 - t291 / 0.2e1) * mrSges(6,2) + t379 + t397;
t41 = t611 * t394;
t399 = -qJD(1) * t21 + qJD(2) * t41;
t395 = m(6) * (-pkin(4) * t591 - qJ(5) * t235);
t378 = ((-t345 + t551) * t235 + (t349 + t458) * t591) * t583;
t13 = -t395 / 0.2e1 + (-(-qJ(5) / 0.2e1 + t561 + t455) * t394 - (t582 - t349 / 0.2e1 - t458 / 0.2e1) * t393) * mrSges(6,2) + t378 + t603 * (t575 + t574);
t381 = -t603 * t551 + (mrSges(6,3) - mrSges(5,2)) * t458;
t143 = -(t345 * t554 + t349 * t365) * t580 - t381;
t371 = (t345 * t89 + t349 * t90 + (t365 * t79 + t554 * t76) * pkin(3)) * t583 - t545 / 0.2e1 + t544 / 0.2e1 - t543 / 0.2e1 - t542 / 0.2e1 + t457 / 0.2e1 + t349 * t259 + t242 * t455 + t243 * t454 + t426 / 0.2e1 + t415 / 0.2e1 + t600 * t420;
t376 = (-pkin(4) * t102 + qJ(5) * t103) * t584 + t525 / 0.2e1 + t524 / 0.2e1 + t523 / 0.2e1 - t522 / 0.2e1 + pkin(4) * t259 - mrSges(6,2) * t503 / 0.2e1;
t9 = t371 + t376;
t386 = t9 * qJD(1) + t13 * qJD(2) - t143 * qJD(3);
t380 = -t357 + ((-qJ(5) - t345) * t369 + t90) * t583;
t27 = -t553 / 0.2e1 + t380;
t325 = m(6) * t345 + mrSges(6,3);
t385 = -qJD(1) * t27 - qJD(3) * t325;
t29 = -t357 + 0.2e1 * (t446 / 0.4e1 - t471 / 0.2e1 + t480 / 0.4e1 - t90 / 0.4e1) * m(6);
t384 = qJD(1) * t29 + qJD(3) * t356 + t608;
t302 = mrSges(6,3) + (qJ(5) + 0.2e1 * t454) * m(6);
t104 = t419 + t552;
t44 = t552 / 0.2e1 + m(6) * t575 + t419;
t28 = t76 * t623 - t357 + t605;
t26 = t605 + t553 / 0.2e1 + t380;
t20 = -mrSges(6,2) * t433 + t513 / 0.2e1 + t379 - t397;
t12 = t519 / 0.2e1 - t521 / 0.2e1 - t520 / 0.2e1 - t518 / 0.2e1 + t395 / 0.2e1 + pkin(4) * t450 + qJ(5) * t449 + (-mrSges(5,1) / 0.2e1 + t622) * t591 + (mrSges(5,2) / 0.2e1 - mrSges(6,3) / 0.2e1) * t235 + (t481 / 0.2e1 - t482 / 0.2e1 + (t365 * t564 + t554 * t566) * pkin(3)) * mrSges(6,2) + t378 + t417;
t8 = t371 - t376 + t418;
t4 = t377 - t373 + mrSges(5,3) * t504 / 0.2e1 + t375;
t1 = (t453 + (pkin(2) * mrSges(4,2) / 0.2e1 - t332 / 0.4e1 - t599 / 0.4e1 + (-Ifges(4,1) / 0.4e1 + t459) * t366) * t366 + (-pkin(2) * mrSges(4,1) / 0.2e1 - t536 / 0.2e1 + t562 - t330 / 0.4e1 + (-Ifges(4,2) / 0.4e1 + t459) * t368 + (t350 * t585 + t576) * pkin(3)) * t368) * t367 + (t483 * pkin(3) + t383) * t585 + (t443 - (-t103 / 0.2e1 + t89 / 0.2e1) * t393) * mrSges(5,3) + t588 + t610;
t14 = [qJD(2) * t5 - qJD(3) * t6 - qJD(4) * t7 + qJD(5) * t25, t1 * qJD(3) + t4 * qJD(4) + t20 * qJD(5) + t507 + (t627 * t617 + t285 * t559 + t154 * t565 + t156 * t566 + t322 * t363 + mrSges(3,2) * t361 + t283 * t556 + (t160 + t158) * t564 + t85 * t510 + t84 * t512 + (t393 * t96 - t394 * t95) * mrSges(5,3) + 0.2e1 * (-t235 * t95 + t327 * t350 + t591 * t96) * t585 + t591 * t241 + t591 * t238 + (Ifges(3,5) + t436 + t332 * t556 + (-mrSges(3,1) + t595) * pkin(6)) * t369 - Ifges(3,6) * t367 + t350 * t173 + t327 * t202 - pkin(2) * t306 - t235 * t244 + t235 * t245 + t132 * t201 + t175 * t172 - t366 * pkin(7) * t324 + t593 * t290 + (t630 * t393 + t631 * t394 + t408) * t558 + t596 * mrSges(4,3) + m(4) * (-pkin(2) * t364 + t596 * pkin(7)) + (t132 * t175 + t235 * t85 + t591 * t84) * t623) * qJD(2), -t506 + t1 * qJD(2) + (t426 + t415 - Ifges(4,5) * t476 - Ifges(4,6) * t474 - t349 * t515 + m(6) * (t102 * t349 + t103 * t345) + t457 + t522 - t525 - t524 - t523 + (-t102 * t554 + t103 * t365) * t581 - t266 * mrSges(4,1) - t265 * mrSges(4,2) + t418) * qJD(3) + t8 * qJD(4) + t26 * qJD(5), -t505 + t4 * qJD(2) + t8 * qJD(3) + (m(6) * (-pkin(4) * t90 + qJ(5) * t89) + t544 - t543 - t542 - t545 + t406 * mrSges(6,2) + t418) * qJD(4) + t28 * qJD(5), qJD(2) * t20 + qJD(3) * t26 + qJD(4) * t28 + t502; qJD(3) * t2 - qJD(4) * t3 + qJD(5) * t21 - t507, qJD(3) * t10 + qJD(4) * t11 - qJD(5) * t41, (-t394 * t463 - t393 * t423 + m(6) * (-t235 * t345 + t349 * t591) + t519 - t521 - t520 - t518 + (-t235 * t365 - t554 * t591) * t581 + t417 + t428 + t595 * pkin(7) + (t481 - t482) * mrSges(6,2)) * qJD(3) + t12 * qJD(4) + t44 * qJD(5) + t401, t12 * qJD(3) + (t404 * mrSges(6,2) + t417 + (-m(6) * pkin(4) - t603) * t591 + (mrSges(5,2) - t356) * t235) * qJD(4) + t104 * qJD(5) + t400, qJD(3) * t44 + qJD(4) * t104 - t399; -qJD(2) * t2 + qJD(4) * t9 + qJD(5) * t27 + t506, qJD(4) * t13 - t401, -qJD(4) * t143 + qJD(5) * t325, ((-pkin(4) * t365 + qJ(5) * t554) * t580 + t381) * qJD(4) + t302 * qJD(5) + t386, qJD(4) * t302 - t385; qJD(2) * t3 - qJD(3) * t9 + qJD(5) * t29 + t505, -qJD(3) * t13 - t400, -t386 + t607, t607, t384; -qJD(2) * t21 - qJD(3) * t27 - qJD(4) * t29 - t502, t399, t385 - t608, -t384, 0;];
Cq = t14;
