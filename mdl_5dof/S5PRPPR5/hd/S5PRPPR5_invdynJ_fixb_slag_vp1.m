% Calculate vector of inverse dynamics joint torques for
% S5PRPPR5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:03
% EndTime: 2019-12-31 17:38:36
% DurationCPUTime: 31.98s
% Computational Cost: add. (8591->659), mult. (24207->944), div. (0->0), fcn. (25686->8), ass. (0->279)
t545 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t323 = sin(qJ(2));
t325 = cos(qJ(2));
t478 = sin(pkin(8));
t479 = cos(pkin(8));
t273 = t323 * t478 + t325 * t479;
t320 = sin(pkin(7));
t216 = t273 * t320;
t274 = t323 * t479 - t325 * t478;
t217 = t274 * t320;
t321 = cos(pkin(7));
t218 = t273 * t321;
t219 = t274 * t321;
t544 = (-Icges(5,5) * t216 - Icges(5,6) * t217 - t321 * t545) * t321 + (Icges(5,5) * t218 + Icges(5,6) * t219 - t320 * t545) * t320;
t250 = t273 * qJD(2);
t177 = t320 * t250;
t518 = qJD(2) * t274;
t178 = t320 * t518;
t179 = t321 * t250;
t180 = t321 * t518;
t543 = (Icges(5,5) * t178 - Icges(5,6) * t177) * t321 + (-Icges(5,5) * t180 + Icges(5,6) * t179) * t320;
t541 = (-Icges(3,6) + Icges(4,6)) * t325 + (-Icges(4,4) - Icges(3,5)) * t323;
t282 = rSges(3,1) * t323 + rSges(3,2) * t325;
t318 = t320 ^ 2;
t319 = t321 ^ 2;
t507 = t318 + t319;
t540 = t282 * t507;
t316 = t325 * pkin(3);
t508 = pkin(2) * t325 + qJ(3) * t323;
t431 = t316 + t508;
t539 = -pkin(4) * t273 + pkin(6) * t274 - t431;
t529 = 0.2e1 * qJD(2);
t528 = 2 * qJDD(2);
t527 = t541 * t320;
t526 = t541 * t321;
t311 = qJD(4) * t321;
t446 = qJD(2) * t323;
t434 = pkin(3) * t446;
t271 = -t320 * t434 + t311;
t444 = qJD(4) * t320;
t272 = -t321 * t434 - t444;
t467 = t320 * t325;
t301 = qJ(3) * t467;
t465 = t321 * t325;
t302 = qJ(3) * t465;
t312 = qJD(3) * t323;
t448 = qJD(2) * t321;
t449 = qJD(2) * t320;
t466 = t321 * t323;
t468 = t320 * t323;
t433 = (-pkin(2) * t468 + t301) * t449 + (-pkin(2) * t466 + t302) * t448 + t312;
t296 = t320 * t312;
t280 = pkin(2) * t323 - qJ(3) * t325;
t365 = qJD(2) * t280;
t173 = -t320 * t365 + t296;
t298 = t321 * t312;
t174 = -t321 * t365 + t298;
t463 = t173 * t320 + t174 * t321;
t524 = t271 * t320 + t272 * t321 - t433 + t463;
t439 = qJD(2) * qJD(3);
t523 = qJDD(3) * t323 + t325 * t439;
t322 = sin(qJ(5));
t324 = cos(qJ(5));
t522 = -rSges(6,1) * t324 + rSges(6,2) * t322;
t516 = (-rSges(5,1) * t178 + rSges(5,2) * t177) * t320 + (-rSges(5,1) * t180 + rSges(5,2) * t179) * t321;
t171 = -qJD(5) * t219 + t449;
t495 = t171 / 0.2e1;
t172 = -qJD(5) * t217 - t448;
t493 = t172 / 0.2e1;
t443 = qJD(5) * t273;
t425 = t443 / 0.2e1;
t509 = rSges(4,1) * t325 + rSges(4,3) * t323;
t510 = t508 + t509;
t506 = g(1) * t321 + g(2) * t320;
t505 = t506 * t323;
t471 = Icges(6,4) * t322;
t398 = Icges(6,1) * t324 - t471;
t115 = Icges(6,5) * t273 + t274 * t398;
t165 = -t216 * t322 + t321 * t324;
t157 = Icges(6,4) * t165;
t167 = -t218 * t322 - t320 * t324;
t158 = Icges(6,4) * t167;
t377 = -t218 * t324 + t320 * t322;
t378 = -t216 * t324 - t321 * t322;
t73 = -Icges(6,1) * t378 - Icges(6,5) * t217 + t157;
t74 = -Icges(6,1) * t377 - Icges(6,5) * t219 + t158;
t332 = t171 * (Icges(6,2) * t377 + t158 + t74) + t172 * (Icges(6,2) * t378 + t157 + t73) + t443 * (t115 + (-Icges(6,2) * t324 - t471) * t274);
t326 = qJD(2) ^ 2;
t473 = Icges(6,4) * t378;
t71 = Icges(6,2) * t165 - Icges(6,6) * t217 - t473;
t404 = -t322 * t71 + t324 * t73;
t69 = -Icges(6,5) * t378 + Icges(6,6) * t165 - Icges(6,3) * t217;
t28 = t273 * t69 + t274 * t404;
t472 = Icges(6,4) * t377;
t72 = Icges(6,2) * t167 - Icges(6,6) * t219 - t472;
t403 = -t322 * t72 + t324 * t74;
t70 = -Icges(6,5) * t377 + Icges(6,6) * t167 - Icges(6,3) * t219;
t29 = t273 * t70 + t274 * t403;
t389 = Icges(6,5) * t324 - Icges(6,6) * t322;
t111 = Icges(6,3) * t273 + t274 * t389;
t470 = Icges(6,4) * t324;
t393 = -Icges(6,2) * t322 + t470;
t113 = Icges(6,6) * t273 + t274 * t393;
t388 = -t113 * t322 + t115 * t324;
t37 = t111 * t273 + t274 * t388;
t502 = t28 * t493 + t29 * t495 + t37 * t425;
t501 = -m(5) - m(6);
t438 = qJDD(2) * t320;
t116 = -qJD(5) * t179 - qJDD(5) * t219 + t438;
t499 = t116 / 0.2e1;
t437 = qJDD(2) * t321;
t117 = -qJD(5) * t177 - qJDD(5) * t217 - t437;
t498 = t117 / 0.2e1;
t156 = -qJD(5) * t518 + qJDD(5) * t273;
t497 = t156 / 0.2e1;
t496 = -t171 / 0.2e1;
t494 = -t172 / 0.2e1;
t492 = -t507 * t446 / 0.2e1;
t489 = -t325 / 0.2e1;
t487 = pkin(3) * t323;
t475 = Icges(5,4) * t217;
t474 = Icges(5,4) * t219;
t462 = -rSges(5,1) * t217 + rSges(5,2) * t216;
t461 = -rSges(5,1) * t219 + rSges(5,2) * t218;
t264 = t508 * t320;
t266 = t508 * t321;
t456 = t264 * t320 + t266 * t321;
t445 = qJD(3) * t325;
t245 = qJD(2) * t508 - t445;
t455 = -qJD(2) * t509 - t245;
t281 = rSges(4,1) * t323 - rSges(4,3) * t325;
t454 = -t280 - t281;
t453 = t523 * t320;
t452 = t523 * t321;
t451 = t296 + t311;
t442 = qJD(5) * t274;
t440 = -m(4) + t501;
t435 = t326 * t316;
t432 = qJDD(4) * t321 + t453;
t426 = -t443 / 0.2e1;
t424 = -t280 - t487;
t422 = -pkin(4) * t217 - pkin(6) * t216;
t421 = -pkin(4) * t219 - pkin(6) * t218;
t416 = qJD(2) * t454;
t415 = t298 - t444;
t275 = pkin(3) * t467 + qJ(4) * t321;
t276 = pkin(3) * t465 - qJ(4) * t320;
t413 = t275 * t320 + t276 * t321 + t456;
t412 = t507 * t487;
t410 = -rSges(5,1) * t274 + rSges(5,2) * t273 + t424;
t408 = -pkin(4) * t274 - pkin(6) * t273 + t424;
t407 = -qJD(2) * t316 - t245;
t285 = rSges(3,1) * t325 - rSges(3,2) * t323;
t402 = -qJDD(4) * t320 + t452;
t136 = rSges(5,1) * t216 + rSges(5,2) * t217 + rSges(5,3) * t321;
t137 = rSges(5,1) * t218 + rSges(5,2) * t219 - rSges(5,3) * t320;
t387 = t136 * t320 + t137 * t321;
t386 = (-Icges(5,5) * t217 + Icges(5,6) * t216) * t321 - (-Icges(5,5) * t219 + Icges(5,6) * t218) * t320;
t146 = pkin(4) * t216 - pkin(6) * t217;
t147 = pkin(4) * t218 - pkin(6) * t219;
t385 = t146 * t320 + t147 * t321;
t209 = -rSges(4,2) * t321 + t320 * t509;
t211 = rSges(4,2) * t320 + t321 * t509;
t380 = t209 * t320 + t211 * t321;
t379 = t507 * t285;
t376 = qJD(2) * t540;
t119 = rSges(6,3) * t273 - t274 * t522;
t375 = -t119 + t408;
t95 = -rSges(6,3) * t216 + t217 * t522;
t96 = -rSges(6,3) * t218 + t219 * t522;
t118 = -rSges(6,3) * t274 - t273 * t522;
t153 = rSges(5,1) * t250 + rSges(5,2) * t518;
t374 = -t153 + t407;
t373 = t264 * t449 + t266 * t448 + qJD(1) - t445;
t372 = qJD(2) * t410;
t371 = qJD(2) * t408;
t370 = -t250 * t322 - t324 * t442;
t369 = t250 * t324 - t322 * t442;
t155 = pkin(4) * t250 - pkin(6) * t518;
t66 = rSges(6,1) * t369 + rSges(6,2) * t370 - rSges(6,3) * t518;
t368 = -t155 + t407 - t66;
t366 = qJD(2) * t281;
t346 = t275 * t449 + t276 * t448 + t373;
t181 = Icges(5,4) * t216;
t130 = Icges(5,2) * t217 + Icges(5,6) * t321 + t181;
t182 = Icges(5,4) * t218;
t131 = Icges(5,2) * t219 - Icges(5,6) * t320 + t182;
t345 = (Icges(5,1) * t219 - t131 - t182) * t320 - (Icges(5,1) * t217 - t130 - t181) * t321;
t132 = Icges(5,1) * t216 + Icges(5,5) * t321 + t475;
t133 = Icges(5,1) * t218 - Icges(5,5) * t320 + t474;
t344 = (-Icges(5,2) * t218 + t133 + t474) * t320 - (-Icges(5,2) * t216 + t132 + t475) * t321;
t341 = (-Icges(6,5) * t322 - Icges(6,6) * t324) * t274 * t443 + t171 * (Icges(6,5) * t167 + Icges(6,6) * t377) + t172 * (Icges(6,5) * t165 + Icges(6,6) * t378);
t340 = qJD(2) * t455 + qJDD(2) * t454;
t339 = -qJDD(3) * t325 + t173 * t449 + t174 * t448 + t264 * t438 + t266 * t437 + t323 * t439 + qJDD(1);
t337 = (-pkin(2) - pkin(3)) * t505;
t331 = t271 * t449 + t272 * t448 + t275 * t438 + t276 * t437 + t339;
t330 = (Icges(6,1) * t167 + t472 - t72) * t171 + (Icges(6,1) * t165 + t473 - t71) * t172 + (-t113 + (-Icges(6,1) * t322 - t470) * t274) * t443;
t329 = -t435 + (-t153 - t245) * qJD(2) + t410 * qJDD(2);
t328 = -t435 + (-t155 - t245) * qJD(2) + t408 * qJDD(2);
t327 = (-Icges(6,3) * t218 - t219 * t389 + t403) * t171 + (-Icges(6,3) * t216 - t217 * t389 + t404) * t172 + (-Icges(6,3) * t274 + t273 * t389 + t388) * t443;
t306 = rSges(4,3) * t465;
t305 = rSges(4,3) * t467;
t299 = t321 * t445;
t297 = t320 * t445;
t265 = t282 * t321;
t263 = t282 * t320;
t236 = t321 * t366;
t234 = t320 * t366;
t154 = (-rSges(6,1) * t322 - rSges(6,2) * t324) * t274;
t149 = t321 * t416 + t298;
t148 = t320 * t416 + t296;
t135 = -pkin(4) * t180 - pkin(6) * t179;
t134 = -pkin(4) * t178 - pkin(6) * t177;
t125 = -Icges(5,1) * t180 + Icges(5,4) * t179;
t124 = -Icges(5,1) * t178 + Icges(5,4) * t177;
t123 = -Icges(5,4) * t180 + Icges(5,2) * t179;
t122 = -Icges(5,4) * t178 + Icges(5,2) * t177;
t114 = -Icges(6,5) * t274 + t273 * t398;
t112 = -Icges(6,6) * t274 + t273 * t393;
t108 = qJD(5) * t167 - t180 * t324;
t107 = qJD(5) * t377 + t180 * t322;
t106 = qJD(5) * t165 - t178 * t324;
t105 = qJD(5) * t378 + t178 * t322;
t104 = rSges(6,1) * t167 + rSges(6,2) * t377;
t103 = rSges(6,1) * t165 + rSges(6,2) * t378;
t94 = -Icges(6,5) * t218 - t219 * t398;
t93 = -Icges(6,5) * t216 - t217 * t398;
t92 = -Icges(6,6) * t218 - t219 * t393;
t91 = -Icges(6,6) * t216 - t217 * t393;
t88 = t321 * t372 + t415;
t87 = t320 * t372 + t451;
t76 = -rSges(6,1) * t377 + rSges(6,2) * t167 - rSges(6,3) * t219;
t75 = -rSges(6,1) * t378 + rSges(6,2) * t165 - rSges(6,3) * t217;
t68 = t321 * t340 + t452;
t67 = t320 * t340 + t453;
t65 = qJD(2) * t380 + t373;
t64 = Icges(6,1) * t369 + Icges(6,4) * t370 - Icges(6,5) * t518;
t63 = Icges(6,4) * t369 + Icges(6,2) * t370 - Icges(6,6) * t518;
t62 = Icges(6,5) * t369 + Icges(6,6) * t370 - Icges(6,3) * t518;
t61 = -qJD(2) * t376 + qJDD(2) * t379 + qJDD(1);
t52 = rSges(6,1) * t108 + rSges(6,2) * t107 - rSges(6,3) * t179;
t51 = rSges(6,1) * t106 + rSges(6,2) * t105 - rSges(6,3) * t177;
t50 = Icges(6,1) * t108 + Icges(6,4) * t107 - Icges(6,5) * t179;
t49 = Icges(6,1) * t106 + Icges(6,4) * t105 - Icges(6,5) * t177;
t48 = Icges(6,4) * t108 + Icges(6,2) * t107 - Icges(6,6) * t179;
t47 = Icges(6,4) * t106 + Icges(6,2) * t105 - Icges(6,6) * t177;
t46 = Icges(6,5) * t108 + Icges(6,6) * t107 - Icges(6,3) * t179;
t45 = Icges(6,5) * t106 + Icges(6,6) * t105 - Icges(6,3) * t177;
t44 = t321 * t329 + t402;
t43 = t320 * t329 + t432;
t38 = qJD(2) * t387 + t346;
t36 = t380 * qJDD(2) + (-t234 * t320 - t236 * t321) * qJD(2) + t339;
t35 = -t111 * t219 + t113 * t167 - t115 * t377;
t34 = -t111 * t217 + t113 * t165 - t115 * t378;
t33 = t119 * t172 + t321 * t371 - t443 * t75 + t415;
t32 = -t119 * t171 + t320 * t371 + t443 * t76 + t451;
t23 = t167 * t72 - t219 * t70 - t377 * t74;
t22 = t167 * t71 - t219 * t69 - t377 * t73;
t21 = t165 * t72 - t217 * t70 - t378 * t74;
t20 = t165 * t71 - t217 * t69 - t378 * t73;
t19 = qJD(2) * t385 + t171 * t75 - t172 * t76 + t346;
t18 = qJD(2) * t516 + qJDD(2) * t387 + t331;
t17 = t117 * t119 - t156 * t75 + t172 * t66 + t321 * t328 - t443 * t51 + t402;
t16 = -t116 * t119 + t156 * t76 - t171 * t66 + t320 * t328 + t443 * t52 + t432;
t15 = t107 * t113 + t108 * t115 - t111 * t179 + t167 * t63 - t219 * t62 - t377 * t64;
t14 = t105 * t113 + t106 * t115 - t111 * t177 + t165 * t63 - t217 * t62 - t378 * t64;
t13 = -t111 * t518 + t273 * t62 + t388 * t250 + (-t322 * t63 + t324 * t64 + (-t113 * t324 - t115 * t322) * qJD(5)) * t274;
t11 = t107 * t72 + t108 * t74 + t167 * t48 - t179 * t70 - t219 * t46 - t377 * t50;
t10 = t107 * t71 + t108 * t73 + t167 * t47 - t179 * t69 - t219 * t45 - t377 * t49;
t9 = t105 * t72 + t106 * t74 + t165 * t48 - t177 * t70 - t217 * t46 - t378 * t50;
t8 = t105 * t71 + t106 * t73 + t165 * t47 - t177 * t69 - t217 * t45 - t378 * t49;
t7 = t171 * t23 + t172 * t22 + t35 * t443;
t6 = t171 * t21 + t172 * t20 + t34 * t443;
t5 = -t518 * t70 + t273 * t46 + t403 * t250 + (-t322 * t48 + t324 * t50 + (-t322 * t74 - t324 * t72) * qJD(5)) * t274;
t4 = -t518 * t69 + t273 * t45 + t404 * t250 + (-t322 * t47 + t324 * t49 + (-t322 * t73 - t324 * t71) * qJD(5)) * t274;
t3 = (t134 * t320 + t135 * t321) * qJD(2) + t385 * qJDD(2) + t116 * t75 - t117 * t76 + t171 * t51 - t172 * t52 + t331;
t2 = t10 * t172 + t11 * t171 + t116 * t23 + t117 * t22 + t15 * t443 + t156 * t35;
t1 = t116 * t21 + t117 * t20 + t14 * t443 + t156 * t34 + t171 * t9 + t172 * t8;
t12 = [m(2) * qJDD(1) + (-m(2) - m(3) + t440) * g(3) + m(3) * t61 + m(4) * t36 + m(5) * t18 + m(6) * t3; ((t165 * t92 - t216 * t70 - t378 * t94) * t171 + (t165 * t91 - t216 * t69 - t378 * t93) * t172 + (-t21 * t218 - t20 * t216 + (-t111 * t216 + t112 * t165 - t114 * t378) * t273 - t34 * t274) * qJD(5) - t327 * t217) * t494 + ((t167 * t92 - t218 * t70 - t377 * t94) * t171 + (t167 * t91 - t218 * t69 - t377 * t93) * t172 + (-t23 * t218 - t22 * t216 + (-t111 * t218 + t112 * t167 - t114 * t377) * t273 - t35 * t274) * qJD(5) - t327 * t219) * t496 + ((-t216 * t28 - t218 * t29) * qJD(5) + ((-t322 * t92 + t324 * t94 - t70) * t171 + (-t322 * t91 + t324 * t93 - t69) * t172 - t37 * qJD(5)) * t274 + ((-t112 * t322 + t114 * t324 - t111) * t442 + t327) * t273) * t426 + (t320 * t5 - t321 * t4) * t425 + (t320 * t9 - t321 * t8) * t493 + (-t10 * t321 + t11 * t320) * t495 + (-t28 * t321 + t29 * t320) * t497 + (-t20 * t321 + t21 * t320) * t498 + (-t22 * t321 + t23 * t320) * t499 + t442 * t502 - (-t216 * t6 - t218 * t7) * qJD(5) / 0.2e1 - (t526 * t318 * qJD(2) - t320 * t448 * t527) * t449 / 0.2e1 + (t527 * qJD(2) * t319 - t321 * t449 * t526) * t448 / 0.2e1 + (-t320 * (-t218 * t345 - t219 * t344 + t320 * t386) / 0.2e1 + t321 * (-t216 * t345 - t217 * t344 - t321 * t386) / 0.2e1) * t326 + (-g(1) * (t302 + t421 + t96) - g(2) * (t301 + t422 + t95) - g(3) * (t118 - t539) - t337 - t33 * (t118 * t172 + t299) - t32 * (-t118 * t171 + t297) - (t33 * (-t119 * t216 - t273 * t95 + t274 * t75) + t32 * (t119 * t218 + t273 * t96 - t274 * t76)) * qJD(5) + t3 * t413 + (t17 * t375 + t33 * t368 + t3 * (t147 + t76)) * t321 + (t16 * t375 + t32 * t368 + t3 * (t146 + t75)) * t320 - (t32 * t320 + t321 * t33) * qJD(2) * t539 + (-t171 * t95 + t172 * t96 - (t216 * t76 - t218 * t75) * qJD(5) - (t320 * t422 + t321 * t421 - t412) * qJD(2) + (t135 + t52) * t321 + (t134 + t51) * t320 + t524) * t19) * m(6) + (-g(1) * (t302 + t461) - g(2) * (t301 + t462) - t337 - t88 * t299 - t87 * t297 + t18 * t413 + (t18 * t137 + t374 * t88 + t410 * t44) * t321 + (t18 * t136 + t374 * t87 + t410 * t43) * t320 + (-(t320 * t462 + t321 * t461 - t412) * qJD(2) + t516 + t524) * t38 + (-g(3) + (t320 * t87 + t321 * t88) * qJD(2)) * (rSges(5,1) * t273 + rSges(5,2) * t274 + t431)) * m(5) + (t36 * t456 + t65 * t463 + (t149 * t455 + t36 * t211 - t65 * t236 + t454 * t68) * t321 + (t148 * t455 + t36 * t209 - t65 * t234 + t454 * t67) * t320 - g(1) * (t302 + t306) - g(2) * (t301 + t305) - g(3) * t510 - (-rSges(4,1) - pkin(2)) * t505 - t149 * t299 - t148 * t297 - t65 * t433 - ((-t149 * t510 + t65 * (-rSges(4,1) * t466 + t306)) * t321 + (-t148 * t510 + t65 * (-rSges(4,1) * t468 + t305)) * t320) * qJD(2)) * m(4) + (g(1) * t265 + g(2) * t263 - g(3) * t285 + t61 * t379 + (qJDD(2) * t282 + t285 * t326) * t540 + (-(-t263 * t320 - t265 * t321) * qJD(2) - t376) * (qJD(2) * t379 + qJD(1))) * m(3) + (t2 + ((-t122 * t219 - t124 * t218 - t130 * t179 + t132 * t180) * t321 + (t123 * t219 + t125 * t218 + t131 * t179 - t133 * t180 - t543) * t320) * t529 + ((-t130 * t219 - t132 * t218) * t321 + (t131 * t219 + t133 * t218 - t544) * t320) * t528) * t320 / 0.2e1 - (t1 + ((t123 * t217 + t125 * t216 + t131 * t177 - t133 * t178) * t320 + (-t122 * t217 - t124 * t216 - t130 * t177 + t132 * t178 + t543) * t321) * t529 + (t320 * (t131 * t217 + t133 * t216) + (-t130 * t217 - t132 * t216 + t544) * t321) * t528) * t321 / 0.2e1; -t440 * g(3) * t325 + 0.2e1 * (t19 * t492 + t3 * t489) * m(6) + 0.2e1 * (t18 * t489 + t38 * t492) * m(5) + 0.2e1 * (t36 * t489 + t492 * t65) * m(4) + (t440 * t506 + m(4) * (qJD(2) * t65 + t320 * t67 + t321 * t68) + m(5) * (qJD(2) * t38 + t320 * t43 + t321 * t44) + m(6) * (qJD(2) * t19 + t16 * t320 + t17 * t321)) * t323; t501 * (-g(1) * t320 + g(2) * t321) + m(5) * (-t320 * t44 + t321 * t43) + m(6) * (t16 * t321 - t17 * t320); -t179 * t7 / 0.2e1 - t219 * t2 / 0.2e1 + (-t217 * t22 - t219 * t23 + t273 * t35) * t499 + (-t10 * t217 - t11 * t219 + t15 * t273 - t177 * t22 - t179 * t23 - t35 * t518) * t495 - t177 * t6 / 0.2e1 - t217 * t1 / 0.2e1 + (-t20 * t217 - t21 * t219 + t273 * t34) * t498 + (t14 * t273 - t177 * t20 - t179 * t21 - t217 * t8 - t219 * t9 - t34 * t518) * t493 - t518 * t502 + t273 * (t116 * t29 + t117 * t28 + t13 * t443 + t156 * t37 + t171 * t5 + t172 * t4) / 0.2e1 + (-t217 * t28 - t219 * t29 + t273 * t37) * t497 + (t13 * t273 - t177 * t28 - t179 * t29 - t217 * t4 - t219 * t5 - t37 * t518) * t425 + (t167 * t332 - t219 * t341 - t330 * t377) * t496 + (t165 * t332 - t217 * t341 - t330 * t378) * t494 + (t341 * t273 + (-t322 * t332 + t330 * t324) * t274) * t426 + (t17 * (-t119 * t217 - t273 * t75) + t16 * (t119 * t219 + t273 * t76) + t3 * (t217 * t76 - t219 * t75) - g(1) * t104 - g(2) * t103 - g(3) * t154 + (t103 * t443 - t119 * t177 - t154 * t172 - t217 * t66 - t273 * t51 + t518 * t75) * t33 + (-t104 * t443 + t119 * t179 + t154 * t171 + t219 * t66 + t273 * t52 - t518 * t76) * t32 + (-t103 * t171 + t104 * t172 + t177 * t76 - t179 * t75 + t217 * t52 - t219 * t51) * t19) * m(6);];
tau = t12;
