% Calculate time derivative of joint inertia matrix for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR9_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:56:48
% EndTime: 2019-03-09 18:57:11
% DurationCPUTime: 10.86s
% Computational Cost: add. (27033->872), mult. (82706->1277), div. (0->0), fcn. (84918->14), ass. (0->369)
t311 = sin(pkin(6));
t458 = 0.2e1 * t311;
t315 = sin(qJ(6));
t319 = cos(qJ(6));
t316 = sin(qJ(5));
t374 = qJD(6) * t316;
t320 = cos(qJ(5));
t378 = qJD(5) * t320;
t453 = t315 * t374 - t319 * t378;
t313 = cos(pkin(7));
t321 = cos(qJ(3));
t394 = t313 * t321;
t298 = pkin(2) * t394;
t310 = sin(pkin(7));
t317 = sin(qJ(3));
t413 = -pkin(10) - qJ(4);
t353 = t413 * t317;
t218 = pkin(3) * t313 + t310 * t353 + t298;
t395 = t313 * t317;
t297 = pkin(2) * t395;
t398 = t310 * t321;
t257 = pkin(10) * t398 + t297;
t233 = qJ(4) * t398 + t257;
t309 = sin(pkin(13));
t312 = cos(pkin(13));
t173 = t309 * t218 + t312 * t233;
t155 = pkin(11) * t313 + t173;
t396 = t312 * t321;
t399 = t310 * t317;
t241 = t309 * t399 - t310 * t396;
t242 = (t309 * t321 + t312 * t317) * t310;
t264 = (-pkin(3) * t321 - pkin(2)) * t310;
t178 = pkin(4) * t241 - pkin(11) * t242 + t264;
t456 = t320 * t155 + t316 * t178;
t457 = qJD(5) * t456;
t421 = t315 / 0.2e1;
t419 = t319 / 0.2e1;
t235 = qJD(3) * t242;
t385 = qJD(3) * t310;
t236 = (-t309 * t317 + t396) * t385;
t184 = Ifges(5,5) * t236 - Ifges(5,6) * t235;
t383 = qJD(3) * t321;
t363 = t310 * t383;
t384 = qJD(3) * t317;
t364 = t310 * t384;
t247 = Ifges(4,5) * t363 - Ifges(4,6) * t364;
t455 = t247 + t184;
t387 = t315 ^ 2 + t319 ^ 2;
t454 = (t398 * t413 - t297) * qJD(3) - qJD(4) * t399;
t373 = qJD(6) * t319;
t328 = -t315 * t378 - t316 * t373;
t303 = -pkin(3) * t312 - pkin(4);
t452 = m(6) * t303 - mrSges(6,1) * t320 + mrSges(6,2) * t316;
t314 = cos(pkin(6));
t322 = cos(qJ(2));
t397 = t311 * t322;
t254 = -t310 * t397 + t314 * t313;
t318 = sin(qJ(2));
t417 = pkin(1) * t314;
t299 = t318 * t417;
t258 = pkin(9) * t397 + t299;
t203 = (t310 * t314 + t313 * t397) * pkin(10) + t258;
t197 = t321 * t203;
t300 = t322 * t417;
t346 = t311 * (-pkin(10) * t313 - pkin(9));
t333 = t318 * t346;
t211 = pkin(2) * t314 + t300 + t333;
t415 = pkin(10) * t310;
t234 = (-pkin(2) * t322 - t318 * t415 - pkin(1)) * t311;
t130 = t211 * t395 + t234 * t399 + t197;
t388 = t321 * t322;
t392 = t317 * t318;
t330 = t313 * t388 - t392;
t205 = t311 * t330 + t314 * t398;
t112 = qJ(4) * t205 + t130;
t129 = -t203 * t317 + t211 * t394 + t234 * t398;
t390 = t318 * t321;
t391 = t317 * t322;
t329 = t313 * t391 + t390;
t206 = t311 * t329 + t314 * t399;
t94 = pkin(3) * t254 - qJ(4) * t206 + t129;
t55 = t312 * t112 + t309 * t94;
t49 = pkin(11) * t254 + t55;
t159 = -t211 * t310 + t313 * t234;
t133 = -pkin(3) * t205 + t159;
t150 = -t312 * t205 + t206 * t309;
t151 = t205 * t309 + t206 * t312;
t80 = pkin(4) * t150 - pkin(11) * t151 + t133;
t412 = t316 * t80 + t320 * t49;
t451 = qJD(5) * t412;
t167 = -t314 * t364 + (-t329 * qJD(3) + (-t313 * t390 - t391) * qJD(2)) * t311;
t168 = t314 * t363 + (t330 * qJD(3) + (-t313 * t392 + t388) * qJD(2)) * t311;
t120 = -t312 * t167 + t168 * t309;
t121 = t167 * t309 + t168 * t312;
t131 = t151 * t316 - t254 * t320;
t386 = qJD(2) * t311;
t365 = t318 * t386;
t347 = t310 * t365;
t78 = -qJD(5) * t131 + t121 * t320 + t316 * t347;
t41 = mrSges(6,1) * t120 - mrSges(6,3) * t78;
t290 = qJD(2) * t300;
t214 = qJD(2) * t333 + t290;
t215 = (t322 * t346 - t299) * qJD(2);
t237 = (pkin(2) * t318 - t322 * t415) * t386;
t84 = -t214 * t317 + t215 * t394 + t237 * t398 + (-t197 + (-t211 * t313 - t234 * t310) * t317) * qJD(3);
t51 = pkin(3) * t347 - qJ(4) * t168 - qJD(4) * t206 + t84;
t362 = t313 * t383;
t83 = -t203 * t384 + t211 * t362 + t321 * t214 + t215 * t395 + t234 * t363 + t237 * t399;
t57 = qJ(4) * t167 + qJD(4) * t205 + t83;
t24 = t309 * t51 + t312 * t57;
t22 = pkin(11) * t347 + t24;
t169 = -t215 * t310 + t313 * t237;
t127 = -pkin(3) * t167 + t169;
t40 = pkin(4) * t120 - pkin(11) * t121 + t127;
t6 = -t22 * t316 + t320 * t40 - t451;
t450 = m(6) * (t6 + t451) + t41;
t212 = t242 * t316 - t320 * t313;
t170 = -qJD(5) * t212 + t236 * t320;
t143 = mrSges(6,1) * t235 - mrSges(6,3) * t170;
t289 = pkin(2) * t362;
t200 = t289 + (qJD(3) * t353 + qJD(4) * t321) * t310;
t147 = t312 * t200 + t309 * t454;
t348 = pkin(3) * t364;
t175 = pkin(4) * t235 - pkin(11) * t236 + t348;
t65 = -t147 * t316 + t175 * t320 - t457;
t449 = m(6) * (t65 + t457) + t143;
t448 = 2 * m(4);
t447 = 2 * m(5);
t446 = 0.2e1 * m(6);
t445 = 2 * m(7);
t444 = -2 * mrSges(3,3);
t443 = -2 * mrSges(4,3);
t442 = -2 * mrSges(5,3);
t132 = t151 * t320 + t254 * t316;
t79 = qJD(5) * t132 + t121 * t316 - t320 * t347;
t27 = Ifges(6,1) * t78 - Ifges(6,4) * t79 + Ifges(6,5) * t120;
t441 = t27 / 0.2e1;
t71 = Ifges(6,1) * t132 - Ifges(6,4) * t131 + Ifges(6,5) * t150;
t440 = t71 / 0.2e1;
t213 = t242 * t320 + t313 * t316;
t171 = qJD(5) * t213 + t236 * t316;
t106 = Ifges(6,1) * t170 - Ifges(6,4) * t171 + Ifges(6,5) * t235;
t439 = t106 / 0.2e1;
t140 = Ifges(6,1) * t213 - Ifges(6,4) * t212 + Ifges(6,5) * t241;
t438 = t140 / 0.2e1;
t406 = Ifges(7,4) * t315;
t282 = Ifges(7,2) * t319 + t406;
t405 = Ifges(7,4) * t319;
t341 = -Ifges(7,2) * t315 + t405;
t191 = -t282 * t374 + (Ifges(7,6) * t316 + t320 * t341) * qJD(5);
t437 = t191 / 0.2e1;
t284 = Ifges(7,1) * t315 + t405;
t342 = Ifges(7,1) * t319 - t406;
t192 = -t284 * t374 + (Ifges(7,5) * t316 + t320 * t342) * qJD(5);
t436 = t192 / 0.2e1;
t435 = t242 / 0.2e1;
t244 = -Ifges(7,6) * t320 + t316 * t341;
t434 = t244 / 0.2e1;
t245 = -Ifges(7,5) * t320 + t316 * t342;
t433 = t245 / 0.2e1;
t305 = Ifges(7,5) * t373;
t375 = qJD(6) * t315;
t267 = -Ifges(7,6) * t375 + t305;
t432 = t267 / 0.2e1;
t306 = Ifges(6,5) * t378;
t380 = qJD(5) * t316;
t431 = -Ifges(6,6) * t380 / 0.2e1 + t306 / 0.2e1;
t269 = t341 * qJD(6);
t430 = t269 / 0.2e1;
t271 = t342 * qJD(6);
t429 = t271 / 0.2e1;
t408 = Ifges(6,4) * t316;
t272 = (Ifges(6,1) * t320 - t408) * qJD(5);
t428 = t272 / 0.2e1;
t427 = Ifges(7,5) * t421 + Ifges(7,6) * t419;
t426 = Ifges(6,5) * t316 / 0.2e1 + Ifges(6,6) * t320 / 0.2e1;
t425 = t282 / 0.2e1;
t424 = t284 / 0.2e1;
t407 = Ifges(6,4) * t320;
t285 = Ifges(6,1) * t316 + t407;
t423 = t285 / 0.2e1;
t422 = -t315 / 0.2e1;
t420 = -t319 / 0.2e1;
t416 = pkin(5) * t316;
t414 = pkin(12) * t320;
t411 = mrSges(7,3) * t316;
t410 = Ifges(4,4) * t317;
t409 = Ifges(4,4) * t321;
t404 = Ifges(7,6) * t315;
t252 = -pkin(9) * t365 + t290;
t403 = t252 * mrSges(3,2);
t253 = t258 * qJD(2);
t402 = t253 * mrSges(3,1);
t113 = Ifges(4,5) * t168 + Ifges(4,6) * t167 + Ifges(4,3) * t347;
t66 = Ifges(5,5) * t121 - Ifges(5,6) * t120 + Ifges(5,3) * t347;
t401 = t113 + t66;
t28 = -t316 * t49 + t320 * t80;
t18 = -pkin(5) * t150 - t28;
t400 = qJD(5) * t18;
t393 = t315 * t320;
t389 = t319 * t320;
t123 = -t155 * t316 + t178 * t320;
t102 = -pkin(5) * t241 - t123;
t382 = qJD(5) * t102;
t381 = qJD(5) * t315;
t379 = qJD(5) * t319;
t260 = -pkin(5) * t320 - pkin(12) * t316 + t303;
t302 = pkin(3) * t309 + pkin(11);
t208 = t260 * t319 - t302 * t393;
t377 = qJD(6) * t208;
t209 = t260 * t315 + t302 * t389;
t376 = qJD(6) * t209;
t372 = t320 * t445;
t87 = -t132 * t315 + t150 * t319;
t32 = qJD(6) * t87 + t120 * t315 + t319 * t78;
t88 = t132 * t319 + t150 * t315;
t33 = -qJD(6) * t88 + t120 * t319 - t315 * t78;
t7 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t79;
t26 = Ifges(6,4) * t78 - Ifges(6,2) * t79 + Ifges(6,6) * t120;
t371 = -t26 / 0.2e1 + t7 / 0.2e1;
t25 = Ifges(6,5) * t78 - Ifges(6,6) * t79 + Ifges(6,3) * t120;
t176 = -t213 * t315 + t241 * t319;
t100 = qJD(6) * t176 + t170 * t319 + t235 * t315;
t177 = t213 * t319 + t241 * t315;
t101 = -qJD(6) * t177 - t170 * t315 + t235 * t319;
t43 = Ifges(7,5) * t100 + Ifges(7,6) * t101 + Ifges(7,3) * t171;
t36 = Ifges(7,5) * t88 + Ifges(7,6) * t87 + Ifges(7,3) * t131;
t70 = Ifges(6,4) * t132 - Ifges(6,2) * t131 + Ifges(6,6) * t150;
t370 = t36 / 0.2e1 - t70 / 0.2e1;
t13 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t4 = -pkin(5) * t120 - t6;
t369 = -m(7) * t4 - t13;
t105 = Ifges(6,4) * t170 - Ifges(6,2) * t171 + Ifges(6,6) * t235;
t368 = -t105 / 0.2e1 + t43 / 0.2e1;
t139 = Ifges(6,4) * t213 - Ifges(6,2) * t212 + Ifges(6,6) * t241;
t95 = Ifges(7,5) * t177 + Ifges(7,6) * t176 + Ifges(7,3) * t212;
t367 = -t139 / 0.2e1 + t95 / 0.2e1;
t104 = Ifges(6,5) * t170 - Ifges(6,6) * t171 + Ifges(6,3) * t235;
t56 = -mrSges(7,1) * t101 + mrSges(7,2) * t100;
t59 = -pkin(5) * t235 - t65;
t366 = -m(7) * t59 - t56;
t360 = t302 * t380;
t190 = -Ifges(7,5) * t453 + Ifges(7,6) * t328 + Ifges(7,3) * t380;
t270 = (-Ifges(6,2) * t316 + t407) * qJD(5);
t355 = -t270 / 0.2e1 + t190 / 0.2e1;
t243 = -Ifges(7,3) * t320 + (Ifges(7,5) * t319 - t404) * t316;
t283 = Ifges(6,2) * t320 + t408;
t354 = -t283 / 0.2e1 + t243 / 0.2e1;
t23 = -t309 * t57 + t312 * t51;
t72 = t120 * mrSges(5,1) + t121 * mrSges(5,2);
t183 = t235 * mrSges(5,1) + t236 * mrSges(5,2);
t54 = -t309 * t112 + t312 * t94;
t146 = t200 * t309 - t312 * t454;
t172 = t218 * t312 - t309 * t233;
t275 = (-t414 + t416) * qJD(5);
t157 = t275 * t315 - t319 * t360 + t377;
t350 = t157 - t377;
t158 = t275 * t319 + t315 * t360 - t376;
t349 = -t158 - t376;
t19 = pkin(12) * t150 + t412;
t48 = -pkin(4) * t254 - t54;
t34 = pkin(5) * t131 - pkin(12) * t132 + t48;
t10 = -t19 * t315 + t319 * t34;
t21 = -pkin(4) * t347 - t23;
t12 = pkin(5) * t79 - pkin(12) * t78 + t21;
t5 = t320 * t22 + t316 * t40 + t80 * t378 - t380 * t49;
t3 = pkin(12) * t120 + t5;
t1 = qJD(6) * t10 + t12 * t315 + t3 * t319;
t11 = t19 * t319 + t315 * t34;
t2 = -qJD(6) * t11 + t12 * t319 - t3 * t315;
t345 = t1 * t319 - t2 * t315;
t344 = t316 * mrSges(6,1) + t320 * mrSges(6,2);
t278 = -mrSges(7,1) * t319 + mrSges(7,2) * t315;
t343 = mrSges(7,1) * t315 + mrSges(7,2) * t319;
t14 = mrSges(7,1) * t79 - mrSges(7,3) * t32;
t15 = -mrSges(7,2) * t79 + mrSges(7,3) * t33;
t340 = -t315 * t14 + t319 * t15;
t64 = t320 * t147 - t155 * t380 + t316 * t175 + t178 * t378;
t58 = pkin(12) * t235 + t64;
t103 = pkin(12) * t241 + t456;
t154 = -pkin(4) * t313 - t172;
t122 = pkin(5) * t212 - pkin(12) * t213 + t154;
t60 = -t103 * t315 + t122 * t319;
t85 = pkin(5) * t171 - pkin(12) * t170 + t146;
t16 = qJD(6) * t60 + t315 * t85 + t319 * t58;
t61 = t103 * t319 + t122 * t315;
t17 = -qJD(6) * t61 - t315 * t58 + t319 * t85;
t339 = t16 * t319 - t17 * t315;
t81 = mrSges(7,1) * t171 - mrSges(7,3) * t100;
t82 = -mrSges(7,2) * t171 + mrSges(7,3) * t101;
t337 = -t315 * t81 + t319 * t82;
t37 = Ifges(7,4) * t88 + Ifges(7,2) * t87 + Ifges(7,6) * t131;
t38 = Ifges(7,1) * t88 + Ifges(7,4) * t87 + Ifges(7,5) * t131;
t332 = t37 * t422 + t38 * t419;
t96 = Ifges(7,4) * t177 + Ifges(7,2) * t176 + Ifges(7,6) * t212;
t97 = Ifges(7,1) * t177 + Ifges(7,4) * t176 + Ifges(7,5) * t212;
t331 = t419 * t97 + t422 * t96;
t42 = -mrSges(6,2) * t120 - mrSges(6,3) * t79;
t46 = -mrSges(7,1) * t87 + mrSges(7,2) * t88;
t90 = mrSges(6,1) * t150 - mrSges(6,3) * t132;
t326 = t42 + m(6) * (-qJD(5) * t28 + t5) + (t46 - t90) * qJD(5);
t128 = -mrSges(7,1) * t176 + mrSges(7,2) * t177;
t144 = -mrSges(6,2) * t235 - mrSges(6,3) * t171;
t180 = mrSges(6,1) * t241 - mrSges(6,3) * t213;
t325 = t144 + m(6) * (-qJD(5) * t123 + t64) + (t128 - t180) * qJD(5);
t324 = -t10 * t373 - t11 * t375 + t345;
t323 = -t373 * t60 - t375 * t61 + t339;
t288 = Ifges(3,5) * t322 * t386;
t274 = -mrSges(7,1) * t320 - t319 * t411;
t273 = mrSges(7,2) * t320 - t315 * t411;
t266 = t344 * qJD(5);
t265 = t343 * qJD(6);
t263 = -mrSges(4,2) * t313 + mrSges(4,3) * t398;
t262 = mrSges(4,1) * t313 - mrSges(4,3) * t399;
t259 = t343 * t316;
t256 = -pkin(9) * t311 * t318 + t300;
t255 = -pkin(10) * t399 + t298;
t251 = t257 * qJD(3);
t250 = -pkin(10) * t364 + t289;
t249 = (Ifges(4,1) * t321 - t410) * t385;
t248 = (-Ifges(4,2) * t317 + t409) * t385;
t246 = (mrSges(4,1) * t317 + mrSges(4,2) * t321) * t385;
t239 = Ifges(4,5) * t313 + (Ifges(4,1) * t317 + t409) * t310;
t238 = Ifges(4,6) * t313 + (Ifges(4,2) * t321 + t410) * t310;
t232 = -mrSges(7,2) * t380 + mrSges(7,3) * t328;
t231 = mrSges(7,1) * t380 + mrSges(7,3) * t453;
t222 = mrSges(5,1) * t313 - mrSges(5,3) * t242;
t221 = -mrSges(5,2) * t313 - mrSges(5,3) * t241;
t202 = mrSges(7,1) * t328 + mrSges(7,2) * t453;
t189 = mrSges(5,1) * t241 + mrSges(5,2) * t242;
t188 = Ifges(5,1) * t242 - Ifges(5,4) * t241 + Ifges(5,5) * t313;
t187 = Ifges(5,4) * t242 - Ifges(5,2) * t241 + Ifges(5,6) * t313;
t186 = Ifges(5,1) * t236 - Ifges(5,4) * t235;
t185 = Ifges(5,4) * t236 - Ifges(5,2) * t235;
t182 = mrSges(4,1) * t254 - mrSges(4,3) * t206;
t181 = -mrSges(4,2) * t254 + mrSges(4,3) * t205;
t179 = -mrSges(6,2) * t241 - mrSges(6,3) * t212;
t153 = mrSges(6,1) * t212 + mrSges(6,2) * t213;
t149 = mrSges(4,1) * t347 - mrSges(4,3) * t168;
t148 = -mrSges(4,2) * t347 + mrSges(4,3) * t167;
t142 = Ifges(4,1) * t206 + Ifges(4,4) * t205 + Ifges(4,5) * t254;
t141 = Ifges(4,4) * t206 + Ifges(4,2) * t205 + Ifges(4,6) * t254;
t138 = Ifges(6,5) * t213 - Ifges(6,6) * t212 + Ifges(6,3) * t241;
t137 = mrSges(7,1) * t212 - mrSges(7,3) * t177;
t136 = -mrSges(7,2) * t212 + mrSges(7,3) * t176;
t135 = mrSges(5,1) * t254 - mrSges(5,3) * t151;
t134 = -mrSges(5,2) * t254 - mrSges(5,3) * t150;
t126 = mrSges(6,1) * t171 + mrSges(6,2) * t170;
t125 = -mrSges(4,1) * t167 + mrSges(4,2) * t168;
t115 = Ifges(4,1) * t168 + Ifges(4,4) * t167 + Ifges(4,5) * t347;
t114 = Ifges(4,4) * t168 + Ifges(4,2) * t167 + Ifges(4,6) * t347;
t109 = mrSges(5,1) * t347 - mrSges(5,3) * t121;
t108 = -mrSges(5,2) * t347 - mrSges(5,3) * t120;
t107 = mrSges(5,1) * t150 + mrSges(5,2) * t151;
t92 = Ifges(5,1) * t151 - Ifges(5,4) * t150 + Ifges(5,5) * t254;
t91 = Ifges(5,4) * t151 - Ifges(5,2) * t150 + Ifges(5,6) * t254;
t89 = -mrSges(6,2) * t150 - mrSges(6,3) * t131;
t86 = mrSges(6,1) * t131 + mrSges(6,2) * t132;
t69 = Ifges(6,5) * t132 - Ifges(6,6) * t131 + Ifges(6,3) * t150;
t68 = Ifges(5,1) * t121 - Ifges(5,4) * t120 + Ifges(5,5) * t347;
t67 = Ifges(5,4) * t121 - Ifges(5,2) * t120 + Ifges(5,6) * t347;
t63 = mrSges(7,1) * t131 - mrSges(7,3) * t88;
t62 = -mrSges(7,2) * t131 + mrSges(7,3) * t87;
t45 = Ifges(7,1) * t100 + Ifges(7,4) * t101 + Ifges(7,5) * t171;
t44 = Ifges(7,4) * t100 + Ifges(7,2) * t101 + Ifges(7,6) * t171;
t35 = mrSges(6,1) * t79 + mrSges(6,2) * t78;
t9 = Ifges(7,1) * t32 + Ifges(7,4) * t33 + Ifges(7,5) * t79;
t8 = Ifges(7,4) * t32 + Ifges(7,2) * t33 + Ifges(7,6) * t79;
t20 = [(t21 * t48 + t28 * t6 + t412 * t5) * t446 + 0.2e1 * t412 * t42 + 0.2e1 * m(3) * (t252 * t258 - t253 * t256) + (t288 - 0.2e1 * t402 - 0.2e1 * t403) * t314 + (t1 * t11 + t10 * t2 + t18 * t4) * t445 + (t127 * t133 + t23 * t54 + t24 * t55) * t447 + (t129 * t84 + t130 * t83 + t159 * t169) * t448 + (-t91 + t69) * t120 + t205 * t114 + 0.2e1 * t169 * (-mrSges(4,1) * t205 + mrSges(4,2) * t206) + t206 * t115 + 0.2e1 * t83 * t181 + 0.2e1 * t84 * t182 + t167 * t141 + t168 * t142 + 0.2e1 * t159 * t125 + t151 * t68 + 0.2e1 * t130 * t148 + 0.2e1 * t129 * t149 + t132 * t27 + 0.2e1 * t133 * t72 + 0.2e1 * t24 * t134 + 0.2e1 * t23 * t135 + 0.2e1 * t127 * t107 + t121 * t92 + (0.2e1 * (t252 * t322 + t253 * t318) * mrSges(3,3) + ((t256 * t444 + Ifges(3,5) * t314 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t322) * t458) * t322 + (t258 * t444 - 0.2e1 * Ifges(3,6) * t314 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t318 + (Ifges(3,1) - Ifges(3,2)) * t322) * t458 + (Ifges(4,5) * t206 + Ifges(5,5) * t151 + Ifges(4,6) * t205 - Ifges(5,6) * t150 + (Ifges(4,3) + Ifges(5,3)) * t254) * t310) * t318) * qJD(2)) * t311 + 0.2e1 * t55 * t108 + 0.2e1 * t54 * t109 + 0.2e1 * t5 * t89 + 0.2e1 * t6 * t90 + 0.2e1 * t21 * t86 + t87 * t8 + t88 * t9 + t78 * t71 + 0.2e1 * t1 * t62 + 0.2e1 * t2 * t63 + 0.2e1 * t48 * t35 + (-t67 + t25) * t150 + 0.2e1 * t4 * t46 + 0.2e1 * t28 * t41 + t33 * t37 + t32 * t38 + 0.2e1 * t18 * t13 + 0.2e1 * t10 * t14 + 0.2e1 * t11 * t15 + (-t26 + t7) * t131 + (t36 - t70) * t79 + t401 * t254; m(6) * (t123 * t6 + t146 * t48 + t154 * t21 + t28 * t65 + t412 * t64 + t456 * t5) + t456 * t42 + t412 * t144 + m(4) * (-pkin(2) * t169 * t310 - t129 * t251 + t130 * t250 + t255 * t84 + t257 * t83) + (t92 / 0.2e1 - t54 * mrSges(5,3)) * t236 + (t69 / 0.2e1 - t91 / 0.2e1 - t55 * mrSges(5,3)) * t235 + (t247 / 0.2e1 + t184 / 0.2e1) * t254 + t78 * t438 + t132 * t439 + t170 * t440 + t213 * t441 + t68 * t435 + (-t185 / 0.2e1 + t104 / 0.2e1) * t150 - t402 - t403 + (-t187 / 0.2e1 + t138 / 0.2e1) * t120 + m(5) * (t127 * t264 + t133 * t348 - t146 * t54 + t147 * t55 + t172 * t23 + t173 * t24) + (t66 / 0.2e1 + t113 / 0.2e1) * t313 + (t25 / 0.2e1 - t67 / 0.2e1) * t241 + t255 * t149 + t257 * t148 + t84 * t262 + t83 * t263 + t264 * t72 + t159 * t246 + t205 * t248 / 0.2e1 + t206 * t249 / 0.2e1 + t250 * t181 - t251 * t182 + t167 * t238 / 0.2e1 + t168 * t239 / 0.2e1 + t24 * t221 + t23 * t222 + t133 * t183 + t151 * t186 / 0.2e1 + t121 * t188 / 0.2e1 + t127 * t189 + t176 * t8 / 0.2e1 + t177 * t9 / 0.2e1 + t5 * t179 + t6 * t180 + t172 * t109 + t173 * t108 + t21 * t153 + t154 * t35 + t147 * t134 + t28 * t143 + t1 * t136 + t2 * t137 + t4 * t128 + t123 * t41 + t48 * t126 + t33 * t96 / 0.2e1 + t32 * t97 / 0.2e1 + t100 * t38 / 0.2e1 + t101 * t37 / 0.2e1 + t102 * t13 + t64 * t89 + t65 * t90 + t87 * t44 / 0.2e1 + t88 * t45 / 0.2e1 + t10 * t81 + t11 * t82 + t16 * t62 + t17 * t63 + t59 * t46 + t60 * t14 + t61 * t15 + t18 * t56 + t288 + t370 * t171 + t371 * t212 - Ifges(3,6) * t365 + t367 * t79 + t368 * t131 + (t169 * (-mrSges(4,1) * t321 + mrSges(4,2) * t317) + t321 * t114 / 0.2e1 + t317 * t115 / 0.2e1 - pkin(2) * t125 + (Ifges(5,5) * t435 - Ifges(5,6) * t241 / 0.2e1 + (Ifges(4,5) * t317 + Ifges(4,6) * t321) * t310 / 0.2e1 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t313) * t365 + ((-t129 * mrSges(4,3) + t142 / 0.2e1) * t321 + (pkin(3) * t107 - t130 * mrSges(4,3) - t141 / 0.2e1) * t317) * qJD(3)) * t310 + m(7) * (t1 * t61 + t10 * t17 + t102 * t4 + t11 * t16 + t18 * t59 + t2 * t60) + (t86 - t135) * t146; (t123 * t65 + t146 * t154 + t456 * t64) * t446 + 0.2e1 * t456 * t144 + (-0.2e1 * pkin(2) * t246 + t321 * t248 + t317 * t249 + ((t255 * t443 + t239) * t321 + (t257 * t443 - t238 + 0.2e1 * (m(5) * t264 + t189) * pkin(3)) * t317) * qJD(3)) * t310 + (t250 * t257 - t251 * t255) * t448 + (t102 * t59 + t16 * t61 + t17 * t60) * t445 + (-t146 * t172 + t147 * t173) * t447 + 0.2e1 * (-t222 + t153) * t146 + t455 * t313 - 0.2e1 * t251 * t262 + 0.2e1 * t250 * t263 + 0.2e1 * t264 * t183 + t242 * t186 + t213 * t106 + 0.2e1 * t147 * t221 + t176 * t44 + t177 * t45 + 0.2e1 * t64 * t179 + 0.2e1 * t65 * t180 + t170 * t140 + 0.2e1 * t154 * t126 + 0.2e1 * t123 * t143 + 0.2e1 * t16 * t136 + 0.2e1 * t17 * t137 + 0.2e1 * t59 * t128 + t100 * t97 + t101 * t96 + 0.2e1 * t102 * t56 + 0.2e1 * t60 * t81 + 0.2e1 * t61 * t82 + (t104 - t185) * t241 + (-t105 + t43) * t212 + (t173 * t442 + t138 - t187) * t235 + (t172 * t442 + t188) * t236 + (t95 - t139) * t171; t401 + t132 * t428 + t150 * t431 + t32 * t433 + t33 * t434 + t88 * t436 + t87 * t437 + t78 * t423 + t120 * t426 + m(7) * (t1 * t209 + t10 * t158 + t11 * t157 + t2 * t208) + t354 * t79 + t355 * t131 + t452 * t21 + (t441 + t9 * t419 + t8 * t422 - t6 * mrSges(6,3) + (t37 * t420 + t38 * t422) * qJD(6) + (-mrSges(6,3) * t412 + t370) * qJD(5) + (-qJD(5) * t89 - t369 - t450) * t302) * t316 + t303 * t35 + t48 * t266 + t1 * t273 + t2 * t274 + t4 * t259 + t10 * t231 + t11 * t232 + t208 * t14 + t209 * t15 - t18 * t202 + t157 * t62 + t158 * t63 - t83 * mrSges(4,2) + t84 * mrSges(4,1) - t24 * mrSges(5,2) + t23 * mrSges(5,1) + (m(5) * (t23 * t312 + t24 * t309) + t312 * t109 + t309 * t108) * pkin(3) + (t5 * mrSges(6,3) + (-t28 * mrSges(6,3) + t332 + t440) * qJD(5) + (m(7) * t400 + t326) * t302 - t371) * t320; (t439 + t45 * t419 + t44 * t422 - t65 * mrSges(6,3) + (t420 * t96 + t422 * t97) * qJD(6) + (-mrSges(6,3) * t456 + t367) * qJD(5) + (-qJD(5) * t179 - t366 - t449) * t302) * t316 + t455 + t241 * t431 + t100 * t433 + t101 * t434 + t177 * t436 + t176 * t437 + t170 * t423 + t235 * t426 + t213 * t428 + t354 * t171 + t355 * t212 + (-mrSges(5,1) + t452) * t146 + t303 * t126 + t154 * t266 + t16 * t273 + t17 * t274 + t59 * t259 - t250 * mrSges(4,2) - t251 * mrSges(4,1) + t60 * t231 + t61 * t232 + t208 * t81 + t209 * t82 - t102 * t202 + t157 * t136 + t158 * t137 - t147 * mrSges(5,2) + m(7) * (t157 * t61 + t158 * t60 + t16 * t209 + t17 * t208) + (m(5) * (-t146 * t312 + t147 * t309) + (-t235 * t309 - t236 * t312) * mrSges(5,3)) * pkin(3) + (t64 * mrSges(6,3) + (-t123 * mrSges(6,3) + t331 + t438) * qJD(5) + (m(7) * t382 + t325) * t302 - t368) * t320; (t157 * t209 + t158 * t208) * t445 + 0.2e1 * t157 * t273 + 0.2e1 * t209 * t232 + 0.2e1 * t158 * t274 + 0.2e1 * t208 * t231 + 0.2e1 * t303 * t266 + (-t190 + t270 + (-t244 * t315 + t245 * t319 + 0.2e1 * t259 * t302 + t285) * qJD(5)) * t320 + (-t315 * t191 + t319 * t192 - 0.2e1 * t302 * t202 + t272 + (-t244 * t319 - t245 * t315) * qJD(6) + (t302 ^ 2 * t372 + t243 - t283) * qJD(5)) * t316; m(5) * t127 + (-t13 + (-t315 * t63 + t319 * t62 + t89) * qJD(5) + m(7) * (-t10 * t381 + t11 * t379 - t4) + t450) * t320 + ((-t315 * t62 - t319 * t63) * qJD(6) + m(7) * (t324 + t400) + t326 + t340) * t316 + t72; m(5) * t348 + (-t56 + (t136 * t319 - t137 * t315 + t179) * qJD(5) + m(7) * (t379 * t61 - t381 * t60 - t59) + t449) * t320 + ((-t136 * t315 - t137 * t319) * qJD(6) + m(7) * (t323 + t382) + t325 + t337) * t316 + t183; t320 * t202 + (m(7) * (t157 * t319 - t158 * t315 - t208 * t373 - t209 * t375) - t273 * t375 + t319 * t232 - t274 * t373 - t315 * t231) * t316 + (m(7) * (-t208 * t393 + t209 * t389 + (t316 ^ 2 - t320 ^ 2) * t302) + t273 * t389 - t274 * t393 + t316 * t259) * qJD(5); (-0.1e1 + t387) * t372 * t380; t8 * t419 + t9 * t421 + t33 * t425 + t32 * t424 + t18 * t265 + t131 * t432 + t87 * t430 + t88 * t429 + t4 * t278 + t79 * t427 + t6 * mrSges(6,1) - t5 * mrSges(6,2) + t332 * qJD(6) + t369 * pkin(5) + ((-t10 * t319 - t11 * t315) * qJD(6) + t345) * mrSges(7,3) + (m(7) * t324 - t373 * t63 - t375 * t62 + t340) * pkin(12) + t25; t44 * t419 + t45 * t421 + t101 * t425 + t100 * t424 + t102 * t265 + t212 * t432 + t176 * t430 + t177 * t429 + t59 * t278 + t171 * t427 - t64 * mrSges(6,2) + t65 * mrSges(6,1) + t331 * qJD(6) + t366 * pkin(5) + ((-t315 * t61 - t319 * t60) * qJD(6) + t339) * mrSges(7,3) + (m(7) * t323 - t136 * t375 - t137 * t373 + t337) * pkin(12) + t104; pkin(5) * t202 + t306 + (-t267 / 0.2e1 + (-m(7) * pkin(5) - mrSges(6,1) + t278) * t302 * qJD(5)) * t320 + (t378 * t424 + t437 + qJD(6) * t433 + t350 * mrSges(7,3) + (m(7) * t350 - qJD(6) * t274 + t232) * pkin(12)) * t319 + (-t282 * t378 / 0.2e1 + t436 - qJD(6) * t244 / 0.2e1 + t349 * mrSges(7,3) + (m(7) * t349 - qJD(6) * t273 - t231) * pkin(12)) * t315 + (t271 * t419 + t269 * t422 + t302 * t265 + (t282 * t420 + t284 * t422) * qJD(6) + (t302 * mrSges(6,2) - Ifges(6,6) + t427) * qJD(5)) * t316; -t320 * t265 + (m(7) * (t387 * t414 - t416) + t316 * t278 + t387 * t320 * mrSges(7,3) - t344) * qJD(5); -0.2e1 * pkin(5) * t265 + t269 * t319 + t271 * t315 + (-t282 * t315 + t284 * t319) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t17 - mrSges(7,2) * t16 + t43; mrSges(7,1) * t158 - mrSges(7,2) * t157 + t190; t202; t305 + (pkin(12) * t278 - t404) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
