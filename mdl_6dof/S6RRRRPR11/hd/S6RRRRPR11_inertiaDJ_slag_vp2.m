% Calculate time derivative of joint inertia matrix for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:12:27
% EndTime: 2019-03-09 23:12:47
% DurationCPUTime: 9.06s
% Computational Cost: add. (18355->832), mult. (48437->1213), div. (0->0), fcn. (47559->12), ass. (0->340)
t327 = sin(qJ(4));
t331 = cos(qJ(4));
t328 = sin(qJ(3));
t364 = qJD(4) * t328;
t332 = cos(qJ(3));
t366 = qJD(3) * t332;
t335 = -t327 * t364 + t331 * t366;
t363 = qJD(4) * t331;
t334 = t327 * t366 + t328 * t363;
t291 = -pkin(3) * t332 - pkin(10) * t328 - pkin(2);
t372 = t331 * t332;
t311 = pkin(9) * t372;
t238 = t327 * t291 + t311;
t431 = qJD(4) * t238;
t322 = sin(pkin(12));
t324 = cos(pkin(12));
t274 = t322 * t331 + t324 * t327;
t326 = sin(qJ(6));
t330 = cos(qJ(6));
t340 = t322 * t327 - t324 * t331;
t199 = -t274 * t326 - t330 * t340;
t409 = t199 / 0.2e1;
t200 = t274 * t330 - t326 * t340;
t408 = t200 / 0.2e1;
t398 = -t340 / 0.2e1;
t397 = t274 / 0.2e1;
t430 = t327 / 0.2e1;
t393 = t331 / 0.2e1;
t325 = cos(pkin(6));
t323 = sin(pkin(6));
t329 = sin(qJ(2));
t376 = t323 * t329;
t266 = t325 * t328 + t332 * t376;
t333 = cos(qJ(2));
t369 = qJD(2) * t323;
t356 = t333 * t369;
t213 = qJD(3) * t266 + t328 * t356;
t265 = -t325 * t332 + t328 * t376;
t214 = -qJD(3) * t265 + t332 * t356;
t375 = t323 * t333;
t336 = -t266 * t331 + t327 * t375;
t368 = qJD(2) * t329;
t357 = t323 * t368;
t135 = qJD(4) * t336 - t214 * t327 + t331 * t357;
t215 = -t266 * t327 - t331 * t375;
t136 = qJD(4) * t215 + t214 * t331 + t327 * t357;
t70 = t135 * t324 - t136 * t322;
t72 = t135 * t322 + t136 * t324;
t31 = Ifges(6,5) * t72 + Ifges(6,6) * t70 + Ifges(6,3) * t213;
t57 = Ifges(5,5) * t136 + Ifges(5,6) * t135 + Ifges(5,3) * t213;
t429 = t31 + t57;
t180 = -t274 * t366 + t340 * t364;
t263 = t274 * qJD(4);
t181 = -t263 * t328 - t340 * t366;
t367 = qJD(3) * t328;
t104 = Ifges(6,5) * t181 + Ifges(6,6) * t180 + Ifges(6,3) * t367;
t189 = Ifges(5,5) * t335 - Ifges(5,6) * t334 + Ifges(5,3) * t367;
t428 = -t104 - t189;
t306 = pkin(8) * t376;
t391 = pkin(1) * t333;
t245 = t306 + (-pkin(2) - t391) * t325;
t157 = t265 * pkin(3) - t266 * pkin(10) + t245;
t270 = t325 * t329 * pkin(1) + pkin(8) * t375;
t246 = pkin(9) * t325 + t270;
t247 = (-pkin(2) * t333 - pkin(9) * t329 - pkin(1)) * t323;
t169 = t332 * t246 + t328 * t247;
t159 = -pkin(10) * t375 + t169;
t95 = t327 * t157 + t331 * t159;
t293 = -mrSges(5,1) * t331 + mrSges(5,2) * t327;
t427 = -m(5) * pkin(3) + t293;
t269 = t325 * t391 - t306;
t426 = 0.2e1 * m(5);
t425 = 2 * m(6);
t424 = 2 * m(7);
t423 = 0.2e1 * pkin(9);
t422 = -2 * mrSges(3,3);
t147 = t215 * t324 + t322 * t336;
t148 = t215 * t322 - t324 * t336;
t89 = t147 * t330 - t148 * t326;
t420 = t89 / 0.2e1;
t90 = t147 * t326 + t148 * t330;
t419 = t90 / 0.2e1;
t264 = t340 * qJD(4);
t132 = qJD(6) * t199 - t263 * t326 - t264 * t330;
t418 = t132 / 0.2e1;
t133 = -qJD(6) * t200 - t263 * t330 + t264 * t326;
t417 = t133 / 0.2e1;
t416 = t135 / 0.2e1;
t139 = Ifges(7,4) * t200 + Ifges(7,2) * t199;
t415 = t139 / 0.2e1;
t140 = Ifges(7,1) * t200 + Ifges(7,4) * t199;
t414 = t140 / 0.2e1;
t413 = t147 / 0.2e1;
t412 = t148 / 0.2e1;
t251 = t274 * t328;
t252 = t340 * t328;
t173 = -t251 * t330 + t252 * t326;
t411 = t173 / 0.2e1;
t174 = -t251 * t326 - t252 * t330;
t410 = t174 / 0.2e1;
t206 = Ifges(6,4) * t274 - Ifges(6,2) * t340;
t407 = t206 / 0.2e1;
t207 = Ifges(6,1) * t274 - Ifges(6,4) * t340;
t406 = t207 / 0.2e1;
t405 = t215 / 0.2e1;
t404 = -t336 / 0.2e1;
t384 = Ifges(5,4) * t327;
t343 = Ifges(5,1) * t331 - t384;
t250 = -Ifges(5,5) * t332 + t328 * t343;
t403 = t250 / 0.2e1;
t402 = -t251 / 0.2e1;
t401 = -t252 / 0.2e1;
t400 = -t263 / 0.2e1;
t399 = -t264 / 0.2e1;
t383 = Ifges(5,4) * t331;
t298 = Ifges(5,1) * t327 + t383;
t396 = t298 / 0.2e1;
t395 = -t327 / 0.2e1;
t394 = -t331 / 0.2e1;
t392 = m(5) * t332;
t390 = pkin(4) * t322;
t389 = pkin(9) * t327;
t388 = pkin(9) * t332;
t321 = t328 * pkin(9);
t387 = -qJ(5) - pkin(10);
t255 = t270 * qJD(2);
t125 = t213 * pkin(3) - t214 * pkin(10) + t255;
t253 = (pkin(2) * t329 - pkin(9) * t333) * t369;
t254 = t269 * qJD(2);
t109 = -t246 * t367 + t247 * t366 + t328 * t253 + t332 * t254;
t99 = pkin(10) * t357 + t109;
t39 = -qJD(4) * t95 + t331 * t125 - t327 * t99;
t26 = pkin(4) * t213 - qJ(5) * t136 + qJD(5) * t336 + t39;
t365 = qJD(4) * t327;
t38 = t327 * t125 + t157 * t363 - t159 * t365 + t331 * t99;
t28 = qJ(5) * t135 + qJD(5) * t215 + t38;
t11 = t322 * t26 + t324 * t28;
t94 = t331 * t157 - t159 * t327;
t63 = pkin(4) * t265 + qJ(5) * t336 + t94;
t83 = qJ(5) * t215 + t95;
t36 = t322 * t63 + t324 * t83;
t386 = Ifges(4,4) * t328;
t385 = Ifges(4,4) * t332;
t382 = Ifges(5,6) * t327;
t381 = pkin(9) * qJD(3);
t380 = t254 * mrSges(3,2);
t379 = t255 * mrSges(3,1);
t378 = t255 * mrSges(4,1);
t377 = t255 * mrSges(4,2);
t374 = t327 * t328;
t373 = t328 * t331;
t74 = Ifges(7,5) * t132 + Ifges(7,6) * t133;
t362 = qJD(5) * t331;
t289 = (pkin(3) * t328 - pkin(10) * t332) * qJD(3);
t370 = t331 * t289 + t367 * t389;
t129 = -t328 * t362 + (pkin(4) * t328 - qJ(5) * t372) * qJD(3) + (-t311 + (qJ(5) * t328 - t291) * t327) * qJD(4) + t370;
t371 = t327 * t289 + t291 * t363;
t141 = (-qJ(5) * qJD(4) - t381) * t373 + (-qJD(5) * t328 + (-pkin(9) * qJD(4) - qJ(5) * qJD(3)) * t332) * t327 + t371;
t71 = t322 * t129 + t324 * t141;
t276 = t331 * t291;
t201 = -qJ(5) * t373 + t276 + (-pkin(4) - t389) * t332;
t217 = -qJ(5) * t374 + t238;
t145 = t322 * t201 + t324 * t217;
t350 = qJD(4) * t387;
t261 = t327 * t350 + t362;
t262 = -qJD(5) * t327 + t331 * t350;
t183 = t324 * t261 + t322 * t262;
t193 = -Ifges(6,5) * t264 - Ifges(6,6) * t263;
t292 = t387 * t327;
t294 = t387 * t331;
t219 = t322 * t292 - t324 * t294;
t290 = pkin(4) * t374 + t321;
t22 = qJD(6) * t89 + t326 * t70 + t330 * t72;
t23 = -qJD(6) * t90 - t326 * t72 + t330 * t70;
t6 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t213;
t91 = qJD(6) * t173 + t180 * t326 + t181 * t330;
t92 = -qJD(6) * t174 + t180 * t330 - t181 * t326;
t43 = Ifges(7,5) * t91 + Ifges(7,6) * t92 + Ifges(7,3) * t367;
t360 = pkin(4) * t365;
t359 = Ifges(4,6) * t375;
t358 = Ifges(4,5) * t214 - Ifges(4,6) * t213 + Ifges(4,3) * t357;
t233 = pkin(4) * t334 + pkin(9) * t366;
t314 = -pkin(4) * t331 - pkin(3);
t124 = -Ifges(5,1) * t336 + Ifges(5,4) * t215 + Ifges(5,5) * t265;
t351 = t124 * t393;
t37 = -t70 * mrSges(6,1) + t72 * mrSges(6,2);
t9 = -t23 * mrSges(7,1) + t22 * mrSges(7,2);
t47 = -t92 * mrSges(7,1) + t91 * mrSges(7,2);
t10 = t324 * t26 - t28 * t322;
t35 = -t322 * t83 + t324 * t63;
t114 = -t180 * mrSges(6,1) + t181 * mrSges(6,2);
t192 = t263 * mrSges(6,1) - t264 * mrSges(6,2);
t73 = -t133 * mrSges(7,1) + t132 * mrSges(7,2);
t69 = t324 * t129 - t141 * t322;
t144 = t324 * t201 - t217 * t322;
t168 = -t328 * t246 + t247 * t332;
t182 = -t261 * t322 + t324 * t262;
t218 = t324 * t292 + t294 * t322;
t166 = (-t331 * t367 - t332 * t365) * pkin(9) + t371;
t237 = -t327 * t388 + t276;
t349 = -qJD(4) * t237 + t166;
t348 = t357 / 0.2e1;
t318 = Ifges(5,5) * t363;
t347 = -Ifges(5,6) * t365 / 0.2e1 + t318 / 0.2e1 + t193 / 0.2e1 + t74 / 0.2e1;
t158 = pkin(3) * t375 - t168;
t346 = Ifges(5,5) * t430 + Ifges(6,5) * t397 + Ifges(7,5) * t408 + Ifges(5,6) * t393 + Ifges(6,6) * t398 + Ifges(7,6) * t409;
t345 = mrSges(5,1) * t327 + mrSges(5,2) * t331;
t313 = pkin(4) * t324 + pkin(5);
t259 = t313 * t330 - t326 * t390;
t243 = t259 * qJD(6);
t260 = t313 * t326 + t330 * t390;
t244 = t260 * qJD(6);
t344 = -t244 * mrSges(7,1) - t243 * mrSges(7,2);
t342 = -Ifges(5,2) * t327 + t383;
t296 = Ifges(5,2) * t331 + t384;
t29 = pkin(5) * t265 - pkin(11) * t148 + t35;
t30 = pkin(11) * t147 + t36;
t12 = t29 * t330 - t30 * t326;
t13 = t29 * t326 + t30 * t330;
t341 = -t327 * t39 + t331 * t38;
t108 = -pkin(5) * t332 + pkin(11) * t252 + t144;
t113 = -pkin(11) * t251 + t145;
t55 = t108 * t330 - t113 * t326;
t56 = t108 * t326 + t113 * t330;
t184 = -pkin(11) * t274 + t218;
t185 = -pkin(11) * t340 + t219;
t111 = t184 * t330 - t185 * t326;
t112 = t184 * t326 + t185 * t330;
t151 = pkin(11) * t264 + t182;
t152 = -pkin(11) * t263 + t183;
t49 = qJD(6) * t111 + t151 * t326 + t152 * t330;
t50 = -qJD(6) * t112 + t151 * t330 - t152 * t326;
t339 = t50 * mrSges(7,1) - t49 * mrSges(7,2) + t74;
t110 = -t246 * t366 - t247 * t367 + t253 * t332 - t328 * t254;
t4 = pkin(5) * t213 - pkin(11) * t72 + t10;
t5 = pkin(11) * t70 + t11;
t2 = qJD(6) * t12 + t326 * t4 + t330 * t5;
t3 = -qJD(6) * t13 - t326 * t5 + t330 * t4;
t338 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t6;
t51 = pkin(5) * t367 - pkin(11) * t181 + t69;
t52 = pkin(11) * t180 + t71;
t15 = qJD(6) * t55 + t326 * t51 + t330 * t52;
t16 = -qJD(6) * t56 - t326 * t52 + t330 * t51;
t337 = t16 * mrSges(7,1) - t15 * mrSges(7,2) + t43;
t117 = -pkin(4) * t215 + t158;
t100 = -pkin(3) * t357 - t110;
t60 = -pkin(4) * t135 + t100;
t319 = Ifges(4,5) * t366;
t301 = Ifges(3,5) * t356;
t299 = Ifges(4,1) * t328 + t385;
t297 = Ifges(4,2) * t332 + t386;
t288 = -mrSges(5,1) * t332 - mrSges(5,3) * t373;
t287 = mrSges(5,2) * t332 - mrSges(5,3) * t374;
t286 = (Ifges(4,1) * t332 - t386) * qJD(3);
t285 = t343 * qJD(4);
t284 = (-Ifges(4,2) * t328 + t385) * qJD(3);
t283 = t342 * qJD(4);
t281 = (mrSges(4,1) * t328 + mrSges(4,2) * t332) * qJD(3);
t280 = t345 * qJD(4);
t271 = t345 * t328;
t249 = -Ifges(5,6) * t332 + t328 * t342;
t248 = -Ifges(5,3) * t332 + (Ifges(5,5) * t331 - t382) * t328;
t239 = pkin(5) * t340 + t314;
t229 = -mrSges(5,2) * t367 - mrSges(5,3) * t334;
t228 = mrSges(5,1) * t367 - mrSges(5,3) * t335;
t227 = pkin(5) * t263 + t360;
t223 = -mrSges(6,1) * t332 + mrSges(6,3) * t252;
t222 = mrSges(6,2) * t332 - mrSges(6,3) * t251;
t221 = -mrSges(4,1) * t375 - t266 * mrSges(4,3);
t220 = mrSges(4,2) * t375 - t265 * mrSges(4,3);
t204 = mrSges(6,1) * t340 + mrSges(6,2) * t274;
t202 = pkin(5) * t251 + t290;
t198 = mrSges(5,1) * t334 + mrSges(5,2) * t335;
t195 = -Ifges(6,1) * t264 - Ifges(6,4) * t263;
t194 = -Ifges(6,4) * t264 - Ifges(6,2) * t263;
t191 = -t298 * t364 + (Ifges(5,5) * t328 + t332 * t343) * qJD(3);
t190 = -t296 * t364 + (Ifges(5,6) * t328 + t332 * t342) * qJD(3);
t188 = mrSges(4,1) * t357 - mrSges(4,3) * t214;
t187 = -mrSges(4,2) * t357 - mrSges(4,3) * t213;
t186 = mrSges(6,1) * t251 - mrSges(6,2) * t252;
t179 = Ifges(4,1) * t266 - Ifges(4,4) * t265 - Ifges(4,5) * t375;
t178 = Ifges(4,4) * t266 - Ifges(4,2) * t265 - t359;
t172 = -Ifges(6,1) * t252 - Ifges(6,4) * t251 - Ifges(6,5) * t332;
t171 = -Ifges(6,4) * t252 - Ifges(6,2) * t251 - Ifges(6,6) * t332;
t170 = -Ifges(6,5) * t252 - Ifges(6,6) * t251 - Ifges(6,3) * t332;
t167 = t370 - t431;
t165 = mrSges(5,1) * t265 + mrSges(5,3) * t336;
t164 = -mrSges(5,2) * t265 + mrSges(5,3) * t215;
t163 = mrSges(6,1) * t367 - mrSges(6,3) * t181;
t162 = -mrSges(6,2) * t367 + mrSges(6,3) * t180;
t161 = -mrSges(7,1) * t332 - mrSges(7,3) * t174;
t160 = mrSges(7,2) * t332 + mrSges(7,3) * t173;
t150 = -mrSges(5,1) * t215 - mrSges(5,2) * t336;
t149 = mrSges(4,1) * t213 + mrSges(4,2) * t214;
t146 = -pkin(5) * t180 + t233;
t143 = Ifges(4,1) * t214 - Ifges(4,4) * t213 + Ifges(4,5) * t357;
t142 = Ifges(4,4) * t214 - Ifges(4,2) * t213 + Ifges(4,6) * t357;
t137 = -mrSges(7,1) * t199 + mrSges(7,2) * t200;
t123 = -Ifges(5,4) * t336 + Ifges(5,2) * t215 + Ifges(5,6) * t265;
t122 = -Ifges(5,5) * t336 + Ifges(5,6) * t215 + Ifges(5,3) * t265;
t116 = mrSges(6,1) * t265 - mrSges(6,3) * t148;
t115 = -mrSges(6,2) * t265 + mrSges(6,3) * t147;
t107 = -mrSges(7,1) * t173 + mrSges(7,2) * t174;
t106 = Ifges(6,1) * t181 + Ifges(6,4) * t180 + Ifges(6,5) * t367;
t105 = Ifges(6,4) * t181 + Ifges(6,2) * t180 + Ifges(6,6) * t367;
t103 = Ifges(7,1) * t174 + Ifges(7,4) * t173 - Ifges(7,5) * t332;
t102 = Ifges(7,4) * t174 + Ifges(7,2) * t173 - Ifges(7,6) * t332;
t101 = Ifges(7,5) * t174 + Ifges(7,6) * t173 - Ifges(7,3) * t332;
t97 = mrSges(5,1) * t213 - mrSges(5,3) * t136;
t96 = -mrSges(5,2) * t213 + mrSges(5,3) * t135;
t93 = -mrSges(6,1) * t147 + mrSges(6,2) * t148;
t85 = -mrSges(7,2) * t367 + mrSges(7,3) * t92;
t84 = mrSges(7,1) * t367 - mrSges(7,3) * t91;
t81 = -pkin(5) * t147 + t117;
t80 = -mrSges(5,1) * t135 + mrSges(5,2) * t136;
t79 = Ifges(6,1) * t148 + Ifges(6,4) * t147 + Ifges(6,5) * t265;
t78 = Ifges(6,4) * t148 + Ifges(6,2) * t147 + Ifges(6,6) * t265;
t77 = Ifges(6,5) * t148 + Ifges(6,6) * t147 + Ifges(6,3) * t265;
t76 = Ifges(7,1) * t132 + Ifges(7,4) * t133;
t75 = Ifges(7,4) * t132 + Ifges(7,2) * t133;
t65 = mrSges(7,1) * t265 - mrSges(7,3) * t90;
t64 = -mrSges(7,2) * t265 + mrSges(7,3) * t89;
t59 = Ifges(5,1) * t136 + Ifges(5,4) * t135 + Ifges(5,5) * t213;
t58 = Ifges(5,4) * t136 + Ifges(5,2) * t135 + Ifges(5,6) * t213;
t54 = mrSges(6,1) * t213 - mrSges(6,3) * t72;
t53 = -mrSges(6,2) * t213 + mrSges(6,3) * t70;
t46 = -mrSges(7,1) * t89 + mrSges(7,2) * t90;
t45 = Ifges(7,1) * t91 + Ifges(7,4) * t92 + Ifges(7,5) * t367;
t44 = Ifges(7,4) * t91 + Ifges(7,2) * t92 + Ifges(7,6) * t367;
t42 = Ifges(7,1) * t90 + Ifges(7,4) * t89 + Ifges(7,5) * t265;
t41 = Ifges(7,4) * t90 + Ifges(7,2) * t89 + Ifges(7,6) * t265;
t40 = Ifges(7,5) * t90 + Ifges(7,6) * t89 + Ifges(7,3) * t265;
t34 = -pkin(5) * t70 + t60;
t33 = Ifges(6,1) * t72 + Ifges(6,4) * t70 + Ifges(6,5) * t213;
t32 = Ifges(6,4) * t72 + Ifges(6,2) * t70 + Ifges(6,6) * t213;
t18 = -mrSges(7,2) * t213 + mrSges(7,3) * t23;
t17 = mrSges(7,1) * t213 - mrSges(7,3) * t22;
t8 = Ifges(7,1) * t22 + Ifges(7,4) * t23 + Ifges(7,5) * t213;
t7 = Ifges(7,4) * t22 + Ifges(7,2) * t23 + Ifges(7,6) * t213;
t1 = [-t336 * t59 + (t12 * t3 + t13 * t2 + t34 * t81) * t424 + (t10 * t35 + t11 * t36 + t117 * t60) * t425 + (t100 * t158 + t38 * t95 + t39 * t94) * t426 + (-t142 + t6 + 0.2e1 * t378 + t429) * t265 + (-t333 * t358 + 0.2e1 * (t254 * t333 + t255 * t329) * mrSges(3,3) + ((t269 * t422 + Ifges(3,5) * t325 + 0.2e1 * (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t333) * t323) * t333 + (t270 * t422 + Ifges(4,5) * t266 - 0.2e1 * Ifges(3,6) * t325 - Ifges(4,6) * t265 + (-0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t329 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3)) * t333) * t323) * t329) * qJD(2)) * t323 + 0.2e1 * t245 * t149 + t215 * t58 + 0.2e1 * t109 * t220 + 0.2e1 * t110 * t221 + t214 * t179 + 0.2e1 * t169 * t187 + 0.2e1 * t168 * t188 + 0.2e1 * t38 * t164 + 0.2e1 * t39 * t165 + 0.2e1 * t158 * t80 + t147 * t32 + t148 * t33 + 0.2e1 * t100 * t150 + t135 * t123 + t136 * t124 + 0.2e1 * t11 * t115 + 0.2e1 * t10 * t116 + 0.2e1 * t117 * t37 + 0.2e1 * t94 * t97 + 0.2e1 * t95 * t96 + t89 * t7 + t90 * t8 + 0.2e1 * t60 * t93 + t70 * t78 + t72 * t79 + 0.2e1 * t81 * t9 + 0.2e1 * t3 * t65 + 0.2e1 * t2 * t64 + 0.2e1 * t36 * t53 + 0.2e1 * t35 * t54 + 0.2e1 * t34 * t46 + t23 * t41 + t22 * t42 + 0.2e1 * t13 * t18 + 0.2e1 * t12 * t17 + (t122 + t77 + t40 - t178) * t213 + 0.2e1 * m(3) * (t254 * t270 - t255 * t269) + 0.2e1 * m(4) * (t109 * t169 + t110 * t168 + t245 * t255) + (t143 + 0.2e1 * t377) * t266 + (t301 - 0.2e1 * t379 - 0.2e1 * t380) * t325; -t379 - t380 + t45 * t419 + t44 * t420 + t33 * t401 + t32 * t402 + t136 * t403 + t191 * t404 + t190 * t405 + t8 * t410 + t7 * t411 + t106 * t412 + t105 * t413 + t249 * t416 + t301 + m(6) * (t10 * t144 + t11 * t145 + t117 * t233 + t290 * t60 + t35 * t69 + t36 * t71) + m(7) * (t12 * t16 + t13 * t15 + t146 * t81 + t2 * t56 + t202 * t34 + t3 * t55) + t245 * t281 + t266 * t286 / 0.2e1 + t38 * t287 + t39 * t288 + t290 * t37 + t214 * t299 / 0.2e1 + t100 * t271 + t238 * t96 + t94 * t228 + t95 * t229 + t233 * t93 + t237 * t97 + t11 * t222 + t10 * t223 + t202 * t9 + t158 * t198 + t180 * t78 / 0.2e1 + t181 * t79 / 0.2e1 + t60 * t186 + t3 * t161 + t36 * t162 + t35 * t163 + t166 * t164 + t167 * t165 + t70 * t171 / 0.2e1 + t72 * t172 / 0.2e1 + t2 * t160 + t146 * t46 - pkin(2) * t149 + t144 * t54 + t145 * t53 + t71 * t115 + t69 * t116 + t117 * t114 + t34 * t107 + t23 * t102 / 0.2e1 + t22 * t103 / 0.2e1 + t91 * t42 / 0.2e1 + t92 * t41 / 0.2e1 + t12 * t84 + t13 * t85 + t81 * t47 + t16 * t65 + t15 * t64 + t56 * t18 + t55 * t17 + (-t284 / 0.2e1 + t189 / 0.2e1 + t104 / 0.2e1 + t43 / 0.2e1) * t265 + (-t297 / 0.2e1 + t248 / 0.2e1 + t170 / 0.2e1 + t101 / 0.2e1) * t213 + ((t179 / 0.2e1 + t123 * t395 - t168 * mrSges(4,3) + t351) * t332 + (t359 / 0.2e1 + t40 / 0.2e1 + t77 / 0.2e1 + t122 / 0.2e1 - t178 / 0.2e1 - t169 * mrSges(4,3)) * t328 + (-t328 * t220 + (t150 - t221) * t332 + m(4) * (-t168 * t332 - t169 * t328) + t158 * t392) * pkin(9)) * qJD(3) + (t377 + t143 / 0.2e1 + t59 * t393 + t58 * t395 - t110 * mrSges(4,3) + Ifges(4,5) * t348 + (t123 * t394 + t124 * t395) * qJD(4) + (t80 - t188) * pkin(9)) * t328 + m(5) * (t100 * t321 + t166 * t95 + t167 * t94 + t237 * t39 + t238 * t38) + m(4) * (-pkin(2) * t255 + t109 * t388 - t110 * t321) + (-t333 * t319 / 0.2e1 - Ifges(3,6) * t368) * t323 + (-t378 + t142 / 0.2e1 - t57 / 0.2e1 - t31 / 0.2e1 - t6 / 0.2e1 + pkin(9) * t187 + t109 * mrSges(4,3) + Ifges(4,6) * t348) * t332; (t146 * t202 + t15 * t56 + t16 * t55) * t424 + (t144 * t69 + t145 * t71 + t233 * t290) * t425 + (t166 * t238 + t167 * t237) * t426 + (t284 - t43 + (-t327 * t249 + t331 * t250 + t271 * t423 + t299) * qJD(3) + t428) * t332 - 0.2e1 * pkin(2) * t281 + 0.2e1 * t166 * t287 + 0.2e1 * t167 * t288 + 0.2e1 * t290 * t114 - t251 * t105 - t252 * t106 + 0.2e1 * t238 * t229 + 0.2e1 * t233 * t186 + 0.2e1 * t237 * t228 + 0.2e1 * t71 * t222 + 0.2e1 * t69 * t223 + 0.2e1 * t202 * t47 + t180 * t171 + t181 * t172 + t173 * t44 + t174 * t45 + 0.2e1 * t16 * t161 + 0.2e1 * t145 * t162 + 0.2e1 * t144 * t163 + 0.2e1 * t15 * t160 + 0.2e1 * t146 * t107 + t92 * t102 + t91 * t103 + 0.2e1 * t55 * t84 + 0.2e1 * t56 * t85 + (t198 * t423 - t327 * t190 + t331 * t191 + t286 + (-t249 * t331 - t250 * t327) * qJD(4) + (0.2e1 * pkin(9) ^ 2 * t392 + t101 + t170 + t248 - t297) * qJD(3)) * t328; t358 + m(7) * (t111 * t3 + t112 * t2 + t12 * t50 + t13 * t49 + t227 * t81 + t239 * t34) + t42 * t418 + t76 * t419 + t75 * t420 + t79 * t399 + t78 * t400 + t285 * t404 + t283 * t405 + t72 * t406 + t70 * t407 + t8 * t408 + t7 * t409 + t195 * t412 + t194 * t413 + t22 * t414 + t23 * t415 + t296 * t416 + t41 * t417 + t58 * t393 + t136 * t396 + t33 * t397 + t32 * t398 + t427 * t100 + t59 * t430 + (-t12 * t132 + t13 * t133 + t199 * t2 - t200 * t3) * mrSges(7,3) + t314 * t37 + t158 * t280 + t239 * t9 + t227 * t46 + t218 * t54 + t219 * t53 + t60 * t204 + t117 * t192 + t182 * t116 + t183 * t115 + t34 * t137 - t109 * mrSges(4,2) + t110 * mrSges(4,1) + t111 * t17 + t112 * t18 - pkin(3) * t80 + t81 * t73 + t50 * t65 + t49 * t64 + ((-t327 * t95 - t331 * t94) * qJD(4) + t341) * mrSges(5,3) + t346 * t213 + t347 * t265 + (-t10 * t274 - t11 * t340 - t263 * t36 + t264 * t35) * mrSges(6,3) + (t351 + (pkin(4) * t93 - t123 / 0.2e1) * t327) * qJD(4) + m(6) * (t10 * t218 + t11 * t219 + t117 * t360 + t182 * t35 + t183 * t36 + t314 * t60) + (-t164 * t365 + m(5) * (-t363 * t94 - t365 * t95 + t341) + t331 * t96 - t327 * t97 - t165 * t363) * pkin(10); (t191 / 0.2e1 - t167 * mrSges(5,3) - t296 * t366 / 0.2e1 + (-t238 * mrSges(5,3) - t249 / 0.2e1 + (m(6) * t290 + t186) * pkin(4)) * qJD(4) + (-qJD(4) * t287 + m(5) * (-t167 - t431) - t228) * pkin(10)) * t327 + m(7) * (t111 * t16 + t112 * t15 + t146 * t239 + t202 * t227 + t49 * t56 + t50 * t55) + m(6) * (t144 * t182 + t145 * t183 + t218 * t69 + t219 * t71 + t233 * t314) + t103 * t418 + t171 * t400 + t195 * t401 + t194 * t402 + t181 * t406 + t180 * t407 + t45 * t408 + t44 * t409 + t76 * t410 + t75 * t411 + t91 * t414 + t92 * t415 + t102 * t417 + t106 * t397 + t105 * t398 + t172 * t399 + t319 + ((-mrSges(4,1) + t427) * t381 - t347) * t332 + (-t132 * t55 + t133 * t56 + t15 * t199 - t16 * t200) * mrSges(7,3) + t314 * t114 + t290 * t192 + t239 * t47 + t227 * t107 + t233 * t204 + t218 * t163 + t219 * t162 + t183 * t222 + t182 * t223 + t202 * t73 - pkin(3) * t198 + t49 * t160 + t50 * t161 + t146 * t137 + t111 * t84 + t112 * t85 + (t285 * t393 + t283 * t395 + pkin(9) * t280 + (t296 * t394 + t298 * t395) * qJD(4) + (pkin(9) * mrSges(4,2) - Ifges(4,6) + t346) * qJD(3)) * t328 + (t190 / 0.2e1 + qJD(4) * t403 + t366 * t396 + t349 * mrSges(5,3) + (m(5) * t349 - qJD(4) * t288 + t229) * pkin(10)) * t331 + (t144 * t264 - t145 * t263 - t274 * t69 - t340 * t71) * mrSges(6,3); -0.2e1 * pkin(3) * t280 + t132 * t140 + t133 * t139 + 0.2e1 * t227 * t137 + 0.2e1 * t314 * t192 - t340 * t194 + t274 * t195 + t199 * t75 + t200 * t76 - t263 * t206 - t264 * t207 + 0.2e1 * t239 * t73 + t331 * t283 + t327 * t285 + (t331 * t298 + (0.2e1 * pkin(4) * t204 - t296) * t327) * qJD(4) + (t182 * t218 + t183 * t219 + t314 * t360) * t425 + (t111 * t50 + t112 * t49 + t227 * t239) * t424 + 0.2e1 * (-t111 * t132 + t112 * t133 + t199 * t49 - t200 * t50) * mrSges(7,3) + 0.2e1 * (-t182 * t274 - t183 * t340 + t218 * t264 - t219 * t263) * mrSges(6,3); t338 + t259 * t17 + t260 * t18 + t243 * t64 - t244 * t65 + m(7) * (-t12 * t244 + t13 * t243 + t2 * t260 + t259 * t3) + (m(6) * (t10 * t324 + t11 * t322) + t322 * t53 + t324 * t54) * pkin(4) - t38 * mrSges(5,2) + t39 * mrSges(5,1) + t10 * mrSges(6,1) - t11 * mrSges(6,2) + t429; t337 + t259 * t84 + t260 * t85 + t243 * t160 - t244 * t161 - t166 * mrSges(5,2) + t167 * mrSges(5,1) + m(7) * (t15 * t260 + t16 * t259 + t243 * t56 - t244 * t55) + (t322 * t162 + m(6) * (t322 * t71 + t324 * t69) + t324 * t163) * pkin(4) - t71 * mrSges(6,2) + t69 * mrSges(6,1) - t428; m(7) * (-t111 * t244 + t112 * t243 + t259 * t50 + t260 * t49) - t183 * mrSges(6,2) + t182 * mrSges(6,1) + t318 + (pkin(10) * t293 - t382) * qJD(4) + (m(6) * (t182 * t324 + t183 * t322) + (-t263 * t322 + t264 * t324) * mrSges(6,3)) * pkin(4) + (-t132 * t259 + t133 * t260 + t199 * t243 + t200 * t244) * mrSges(7,3) + t339 + t193; 0.2e1 * m(7) * (t243 * t260 - t244 * t259) + 0.2e1 * t344; m(6) * t60 + m(7) * t34 + t37 + t9; m(6) * t233 + m(7) * t146 + t114 + t47; m(6) * t360 + m(7) * t227 + t192 + t73; 0; 0; t338; t337; t339; t344; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
