% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2018-11-23 15:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:27:32
% EndTime: 2018-11-23 15:27:42
% DurationCPUTime: 10.18s
% Computational Cost: add. (7971->553), mult. (20020->729), div. (0->0), fcn. (14503->10), ass. (0->265)
t436 = Ifges(6,4) + Ifges(7,4);
t437 = Ifges(6,1) + Ifges(7,1);
t418 = Ifges(6,5) + Ifges(7,5);
t435 = Ifges(6,2) + Ifges(7,2);
t417 = Ifges(7,6) + Ifges(6,6);
t247 = sin(qJ(4));
t248 = sin(qJ(3));
t251 = cos(qJ(4));
t252 = cos(qJ(3));
t214 = t247 * t248 - t251 * t252;
t243 = qJD(3) + qJD(4);
t176 = t243 * t214;
t215 = t247 * t252 + t248 * t251;
t177 = t243 * t215;
t249 = sin(qJ(2));
t244 = sin(pkin(6));
t333 = qJD(1) * t244;
t313 = t249 * t333;
t327 = qJD(3) * t248;
t323 = pkin(3) * t327;
t440 = pkin(4) * t177 + pkin(10) * t176 - t313 + t323;
t389 = -pkin(9) - pkin(8);
t314 = qJD(3) * t389;
t219 = t248 * t314;
t220 = t252 * t314;
t253 = cos(qJ(2));
t312 = t253 * t333;
t229 = t389 * t248;
t230 = t389 * t252;
t394 = t251 * t229 + t230 * t247;
t427 = qJD(4) * t394 + t214 * t312 + t219 * t251 + t220 * t247;
t250 = cos(qJ(5));
t439 = t436 * t250;
t246 = sin(qJ(5));
t438 = t436 * t246;
t407 = -t246 * t417 + t250 * t418;
t406 = -t246 * t435 + t439;
t403 = t250 * t437 - t438;
t434 = -t246 * t427 + t250 * t440;
t239 = -pkin(3) * t252 - pkin(2);
t167 = pkin(4) * t214 - pkin(10) * t215 + t239;
t325 = qJD(5) * t250;
t433 = t167 * t325 + t246 * t440 + t250 * t427;
t432 = t243 * Ifges(5,6) / 0.2e1;
t190 = t229 * t247 - t230 * t251;
t180 = t250 * t190;
t270 = qJ(6) * t176 - qJD(6) * t215;
t431 = pkin(5) * t177 + t270 * t250 + (-t180 + (qJ(6) * t215 - t167) * t246) * qJD(5) + t434;
t309 = t215 * t325;
t430 = -qJ(6) * t309 + (-qJD(5) * t190 + t270) * t246 + t433;
t326 = qJD(5) * t246;
t429 = -t190 * t326 + t433;
t105 = t246 * t167 + t180;
t428 = -qJD(5) * t105 + t434;
t110 = qJD(4) * t190 + t219 * t247 - t251 * t220;
t185 = t215 * t312;
t426 = t110 - t185;
t208 = t214 * qJD(2);
t209 = t215 * qJD(2);
t165 = mrSges(5,1) * t208 + mrSges(5,2) * t209;
t198 = qJD(2) * t239 - t312;
t425 = -m(5) * t198 - t165;
t187 = -t209 * t246 + t243 * t250;
t424 = t436 * t187;
t345 = t208 * t246;
t423 = -qJ(6) * t345 + t250 * qJD(6);
t188 = t209 * t250 + t243 * t246;
t422 = t436 * t188;
t222 = qJD(2) * pkin(8) + t313;
t299 = pkin(9) * qJD(2) + t222;
t245 = cos(pkin(6));
t332 = qJD(1) * t248;
t308 = t245 * t332;
t182 = t252 * t299 + t308;
t169 = t247 * t182;
t337 = t245 * t252;
t233 = qJD(1) * t337;
t274 = t299 * t248;
t181 = t233 - t274;
t171 = qJD(3) * pkin(3) + t181;
t102 = t171 * t251 - t169;
t352 = t243 * Ifges(5,5);
t421 = t198 * mrSges(5,2) - t102 * mrSges(5,3) + t352 / 0.2e1;
t353 = t209 * Ifges(5,4);
t420 = t432 + t353 / 0.2e1 - t208 * Ifges(5,2) / 0.2e1;
t162 = t176 * qJD(2);
t100 = -qJD(5) * t188 + t162 * t246;
t163 = t177 * qJD(2);
t99 = qJD(5) * t187 - t162 * t250;
t416 = t435 * t100 + t417 * t163 + t436 * t99;
t415 = t436 * t100 + t418 * t163 + t437 * t99;
t202 = qJD(5) + t208;
t398 = t435 * t187 + t417 * t202 + t422;
t397 = t437 * t188 + t418 * t202 + t424;
t170 = t251 * t182;
t106 = t181 * t247 + t170;
t195 = pkin(5) * t345;
t321 = pkin(5) * t326;
t362 = pkin(3) * qJD(4);
t413 = t247 * t362 - t106 + t195 + t321;
t242 = t250 * qJ(6);
t292 = t209 * pkin(5) + t208 * t242;
t236 = pkin(3) * t247 + pkin(10);
t336 = -qJ(6) - t236;
t297 = qJD(5) * t336;
t322 = t251 * t362;
t107 = t181 * t251 - t169;
t166 = pkin(4) * t209 + pkin(10) * t208;
t330 = qJD(2) * t248;
t140 = pkin(3) * t330 + t166;
t40 = -t107 * t246 + t250 * t140;
t412 = -t292 - t40 + (-qJD(6) - t322) * t246 + t250 * t297;
t103 = t171 * t247 + t170;
t331 = qJD(2) * t244;
t303 = qJD(1) * t331;
t296 = t253 * t303;
t334 = qJD(3) * t233 + t252 * t296;
t129 = -qJD(3) * t274 + t334;
t393 = -t182 * qJD(3) - t248 * t296;
t28 = qJD(4) * t103 + t129 * t247 - t251 * t393;
t34 = -mrSges(6,1) * t100 + mrSges(6,2) * t99;
t411 = m(6) * t28 + t34;
t41 = t250 * t107 + t246 * t140;
t410 = t246 * t297 + t250 * t322 - t41 + t423;
t288 = mrSges(7,1) * t246 + mrSges(7,2) * t250;
t290 = mrSges(6,1) * t246 + mrSges(6,2) * t250;
t90 = -pkin(4) * t243 - t102;
t51 = -pkin(5) * t187 + qJD(6) + t90;
t409 = t51 * t288 + t90 * t290;
t408 = t418 * t246 + t417 * t250;
t405 = t435 * t250 + t438;
t404 = t437 * t246 + t439;
t382 = t188 / 0.2e1;
t402 = t409 + t406 * t187 / 0.2e1 + t403 * t382 + t407 * t202 / 0.2e1;
t121 = pkin(4) * t208 - pkin(10) * t209 + t198;
t91 = pkin(10) * t243 + t103;
t36 = t250 * t121 - t246 * t91;
t29 = -qJ(6) * t188 + t36;
t21 = pkin(5) * t202 + t29;
t37 = t121 * t246 + t250 * t91;
t30 = qJ(6) * t187 + t37;
t401 = -t198 * mrSges(5,1) - t36 * mrSges(6,1) - t21 * mrSges(7,1) + t37 * mrSges(6,2) + t30 * mrSges(7,2) + t420;
t400 = -Ifges(4,1) / 0.2e1;
t328 = qJD(2) * t252;
t399 = -Ifges(4,4) * t328 / 0.2e1;
t369 = -qJ(6) - pkin(10);
t300 = qJD(5) * t369;
t43 = -t102 * t246 + t250 * t166;
t396 = -qJD(6) * t246 + t250 * t300 - t292 - t43;
t44 = t250 * t102 + t246 * t166;
t395 = t246 * t300 + t423 - t44;
t392 = -t246 * t36 + t250 * t37;
t390 = t99 / 0.2e1;
t387 = t100 / 0.2e1;
t386 = t163 / 0.2e1;
t385 = -t187 / 0.2e1;
t383 = -t188 / 0.2e1;
t381 = -t202 / 0.2e1;
t379 = t208 / 0.2e1;
t377 = -t246 / 0.2e1;
t374 = t250 / 0.2e1;
t372 = pkin(3) * t251;
t371 = pkin(10) * t250;
t27 = qJD(4) * t102 + t251 * t129 + t247 * t393;
t203 = qJD(2) * t323 + t249 * t303;
t61 = pkin(4) * t163 + pkin(10) * t162 + t203;
t5 = t121 * t325 + t246 * t61 + t250 * t27 - t326 * t91;
t370 = t250 * t5;
t368 = mrSges(5,3) * t209;
t367 = Ifges(4,4) * t248;
t197 = Ifges(5,4) * t208;
t360 = t103 * mrSges(5,3);
t339 = t244 * t249;
t204 = -t248 * t339 + t337;
t205 = t245 * t248 + t252 * t339;
t272 = t251 * t204 - t205 * t247;
t359 = t272 * t28;
t356 = t394 * t28;
t354 = t209 * Ifges(5,1);
t348 = Ifges(4,5) * qJD(3);
t347 = Ifges(4,6) * qJD(3);
t344 = t208 * t250;
t343 = t215 * t246;
t342 = t222 * t252;
t340 = t236 * t250;
t338 = t244 * t253;
t335 = -mrSges(5,1) * t243 - mrSges(6,1) * t187 + mrSges(6,2) * t188 + t368;
t329 = qJD(2) * t249;
t324 = qJD(2) * qJD(3);
t320 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t319 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t318 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t113 = -mrSges(7,1) * t187 + mrSges(7,2) * t188;
t315 = t113 + t335;
t238 = -pkin(5) * t250 - pkin(4);
t311 = t244 * t329;
t310 = t253 * t331;
t307 = t348 / 0.2e1;
t306 = -t347 / 0.2e1;
t33 = -t100 * mrSges(7,1) + t99 * mrSges(7,2);
t104 = t250 * t167 - t190 * t246;
t6 = -qJD(5) * t37 - t246 * t27 + t250 * t61;
t1 = pkin(5) * t163 - qJ(6) * t99 - qJD(6) * t188 + t6;
t295 = -t6 * mrSges(6,3) - t1 * mrSges(7,3);
t294 = -t36 * mrSges(6,3) - t21 * mrSges(7,3);
t293 = -t37 * mrSges(6,3) - t30 * mrSges(7,3);
t291 = mrSges(6,1) * t250 - mrSges(6,2) * t246;
t289 = mrSges(7,1) * t250 - mrSges(7,2) * t246;
t275 = -t246 * t37 - t250 * t36;
t150 = -t222 * t327 + t334;
t151 = -qJD(3) * t342 + (-qJD(3) * t245 - t310) * t332;
t273 = t150 * t252 - t151 * t248;
t146 = t204 * t247 + t205 * t251;
t119 = -t146 * t246 - t250 * t338;
t269 = -t146 * t250 + t246 * t338;
t191 = -t222 * t248 + t233;
t223 = -qJD(2) * pkin(2) - t312;
t266 = t191 * mrSges(4,3) + t330 * t400 + t399 - t348 / 0.2e1 - t223 * mrSges(4,2);
t192 = t308 + t342;
t265 = t192 * mrSges(4,3) + t347 / 0.2e1 + (t252 * Ifges(4,2) + t367) * qJD(2) / 0.2e1 - t223 * mrSges(4,1);
t3 = qJ(6) * t100 + qJD(6) * t187 + t5;
t258 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t256 = m(6) * (qJD(5) * t275 - t6 * t246 + t370);
t10 = -pkin(5) * t100 + t28;
t142 = -t197 + t352 + t354;
t82 = t188 * Ifges(7,5) + t187 * Ifges(7,6) + t202 * Ifges(7,3);
t83 = t188 * Ifges(6,5) + t187 * Ifges(6,6) + t202 * Ifges(6,3);
t255 = -t27 * mrSges(5,2) - Ifges(5,5) * t162 - Ifges(5,6) * t163 - t10 * t289 + t404 * t390 + t405 * t387 + t408 * t386 + (-t197 + t142) * t379 + t415 * t246 / 0.2e1 + t416 * t374 + (-mrSges(5,1) - t291) * t28 + (-t21 * t344 + t250 * t3 - t30 * t345) * mrSges(7,3) + (-t344 * t36 - t345 * t37 + t370) * mrSges(6,3) + (-Ifges(5,2) * t379 + t432 + t417 * t385 + t418 * t383 + (Ifges(7,3) + Ifges(6,3)) * t381 + t401) * t209 + (-t407 * t381 - t403 * t383 - t406 * t385 + t409 + t421) * t208 - (-Ifges(5,1) * t208 - t353 + t82 + t83) * t209 / 0.2e1 + t402 * qJD(5) + (-t345 / 0.2e1 - t326 / 0.2e1) * t398 + (t344 / 0.2e1 + t325 / 0.2e1) * t397;
t254 = qJD(2) ^ 2;
t228 = t242 + t371;
t227 = t369 * t246;
t226 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t328;
t225 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t330;
t224 = t238 - t372;
t213 = t242 + t340;
t212 = t336 * t246;
t211 = (mrSges(4,1) * t248 + mrSges(4,2) * t252) * t324;
t193 = -mrSges(5,2) * t243 - mrSges(5,3) * t208;
t179 = -qJD(3) * t205 - t248 * t310;
t178 = qJD(3) * t204 + t252 * t310;
t159 = Ifges(6,3) * t163;
t158 = Ifges(7,3) * t163;
t132 = pkin(5) * t343 - t394;
t126 = mrSges(6,1) * t202 - mrSges(6,3) * t188;
t125 = mrSges(7,1) * t202 - mrSges(7,3) * t188;
t124 = -mrSges(6,2) * t202 + mrSges(6,3) * t187;
t123 = -mrSges(7,2) * t202 + mrSges(7,3) * t187;
t97 = Ifges(6,5) * t99;
t96 = Ifges(7,5) * t99;
t95 = Ifges(6,6) * t100;
t94 = Ifges(7,6) * t100;
t81 = mrSges(5,1) * t163 - mrSges(5,2) * t162;
t70 = -qJ(6) * t343 + t105;
t62 = t103 - t195;
t52 = pkin(5) * t214 - t215 * t242 + t104;
t50 = qJD(4) * t146 + t178 * t247 - t251 * t179;
t49 = qJD(4) * t272 + t178 * t251 + t179 * t247;
t48 = -mrSges(6,2) * t163 + mrSges(6,3) * t100;
t47 = -mrSges(7,2) * t163 + mrSges(7,3) * t100;
t46 = mrSges(6,1) * t163 - mrSges(6,3) * t99;
t45 = mrSges(7,1) * t163 - mrSges(7,3) * t99;
t42 = (-t246 * t176 + t309) * pkin(5) + t110;
t24 = qJD(5) * t269 - t246 * t49 + t250 * t311;
t23 = qJD(5) * t119 + t246 * t311 + t250 * t49;
t2 = [-t146 * t163 * mrSges(5,3) + t178 * t226 + t179 * t225 + t49 * t193 + (t125 + t126) * t24 + (t123 + t124) * t23 - (t47 + t48) * t269 + (t45 + t46) * t119 + (-t204 * t252 - t205 * t248) * mrSges(4,3) * t324 + t315 * t50 - (-t162 * mrSges(5,3) + t33 + t34) * t272 + ((-mrSges(3,2) * t254 - t211 - t81) * t253 + (-mrSges(3,1) * t254 + (t165 + qJD(2) * (-mrSges(4,1) * t252 + mrSges(4,2) * t248)) * qJD(2)) * t249) * t244 + m(4) * (t150 * t205 + t151 * t204 + t178 * t192 + t179 * t191 + (t223 - t312) * t311) + m(5) * (-t102 * t50 + t103 * t49 - t359 + t146 * t27 + (t198 * t329 - t203 * t253) * t244) + m(6) * (t119 * t6 + t23 * t37 + t24 * t36 - t269 * t5 + t50 * t90 - t359) + m(7) * (t1 * t119 - t10 * t272 + t21 * t24 + t23 * t30 - t269 * t3 + t50 * t51); (t10 * t288 + t203 * mrSges(5,2) - Ifges(5,1) * t162 - Ifges(5,4) * t163 + (mrSges(5,3) + t290) * t28 + (-t1 * t250 - t246 * t3) * mrSges(7,3) + (-t246 * t5 - t250 * t6) * mrSges(6,3) + (t51 * t289 + t90 * t291 + (t21 * t246 - t250 * t30) * mrSges(7,3) - t392 * mrSges(6,3) + t405 * t385 + t404 * t383 + t408 * t381 - t398 * t250 / 0.2e1) * qJD(5) + t403 * t390 + t406 * t387 + t407 * t386 + (qJD(5) * t397 + t416) * t377 + t415 * t374) * t215 + (t162 * t394 - t163 * t190) * mrSges(5,3) + (-t27 * mrSges(5,3) + t96 / 0.2e1 + t94 / 0.2e1 + t158 / 0.2e1 + t97 / 0.2e1 + t95 / 0.2e1 + t159 / 0.2e1 + Ifges(5,4) * t162 + t203 * mrSges(5,1) + t320 * t99 + t319 * t100 + (Ifges(5,2) + t318) * t163 + t258) * t214 + t52 * t45 + (-t360 + t82 / 0.2e1 + t83 / 0.2e1 + t318 * t202 + t320 * t188 + t319 * t187 - t401 - t420) * t177 + t428 * t126 + (t104 * t6 + t105 * t5 + t36 * t428 + t37 * t429 + t426 * t90 - t356) * m(6) + t429 * t124 + ((t225 * t248 - t226 * t252) * t253 - (pkin(2) * t329 + (-t191 * t248 + t192 * t252) * t253) * m(4) + (-t223 * m(4) + t425) * t249) * t333 + (-t102 * t426 + t103 * t427 + t190 * t27 + t198 * t323 + t203 * t239 - t356) * m(5) - (t142 / 0.2e1 + t354 / 0.2e1 + t402 + (-t21 * t250 - t246 * t30) * mrSges(7,3) + t275 * mrSges(6,3) + t397 * t374 + t398 * t377 - t197 / 0.2e1 + t421) * t176 + t430 * t123 + (t1 * t52 + t10 * t132 + t3 * t70 + (-t185 + t42) * t51 + t430 * t30 + t431 * t21) * m(7) + t431 * t125 + t273 * mrSges(4,3) + (0.3e1 / 0.2e1 * t252 ^ 2 - 0.3e1 / 0.2e1 * t248 ^ 2) * Ifges(4,4) * t324 - t394 * t34 + m(4) * ((-t191 * t252 - t192 * t248) * qJD(3) + t273) * pkin(8) + t427 * t193 + t335 * t110 + ((-pkin(8) * t225 - t266 + t307) * t252 + (-pkin(8) * t226 + pkin(3) * t165 + t306 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t328 - t265) * t248) * qJD(3) + t239 * t81 + t42 * t113 + t132 * t33 - pkin(2) * t211 + t70 * t47 + t104 * t46 + t105 * t48 - t315 * t185; (m(5) * (t247 * t27 - t251 * t28) + (t162 * t251 - t163 * t247) * mrSges(5,3) + ((-m(5) * t102 + m(6) * t90 + t335) * t247 + (m(5) * t103 + m(6) * t392 + t124 * t250 - t126 * t246 + t193) * t251) * qJD(4)) * pkin(3) + t255 - m(5) * (-t102 * t106 + t103 * t107) + t209 * t360 + t48 * t340 + ((t307 + t399 + t266) * t252 + (t306 + (t367 / 0.2e1 + (Ifges(4,2) / 0.2e1 + t400) * t252) * qJD(2) + t425 * pkin(3) + t265) * t248) * qJD(2) + t410 * t123 + t411 * (-pkin(4) - t372) + t412 * t125 + (t1 * t212 + t10 * t224 + t21 * t412 + t213 * t3 + t30 * t410 + t413 * t51) * m(7) + t413 * t113 - m(6) * (t106 * t90 + t36 * t40 + t37 * t41) - t335 * t106 + t236 * t256 + t224 * t33 + t192 * t225 - t191 * t226 - t41 * t124 - t40 * t126 - t150 * mrSges(4,2) + t151 * mrSges(4,1) + t212 * t45 + t213 * t47 + ((-t126 * t236 + t294) * t250 + (-t124 * t236 + t293) * t246) * qJD(5) + (-t236 * t46 + t295) * t246 - t107 * t193; -m(6) * (t103 * t90 + t36 * t43 + t37 * t44) + t396 * t125 + t395 * t123 + t255 + pkin(10) * t256 + (-t335 + t368) * t103 + ((-pkin(10) * t126 + t294) * t250 + (pkin(5) * t113 - pkin(10) * t124 + t293) * t246) * qJD(5) + (-pkin(10) * t46 + t295) * t246 + t48 * t371 + t227 * t45 + t228 * t47 + t238 * t33 - t62 * t113 - t44 * t124 - t43 * t126 - t102 * t193 - t411 * pkin(4) + (t1 * t227 + t10 * t238 + t228 * t3 + (t321 - t62) * t51 + t395 * t30 + t396 * t21) * m(7); t96 + t158 + t159 + (t187 * t21 + t188 * t30) * mrSges(7,3) + t94 + t97 + t258 + (-t113 * t188 + t45) * pkin(5) + t95 + (t187 * t36 + t188 * t37) * mrSges(6,3) + (-(-t21 + t29) * t30 + (-t188 * t51 + t1) * pkin(5)) * m(7) - t29 * t123 - t36 * t124 + t30 * t125 + t37 * t126 - t51 * (mrSges(7,1) * t188 + mrSges(7,2) * t187) - t90 * (mrSges(6,1) * t188 + mrSges(6,2) * t187) + (t437 * t187 - t422) * t383 + t398 * t382 + (t187 * t418 - t188 * t417) * t381 + (-t435 * t188 + t397 + t424) * t385; -t187 * t123 + t188 * t125 + 0.2e1 * (t10 / 0.2e1 + t30 * t385 + t21 * t382) * m(7) + t33;];
tauc  = t2(:);
