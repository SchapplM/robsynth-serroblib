% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR12_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:06
% EndTime: 2019-12-31 20:29:16
% DurationCPUTime: 4.67s
% Computational Cost: add. (7725->388), mult. (15094->513), div. (0->0), fcn. (14572->6), ass. (0->218)
t245 = cos(qJ(5));
t353 = t245 * mrSges(6,1);
t242 = sin(qJ(5));
t371 = mrSges(6,2) * t242;
t287 = t353 - t371;
t442 = -t287 - mrSges(5,1);
t244 = sin(qJ(2));
t247 = cos(qJ(2));
t404 = -pkin(2) * t247 - qJ(3) * t244;
t213 = -pkin(1) + t404;
t186 = pkin(3) * t247 - t213;
t243 = sin(qJ(4));
t246 = cos(qJ(4));
t204 = t243 * t244 + t246 * t247;
t205 = -t243 * t247 + t244 * t246;
t441 = m(5) * t186 + mrSges(5,1) * t204 + mrSges(5,2) * t205;
t427 = pkin(6) - pkin(7);
t297 = t427 * t244;
t316 = t427 * t247;
t157 = t243 * t297 + t246 * t316;
t399 = t243 * t316 - t246 * t297;
t440 = t399 * mrSges(5,2) + t157 * t442;
t376 = -t242 / 0.2e1;
t374 = t245 / 0.2e1;
t339 = t399 * t243;
t340 = t157 * t246;
t438 = t340 + t339;
t365 = Ifges(6,2) * t245;
t368 = Ifges(6,4) * t242;
t219 = t365 + t368;
t367 = Ifges(6,4) * t245;
t369 = Ifges(6,1) * t242;
t221 = t367 + t369;
t424 = -t245 / 0.2e1;
t425 = t242 / 0.2e1;
t271 = t219 * t425 + t221 * t424;
t267 = -Ifges(5,5) + t271;
t436 = t267 * t204 + t440;
t352 = t245 * mrSges(6,2);
t357 = t242 * mrSges(6,1);
t216 = t352 + t357;
t117 = t216 * t205;
t337 = t204 * t242;
t121 = -mrSges(6,2) * t205 + mrSges(6,3) * t337;
t336 = t204 * t245;
t124 = mrSges(6,1) * t205 + mrSges(6,3) * t336;
t233 = Ifges(6,5) * t245;
t362 = Ifges(6,6) * t242;
t402 = -t362 + t233;
t406 = t216 * t204;
t426 = t204 / 0.2e1;
t89 = pkin(4) * t204 - pkin(8) * t205 + t186;
t45 = -t157 * t242 + t245 * t89;
t46 = t157 * t245 + t242 * t89;
t220 = -Ifges(6,2) * t242 + t367;
t74 = -Ifges(6,6) * t205 + t204 * t220;
t222 = Ifges(6,1) * t245 - t368;
t77 = -Ifges(6,5) * t205 + t204 * t222;
t435 = -(t77 * t374 + t74 * t376 - t186 * mrSges(5,1) - (Ifges(5,2) + Ifges(6,3)) * t426 - (-Ifges(5,4) + t402 / 0.2e1) * t205) * t205 + t157 * t117 + t46 * t121 + t45 * t124 - t399 * t406;
t248 = -pkin(2) - pkin(3);
t211 = -qJ(3) * t243 + t246 * t248;
t209 = pkin(4) - t211;
t433 = t157 * t209;
t431 = t157 * t399;
t239 = t242 ^ 2;
t240 = t245 ^ 2;
t322 = t239 + t240;
t334 = t205 * t242;
t317 = mrSges(6,3) * t334;
t122 = -mrSges(6,2) * t204 - t317;
t333 = t205 * t245;
t125 = mrSges(6,1) * t204 - mrSges(6,3) * t333;
t414 = t122 * t374 + t125 * t376;
t423 = t287 / 0.2e1;
t212 = qJ(3) * t246 + t243 * t248;
t210 = -pkin(8) + t212;
t420 = t210 * t121;
t419 = t210 * t124;
t418 = t212 * t399;
t417 = t242 * t399;
t202 = t243 * t287;
t416 = t245 * t399;
t294 = t322 * t246;
t325 = t245 * t125;
t331 = t242 * t122;
t272 = t331 / 0.2e1 + t325 / 0.2e1;
t412 = Ifges(5,4) * t426;
t409 = -t205 / 0.2e1;
t408 = Ifges(3,4) - Ifges(4,5);
t288 = mrSges(6,3) * (t240 / 0.2e1 + t239 / 0.2e1);
t289 = t246 * mrSges(5,2) - mrSges(6,3) * t294;
t405 = -t243 * mrSges(5,1) - t202 - t289;
t327 = t245 * t121;
t330 = t242 * t124;
t403 = t327 / 0.2e1 - t330 / 0.2e1;
t270 = -t121 * t374 - t124 * t376;
t79 = t204 * Ifges(6,5) + t205 * t222;
t348 = t245 * t79;
t76 = t204 * Ifges(6,6) + t205 * t220;
t355 = t242 * t76;
t389 = Ifges(6,3) / 0.2e1;
t400 = (t389 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t205 + Ifges(5,1) * t409 + t412 - t186 * mrSges(5,2) + t355 / 0.2e1 - t348 / 0.2e1;
t398 = mrSges(6,3) * t322 * t409 - t272;
t397 = m(5) / 0.2e1;
t396 = -m(6) / 0.2e1;
t395 = m(6) / 0.2e1;
t394 = -pkin(4) / 0.2e1;
t393 = mrSges(6,1) / 0.2e1;
t392 = -mrSges(6,2) / 0.2e1;
t390 = -Ifges(6,3) / 0.2e1;
t135 = pkin(4) * t205 + pkin(8) * t204;
t59 = t135 * t245 + t417;
t388 = -t59 / 0.2e1;
t60 = t135 * t242 - t416;
t387 = t60 / 0.2e1;
t386 = t74 / 0.4e1;
t385 = t77 / 0.4e1;
t384 = -t406 / 0.2e1;
t383 = -t117 / 0.2e1;
t382 = -t209 / 0.2e1;
t381 = t209 / 0.2e1;
t380 = -t211 / 0.2e1;
t379 = t216 / 0.2e1;
t218 = Ifges(6,5) * t242 + Ifges(6,6) * t245;
t378 = t218 / 0.2e1;
t377 = t221 / 0.2e1;
t375 = t243 / 0.2e1;
t373 = m(6) * (-0.1e1 + t322) * t246 * t243;
t372 = pkin(4) * t216;
t364 = Ifges(5,6) * t205;
t231 = t247 * qJ(3);
t201 = t244 * t248 + t231;
t90 = -t135 + t201;
t47 = t245 * t90 - t417;
t356 = t242 * t47;
t354 = t242 * t77;
t351 = t245 * t46;
t48 = t242 * t90 + t416;
t350 = t245 * t48;
t349 = t245 * t74;
t314 = t233 / 0.2e1;
t275 = t314 - t362 / 0.2e1;
t215 = -t247 * mrSges(4,1) - t244 * mrSges(4,3);
t298 = m(4) * t213 + t215;
t3 = t48 * t122 + t47 * t125 + t298 * (pkin(2) * t244 - t231) + m(6) * (t45 * t47 + t46 * t48 - t431) + (-pkin(1) * mrSges(3,2) - t213 * mrSges(4,3) + t408 * t247) * t247 + ((-Ifges(5,4) / 0.2e1 + t275) * t204 - t400) * t204 + (-pkin(1) * mrSges(3,1) + t213 * mrSges(4,1) + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t247 - t408 * t244) * t244 - t435 + t441 * t201;
t347 = t3 * qJD(1);
t346 = t59 * t242;
t345 = t60 * t245;
t269 = t287 * t205;
t7 = t46 * t125 - t399 * t269 + (mrSges(6,3) * t351 + t204 * t378 - t271 * t205 + t76 * t374 + t425 * t79) * t205 + (-t317 - t122) * t45;
t344 = t7 * qJD(1);
t18 = (t325 + m(6) * (t242 * t46 + t245 * t45) + t331 - t298 + t441) * t244;
t343 = qJD(1) * t18;
t335 = t205 * t218;
t332 = t209 * t216;
t323 = t246 * t216;
t321 = mrSges(6,3) * t346;
t320 = mrSges(6,3) * t345;
t319 = t117 * t375 + t246 * t414;
t318 = qJD(3) * t373;
t301 = -t219 / 0.2e1 + t222 / 0.2e1;
t300 = t220 / 0.2e1 + t377;
t295 = t322 * t211;
t286 = -t242 * t45 + t351;
t285 = t350 - t356;
t284 = t345 - t346;
t250 = (t210 * t284 + t211 * t286 + t418 + t433) * t396 - t406 * t382 + t212 * t383;
t255 = (pkin(4) * t157 + pkin(8) * t285) * t395 + pkin(4) * t384;
t2 = (-pkin(8) * t121 / 0.2e1 + t122 * t380 - t420 / 0.2e1 - t74 / 0.4e1 + t386) * t245 + (pkin(8) * t124 / 0.2e1 + t211 * t125 / 0.2e1 + t419 / 0.2e1 - t77 / 0.4e1 + t385) * t242 + ((t387 + t48 / 0.2e1) * t245 + (t388 - t47 / 0.2e1) * t242) * mrSges(6,3) + t250 + t255;
t257 = -t211 * mrSges(5,2) + mrSges(6,3) * t295 + t212 * t442;
t21 = -m(6) * (t209 * t212 + t210 * t295) + t257;
t283 = -t2 * qJD(1) - t21 * qJD(2);
t258 = t301 * t242 + t300 * t245;
t44 = t258 - t332;
t262 = t399 * t379 - t355 / 0.4e1 + t348 / 0.4e1;
t253 = (t314 - 0.3e1 / 0.4e1 * t362 + t233 / 0.4e1) * t204 + t262;
t261 = -t220 / 0.4e1 - t221 / 0.4e1 - t369 / 0.4e1 - t367 / 0.2e1;
t273 = -t219 / 0.4e1 + t222 / 0.4e1 - t365 / 0.4e1;
t277 = t392 * t48 + t393 * t47;
t5 = t272 * t210 + (t390 + t210 * t288 + (mrSges(6,1) * t382 + t273) * t245 + (mrSges(6,2) * t381 + t261) * t242) * t205 + t253 + t277;
t282 = -qJD(1) * t5 + qJD(2) * t44;
t13 = ((-t157 + t286) * t395 + t406 / 0.2e1) * t246 + ((t399 + t284) * t395 + t403) * t243 + t319;
t4 = t59 * t125 + m(6) * (t45 * t59 + t46 * t60 + t431) + t60 * t122 + (-t275 * t204 + t400 + t412) * t204 + t435;
t281 = t4 * qJD(1) + t13 * qJD(3);
t251 = (t243 * t285 + t340) * t395 + t438 * t397;
t252 = -m(5) * t438 / 0.2e1 + (t246 * t286 + t339) * t396;
t11 = (t384 - t414) * t246 + (t383 + t270) * t243 + t251 + t252;
t32 = -mrSges(4,3) - m(4) * qJ(3) - m(6) * (t209 * t243 + t210 * t294) - m(5) * (-t211 * t243 + t212 * t246) + t405;
t280 = qJD(1) * t11 + qJD(2) * t32;
t274 = -t352 / 0.2e1 - t357 / 0.2e1;
t129 = (-t216 / 0.2e1 + t274) * t246;
t16 = (t125 * t375 + t244 * t393) * t245 + (t122 * t375 + t244 * t392) * t242 + (t243 * t288 + t246 * t423) * t205;
t279 = -t16 * qJD(1) - t129 * qJD(2);
t276 = mrSges(6,1) * t388 + mrSges(6,2) * t387;
t266 = t272 * pkin(8);
t265 = t274 * t246;
t264 = m(6) * (-pkin(4) * t243 + pkin(8) * t294);
t254 = t202 / 0.2e1 + ((t210 * t322 - t212) * t246 + (t209 + t295) * t243) * t395;
t20 = (mrSges(5,1) + t423) * t243 - t264 / 0.2e1 + t254 + t289;
t263 = qJD(1) * t13 + qJD(2) * t20 + t318;
t128 = (t379 + t274) * t246;
t27 = (t382 + t394) * t216 + (mrSges(6,2) * t380 + t300) * t245 + (mrSges(6,1) * t380 + t301) * t242;
t53 = t258 - t372;
t249 = -pkin(8) * t288 + (mrSges(6,1) * t394 + t273) * t245 + (pkin(4) * mrSges(6,2) / 0.2e1 + t261) * t242;
t9 = -t266 + (t390 + t249) * t205 + t253 + t276;
t259 = t9 * qJD(1) - t27 * qJD(2) - t128 * qJD(3) + t53 * qJD(4);
t256 = t205 * t389 + t262;
t131 = t323 / 0.2e1 + t265;
t130 = -t323 / 0.2e1 + t265;
t28 = t332 / 0.2e1 + t220 * t424 + t222 * t376 + t372 / 0.2e1 + t274 * t211 + t271;
t25 = t264 / 0.2e1 - t202 / 0.2e1 + t254;
t17 = -t246 * t269 / 0.2e1 + (t353 / 0.2e1 - t371 / 0.2e1) * t244 + t398 * t243;
t12 = t13 * qJD(4);
t10 = (m(4) * pkin(6) + mrSges(4,2)) * t247 + (t205 * mrSges(5,3) + t270) * t243 + t251 - t252 + t319 + (-mrSges(5,3) * t204 + t384) * t246;
t8 = t204 * t402 / 0.4e1 - Ifges(6,5) * t336 / 0.2e1 + Ifges(6,6) * t337 / 0.2e1 - t266 + t249 * t205 + t256 - t276;
t6 = t269 * t381 - t256 + t277 + (t377 + t220 / 0.4e1) * t334 + (t219 / 0.2e1 - t222 / 0.4e1) * t333 + (-t402 / 0.4e1 + t275) * t204 + t398 * t210;
t1 = t255 + t321 / 0.2e1 - t320 / 0.2e1 + t364 / 0.2e1 + t354 / 0.4e1 + t349 / 0.4e1 - t335 / 0.4e1 + (t350 / 0.2e1 - t356 / 0.2e1) * mrSges(6,3) + t242 * t385 + t245 * t386 + (Ifges(5,6) / 0.2e1 - t218 / 0.4e1) * t205 - t250 + t414 * t211 + t403 * t210 + t270 * pkin(8) - t436;
t14 = [qJD(2) * t3 + qJD(3) * t18 + qJD(4) * t4 - qJD(5) * t7, t10 * qJD(3) + t1 * qJD(4) + t6 * qJD(5) + t347 + (t209 * t406 + (-t420 - t48 * mrSges(6,3) - t74 / 0.2e1) * t245 + (t419 + t47 * mrSges(6,3) - t77 / 0.2e1) * t242 + 0.2e1 * (t210 * t285 - t433) * t395 + 0.2e1 * (t157 * t211 + t418) * t397 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t247 + (-mrSges(4,2) * qJ(3) - Ifges(3,6) + Ifges(4,6)) * t244 + (t212 * mrSges(5,3) - Ifges(5,6) + t378) * t205 + (-t211 * mrSges(5,3) + t267) * t204 + (m(4) * t404 - t247 * mrSges(3,1) + t244 * mrSges(3,2) + t215) * pkin(6) + t440) * qJD(2), qJD(2) * t10 + qJD(5) * t17 + t12 + t343, t1 * qJD(2) + t8 * qJD(5) + t281 + (-t321 - t354 / 0.2e1 - t349 / 0.2e1 + t335 / 0.2e1 + t320 - t364 + (-m(6) * t157 + t406) * pkin(4) + (m(6) * t284 + t327 - t330) * pkin(8) + t436) * qJD(4), -t344 + t6 * qJD(2) + t17 * qJD(3) + t8 * qJD(4) + (-mrSges(6,1) * t46 - mrSges(6,2) * t45 - t335) * qJD(5); -qJD(3) * t11 - qJD(4) * t2 - qJD(5) * t5 - t347, -qJD(3) * t32 - qJD(4) * t21 + qJD(5) * t44, qJD(4) * t25 + qJD(5) * t131 - t280 + t318, t25 * qJD(3) + (m(6) * (-pkin(4) * t212 + pkin(8) * t295) + t257) * qJD(4) + t28 * qJD(5) + t283, t131 * qJD(3) + t28 * qJD(4) + (-t210 * t287 - t402) * qJD(5) + t282; qJD(2) * t11 - qJD(5) * t16 + t12 - t343, qJD(4) * t20 - qJD(5) * t129 + t280, qJD(4) * t373, (t264 + t405) * qJD(4) + t130 * qJD(5) + t263, t130 * qJD(4) - qJD(5) * t202 + t279; qJD(2) * t2 + qJD(5) * t9 - t281, -qJD(3) * t20 - qJD(5) * t27 - t283, -qJD(5) * t128 - t263, t53 * qJD(5), (-pkin(8) * t287 + t402) * qJD(5) + t259; qJD(2) * t5 + qJD(3) * t16 - qJD(4) * t9 + t344, qJD(3) * t129 + qJD(4) * t27 - t282, qJD(4) * t128 - t279, -t259, 0;];
Cq = t14;
