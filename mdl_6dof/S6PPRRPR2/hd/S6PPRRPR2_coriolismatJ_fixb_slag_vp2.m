% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PPRRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:49:16
% EndTime: 2019-03-08 18:49:24
% DurationCPUTime: 4.60s
% Computational Cost: add. (7325->407), mult. (20500->621), div. (0->0), fcn. (22586->12), ass. (0->240)
t245 = sin(pkin(7));
t386 = cos(qJ(3));
t329 = t245 * t386;
t431 = -t329 / 0.2e1;
t247 = sin(qJ(4));
t250 = cos(qJ(4));
t416 = t247 ^ 2 + t250 ^ 2;
t430 = pkin(9) * t416;
t373 = t247 * mrSges(6,2);
t217 = -t250 * mrSges(6,3) - t373;
t357 = qJ(5) * t250;
t216 = t247 * pkin(4) - t357;
t422 = m(6) * t216;
t429 = t217 + t422;
t428 = t422 / 0.2e1;
t248 = sin(qJ(3));
t360 = sin(pkin(6));
t317 = sin(pkin(12)) * t360;
t302 = cos(pkin(12)) * t360;
t361 = cos(pkin(7));
t362 = cos(pkin(6));
t415 = t362 * t245 + t361 * t302;
t138 = t248 * t415 + t386 * t317;
t258 = -t245 * t302 + t361 * t362;
t104 = t138 * t250 + t247 * t258;
t384 = m(7) * t104;
t350 = t245 * t248;
t188 = t247 * t361 + t250 * t350;
t427 = m(7) * t188;
t103 = t138 * t247 - t250 * t258;
t249 = cos(qJ(6));
t137 = t248 * t317 - t386 * t415;
t246 = sin(qJ(6));
t66 = t103 * t249 - t137 * t246;
t367 = t249 * t66;
t67 = t103 * t246 + t137 * t249;
t377 = t246 * t67;
t426 = -t103 + t367 + t377;
t187 = t247 * t350 - t250 * t361;
t327 = t249 * t386;
t144 = -t246 * t187 + t245 * t327;
t352 = t144 * t246;
t143 = t249 * t187 + t246 * t329;
t353 = t143 * t249;
t425 = -t187 - t352 + t353;
t197 = pkin(10) * t247 + t216;
t401 = pkin(5) + pkin(9);
t225 = t401 * t250;
t141 = -t197 * t246 + t225 * t249;
t142 = t197 * t249 + t225 * t246;
t218 = t247 * mrSges(5,1) + t250 * mrSges(5,2);
t251 = -pkin(4) - pkin(10);
t358 = qJ(5) * t247;
t194 = t250 * t251 - pkin(3) - t358;
t224 = t401 * t247;
t139 = -t194 * t246 + t224 * t249;
t140 = t194 * t249 + t224 * t246;
t369 = t249 * mrSges(7,1);
t378 = t246 * mrSges(7,2);
t212 = t369 - t378;
t189 = t212 * t247;
t345 = t246 * t250;
t374 = t247 * mrSges(7,1);
t200 = mrSges(7,3) * t345 + t374;
t340 = t249 * t250;
t372 = t247 * mrSges(7,2);
t202 = -mrSges(7,3) * t340 - t372;
t387 = t249 / 0.2e1;
t389 = t246 / 0.2e1;
t405 = m(7) / 0.2e1;
t256 = (t139 * t249 + t140 * t246 - t224) * t405 + t200 * t387 + t202 * t389 - t189 / 0.2e1;
t344 = t247 * t249;
t363 = t250 * mrSges(7,2);
t201 = mrSges(7,3) * t344 - t363;
t394 = -t201 / 0.2e1;
t346 = t246 * t247;
t365 = t250 * mrSges(7,1);
t199 = -mrSges(7,3) * t346 + t365;
t397 = t199 / 0.2e1;
t190 = t212 * t250;
t399 = -t190 / 0.2e1;
t424 = t143 * t397 + t144 * t394 + t187 * t399 + (t141 * t143 - t142 * t144 - t187 * t225) * t405 + (t218 + t429) * t431 + t256 * t188;
t423 = t416 * (mrSges(6,1) + mrSges(5,3)) - mrSges(4,2);
t368 = t249 * mrSges(7,2);
t379 = t246 * mrSges(7,1);
t215 = t368 + t379;
t420 = mrSges(6,3) + t215;
t388 = -t249 / 0.2e1;
t390 = -t246 / 0.2e1;
t419 = Ifges(7,5) * t390 + Ifges(7,6) * t388 + Ifges(5,4) + Ifges(6,6);
t213 = t250 * mrSges(6,2) - t247 * mrSges(6,3);
t417 = -t250 * mrSges(5,1) + t247 * mrSges(5,2) + t213;
t240 = t246 ^ 2;
t242 = t249 ^ 2;
t336 = t240 + t242;
t282 = -t378 / 0.2e1 + t369 / 0.2e1;
t414 = t282 + t212 / 0.2e1;
t343 = t249 * t141;
t349 = t246 * t142;
t291 = t343 + t349;
t413 = m(6) / 0.4e1 + m(5) / 0.4e1;
t411 = 2 * qJD(3);
t410 = 2 * qJD(4);
t409 = m(5) / 0.2e1;
t407 = m(6) / 0.2e1;
t404 = m(7) / 0.4e1;
t403 = mrSges(7,1) / 0.2e1;
t402 = -mrSges(7,2) / 0.2e1;
t191 = t250 * t215;
t398 = -t191 / 0.2e1;
t396 = -t200 / 0.2e1;
t395 = t200 / 0.2e1;
t393 = t202 / 0.2e1;
t391 = t225 / 0.2e1;
t382 = mrSges(7,3) * t250;
t381 = Ifges(7,4) * t246;
t380 = Ifges(7,4) * t249;
t81 = -t137 * t346 + t138 * t249;
t376 = t246 * t81;
t371 = t247 * Ifges(7,5);
t370 = t247 * Ifges(7,6);
t80 = -t137 * t344 - t138 * t246;
t366 = t249 * t80;
t359 = qJ(5) * t187;
t356 = t103 * qJ(5);
t355 = t137 * t247;
t354 = t137 * t250;
t351 = t188 * t250;
t348 = t246 * t200;
t347 = t246 * t201;
t342 = t249 * t199;
t341 = t249 * t202;
t319 = t240 / 0.2e1 + t242 / 0.2e1;
t50 = (t372 / 0.2e1 - t202 / 0.2e1) * t249 + (t374 / 0.2e1 + t395) * t246 - t319 * t382;
t339 = t50 * qJD(3);
t338 = t137 * t430;
t337 = t329 * t430;
t333 = -mrSges(4,1) + t417;
t331 = mrSges(7,3) * t390;
t330 = mrSges(7,3) * t388;
t328 = t247 * t386;
t326 = t250 * t386;
t325 = -t355 / 0.2e1;
t220 = -Ifges(7,2) * t246 + t380;
t301 = Ifges(7,1) * t246 + t380;
t321 = t301 / 0.4e1 + t220 / 0.4e1;
t222 = Ifges(7,1) * t249 - t381;
t300 = Ifges(7,2) * t249 + t381;
t320 = -t222 / 0.4e1 + t300 / 0.4e1;
t318 = m(7) * t336;
t315 = t336 * t251;
t313 = t245 ^ 2 * t386 * t248;
t312 = qJ(5) * t326;
t311 = t104 * t326;
t310 = t187 * t328;
t309 = t188 * t326;
t304 = t240 / 0.4e1 + t242 / 0.4e1 - 0.1e1 / 0.4e1;
t303 = qJD(4) * (mrSges(5,2) - t420);
t299 = -pkin(4) * t250 - t358;
t160 = (-t246 * t248 + t247 * t327) * t245;
t161 = (t246 * t328 + t248 * t249) * t245;
t4 = (t143 * t80 - t144 * t81 + t160 * t66 + t161 * t67) * t405 + 0.2e1 * (-m(7) * t351 / 0.4e1 + t413 * (-t187 * t247 - t351)) * t137 + 0.2e1 * (t311 * t404 + t413 * (t103 * t328 + t137 * t248 - t138 * t386 + t311)) * t245;
t68 = t104 * t354;
t94 = t137 * t138;
t8 = m(6) * (t94 + (-t103 * t247 - t104 * t250) * t137) + m(7) * (t66 * t80 + t67 * t81 - t68) + m(5) * (-t103 * t355 - t68 + t94);
t298 = t8 * qJD(1) + t4 * qJD(2);
t296 = t366 + t376;
t11 = t426 * t384;
t7 = (t104 * t425 + t188 * t426) * t405;
t294 = t11 * qJD(1) + t7 * qJD(2);
t33 = t425 * t427;
t293 = t7 * qJD(1) + t33 * qJD(2);
t149 = t245 * t309;
t34 = m(7) * (t143 * t160 - t144 * t161 + t149) + m(5) * (t245 * t310 + t149 - t313) + m(6) * (-t313 + (t310 + t309) * t245);
t292 = t4 * qJD(1) + t34 * qJD(2);
t289 = t249 * t160 + t246 * t161;
t286 = t319 * t251 * mrSges(7,3);
t285 = t402 * t81 + t403 * t80;
t284 = -t141 * mrSges(7,1) / 0.2e1 + t142 * mrSges(7,2) / 0.2e1;
t283 = t160 * t403 + t161 * t402;
t281 = t428 + t218 / 0.2e1 + t217 / 0.2e1;
t279 = m(7) * t296;
t277 = qJ(5) * t398 + t212 * t391;
t172 = -t250 * t300 + t370;
t193 = t250 * t222;
t276 = -t172 / 0.4e1 - t193 / 0.4e1 + t251 * t393;
t174 = -t250 * t301 + t371;
t192 = t250 * t220;
t275 = -t174 / 0.4e1 + t192 / 0.4e1 + t251 * t396;
t274 = m(7) * t289;
t271 = (t246 * t66 - t249 * t67) * t405;
t270 = -t279 / 0.2e1;
t269 = qJD(4) * (-mrSges(7,3) * t336 - mrSges(5,1) + mrSges(6,2));
t171 = Ifges(7,6) * t250 + t247 * t300;
t173 = Ifges(7,5) * t250 + t247 * t301;
t209 = -pkin(3) + t299;
t12 = m(7) * (t139 * t141 + t140 * t142 - t224 * t225) + t216 * t213 - pkin(3) * t218 - t224 * t190 - t225 * t189 + t139 * t199 + t141 * t200 + t140 * t201 + t142 * t202 + t429 * t209 + (t172 * t387 + t174 * t389 - t419 * t247) * t247 + (t171 * t388 + t173 * t390 + (Ifges(7,3) + Ifges(5,1) + Ifges(6,2) - Ifges(5,2) - Ifges(6,3)) * t247 + t419 * t250) * t250;
t253 = t256 * t104 + (-t103 * t225 + t141 * t66 + t142 * t67) * t405 + t103 * t399 + t66 * t397 + t67 * t201 / 0.2e1;
t3 = (t376 / 0.2e1 + t366 / 0.2e1) * mrSges(7,3) + t251 * t270 + ((-mrSges(5,1) / 0.2e1 + mrSges(6,2) / 0.2e1) * t247 + (t215 / 0.2e1 - mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1) * t250 - t422 / 0.2e1 + t357 * t405 + t281) * t137 + t253;
t252 = (t245 * t312 + t251 * t289) * t405 + (-pkin(4) * t328 + t312) * t245 * t407 + t160 * t330 + t161 * t331 + t218 * t431 + (t250 * t420 + t373) * t329 / 0.2e1;
t9 = -t252 + t424;
t267 = t3 * qJD(1) + t9 * qJD(2) + t12 * qJD(3);
t266 = (t143 * t246 + t144 * t249) * t405;
t265 = -t212 / 0.2e1 + t282;
t254 = (-t352 / 0.2e1 + t353 / 0.2e1) * t382 + t143 * t393 + t144 * t395 + t188 * t398;
t22 = t254 - t283;
t234 = Ifges(7,6) * t345;
t26 = -t140 * t200 + t247 * t234 / 0.2e1 + t139 * t202 - t225 * t191 + ((t140 * mrSges(7,3) + t193 / 0.2e1 + t172 / 0.2e1) * t246 + (t139 * mrSges(7,3) - t371 / 0.2e1 - t174 / 0.2e1 + t192 / 0.2e1) * t249) * t250;
t255 = (t377 / 0.2e1 + t367 / 0.2e1) * t382 + t104 * t398 + t66 * t393 + t67 * t396;
t5 = t255 - t285;
t264 = t5 * qJD(1) + t22 * qJD(2) + t26 * qJD(3);
t25 = t247 * t271 + t270;
t48 = (t348 - t341 + m(7) * (t139 * t246 - t140 * t249) - m(6) * t209 - t213) * t247;
t53 = -t274 / 0.2e1 + t247 * t266;
t263 = qJD(1) * t25 + qJD(2) * t53 + qJD(3) * t48;
t13 = t265 * t104;
t21 = (-Ifges(7,3) / 0.2e1 + t286) * t250 + (-0.3e1 / 0.4e1 * t370 + t320 * t250 + t276) * t249 + (-0.3e1 / 0.4e1 * t371 + t321 * t250 + t275) * t246 + t277 + t284;
t37 = t265 * t188;
t91 = qJ(5) * t212 + (-t301 / 0.2e1 - t220 / 0.2e1) * t249 + (-t222 / 0.2e1 + t300 / 0.2e1) * t246;
t260 = t13 * qJD(1) + t37 * qJD(2) - t21 * qJD(3) - t91 * qJD(4);
t186 = (m(6) + m(7)) * qJ(5) + t420;
t44 = 0.2e1 * t304 * t384;
t52 = (t365 / 0.2e1 - t199 / 0.2e1) * t249 + (-t363 / 0.2e1 + t394) * t246 + 0.2e1 * (t225 / 0.4e1 - t343 / 0.4e1 - t349 / 0.4e1) * m(7);
t93 = -0.2e1 * t304 * t427;
t259 = qJD(1) * t44 - qJD(2) * t93 - qJD(3) * t52 - qJD(4) * t186;
t69 = (t318 / 0.2e1 + 0.2e1 * t404 + 0.2e1 * t407) * t188;
t51 = t341 / 0.2e1 - t348 / 0.2e1 + (t379 / 0.2e1 + t368 / 0.2e1) * t247 + t336 * t382 / 0.2e1;
t49 = t274 / 0.2e1 + (m(6) * t329 + t266) * t247;
t47 = m(7) * t391 + t347 / 0.2e1 + t342 / 0.2e1 + t291 * t405 + (m(6) * pkin(9) + mrSges(6,1) + t282) * t250;
t38 = t414 * t188;
t36 = t384 / 0.2e1 + 0.2e1 * (t407 + t318 / 0.4e1) * t104;
t24 = m(6) * t325 + t279 / 0.2e1 + (-m(6) * t137 / 0.2e1 + t271) * t247;
t23 = t254 + t283;
t20 = Ifges(7,5) * t346 / 0.2e1 + Ifges(7,6) * t344 / 0.2e1 + (-t370 / 0.4e1 + t276) * t249 + (-t371 / 0.4e1 + t275) * t246 + t277 - t284 + (Ifges(7,3) / 0.2e1 + t321 * t246 + t320 * t249 + t286) * t250;
t14 = t414 * t104;
t10 = t252 + t424;
t6 = t255 + t285;
t2 = t296 * t251 * t405 + t81 * t331 + t80 * t330 + mrSges(5,1) * t355 / 0.2e1 + mrSges(6,2) * t325 + t253 + (t428 + t281) * t137 + (-qJ(5) * t405 + mrSges(5,2) / 0.2e1 - t420 / 0.2e1) * t354;
t1 = t4 * qJD(3) + t7 * qJD(4);
t15 = [t8 * qJD(3) + t11 * qJD(4), t1, t2 * qJD(4) + t24 * qJD(5) + t6 * qJD(6) + ((t138 * t209 - t338) * t407 + (t139 * t80 + t140 * t81 - t225 * t354) * t405 + (-pkin(3) * t138 - t338) * t409) * t411 + t298 + (t80 * t200 + t81 * t202 + t333 * t138 + (-t250 * t190 - t423) * t137) * qJD(3), t2 * qJD(3) + t36 * qJD(5) + t14 * qJD(6) + t103 * t303 + t104 * t269 + ((-pkin(4) * t104 - t356) * t407 + (t104 * t315 - t356) * t405) * t410 + t294, qJD(3) * t24 + qJD(4) * t36, t6 * qJD(3) + t14 * qJD(4) + (-mrSges(7,1) * t67 - mrSges(7,2) * t66) * qJD(6); t1, t34 * qJD(3) + t33 * qJD(4), t10 * qJD(4) + t49 * qJD(5) + t23 * qJD(6) + ((t225 * t245 * t326 + t139 * t160 + t140 * t161) * t405 + (-pkin(3) * t350 + t337) * t409 + (t209 * t350 + t337) * t407) * t411 + t292 + (t160 * t200 + t161 * t202 + (t190 * t326 + t333 * t248 + t386 * t423) * t245) * qJD(3), t10 * qJD(3) + t69 * qJD(5) + t38 * qJD(6) + t187 * t303 + t188 * t269 + ((t188 * t315 - t359) * t405 + (-pkin(4) * t188 - t359) * t407) * t410 + t293, qJD(3) * t49 + qJD(4) * t69, t23 * qJD(3) + t38 * qJD(4) + (mrSges(7,1) * t144 - mrSges(7,2) * t143) * qJD(6); qJD(4) * t3 + qJD(5) * t25 + qJD(6) * t5 - t298, qJD(4) * t9 + qJD(5) * t53 + qJD(6) * t22 - t292, qJD(4) * t12 + qJD(5) * t48 + qJD(6) * t26, t47 * qJD(5) + t20 * qJD(6) + t267 + (t173 * t387 + t171 * t390 - t224 * t215 + (-pkin(4) * mrSges(6,1) + Ifges(7,5) * t387 + Ifges(7,6) * t390 - Ifges(6,4) + Ifges(5,5)) * t250 + (t220 * t387 + t222 * t389 + Ifges(6,5) - Ifges(5,6)) * t247 + (m(6) * t299 + t417) * pkin(9) + (-m(7) * t224 - t247 * mrSges(6,1) - t189) * qJ(5) + (m(7) * t291 + t342 + t347) * t251 - t291 * mrSges(7,3)) * qJD(4), qJD(4) * t47 + qJD(6) * t51 + t263, t20 * qJD(4) + t51 * qJD(5) + (-mrSges(7,1) * t140 - mrSges(7,2) * t139 - Ifges(7,5) * t340 + t234) * qJD(6) + t264; -qJD(3) * t3 - qJD(5) * t44 - qJD(6) * t13 - t294, -qJD(3) * t9 + qJD(5) * t93 - qJD(6) * t37 - t293, qJD(5) * t52 + qJD(6) * t21 - t267, qJD(5) * t186 + qJD(6) * t91, -t259 ((-mrSges(7,2) * t251 - Ifges(7,6)) * t249 + (-mrSges(7,1) * t251 - Ifges(7,5)) * t246) * qJD(6) - t260; -qJD(3) * t25 + qJD(4) * t44, -qJD(3) * t53 - qJD(4) * t93, -qJD(4) * t52 - qJD(6) * t50 - t263, t259, 0, -t215 * qJD(6) - t339; -t5 * qJD(3) + t13 * qJD(4), -t22 * qJD(3) + t37 * qJD(4), -qJD(4) * t21 + qJD(5) * t50 - t264, t260, t339, 0;];
Cq  = t15;
