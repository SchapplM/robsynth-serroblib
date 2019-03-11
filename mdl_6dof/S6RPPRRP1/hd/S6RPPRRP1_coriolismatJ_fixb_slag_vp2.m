% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:41
% EndTime: 2019-03-09 01:57:49
% DurationCPUTime: 5.25s
% Computational Cost: add. (9959->362), mult. (19519->497), div. (0->0), fcn. (20327->8), ass. (0->211)
t235 = sin(qJ(5));
t360 = Ifges(6,4) + Ifges(7,4);
t425 = t360 * t235;
t426 = Ifges(6,1) + Ifges(7,1);
t421 = Ifges(6,2) + Ifges(7,2);
t230 = t235 ^ 2;
t237 = cos(qJ(5));
t231 = t237 ^ 2;
t333 = t231 + t230;
t412 = mrSges(6,3) + mrSges(7,3);
t424 = t333 * t412;
t354 = Ifges(7,4) * t237;
t356 = Ifges(6,4) * t237;
t423 = -t235 * t421 + t354 + t356;
t422 = t237 * t426 - t425;
t411 = Ifges(6,5) + Ifges(7,5);
t410 = Ifges(7,6) + Ifges(6,6);
t379 = -t235 / 0.2e1;
t232 = sin(pkin(10));
t233 = cos(pkin(10));
t236 = sin(qJ(4));
t373 = cos(qJ(4));
t218 = t232 * t236 - t233 * t373;
t220 = t232 * t373 + t236 * t233;
t342 = t220 * t235;
t323 = mrSges(7,3) * t342;
t269 = -t218 * mrSges(7,2) - t323;
t420 = t269 / 0.2e1;
t325 = mrSges(6,3) * t342;
t270 = -t218 * mrSges(6,2) - t325;
t341 = t220 * t237;
t324 = mrSges(6,3) * t341;
t272 = t218 * mrSges(6,1) - t324;
t322 = mrSges(7,3) * t341;
t271 = t218 * mrSges(7,1) - t322;
t335 = t237 * t420 + t271 * t379;
t375 = t237 / 0.2e1;
t377 = t235 / 0.2e1;
t275 = -cos(pkin(9)) * pkin(1) - pkin(3) * t233 - pkin(2);
t368 = pkin(8) * t220;
t129 = pkin(4) * t218 + t275 - t368;
t299 = sin(pkin(9)) * pkin(1) + qJ(3);
t274 = pkin(7) + t299;
t205 = t274 * t233;
t260 = t232 * t274;
t137 = t205 * t373 - t236 * t260;
t58 = t237 * t129 - t235 * t137;
t259 = -qJ(6) * t341 + t58;
t41 = t218 * pkin(5) + t259;
t251 = t259 - t41;
t413 = m(7) * t251;
t419 = t270 * t375 + t272 * t379 + t377 * t413 + t335;
t418 = t423 * t220;
t417 = t422 * t220;
t415 = t218 / 0.2e1;
t380 = t220 / 0.2e1;
t361 = mrSges(6,2) + mrSges(7,2);
t409 = Ifges(6,3) + Ifges(7,3);
t297 = mrSges(6,1) * t237 - t235 * mrSges(6,2);
t408 = -mrSges(5,1) - t297;
t386 = m(7) * pkin(5);
t318 = mrSges(7,1) + t386;
t350 = t237 * mrSges(7,2);
t352 = t235 * mrSges(7,1);
t268 = t350 / 0.2e1 + t352 / 0.2e1;
t407 = t268 * t218;
t295 = mrSges(7,1) * t237 - t235 * mrSges(7,2);
t266 = t295 * t220;
t151 = t333 * t220;
t406 = t360 * t237;
t358 = -qJ(6) - pkin(8);
t222 = t358 * t235;
t225 = t358 * t237;
t405 = -t222 * t235 - t225 * t237;
t403 = t410 * t218 + t418;
t402 = t411 * t218 + t417;
t401 = t235 * t270 + t237 * t272;
t400 = -t410 * t235 + t411 * t237;
t399 = t421 * t237 + t425;
t290 = Ifges(7,1) * t235 + t354;
t292 = Ifges(6,1) * t235 + t356;
t398 = t290 + t292;
t294 = t350 + t352;
t149 = t294 * t220;
t136 = t205 * t236 + t373 * t260;
t83 = pkin(5) * t342 + t136;
t396 = -m(7) * t83 - t149;
t319 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t338 = t235 * t218;
t369 = pkin(8) * t218;
t370 = pkin(4) * t220;
t159 = t369 + t370;
t64 = -t237 * t136 + t235 * t159;
t51 = qJ(6) * t338 + t64;
t395 = -t64 * mrSges(6,3) - t51 * mrSges(7,3) + t423 * t415 + (t319 - t410 / 0.2e1) * t220;
t343 = t218 * t237;
t156 = t220 * mrSges(7,1) + mrSges(7,3) * t343;
t363 = t64 * mrSges(6,2);
t63 = t235 * t136 + t237 * t159;
t364 = t63 * mrSges(6,1);
t365 = t51 * mrSges(7,2);
t47 = t220 * pkin(5) + qJ(6) * t343 + t63;
t366 = t47 * mrSges(7,1);
t388 = m(7) / 0.2e1;
t394 = t366 / 0.2e1 - t365 / 0.2e1 + t364 / 0.2e1 - t363 / 0.2e1 + (t47 * t388 + t156 / 0.2e1) * pkin(5);
t393 = 0.2e1 * m(7);
t392 = t220 ^ 2;
t391 = 2 * qJD(4);
t390 = m(6) / 0.2e1;
t389 = m(6) / 0.4e1;
t385 = t149 / 0.2e1;
t381 = -t220 / 0.4e1;
t372 = m(7) * t151;
t367 = t235 * pkin(5);
t362 = mrSges(6,1) + mrSges(7,1);
t353 = t235 * mrSges(6,1);
t351 = t237 * mrSges(6,2);
t349 = t237 * t41;
t267 = t297 * t220;
t254 = -t267 / 0.2e1;
t311 = t341 / 0.2e1;
t313 = t342 / 0.2e1;
t328 = t392 / 0.2e1;
t7 = t269 * t313 + t271 * t311 + t401 * t380 + (t254 - t266 / 0.2e1) * t218 + t328 * t424;
t348 = t7 * qJD(1);
t206 = t218 * mrSges(5,2);
t227 = -pkin(5) * t237 - pkin(4);
t339 = t227 * t220;
t273 = (-t333 * t369 - t370) * t390 + (-t218 * t405 + t339) * t388;
t240 = t273 + t412 * t218 * (-t231 / 0.2e1 - t230 / 0.2e1);
t241 = -m(6) * (t235 * t64 + t237 * t63) / 0.2e1 - m(7) * (t235 * t51 + t237 * t47) / 0.2e1;
t305 = -t297 / 0.2e1 - t295 / 0.2e1;
t157 = t220 * mrSges(6,1) + mrSges(6,3) * t343;
t306 = -t157 / 0.2e1 - t156 / 0.2e1;
t154 = -mrSges(7,2) * t220 + mrSges(7,3) * t338;
t155 = -mrSges(6,2) * t220 + mrSges(6,3) * t338;
t307 = t155 / 0.2e1 + t154 / 0.2e1;
t11 = t206 + t306 * t237 - t307 * t235 + (-mrSges(5,1) + t305) * t220 + t240 + t241;
t347 = qJD(1) * t11;
t345 = t137 * t237;
t48 = t345 + (-qJ(6) * t220 + t129) * t235;
t21 = (m(7) * (-t235 * t48 - t349) - t295 * t218 + mrSges(7,3) * t151) * t220;
t346 = qJD(1) * t21;
t344 = t218 * t220;
t301 = (m(7) / 0.4e1 + t389) * (0.1e1 - t333) * t344;
t25 = 0.2e1 * t301;
t337 = t25 * qJD(1);
t76 = (-t151 / 0.4e1 + t381) * t393;
t336 = t76 * qJD(1);
t308 = t338 / 0.2e1;
t314 = t343 / 0.2e1;
t334 = mrSges(7,1) * t308 + mrSges(7,2) * t314;
t332 = qJD(4) * t218;
t331 = qJD(5) * t220;
t330 = qJD(5) * t235;
t329 = qJD(5) * t237;
t327 = pkin(5) * t341;
t320 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t317 = t351 / 0.2e1;
t316 = -t344 / 0.2e1;
t315 = -t343 / 0.2e1;
t312 = -t341 / 0.2e1;
t302 = mrSges(7,3) * pkin(5) - t411;
t300 = t308 * t386;
t84 = -pkin(5) * t338 + t137;
t298 = t84 * t388 - mrSges(7,1) * t338 / 0.2e1 + mrSges(7,2) * t315;
t296 = t351 + t353;
t147 = t218 * t294;
t148 = t218 * t296;
t150 = t296 * t220;
t244 = -t63 * mrSges(6,3) - t47 * mrSges(7,3) + t320 * t220 - t422 * t218 / 0.2e1 + t411 * t380;
t59 = t235 * t129 + t345;
t1 = -t83 * t147 - t136 * t148 + t84 * t149 + t137 * t150 + t48 * t154 + t59 * t155 + t41 * t156 + t58 * t157 - t275 * t206 + m(6) * (t136 * t137 + t58 * t63 + t59 * t64) + m(7) * (t41 * t47 + t48 * t51 + t83 * t84) + (t275 * mrSges(5,1) - Ifges(5,4) * t220 + t395 * t235 + t244 * t237) * t220 + (-t365 - t363 + t366 + t364 + (-Ifges(5,1) + Ifges(5,2) + (-Ifges(6,1) / 0.2e1 - Ifges(7,1) / 0.2e1) * t231 + ((-Ifges(6,2) / 0.2e1 - Ifges(7,2) / 0.2e1) * t235 + t406) * t235 + t409) * t220 + (Ifges(5,4) - t400) * t218) * t218;
t279 = -t235 * t63 + t237 * t64;
t280 = t235 * t58 - t237 * t59;
t281 = t235 * t41 - t237 * t48;
t6 = (t150 / 0.2e1 + t385 + t307 * t237 + t306 * t235 + (-t235 * t47 + t237 * t51 + t83) * t388 + (t136 + t279) * t390) * t220 + (-t148 / 0.2e1 - t147 / 0.2e1 + (t317 + t353 / 0.2e1 + t268) * t218 + (t281 + t84) * t388 + (t137 + t280) * t390) * t218;
t283 = t1 * qJD(1) + t6 * qJD(2);
t5 = -t136 * t267 - t259 * t269 - t83 * t266 - t41 * t323 + (t324 + t272) * t59 + (-t270 - t325) * t58 + t396 * t327 + (t271 + t322 - t413) * t48 + t399 * t392 * t379 + t402 * t313 + t403 * t311 - t411 * t316 * t235 + (-t316 * t410 + t328 * t398) * t237;
t282 = -t5 * qJD(1) - t7 * qJD(2);
t24 = 0.4e1 * t301;
t278 = t6 * qJD(1) + t24 * qJD(2);
t8 = (mrSges(5,3) * t220 + t150 + 0.4e1 * (t389 + m(5) / 0.4e1) * t136 - t396) * t220 + ((t235 * t362 + t237 * t361 + mrSges(5,3)) * t218 + m(7) * t281 + m(6) * t280 - m(5) * t137) * t218 + (m(4) * t299 + mrSges(4,3)) * (t232 ^ 2 + t233 ^ 2);
t277 = qJD(1) * t8 + qJD(2) * t25;
t208 = t235 * t318 + t350;
t200 = m(7) * t327;
t87 = -t200 - t266;
t276 = qJD(1) * t87 - qJD(4) * t208;
t9 = mrSges(6,1) * t308 + mrSges(6,2) * t314 + t300 + t334 - t419;
t261 = t9 * qJD(1);
t130 = m(7) * t405 + t333 * mrSges(7,3);
t242 = m(7) * ((-t222 * t237 + t225 * t235) * t220 - t281);
t18 = -t242 / 0.2e1 + t407 + t298;
t78 = (t151 / 0.4e1 + t381) * t393;
t252 = -qJD(1) * t18 + qJD(2) * t78 + qJD(4) * t130;
t26 = -t227 * t294 + t295 * t367 + (pkin(4) * mrSges(6,2) - t406) * t237 + (-t227 * t386 + pkin(4) * mrSges(6,1) + t425 + (t421 - t426) * t237) * t235;
t27 = t334 - t407;
t238 = pkin(4) * t254 + t225 * t271 / 0.2e1 + (-t251 * t225 + (t83 * t235 + t237 * t339) * pkin(5)) * t388 + t367 * t385 - pkin(5) * t295 * t311 + t222 * t420 + t227 * t266 / 0.2e1 + t83 * t294 / 0.2e1 + t136 * t296 / 0.2e1 - t401 * pkin(8) / 0.2e1 + t400 * t218 / 0.4e1 + (-t292 / 0.2e1 - t290 / 0.2e1) * t342 + t399 * t312 - t333 * mrSges(6,3) * t368 / 0.2e1 - (t403 + t418) * t235 / 0.4e1 + (t402 + t417) * t237 / 0.4e1 + (t259 * t375 - t225 * t312 + t222 * t313 - t349 / 0.2e1) * mrSges(7,3);
t3 = t308 * t410 + t315 * t411 + t380 * t409 - t238 + t394;
t245 = t3 * qJD(1) + t27 * qJD(2) + t26 * qJD(4);
t243 = (t317 + (mrSges(6,1) / 0.2e1 + t386 / 0.2e1) * t235) * t218 + t334;
t201 = m(7) * t380;
t77 = t372 / 0.2e1 + t201;
t75 = -t372 / 0.2e1 + t201;
t28 = t243 + t300 + (t294 + t296) * t415;
t19 = t242 / 0.2e1 + t298 + t335;
t12 = t305 * t220 + t240 - t241 + (t154 + t155) * t377 + (t156 + t157) * t375;
t10 = t243 + t419;
t4 = t238 + (-t235 * t319 - t237 * t320) * t218 + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t220 + t394;
t2 = qJD(3) * t25 + qJD(4) * t6 - qJD(5) * t7;
t13 = [qJD(3) * t8 + qJD(4) * t1 - qJD(5) * t5 + qJD(6) * t21, t2, qJD(4) * t12 + qJD(5) * t10 + qJD(6) * t75 + t277, t12 * qJD(3) + t4 * qJD(5) + t19 * qJD(6) + ((-pkin(4) * t137 + pkin(8) * t279) * t390 + (t222 * t47 - t225 * t51 + t227 * t84) * t388) * t391 + (-Ifges(5,5) + t399 * t377 - t398 * t237 / 0.2e1) * t332 + t283 + (t136 * mrSges(5,2) - Ifges(5,6) * t220 + pkin(4) * t148 - t227 * t147 - t225 * t154 + t222 * t156 - t84 * t295 + (pkin(8) * t155 - t395) * t237 + (-pkin(8) * t157 + t244) * t235 + t408 * t137) * qJD(4), t10 * qJD(3) + t4 * qJD(4) + (-t59 * mrSges(6,1) - t318 * t48 - t361 * t58) * qJD(5) + ((mrSges(7,2) * qJ(6) - t410) * t237 + t302 * t235) * t331 + t282, qJD(3) * t75 + qJD(4) * t19 + t346; t2, t24 * qJD(4), t337, t28 * qJD(5) + t77 * qJD(6) + t273 * t391 + t278 - t332 * t424 + (t206 + (-t295 + t408) * t220) * qJD(4), -t348 + t28 * qJD(4) - t200 * qJD(5) + (t235 * t361 - t237 * t362) * t331, t77 * qJD(4); -qJD(4) * t11 - qJD(5) * t9 + qJD(6) * t76 - t277, -t337, 0, -t347, -t361 * t329 + (-mrSges(6,1) - t318) * t330 - t261, t336; qJD(3) * t11 - qJD(5) * t3 - qJD(6) * t18 - t283, -qJD(5) * t27 + qJD(6) * t78 - t278, t347, -qJD(5) * t26 + qJD(6) * t130 (-t222 * mrSges(7,2) + t225 * t318) * qJD(5) + (pkin(8) * mrSges(6,2) - t410) * t330 + (-pkin(8) * mrSges(6,1) - t302) * t329 - t245, t252; qJD(3) * t9 + qJD(4) * t3 + qJD(6) * t87 - t282, qJD(4) * t27 + t348, t261, -t208 * qJD(6) + t245, 0, t276; -qJD(3) * t76 + qJD(4) * t18 - qJD(5) * t87 - t346, -t78 * qJD(4), -t336, qJD(5) * t208 - t252, -t276, 0;];
Cq  = t13;
