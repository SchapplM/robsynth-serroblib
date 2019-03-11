% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRP8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:15:29
% EndTime: 2019-03-09 02:15:36
% DurationCPUTime: 4.92s
% Computational Cost: add. (10007->452), mult. (18663->589), div. (0->0), fcn. (19352->6), ass. (0->234)
t259 = sin(qJ(5));
t254 = t259 ^ 2;
t261 = cos(qJ(5));
t255 = t261 ^ 2;
t323 = t254 + t255;
t414 = Ifges(7,4) + Ifges(6,5);
t415 = mrSges(7,2) + mrSges(6,3);
t281 = t415 * (t255 / 0.2e1 + t254 / 0.2e1);
t256 = sin(pkin(9));
t257 = cos(pkin(9));
t260 = sin(qJ(4));
t383 = cos(qJ(4));
t215 = t256 * t260 - t257 * t383;
t216 = t256 * t383 + t260 * t257;
t247 = Ifges(7,5) * t259;
t234 = Ifges(7,1) * t261 + t247;
t100 = t216 * Ifges(7,4) - t215 * t234;
t364 = Ifges(6,4) * t259;
t236 = Ifges(6,1) * t261 - t364;
t102 = t216 * Ifges(6,5) - t215 * t236;
t318 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t303 = t318 * t216;
t378 = pkin(5) * t216;
t237 = t256 * pkin(3) + qJ(2);
t379 = pkin(4) * t216;
t133 = pkin(8) * t215 + t237 + t379;
t258 = -pkin(1) - qJ(3);
t371 = -pkin(7) + t258;
t221 = t371 * t256;
t305 = t371 * t257;
t153 = t221 * t383 + t260 * t305;
t59 = t133 * t261 - t153 * t259;
t45 = -t59 - t378;
t421 = t303 - t59 * mrSges(6,3) + t100 / 0.2e1 + t102 / 0.2e1 + t45 * mrSges(7,2);
t420 = t323 * pkin(8);
t392 = -t215 / 0.2e1;
t419 = t215 / 0.2e1;
t418 = -t216 / 0.2e1;
t416 = t236 / 0.2e1;
t372 = mrSges(6,1) + mrSges(7,1);
t413 = Ifges(7,2) + Ifges(6,3);
t337 = t215 * t259;
t144 = -t216 * mrSges(6,2) + mrSges(6,3) * t337;
t358 = t216 * mrSges(7,3);
t147 = mrSges(7,2) * t337 + t358;
t329 = t144 + t147;
t336 = t215 * t261;
t145 = t216 * mrSges(6,1) + mrSges(6,3) * t336;
t146 = -t216 * mrSges(7,1) - mrSges(7,2) * t336;
t328 = t145 - t146;
t351 = t261 * mrSges(7,1);
t353 = t259 * mrSges(7,3);
t224 = -t351 - t353;
t352 = t261 * mrSges(6,1);
t354 = t259 * mrSges(6,2);
t225 = -t352 + t354;
t412 = t224 + t225;
t229 = -Ifges(7,3) * t261 + t247;
t411 = Ifges(7,6) * t259 + t261 * t414;
t250 = Ifges(6,4) * t261;
t235 = Ifges(6,1) * t259 + t250;
t300 = pkin(5) * t261 + qJ(6) * t259;
t301 = Ifges(6,2) * t259 - t250;
t167 = t215 ^ 2;
t408 = t216 ^ 2;
t410 = t167 / 0.2e1 + t408 / 0.2e1;
t244 = m(7) * qJ(6) + mrSges(7,3);
t202 = Ifges(7,5) * t336;
t140 = Ifges(7,1) * t337 - t202;
t141 = t215 * t235;
t363 = Ifges(7,5) * t261;
t230 = Ifges(7,3) * t259 + t363;
t233 = Ifges(7,1) * t259 - t363;
t335 = t216 * qJ(6);
t330 = t261 * t153;
t331 = t259 * t133;
t60 = t330 + t331;
t44 = t60 + t335;
t370 = t44 - t60;
t404 = -m(7) / 0.2e1;
t96 = Ifges(7,6) * t216 - Ifges(7,3) * t337 - t202;
t357 = t216 * Ifges(6,6);
t98 = t215 * t301 + t357;
t409 = (t235 / 0.4e1 + t233 / 0.4e1 - t301 / 0.4e1 - t230 / 0.4e1) * t215 + (t370 * t404 - t147 / 0.2e1 - t144 / 0.2e1) * pkin(8) + t140 / 0.4e1 + t141 / 0.4e1 + t96 / 0.4e1 - t98 / 0.4e1;
t407 = 2 * qJD(4);
t406 = -m(6) / 0.2e1;
t405 = m(6) / 0.2e1;
t403 = m(7) / 0.2e1;
t402 = mrSges(7,1) / 0.2e1;
t401 = mrSges(6,2) / 0.2e1;
t400 = -mrSges(7,3) / 0.2e1;
t399 = Ifges(7,6) / 0.2e1;
t398 = t60 / 0.2e1;
t152 = t221 * t260 - t383 * t305;
t376 = pkin(8) * t216;
t380 = pkin(4) * t215;
t154 = t376 - t380;
t63 = t152 * t259 + t154 * t261;
t52 = pkin(5) * t215 - t63;
t397 = m(7) * t52;
t396 = qJ(6) / 0.2e1;
t134 = t300 * t215;
t395 = -t134 / 0.2e1;
t135 = t215 * t224;
t394 = t135 / 0.2e1;
t349 = t261 * mrSges(7,3);
t355 = t259 * mrSges(7,1);
t227 = -t349 + t355;
t137 = t215 * t227;
t393 = t137 / 0.2e1;
t390 = -t225 / 0.2e1;
t389 = t227 / 0.2e1;
t350 = t261 * mrSges(6,2);
t356 = t259 * mrSges(6,1);
t228 = t350 + t356;
t388 = t228 / 0.2e1;
t387 = -t259 / 0.2e1;
t386 = t259 / 0.2e1;
t385 = -t261 / 0.2e1;
t384 = t261 / 0.2e1;
t345 = qJ(6) * t261;
t373 = t259 * pkin(5);
t299 = -t345 + t373;
t382 = m(7) * t299;
t381 = m(7) * t261;
t369 = t45 + t59;
t368 = m(7) * qJD(5);
t367 = m(7) * qJD(6);
t293 = t259 * t59 - t261 * t60;
t297 = -t259 * t45 - t261 * t44;
t324 = t256 ^ 2 + t257 ^ 2;
t306 = m(4) * t324;
t339 = t152 * t215;
t65 = -t215 * t299 + t152;
t9 = t324 * mrSges(4,3) + (mrSges(5,3) * t216 + t259 * t328 - t261 * t329) * t216 + (t137 + (mrSges(5,3) + t228) * t215) * t215 + m(5) * (-t153 * t216 - t339) + m(7) * (-t215 * t65 + t216 * t297) + m(6) * (t216 * t293 - t339) - t258 * t306;
t360 = qJD(1) * t9;
t359 = t215 * mrSges(7,1);
t136 = t215 * t225;
t138 = t215 * t229;
t231 = Ifges(6,2) * t261 + t364;
t139 = t215 * t231;
t201 = Ifges(7,6) * t336;
t284 = t44 * mrSges(7,2) + t60 * mrSges(6,3) - t96 / 0.2e1 + t98 / 0.2e1;
t313 = t357 / 0.2e1;
t4 = t201 * t418 + t152 * t136 + t65 * t135 + ((t313 - t140 / 0.2e1 - t141 / 0.2e1 + t284) * t261 + (t139 / 0.2e1 - t138 / 0.2e1 + t421) * t259) * t215 - (m(7) * t65 - t137) * t134 + (m(7) * t45 - t328) * t60 + (m(7) * t44 + t329) * t59;
t348 = t4 * qJD(1);
t294 = t259 * t60 + t261 * t59;
t296 = t259 * t44 - t261 * t45;
t265 = (m(7) * t395 + t136 / 0.2e1 + t394) * t215 + ((t294 - t296) * t403 + t145 * t385 + t146 * t384 + t281 * t215 + t329 * t387) * t216;
t270 = t300 * t403 - t354 / 0.2e1 + t353 / 0.2e1 + t352 / 0.2e1 + t351 / 0.2e1;
t7 = t265 - t270;
t347 = t7 * qJD(1);
t212 = t216 * mrSges(5,1);
t15 = t256 * mrSges(4,1) + t257 * mrSges(4,2) - t215 * mrSges(5,2) + mrSges(3,3) + t212 + t328 * t261 + t329 * t259 + (m(4) + m(3)) * qJ(2) + m(5) * t237 + m(7) * t296 + m(6) * t294;
t344 = qJD(1) * t15;
t22 = m(7) * (t216 * t44 + t336 * t65) + t216 * t147 - t137 * t336;
t343 = qJD(1) * t22;
t69 = (0.1e1 / 0.2e1 + t410) * t381;
t342 = qJD(1) * t69;
t278 = t145 * t386 + t146 * t387 + t329 * t385;
t267 = (t259 * t369 + t261 * t370) * t404 + t278;
t269 = t350 / 0.2e1 + t356 / 0.2e1 + t355 / 0.2e1 - t349 / 0.2e1 + t299 * t403;
t268 = t269 * t216;
t10 = t268 + t267;
t341 = t10 * qJD(1);
t211 = t216 * mrSges(5,2);
t334 = t216 * t259;
t287 = mrSges(7,2) * t334 - t215 * mrSges(7,3);
t333 = t216 * t261;
t319 = mrSges(7,2) * t333;
t288 = -t319 + t359;
t64 = -t261 * t152 + t259 * t154;
t51 = -qJ(6) * t215 + t64;
t264 = (t259 * t64 + t261 * t63) * t405 + (t259 * t51 - t261 * t52) * t403 + (-t215 * mrSges(6,1) + mrSges(6,3) * t333) * t384 + t288 * t385 + (t215 * mrSges(6,2) + mrSges(6,3) * t334 + t287) * t386;
t327 = t420 * t216;
t222 = -pkin(4) - t300;
t338 = t215 * t222;
t277 = (-t327 + t380) * t405 + (-t327 - t338) * t403;
t12 = -t135 / 0.2e1 + t211 - t264 + t277 + t323 * t216 * (-mrSges(6,3) / 0.2e1 - mrSges(7,2) / 0.2e1) + (t390 + mrSges(5,1)) * t215;
t340 = t12 * qJD(1);
t320 = m(7) / 0.4e1 + m(6) / 0.4e1;
t276 = -t306 / 0.2e1 + m(5) * (-t167 - t408) / 0.2e1 + 0.2e1 * t320 * (-t323 * t408 - t167);
t280 = -m(4) / 0.2e1 - m(5) / 0.2e1 + (t406 + t404) * t323;
t24 = t276 + t280;
t332 = t24 * qJD(1);
t326 = t420 * t215;
t322 = qJD(4) * t215;
t130 = m(7) * t334;
t321 = t130 * qJD(1);
t317 = t399 - Ifges(6,6) / 0.2e1;
t316 = t59 / 0.2e1 + t45 / 0.2e1;
t315 = t398 - t44 / 0.2e1;
t314 = Ifges(6,6) * t392 + Ifges(7,6) * t419 + (t230 + t301) * t216 / 0.2e1;
t312 = t216 * t416 - t234 * t418 + t414 * t419;
t310 = t224 * t384;
t309 = -t229 / 0.2e1 + t231 / 0.2e1;
t308 = t233 / 0.2e1 + t235 / 0.2e1;
t66 = -t216 * t299 + t153;
t1 = -t66 * t137 + t64 * t144 + t63 * t145 + t52 * t146 + t51 * t147 - t237 * t211 + m(6) * (t152 * t153 + t59 * t63 + t60 * t64) + m(7) * (t44 * t51 + t45 * t52 + t65 * t66) + (-t237 * mrSges(5,1) - t59 * mrSges(6,1) + t45 * mrSges(7,1) + t60 * mrSges(6,2) - t44 * mrSges(7,3) - Ifges(5,4) * t215 + (-t153 * mrSges(6,2) + t215 * t318 + t312) * t261 + (-t153 * mrSges(6,1) + t215 * t317 + t314) * t259) * t215 + (Ifges(5,4) * t216 + (-Ifges(5,2) + Ifges(5,1) - t413) * t215 + (-t152 * mrSges(6,2) + t65 * mrSges(7,3) - t421) * t261 + (-t152 * mrSges(6,1) - t65 * mrSges(7,1) - t216 * t317 + t284) * t259) * t216;
t292 = -t259 * t63 + t261 * t64;
t295 = t259 * t52 + t261 * t51;
t6 = (t393 + (t295 + t65) * t404 + (t152 + t292) * t406) * t216 + (t216 * t388 + (t297 + t66) * t404 + (t153 + t293) * t406 - t278) * t215;
t298 = t1 * qJD(1) - t6 * qJD(2);
t25 = 0.4e1 * t320 * (0.1e1 - t323) * t216 * t215;
t291 = -t6 * qJD(1) + t25 * qJD(2);
t157 = (m(7) * t222 + t224) * t259;
t275 = (-t259 * t65 + (t338 + t376) * t261) * t403 - t137 * t387;
t21 = t319 + (t310 - mrSges(7,1) / 0.2e1) * t215 - t397 / 0.2e1 + t275;
t290 = -qJD(1) * t21 + qJD(4) * t157;
t27 = t358 + 0.2e1 * (t335 / 0.2e1 + t331 / 0.4e1 + t330 / 0.4e1 - t60 / 0.4e1) * m(7);
t289 = qJD(1) * t27 + qJD(5) * t244;
t285 = -t236 / 0.4e1 - t234 / 0.4e1 + t231 / 0.4e1 - t229 / 0.4e1;
t28 = (-t228 / 0.2e1 - t227 / 0.2e1 + (t401 + t400) * t261 + (mrSges(6,1) / 0.2e1 + t402) * t259 + (t373 / 0.2e1 - t345 / 0.2e1 - t299 / 0.2e1) * m(7)) * t215;
t266 = (-t134 * t222 + t299 * t65) * t403 - pkin(4) * t136 / 0.2e1 + t224 * t395 + t152 * t388 + t222 * t394 - t299 * t393 + t65 * t389 + t411 * t216 / 0.4e1;
t271 = (-pkin(5) * t52 + qJ(6) * t51) * t404 + t51 * t400 + t52 * t402 - t63 * mrSges(6,1) / 0.2e1 + t64 * t401;
t273 = (t369 * t403 - t145 / 0.2e1 + t146 / 0.2e1) * pkin(8) + t100 / 0.4e1 + t102 / 0.4e1 - t138 / 0.4e1 + t139 / 0.4e1;
t279 = t281 * pkin(8);
t3 = (t303 + (-t378 / 0.2e1 + t316) * mrSges(7,2) + t285 * t215 + t273) * t261 + ((-0.3e1 / 0.4e1 * Ifges(6,6) + t399) * t216 + (-t335 / 0.2e1 + t315) * mrSges(7,2) + t409) * t259 + (mrSges(7,3) * t396 + pkin(5) * t402 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + t279) * t215 + t266 + t271;
t32 = -pkin(4) * t228 + t299 * t224 + (t227 + t382) * t222 + (-t301 / 0.2e1 - t230 / 0.2e1 + t308) * t261 + (t234 / 0.2e1 + t416 - t309) * t259;
t283 = t3 * qJD(1) - t28 * qJD(2) + t32 * qJD(4);
t272 = (-m(7) * t300 + t412) * qJD(5);
t223 = (m(7) * pkin(8) + mrSges(7,2)) * t261;
t127 = m(7) * t337;
t70 = (-0.1e1 / 0.2e1 + t410) * t381;
t29 = -t382 * t392 + (t269 + t388 + t389) * t215;
t26 = (t60 + 0.2e1 * t335) * t403 + m(7) * t398 + t147;
t23 = t276 - t280;
t20 = t215 * t310 + t397 / 0.2e1 + t359 / 0.2e1 + t275;
t13 = (t390 - t224 / 0.2e1) * t215 + t264 + t277 - t281 * t216;
t11 = t268 - t267;
t8 = t265 + t270;
t5 = t6 * qJD(4);
t2 = (mrSges(7,2) * t316 + t273) * t261 - Ifges(7,6) * t334 / 0.2e1 + t266 - pkin(5) * t288 / 0.2e1 + t287 * t396 + (t285 * t261 + t279) * t215 - t271 + t413 * t392 - t414 * t333 / 0.2e1 + (-t357 / 0.4e1 + t315 * mrSges(7,2) + t313 + t409) * t259;
t14 = [qJD(2) * t15 + qJD(3) * t9 + qJD(4) * t1 + qJD(5) * t4 + qJD(6) * t22, qJD(3) * t23 + qJD(5) * t8 + qJD(6) * t70 + t344 - t5, qJD(2) * t23 + qJD(4) * t13 + qJD(5) * t11 + t360, t13 * qJD(3) + t2 * qJD(5) + t20 * qJD(6) + ((pkin(8) * t295 + t222 * t66) * t403 + (-pkin(4) * t153 + pkin(8) * t292) * t405) * t407 + (Ifges(5,6) + t317 * t261 - t318 * t259 + ((mrSges(6,2) - mrSges(7,3)) * t261 + t372 * t259) * pkin(8)) * t322 + t298 + (t152 * mrSges(5,2) + t66 * t224 + (t51 * mrSges(7,2) + t64 * mrSges(6,3) + t314) * t261 + (t52 * mrSges(7,2) - t63 * mrSges(6,3) - t312) * t259 + (-Ifges(5,5) + (pkin(4) * mrSges(6,2) + t222 * mrSges(7,3) - t308) * t261 + (pkin(4) * mrSges(6,1) - t222 * mrSges(7,1) + t309) * t259) * t216 + (t225 - mrSges(5,1)) * t153) * qJD(4), t8 * qJD(2) + t11 * qJD(3) + t2 * qJD(4) + t26 * qJD(6) + t348 + (-t201 + (-m(7) * pkin(5) - t372) * t60 + (-mrSges(6,2) + t244) * t59 + ((qJ(6) * mrSges(7,2) + Ifges(6,6)) * t261 + (-pkin(5) * mrSges(7,2) + t414) * t259) * t215) * qJD(5), qJD(2) * t70 + qJD(4) * t20 + qJD(5) * t26 + t343; qJD(3) * t24 + qJD(5) * t7 + qJD(6) * t69 - t344 - t5, t25 * qJD(4), t332 (t216 * t412 - t212) * qJD(4) + t29 * qJD(5) - t127 * qJD(6) + ((-t326 - t379) * t405 + (t216 * t222 - t326) * t403) * t407 + (-t323 * t415 + mrSges(5,2)) * t322 + t291, t347 + t29 * qJD(4) + (t261 * t367 + t272) * t216, -qJD(4) * t127 + t333 * t368 + t342; -qJD(2) * t24 - qJD(4) * t12 - qJD(5) * t10 + qJD(6) * t130 - t360, -t332, 0, -t340, -t341 + (t349 - t350 - t382) * qJD(5) + (-qJD(5) * t372 + t367) * t259, t259 * t368 + t321; qJD(3) * t12 + qJD(5) * t3 + qJD(6) * t21 - t298, -qJD(5) * t28 - t291, t340, qJD(5) * t32 - qJD(6) * t157 (-mrSges(7,2) * t300 - Ifges(6,6) * t259 + t411) * qJD(5) + t223 * qJD(6) + pkin(8) * t272 + t283, qJD(5) * t223 - t290; -qJD(2) * t7 + qJD(3) * t10 - qJD(4) * t3 + qJD(6) * t27 - t348, qJD(4) * t28 - t347, t341, -t283, t244 * qJD(6), t289; -qJD(2) * t69 - qJD(3) * t130 - qJD(4) * t21 - qJD(5) * t27 - t343, -t342, -t321, t290, -t289, 0;];
Cq  = t14;
