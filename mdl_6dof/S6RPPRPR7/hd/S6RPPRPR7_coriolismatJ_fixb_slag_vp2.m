% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:52:57
% EndTime: 2019-03-09 01:53:06
% DurationCPUTime: 5.01s
% Computational Cost: add. (14174->430), mult. (27009->603), div. (0->0), fcn. (30608->8), ass. (0->241)
t238 = sin(pkin(9));
t240 = cos(pkin(9));
t242 = sin(qJ(4));
t363 = cos(qJ(4));
t219 = t238 * t242 - t240 * t363;
t237 = sin(pkin(10));
t239 = cos(pkin(10));
t241 = sin(qJ(6));
t243 = cos(qJ(6));
t394 = -t237 * t241 + t239 * t243;
t403 = t394 * t219;
t375 = -t403 / 0.2e1;
t220 = t238 * t363 + t242 * t240;
t402 = t394 * t220;
t377 = -t402 / 0.2e1;
t222 = t237 * t243 + t239 * t241;
t339 = t222 * mrSges(7,3);
t405 = t339 * t375;
t396 = t220 * t222;
t404 = t222 * t402 - t394 * t396;
t359 = pkin(8) + qJ(5);
t226 = t359 * t237;
t228 = t359 * t239;
t194 = -t226 * t243 - t228 * t241;
t195 = -t226 * t241 + t228 * t243;
t260 = t222 * t219;
t319 = t219 * t237;
t181 = -t220 * mrSges(6,2) + mrSges(6,3) * t319;
t310 = t239 * t181;
t318 = t219 * t239;
t183 = t220 * mrSges(6,1) + mrSges(6,3) * t318;
t314 = t237 * t183;
t267 = t314 / 0.2e1 - t310 / 0.2e1;
t114 = mrSges(7,1) * t220 + mrSges(7,3) * t403;
t368 = t222 / 0.2e1;
t346 = t260 * mrSges(7,3);
t112 = -mrSges(7,2) * t220 + t346;
t400 = -t112 / 0.2e1;
t268 = t114 * t368 + t394 * t400;
t229 = t238 * pkin(3) + qJ(2);
t331 = qJ(5) * t219;
t361 = pkin(4) * t220;
t174 = t229 + t331 + t361;
t360 = pkin(1) + qJ(3);
t299 = -pkin(7) - t360;
t225 = t299 * t238;
t283 = t299 * t240;
t186 = t225 * t363 + t242 * t283;
t95 = t239 * t174 - t186 * t237;
t96 = t237 * t174 + t239 * t186;
t275 = t237 * t95 - t239 * t96;
t65 = pkin(5) * t220 + pkin(8) * t318 + t95;
t83 = pkin(8) * t319 + t96;
t36 = -t241 * t83 + t243 * t65;
t37 = t241 * t65 + t243 * t83;
t388 = -m(7) / 0.2e1;
t389 = m(6) / 0.2e1;
t401 = t275 * t389 + (t194 * t403 + t195 * t260 - t222 * t36 + t37 * t394) * t388 + t267 + t268;
t384 = mrSges(7,3) / 0.2e1;
t379 = t260 / 0.2e1;
t233 = t237 ^ 2;
t235 = t239 ^ 2;
t302 = t233 + t235;
t398 = mrSges(6,3) * t302;
t397 = t394 * t260;
t301 = t238 ^ 2 + t240 ^ 2;
t216 = Ifges(7,4) * t394;
t191 = Ifges(7,1) * t222 + t216;
t395 = -Ifges(7,2) * t222 + t191 + t216;
t185 = t225 * t242 - t363 * t283;
t362 = pkin(4) * t219;
t188 = qJ(5) * t220 - t362;
t107 = t185 * t237 + t239 * t188;
t108 = -t239 * t185 + t237 * t188;
t274 = -t107 * t237 + t108 * t239;
t345 = t402 * mrSges(7,2);
t349 = t396 * mrSges(7,1);
t305 = t349 / 0.2e1 + t345 / 0.2e1;
t393 = t220 ^ 2;
t392 = 0.2e1 * t220;
t391 = 2 * qJD(4);
t390 = -m(4) / 0.2e1;
t387 = m(7) / 0.2e1;
t386 = mrSges(7,1) / 0.2e1;
t385 = -mrSges(7,2) / 0.2e1;
t382 = m(7) * (t260 * t402 - t396 * t403);
t381 = -t396 / 0.2e1;
t380 = t396 / 0.2e1;
t376 = t403 / 0.2e1;
t374 = t195 / 0.2e1;
t373 = -t394 / 0.4e1;
t372 = t394 / 0.2e1;
t371 = -t219 / 0.2e1;
t370 = t219 / 0.2e1;
t369 = -t222 / 0.2e1;
t367 = t222 / 0.4e1;
t366 = t237 / 0.2e1;
t365 = -t239 / 0.2e1;
t364 = t239 / 0.2e1;
t316 = t222 * t403;
t64 = -t316 + t397;
t294 = t64 * t387;
t358 = qJD(4) * t294;
t357 = m(7) * qJD(3);
t356 = Ifges(6,4) * t237;
t355 = Ifges(6,4) * t239;
t354 = Ifges(7,4) * t403;
t353 = Ifges(7,4) * t222;
t352 = Ifges(6,5) * t239;
t351 = Ifges(6,6) * t237;
t111 = mrSges(7,2) * t219 + mrSges(7,3) * t396;
t113 = -mrSges(7,1) * t219 + mrSges(7,3) * t402;
t337 = t237 * Ifges(6,2);
t132 = -Ifges(6,6) * t219 + (t337 - t355) * t220;
t133 = -Ifges(6,5) * t219 + (-Ifges(6,1) * t239 + t356) * t220;
t134 = -pkin(5) * t319 + t185;
t313 = t237 * t220;
t135 = -pkin(5) * t313 + t186;
t336 = t239 * mrSges(6,2);
t338 = t237 * mrSges(6,1);
t280 = -t336 - t338;
t172 = t280 * t220;
t173 = t280 * t219;
t180 = mrSges(6,2) * t219 + mrSges(6,3) * t313;
t309 = t239 * t220;
t182 = -mrSges(6,1) * t219 + mrSges(6,3) * t309;
t213 = t220 * mrSges(5,2);
t271 = Ifges(7,5) * t377 + Ifges(7,6) * t380;
t66 = -pkin(5) * t219 + pkin(8) * t309 + t107;
t86 = pkin(8) * t313 + t108;
t40 = -t241 * t86 + t243 * t66;
t41 = t241 * t66 + t243 * t86;
t78 = -Ifges(7,4) * t402 + Ifges(7,2) * t396 - Ifges(7,6) * t219;
t79 = Ifges(7,2) * t260 + Ifges(7,6) * t220 - t354;
t80 = -Ifges(7,1) * t402 + Ifges(7,4) * t396 - Ifges(7,5) * t219;
t155 = Ifges(7,4) * t260;
t81 = -Ifges(7,1) * t403 + Ifges(7,5) * t220 + t155;
t90 = -t345 - t349;
t343 = t403 * mrSges(7,2);
t348 = t260 * mrSges(7,1);
t91 = -t343 - t348;
t1 = -t229 * t213 + t107 * t183 + t185 * t172 + t186 * t173 + t96 * t180 + t108 * t181 + t95 * t182 + t81 * t377 + t80 * t375 + t79 * t380 + t78 * t379 + t134 * t90 + t135 * t91 + t40 * t114 + t37 * t111 + t41 * t112 + t36 * t113 + m(6) * (t107 * t95 + t108 * t96 + t185 * t186) + m(7) * (t134 * t135 + t36 * t40 + t37 * t41) + (t133 * t365 + t132 * t366 - t229 * mrSges(5,1) + Ifges(7,5) * t376 - Ifges(7,6) * t260 / 0.2e1 + (t352 / 0.2e1 - t351 / 0.2e1 - Ifges(5,4)) * t219) * t219 + ((t235 * Ifges(6,1) / 0.2e1 - Ifges(6,3) - Ifges(7,3) + Ifges(5,1) - Ifges(5,2) + (-t355 + t337 / 0.2e1) * t237) * t219 + t271 + (Ifges(5,4) + t351 - t352) * t220) * t220;
t350 = t1 * qJD(1);
t347 = t260 * mrSges(7,2);
t344 = t403 * mrSges(7,1);
t342 = t394 * mrSges(7,1);
t341 = t394 * mrSges(7,3);
t340 = t222 * mrSges(7,2);
t276 = -t260 * t36 + t37 * t403;
t279 = Ifges(7,1) * t260 + t354;
t286 = Ifges(7,2) * t403 + t155;
t303 = Ifges(7,5) * t260 + Ifges(7,6) * t403;
t89 = -t344 + t347;
t6 = -t37 * t114 + t36 * t112 + t279 * t375 + t79 * t376 + t134 * t89 + t220 * t303 / 0.2e1 + t276 * mrSges(7,3) + (t81 + t286) * t379;
t335 = t6 * qJD(1);
t325 = t402 * t403;
t327 = t396 * t260;
t252 = (t327 / 0.2e1 + t325 / 0.2e1) * mrSges(7,3) + t112 * t381 + t114 * t377 + t89 * t370;
t269 = t342 / 0.2e1 - t340 / 0.2e1;
t7 = t252 - t269;
t334 = t7 * qJD(1);
t324 = t185 * t219;
t9 = -t402 * t112 + t396 * t114 + (mrSges(5,3) * t220 - t310 + t314) * t220 + (mrSges(5,3) * t219 - t173 - t91) * t219 + m(7) * (-t134 * t219 + t36 * t396 - t37 * t402) + m(6) * (t220 * t275 - t324) + m(5) * (-t186 * t220 - t324) + (m(4) * t360 + mrSges(4,3)) * t301;
t333 = t9 * qJD(1);
t197 = t219 ^ 2;
t251 = (-t396 ^ 2 - t402 ^ 2 - t197) * t387 + (-t302 * t393 - t197) * t389 + m(5) * (-t197 - t393) / 0.2e1 + t301 * t390;
t287 = m(6) * t302;
t288 = t222 ^ 2 + t394 ^ 2;
t255 = t390 - m(5) / 0.2e1 - t287 / 0.2e1 + t288 * t388;
t29 = t251 + t255;
t330 = qJD(1) * t29;
t323 = t195 * t403;
t214 = t220 * mrSges(5,1);
t256 = m(6) * (t237 * t96 + t239 * t95) + t237 * t181 + t239 * t183;
t21 = t238 * mrSges(4,1) + t240 * mrSges(4,2) - t219 * mrSges(5,2) + t222 * t112 + t394 * t114 + mrSges(3,3) + t214 + (m(4) + m(3)) * qJ(2) + m(7) * (t222 * t37 + t36 * t394) + m(5) * t229 + t256;
t322 = t21 * qJD(1);
t315 = t237 * t182;
t311 = t239 * t180;
t266 = m(7) * (t222 * t260 + t394 * t403);
t282 = -t287 / 0.4e1;
t298 = m(7) / 0.4e1 + m(6) / 0.4e1;
t44 = -t266 / 0.2e1 + 0.2e1 * (t282 - t298) * t219;
t307 = t44 * qJD(1);
t187 = t222 * mrSges(7,1) + mrSges(7,2) * t394;
t300 = t187 * qJD(6);
t296 = -t382 / 0.2e1;
t295 = t382 / 0.2e1;
t293 = Ifges(7,2) / 0.4e1 - Ifges(7,1) / 0.4e1;
t292 = -t346 / 0.2e1;
t290 = t187 * t370;
t189 = t340 - t342;
t227 = -mrSges(6,1) * t239 + mrSges(6,2) * t237;
t289 = -t189 / 0.2e1 - t227 / 0.2e1;
t284 = t302 * t220;
t281 = t287 / 0.2e1;
t278 = Ifges(7,1) * t394 - t353;
t277 = Ifges(7,5) * t394 - Ifges(7,6) * t222;
t49 = 0.2e1 * mrSges(7,1) * t375 + 0.2e1 * t379 * mrSges(7,2);
t272 = qJD(1) * t49 + qJD(4) * t187;
t270 = t348 / 0.2e1 + t343 / 0.2e1;
t265 = m(7) * (t222 * t396 + t394 * t402);
t264 = m(7) * (t194 * t394 + t195 * t222);
t231 = -pkin(5) * t239 - pkin(4);
t247 = (t369 * t396 - t372 * t402) * mrSges(7,3) + (-t235 / 0.2e1 - t233 / 0.2e1) * t220 * mrSges(6,3) + (-qJ(5) * t284 + t362) * t389 + (t194 * t396 - t195 * t402 - t231 * t219) * t387;
t248 = (t239 * t107 + t237 * t108) * t389 + (t222 * t41 + t394 * t40) * t387 + t113 * t372 + t111 * t368 + t180 * t366 + t182 * t364;
t12 = t213 + (mrSges(5,1) + t289) * t219 + t247 - t248;
t263 = -t12 * qJD(1) + qJD(2) * t294;
t14 = (t397 / 0.2e1 - t316 / 0.2e1) * mrSges(7,3) + t268 + t305;
t262 = t14 * qJD(1);
t18 = t403 * t114 + t260 * t112 + m(7) * (t260 * t37 + t36 * t403) + t256 * t219;
t261 = -t18 * qJD(1) + qJD(2) * t296;
t253 = (-t336 / 0.2e1 - t338 / 0.2e1) * t220 + t186 * t389 + t135 * t387 - t305;
t10 = (-t397 / 0.2e1 + t403 * t368) * mrSges(7,3) + t253 + t401;
t43 = -t265 / 0.2e1 + (t282 + t298) * t392;
t57 = m(7) * (-t194 * t222 + t195 * t394) + mrSges(7,3) * t288 + qJ(5) * t287 + t398;
t259 = -t10 * qJD(1) - t43 * qJD(2) + t57 * qJD(4);
t196 = t219 * t220;
t30 = m(7) * (t196 - t325 - t327) + m(6) * (-t219 * t284 + t196);
t246 = (t90 / 0.2e1 + t172 / 0.2e1 + t267) * t219 + (t91 / 0.2e1 + t311 / 0.2e1 - t315 / 0.2e1 + t173 / 0.2e1) * t220 + ((t185 + t274) * t220 + (t186 + t275) * t219) * t389 + (t134 * t220 + t219 * t135 - t396 * t40 + t402 * t41 - t276) * t387 + t113 * t381 + t114 * t379 + t402 * t111 / 0.2e1 + t112 * t375;
t5 = -t264 / 0.2e1 + t246;
t258 = -t30 * qJD(2) - t5 * qJD(1) - t64 * t357 / 0.2e1;
t190 = Ifges(7,2) * t394 + t353;
t249 = -t134 * t187 / 0.2e1 + t194 * t400 + t114 * t374 + t81 * t373 - t220 * t277 / 0.4e1 + t79 * t367 - t231 * t89 / 0.2e1;
t254 = Ifges(7,3) * t371 + t385 * t41 + t386 * t40 + t271;
t2 = t155 * t373 + (t194 * t384 - t191 / 0.4e1 - t216 / 0.4e1 + t293 * t222) * t260 - (t353 / 0.2e1 + t190 / 0.4e1 + mrSges(7,3) * t374 + t293 * t394) * t403 + t249 + t254;
t23 = t290 - t270;
t27 = t231 * t187 + t190 * t369 + t278 * t368 + t372 * t395;
t257 = -t2 * qJD(1) + t23 * qJD(2) + t27 * qJD(4);
t51 = qJD(5) * t295;
t50 = t347 / 0.2e1 - t344 / 0.2e1 + t260 * t385 + t403 * t386;
t46 = t266 / 0.2e1 + (t281 - 0.2e1 * t298) * t219;
t45 = t265 / 0.2e1 + t220 * t281 + t298 * t392;
t28 = t251 - t255;
t24 = t290 + t270;
t15 = t292 * t394 - t268 + t305 - t405;
t13 = t219 * t289 + t247 + t248;
t11 = t384 * t397 + t253 - t401 + t405;
t8 = t252 + t269;
t4 = t264 / 0.2e1 + t246;
t3 = t194 * t292 + t323 * t384 - t249 + t279 * t367 + t394 * t286 / 0.4e1 + t254 + t395 * t260 / 0.4e1 - (t278 / 0.4e1 - t190 / 0.4e1) * t403;
t16 = [qJD(2) * t21 + qJD(3) * t9 + qJD(4) * t1 + qJD(5) * t18 + qJD(6) * t6, m(7) * qJD(2) * t404 + t28 * qJD(3) + t4 * qJD(4) + t8 * qJD(6) + t322 + t51, t28 * qJD(2) + t13 * qJD(4) + t46 * qJD(5) + t15 * qJD(6) - t357 * t404 + t333, t350 + t4 * qJD(2) + t13 * qJD(3) + t11 * qJD(5) + t3 * qJD(6) + ((-pkin(4) * t186 + qJ(5) * t274) * t389 + (t135 * t231 + t194 * t40 + t195 * t41) * t387) * t391 + (-t40 * t339 + t41 * t341 + t132 * t364 + t133 * t366 + t231 * t90 + t80 * t368 + t78 * t372 + Ifges(5,6) * t219 + t135 * t189 + t190 * t380 + t191 * t377 + t194 * t113 + t195 * t111 + t185 * mrSges(5,2) - pkin(4) * t172 + ((Ifges(6,1) * t237 + t355) * t365 + (Ifges(6,2) * t239 + t356) * t366 - Ifges(5,5)) * t220 + (Ifges(6,5) * t237 + Ifges(7,5) * t222 + Ifges(6,6) * t239 + Ifges(7,6) * t394) * t371 + (t227 - mrSges(5,1)) * t186 + (-t315 + t311) * qJ(5) + t274 * mrSges(6,3)) * qJD(4), t46 * qJD(3) + t11 * qJD(4) + t50 * qJD(6) - t261, t335 + t8 * qJD(2) + t15 * qJD(3) + t3 * qJD(4) + t50 * qJD(5) + (-mrSges(7,1) * t37 - mrSges(7,2) * t36 + t303) * qJD(6); qJD(3) * t29 + qJD(4) * t5 + qJD(6) * t7 - t322 + t51, t30 * qJD(4), t330 + t358, t45 * qJD(5) + t24 * qJD(6) + ((t194 * t260 + t231 * t220 - t323) * t387 + (-t302 * t331 - t361) * t389) * t391 - t258 + (-t260 * t339 - t403 * t341 - t214 + (mrSges(5,2) - t398) * t219 + (t189 + t227) * t220) * qJD(4), qJD(1) * t295 + t45 * qJD(4), t334 + t24 * qJD(4) + (-mrSges(7,1) * t402 + mrSges(7,2) * t396) * qJD(6); -qJD(2) * t29 - qJD(4) * t12 - qJD(5) * t44 - qJD(6) * t14 - t333, -t330 + t358, 0, t263, -t307, -t262 - t300; -qJD(2) * t5 + qJD(3) * t12 - qJD(5) * t10 - qJD(6) * t2 - t350, -t43 * qJD(5) + t23 * qJD(6) + t258, -t263, qJD(5) * t57 + qJD(6) * t27, t259 (-mrSges(7,1) * t195 - mrSges(7,2) * t194 + t277) * qJD(6) + t257; t44 * qJD(3) + t10 * qJD(4) + t49 * qJD(6) + t261, qJD(1) * t296 + t43 * qJD(4), t307, -t259 + t300, 0, t272; -qJD(2) * t7 + qJD(3) * t14 + qJD(4) * t2 - qJD(5) * t49 - t335, -t23 * qJD(4) - t334, t262, -t187 * qJD(5) - t257, -t272, 0;];
Cq  = t16;
