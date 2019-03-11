% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:26
% EndTime: 2019-03-09 03:01:34
% DurationCPUTime: 4.97s
% Computational Cost: add. (10435->426), mult. (20134->578), div. (0->0), fcn. (20914->8), ass. (0->223)
t262 = cos(qJ(5));
t367 = Ifges(6,6) + Ifges(7,6);
t418 = t367 * t262;
t408 = Ifges(6,5) + Ifges(7,5);
t373 = sin(qJ(3));
t324 = t373 * pkin(3);
t417 = m(5) * t324;
t344 = sin(pkin(10));
t345 = cos(pkin(10));
t374 = cos(qJ(3));
t237 = -t344 * t374 - t345 * t373;
t371 = m(7) * t237;
t212 = -t371 / 0.2e1;
t261 = sin(qJ(5));
t258 = t261 ^ 2;
t259 = t262 ^ 2;
t328 = t258 + t259;
t296 = t328 * t371;
t413 = t296 / 0.2e1 + t212;
t416 = qJD(6) * t413;
t356 = t261 * mrSges(6,2);
t241 = -mrSges(6,1) * t262 + t356;
t415 = mrSges(5,1) - t241;
t414 = t261 * t408 + t418;
t412 = m(5) * pkin(3);
t411 = mrSges(7,3) / 0.2e1;
t235 = t344 * t373 - t345 * t374;
t337 = t237 * t261;
t322 = mrSges(7,3) * t337;
t169 = -t235 * mrSges(7,2) + t322;
t410 = -t169 / 0.2e1;
t407 = Ifges(7,3) + Ifges(6,3);
t395 = m(7) * pkin(5);
t317 = mrSges(7,1) + t395;
t311 = t344 * pkin(3);
t250 = t311 + pkin(8);
t332 = qJ(6) + t250;
t221 = t332 * t261;
t223 = t332 * t262;
t406 = t221 * t261 + t223 * t262;
t362 = Ifges(7,4) * t261;
t242 = t262 * Ifges(7,2) + t362;
t363 = Ifges(6,4) * t261;
t243 = t262 * Ifges(6,2) + t363;
t405 = t242 + t243;
t255 = Ifges(7,4) * t262;
t244 = t261 * Ifges(7,1) + t255;
t256 = Ifges(6,4) * t262;
t245 = t261 * Ifges(6,1) + t256;
t404 = t245 + t244;
t251 = sin(pkin(9)) * pkin(1) + pkin(7);
t312 = t373 * t251;
t222 = -qJ(4) * t373 - t312;
t314 = t374 * t251;
t224 = qJ(4) * t374 + t314;
t403 = t344 * t222 + t345 * t224;
t154 = -t345 * t222 + t224 * t344;
t164 = -t237 * pkin(4) + t235 * pkin(8) + t324;
t67 = t154 * t261 + t262 * t164;
t68 = -t154 * t262 + t261 * t164;
t282 = -t261 * t67 + t262 * t68;
t338 = t235 * t262;
t44 = -t237 * pkin(5) + qJ(6) * t338 + t67;
t339 = t235 * t261;
t52 = qJ(6) * t339 + t68;
t402 = -t261 * t44 + t262 * t52;
t401 = t235 * (-mrSges(6,3) - mrSges(7,3));
t399 = m(6) / 0.4e1;
t398 = -m(7) / 0.2e1;
t397 = m(7) / 0.2e1;
t396 = m(7) / 0.4e1;
t394 = -mrSges(6,1) / 0.2e1;
t393 = -mrSges(7,1) / 0.2e1;
t392 = mrSges(6,2) / 0.2e1;
t391 = mrSges(7,2) / 0.2e1;
t390 = -t44 / 0.2e1;
t88 = -pkin(5) * t337 + t154;
t389 = t88 / 0.2e1;
t87 = -pkin(5) * t339 + t403;
t388 = m(7) * t87;
t387 = m(7) * t88;
t386 = t154 / 0.2e1;
t171 = -t237 * mrSges(7,1) + mrSges(7,3) * t338;
t385 = t171 / 0.2e1;
t336 = t237 * t262;
t360 = t235 * mrSges(6,1);
t174 = mrSges(6,3) * t336 + t360;
t384 = t174 / 0.2e1;
t383 = -t235 / 0.2e1;
t381 = -t237 / 0.2e1;
t240 = -mrSges(7,1) * t262 + t261 * mrSges(7,2);
t380 = -t240 / 0.2e1;
t379 = -t250 / 0.2e1;
t378 = -t261 / 0.2e1;
t377 = t261 / 0.2e1;
t375 = t262 / 0.2e1;
t313 = t345 * pkin(3);
t252 = -t313 - pkin(4);
t368 = t262 * pkin(5);
t238 = t252 - t368;
t370 = m(7) * t238;
t369 = t235 * pkin(5);
t315 = -cos(pkin(9)) * pkin(1) - pkin(2);
t239 = -pkin(3) * t374 + t315;
t140 = t235 * pkin(4) + t237 * pkin(8) + t239;
t62 = t262 * t140 - t261 * t403;
t48 = qJ(6) * t336 + t62;
t41 = t48 + t369;
t366 = -t41 + t48;
t365 = mrSges(5,3) * t235;
t364 = mrSges(5,3) * t237;
t219 = t235 * mrSges(5,2);
t272 = (-t235 * t406 - t238 * t237) * t397 + (-t235 * t344 + t237 * t345) * t412 / 0.2e1 + m(6) * (-t235 * t250 * t328 - t237 * t252) / 0.2e1;
t302 = t259 / 0.2e1 + t258 / 0.2e1;
t265 = t302 * t401 + t272;
t266 = -m(6) * (t261 * t68 + t262 * t67) / 0.2e1 + (t261 * t52 + t262 * t44) * t398 - t417 / 0.2e1;
t304 = t380 - t241 / 0.2e1;
t172 = -t237 * mrSges(6,1) + mrSges(6,3) * t338;
t306 = t385 + t172 / 0.2e1;
t167 = mrSges(7,2) * t237 + mrSges(7,3) * t339;
t168 = mrSges(6,2) * t237 + mrSges(6,3) * t339;
t308 = -t167 / 0.2e1 - t168 / 0.2e1;
t9 = t219 - t306 * t262 + t308 * t261 + (mrSges(5,1) + t304) * t237 + t265 + t266;
t361 = qJD(1) * t9;
t359 = t235 * mrSges(7,1);
t358 = t235 * mrSges(6,2);
t357 = t261 * mrSges(7,1);
t351 = t262 * mrSges(7,2);
t63 = t261 * t140 + t262 * t403;
t49 = qJ(6) * t337 + t63;
t350 = t262 * t49;
t348 = t262 * t63;
t215 = mrSges(7,2) * t337;
t295 = t302 * mrSges(6,3);
t173 = mrSges(7,3) * t336 + t359;
t305 = t173 / 0.2e1 + t384;
t323 = mrSges(6,3) * t337;
t170 = t323 - t358;
t307 = t410 - t170 / 0.2e1;
t7 = t215 * t383 + ((mrSges(7,3) * t302 + t295) * t237 + (-t358 / 0.2e1 + t307) * t261 + (t360 / 0.2e1 + t359 / 0.2e1 + (-t366 - t369) * t398 - t305) * t262) * t237;
t346 = t7 * qJD(1);
t21 = (m(7) * (t261 * t49 + t41 * t262) + t261 * t169 + t262 * t173) * t237;
t343 = qJD(1) * t21;
t341 = t154 * t237;
t335 = t250 * t262;
t299 = (t396 + t399) * (-0.1e1 + t328) * t237 * t235;
t26 = 0.2e1 * t299;
t334 = t26 * qJD(1);
t80 = 0.2e1 * (t258 / 0.4e1 + t259 / 0.4e1 + 0.1e1 / 0.4e1) * t371;
t333 = t80 * qJD(1);
t331 = t169 * t375 + t173 * t378;
t327 = qJD(5) * t237;
t326 = qJD(5) * t261;
t325 = t395 / 0.2e1;
t320 = t391 + t392;
t319 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t318 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t316 = m(7) * t366;
t310 = t339 / 0.2e1;
t309 = -t338 / 0.2e1;
t253 = Ifges(7,5) * t262;
t254 = Ifges(6,5) * t262;
t303 = t254 / 0.4e1 + t253 / 0.4e1;
t300 = mrSges(6,1) + t317;
t298 = t316 / 0.2e1;
t294 = t261 * mrSges(6,1) + t262 * mrSges(6,2);
t293 = t351 + t357;
t292 = -Ifges(6,1) * t262 + t363;
t291 = -Ifges(7,1) * t262 + t362;
t290 = Ifges(6,2) * t261 - t256;
t289 = Ifges(7,2) * t261 - t255;
t115 = -Ifges(7,6) * t237 + t235 * t289;
t116 = Ifges(7,6) * t235 + t237 * t289;
t117 = -Ifges(6,6) * t237 + t235 * t290;
t118 = Ifges(6,6) * t235 + t237 * t290;
t119 = -Ifges(7,5) * t237 + t235 * t291;
t120 = Ifges(7,5) * t235 + t237 * t291;
t121 = -Ifges(6,5) * t237 + t235 * t292;
t122 = Ifges(6,5) * t235 + t237 * t292;
t160 = t235 * t293;
t161 = t235 * t294;
t162 = t293 * t237;
t163 = t294 * t237;
t274 = mrSges(4,1) * t373 + mrSges(4,2) * t374;
t1 = t315 * t274 - t154 * t161 - t87 * t162 - t403 * t163 + t49 * t167 + t63 * t168 + t52 * t169 + t68 * t170 + t41 * t171 + t62 * t172 + t44 * t173 + t67 * t174 - t88 * t160 + m(6) * (t154 * t403 + t62 * t67 + t63 * t68) + m(7) * (t41 * t44 + t49 * t52 + t87 * t88) + (mrSges(5,1) * t324 + Ifges(5,4) * t235 + (-t120 / 0.2e1 - t122 / 0.2e1 - t319 * t235) * t262 + (t116 / 0.2e1 + t118 / 0.2e1 + t318 * t235) * t261) * t235 + (-mrSges(5,2) * t324 - Ifges(5,4) * t237 + (-t119 / 0.2e1 - t121 / 0.2e1 + t319 * t237) * t262 + (t115 / 0.2e1 + t117 / 0.2e1 - t318 * t237) * t261 + (Ifges(5,1) - Ifges(5,2) - t407) * t235) * t237 + (-Ifges(4,2) + Ifges(4,1)) * t374 * t373 + (-t373 ^ 2 + t374 ^ 2) * Ifges(4,4) + (-t237 * mrSges(5,1) - t219 + t417) * t239;
t283 = t261 * t62 - t348;
t284 = t261 * t41 - t350;
t6 = (t163 / 0.2e1 + t162 / 0.2e1 + t308 * t262 + t306 * t261) * t237 + (-t161 / 0.2e1 - t160 / 0.2e1 + t307 * t262 + t305 * t261) * t235 + 0.2e1 * ((-t88 - t402) * t396 + (-t154 - t282) * t399) * t237 + 0.2e1 * ((t284 + t87) * t396 + (t403 + t283) * t399) * t235;
t286 = t1 * qJD(1) + t6 * qJD(2);
t5 = -t48 * t169 + t41 * t322 - pkin(5) * t162 * t336 - t88 * (-mrSges(7,1) * t336 + t215) - t241 * t341 + t63 * t174 + (t404 * t375 * t237 - mrSges(6,3) * t348 - mrSges(7,3) * t350 + t368 * t387 + t414 * t383 + (t237 * t405 + t120 + t122) * t378 - (t116 + t118) * t262 / 0.2e1) * t237 + (-t170 + t323) * t62 + (t173 - t316) * t49;
t285 = -t5 * qJD(1) - t7 * qJD(2);
t25 = 0.4e1 * t299;
t281 = t6 * qJD(1) + t25 * qJD(2);
t8 = (t162 + t163 + t364) * t237 + (t365 + (-t169 - t170) * t262 + (t173 + t174) * t261) * t235 + m(7) * (t235 * t284 - t88 * t237) + m(6) * (t235 * t283 - t341) + m(5) * (-t235 * t403 - t341);
t280 = qJD(1) * t8 + qJD(2) * t26;
t225 = t261 * t317 + t351;
t89 = -t317 * t336 + t215;
t279 = qJD(1) * t89 + qJD(3) * t225;
t278 = mrSges(7,1) / 0.2e1 + mrSges(6,1) / 0.2e1 + t325;
t277 = t298 - t173 / 0.2e1;
t276 = -Ifges(6,2) / 0.4e1 - Ifges(7,2) / 0.4e1 + Ifges(6,1) / 0.4e1 + Ifges(7,1) / 0.4e1;
t11 = (t235 * t320 + t307) * t262 + (t235 * t278 - t277 + t384) * t261;
t275 = t11 * qJD(1);
t123 = m(7) * t406 + t328 * mrSges(7,3);
t268 = ((-t221 * t262 + t223 * t261) * t237 - t284) * t397 + t331;
t19 = (t351 / 0.2e1 + t357 / 0.2e1) * t235 - t388 / 0.2e1 + t268;
t273 = qJD(1) * t19 - qJD(2) * t413 + qJD(3) * t123;
t271 = mrSges(7,1) * t390 + t391 * t52 + t392 * t68 + t394 * t67;
t24 = -t238 * t293 - t252 * t294 + (-t256 / 0.2e1 - t244 / 0.2e1 - t255 / 0.2e1 - t245 / 0.2e1) * t262 + (t243 / 0.2e1 + t242 / 0.2e1 + (Ifges(7,4) / 0.2e1 + Ifges(6,4) / 0.2e1) * t261 + (Ifges(6,2) / 0.2e1 - Ifges(7,1) / 0.2e1 + Ifges(7,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t262 + (-t240 - t370) * pkin(5)) * t261;
t263 = t250 * t295 + (t245 / 0.4e1 + t244 / 0.4e1 + t256 / 0.4e1 + t255 / 0.4e1 + t221 * t411 + t252 * t392 + t276 * t261) * t261 + (t243 / 0.4e1 + t242 / 0.4e1 + t223 * t411 + t238 * t393 + t252 * t394 - t276 * t262 + (0.3e1 / 0.4e1 * Ifges(6,4) + 0.3e1 / 0.4e1 * Ifges(7,4)) * t261 + (-t370 / 0.2e1 + t380) * pkin(5)) * t262;
t264 = t277 * t223 + (t174 * t379 + mrSges(6,2) * t386 + t122 / 0.4e1 + t120 / 0.4e1 + mrSges(7,2) * t389 + (t48 / 0.2e1 - t41 / 0.2e1) * mrSges(7,3)) * t262 + t221 * t410 + t238 * t215 / 0.2e1;
t267 = t170 * t379 + mrSges(6,1) * t386 + mrSges(7,1) * t389 - t118 / 0.4e1 - t116 / 0.4e1 + (t387 / 0.2e1 - t162 / 0.2e1) * pkin(5);
t3 = (m(7) * t390 - t171 / 0.2e1) * pkin(5) + t267 * t261 + (t319 * t262 + (-0.3e1 / 0.4e1 * Ifges(6,6) - 0.3e1 / 0.4e1 * Ifges(7,6)) * t261 + t303) * t235 + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + t263) * t237 + t264 + t271;
t270 = -t3 * qJD(1) + t24 * qJD(3);
t269 = (t261 * t278 + t262 * t320) * t235;
t82 = t212 - t296 / 0.2e1;
t28 = t310 * t395 + t269 + (t293 + t294) * t235 / 0.2e1;
t18 = t388 / 0.2e1 + mrSges(7,2) * t309 + t339 * t393 + t268;
t12 = t170 * t375 + t174 * t378 + t261 * t298 + t269 + t331;
t10 = t304 * t237 + t265 - t266 + (t167 + t168) * t377 + (t171 + t172) * t375;
t4 = t44 * t325 + pkin(5) * t385 + t303 * t235 + ((-Ifges(6,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t235 + t267) * t261 + t263 * t237 + t264 - t271 + t407 * t381 + t367 * t310 + t408 * t309;
t2 = qJD(3) * t6 + qJD(4) * t26 - qJD(5) * t7;
t13 = [qJD(3) * t1 + qJD(4) * t8 - qJD(5) * t5 + qJD(6) * t21, t2 ((-t344 * t412 + mrSges(5,2)) * t154 + t311 * t364 + t313 * t365 + mrSges(4,2) * t312 + t168 * t335 + (t121 + t119) * t377 + (t117 + t115) * t375 + m(7) * (-t221 * t44 + t223 * t52 + t238 * t87) + (m(6) * t282 - t261 * t172) * t250 + t404 * t309 + t405 * t310 + t402 * mrSges(7,3) + t282 * mrSges(6,3) + (m(6) * t252 - t345 * t412 - t415) * t403 + t414 * t381 - t252 * t161 - t238 * t160 + t87 * t240 - Ifges(5,5) * t235 + Ifges(5,6) * t237 - t221 * t171 + t223 * t167 - Ifges(4,6) * t373 + Ifges(4,5) * t374 - mrSges(4,1) * t314) * qJD(3) + t10 * qJD(4) + t4 * qJD(5) + t18 * qJD(6) + t286, qJD(3) * t10 + qJD(5) * t12 + t280 + t416, t4 * qJD(3) + t12 * qJD(4) + (-t63 * mrSges(6,1) - t62 * mrSges(6,2) - t48 * mrSges(7,2) - t317 * t49) * qJD(5) + (t418 + (-mrSges(7,3) * pkin(5) + t408) * t261) * t327 + t285, qJD(3) * t18 + qJD(4) * t413 + t343; t2, t25 * qJD(3), t28 * qJD(5) + t82 * qJD(6) + t281 + (t219 - t274 + (-t240 + t415) * t237 + 0.2e1 * t272 + t328 * t401) * qJD(3), t334, -t346 + t28 * qJD(3) - t215 * qJD(5) + (t262 * t300 - t356) * t327, t82 * qJD(3); qJD(4) * t9 + qJD(5) * t3 + qJD(6) * t19 - t286, -t281 - t416, -qJD(5) * t24 + qJD(6) * t123, t361 (-mrSges(6,1) * t335 + t221 * mrSges(7,2) - mrSges(7,3) * t368 - t223 * t317 + t253 + t254) * qJD(5) + (mrSges(6,2) * t250 - t367) * t326 - t270, t273; -qJD(3) * t9 - qJD(5) * t11 + qJD(6) * t80 - t280, -t334, -t361, 0 (-mrSges(6,2) - mrSges(7,2)) * qJD(5) * t262 - t300 * t326 - t275, t333; -qJD(3) * t3 + qJD(4) * t11 - qJD(6) * t89 - t285, t346, -t225 * qJD(6) + t270, t275, 0, -t279; -qJD(3) * t19 - qJD(4) * t80 + qJD(5) * t89 - t343, t413 * qJD(3), qJD(5) * t225 - t273, -t333, t279, 0;];
Cq  = t13;
